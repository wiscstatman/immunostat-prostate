library(tidyverse) # make sure you have the latest tidyverse !
library(lme4) # linear mixed effects model
library(fdrtool)

####################################################################################### 
#                           Some Common Variables and Functions                       #
####################################################################################### 

# specified color and shape scheme 
pal <- c("navy", "cornflowerblue", "turquoise1", "orchid1", "darkorange1", "firebrick1")
names(pal) <- c("normal",  "new_dx", "nmCSPC", "mCSPC", "nmCRPC", "mCRPC")
shp <- c(8, 15, 16, 3, 17, 18)
names(shp) <- names(pal)

# function to plot PCA loadings
PCload.func <- function(cols, shapes, U, D, x, y, pca.vec, title){
  Z <- U %*% diag(D)
  plot(Z[,x], Z[,y], col = cols, pch = shapes, main = title,
       xlab = paste0("PC",x," (", round(pca.vec[x]*100, 1) ,"% variance explained)"),
       ylab = paste0("PC",y," (", round(pca.vec[y]*100, 1) ,"% variance explained)") 
  )
}
## cols = color vector. 
## U = matrix U from SVD
## D = (vector) the diagonals of D matrix from SVD
## (x,y) = which PC loadings to plot
## pca.vec = pca % variance explained 

# function to count peptides at different FDR thresholds
count.func <- function(pval.vec, thresh.vec){
  counter <- NULL
  for (i in thresh.vec ){
    counter <- c(counter, length(pval.vec[pval.vec <= i]))
  }
  countab <- rbind(c("FDR threshold", thresh.vec),
                   c("Peptide counts", counter))
  return(countab)
}


####################################################################################### 
#                                      Data Processing                                #
####################################################################################### 

array_id_key = read_csv("sample_key_project1.csv") %>%
  janitor::clean_names() %>%  # unify all column names to snake_case
  rename(stage = condition,
         array_id = "file_name") %>%
  mutate(id = str_replace_all(id, " ", ""),
         id = str_to_lower(id),
         id = str_replace_all(id, "/", ""),
         stage = as_factor(stage),
         stage = fct_recode(stage,  # rename factor levels
                            "normal" = "Normal_male_controls",
                            "new_dx" = "Newly_diagnosed",
                            "nmCSPC" = "PSA-recurrent_nonMet",
                            "mCSPC" = "Met",
                            "nmCRPC" = "Castration-resistent_nonMet",
                            "mCRPC" = "Castration-resistent_Met",
                            "binding_buffer" = "Binding buffer alone"))

# drop binding buffer
length(array_id_key$stage[array_id_key$stage=="binding_buffer"]) # 1 binding_buffer
array_id_key <- array_id_key[!(array_id_key$stage=="binding_buffer"),]
array_id_key$stage <- factor(array_id_key$stage) # remove binding_buffer level

# remove patients whose rep == 1
sample_key = array_id_key %>%
  group_by(id, stage) %>%
  summarize(n = n()) %>%  # assign (id,stage) group counts to the variable "n"
  ungroup() %>% 
  filter(n >= 2) %>%
  select(-n) # drop the column "n"
array_id_key = array_id_key %>%
  filter(id %in% sample_key$id)

# check patient counts (NOT distinct patients)
array_id_key %>%
  group_by(id, stage) %>%
  summarize() %>%
  group_by(stage) %>%
  tally()

#-----------------------------------------------------------------------------------------------
# Ensuring only distinct patients

# 11 patients with replicates in 2 stages
ids_to_remove = c("pdv004", #1
                  "arv003", #2
                  "pap092", #3
                  "pap032", #4
                  "adt034", #5
                  "adt181", #6
                  "adt223", #7
                  "pdv008", #8
                  "pap123", #9
                  "pap067", #10
                  "adt143") #11

# drop patients' earlier records
array_id_key = array_id_key %>%
  filter(!(id %in% ids_to_remove))

# check patient counts (finally distinct patients)
distinct_patients_with_reps <- array_id_key %>%
  group_by(id, stage) %>%
  summarize() %>%
  group_by(stage) %>%
  tally()
distinct_patients_with_reps

# how many distinct patients with replicates
n <- sum(distinct_patients_with_reps$n); n

#-----------------------------------------------------------------------------------------------
# Now read in fluorescence intensity

raw_data = read_csv("raw_data_complete.csv")

# some checks
sum(as.numeric(raw_data$X == raw_data$COL_NUM)) == nrow(raw_data) # X = COL_NUM
sum(as.numeric(raw_data$Y == raw_data$ROW_NUM)) == nrow(raw_data) # Y = ROW_NUM
sum(as.numeric(raw_data$MATCH_INDEX == raw_data$FEATURE_ID)) == nrow(raw_data) # MATCH_INDEX = FEATURE_ID
unique(raw_data$DESIGN_NOTE) # only NA
unique(raw_data$SELECTION_CRITERIA) # only NA
unique(raw_data$MISMATCH) # only 0
unique(raw_data$PROBE_CLASS) # only NA
# we can drop X, Y, MATCH_INDEX, DESIGN_NOTE, SELECTION_CRITERIA, MISMATCH, PROBE_CLASS

raw_data = raw_data %>%
  select(PROBE_DESIGN_ID:Y, any_of(array_id_key$array_id)) %>% # drop patients with records at different stages 
  select( -c(X, Y, MATCH_INDEX, DESIGN_NOTE, SELECTION_CRITERIA, MISMATCH, PROBE_CLASS) ) # drop unhelpful columns

colnames(raw_data)[1:15] # check


# arrange array_id_key according to order in raw_data
iii = match( colnames(raw_data %>% select(contains("dat"))), array_id_key$array_id )
array_id_key <- array_id_key[iii,]
sum(as.numeric(colnames(raw_data %>% select(contains("dat"))) == array_id_key$array_id)) == nrow(array_id_key) 


#-----------------------------------------------------------------------------------------------
# take log2 transformation

raw_data <- raw_data %>%
  mutate_at(vars(matches("dat")), log2)


####################################################################################### 
#      Evaluate Reproducibility of Replicates via Linear Mixed Effects Model          #         
#######################################################################################

numcol <- ncol(raw_data) ; numrep <- nrow(array_id_key)
lmer_result <- matrix(NA, nrow = nrow(raw_data), ncol = 4)
colnames(lmer_result) <- c("variance_id", "variance_residual", "lrstat", "singularity")

for(i in 1:nrow(raw_data)){
  y <- as.numeric(raw_data[i, (numcol - numrep + 1) : numcol])
  fit1 <- lmer(y ~ stage + (1|id), data = array_id_key)
  fit2 <- lmer(y ~ 1 + (1|id), data = array_id_key)
  lmer_result[i,] <- c(
    as.data.frame(VarCorr(fit1))$'vcov',
    as.numeric(-2*(logLik(fit2, REML=T) - logLik(fit1, REML=T))),
    ( isSingular(fit1) | isSingular(fit2) )
  )
  print(i)
}
save(lmer_result, file = "lmer_result.RData")


# check how many singular fits
sum(as.numeric(lmer_result[,'singularity']))
max(as.numeric(lmer_result[,'variance_id'][ lmer_result[,'singularity']==T ]))
min(as.numeric(lmer_result[,'variance_residual'] ))

# get estimated proportion of variances
lmer_var_ratio <- lmer_result[,'variance_id'] / ( lmer_result[,'variance_id'] + lmer_result[,'variance_residual']  )
hist(lmer_var_ratio, breaks = 100, xlab = "estimated proportion of variances",
     main = "Histogram of peptide-level proportion of random-effect variance to total variance")


# get p-val, BH-FDR, q-val

df <- ( array_id_key %>% group_by(stage) %>% tally() %>% nrow() ) - 1

lmer_anova_pval <- 1 - pchisq(lmer_result[,'lrstat'], df)
# lmer_anova_pval <- 1 - pchisq(lmer_result[,'lrstat'][ lmer_result[,'singularity']==F ], df) # exclude singular peptides

lmer_anova_BHpval <- p.adjust(lmer_anova_pval, method = "BH")
lmer_anova_qval <- fdrtool(lmer_anova_pval, statistic = "pvalue", verbose = F, plot  = F)$qval
lmer_qval_eta0 <- unname(fdrtool(lmer_anova_pval, statistic = "pvalue", verbose = F, plot  = F)$param[,"eta0"])


# tabulate BH-FDR and q-val
count.func(lmer_anova_BHpval, seq(0.01, 0.1, by = 0.01))
count.func(lmer_anova_qval, seq(0.01, 0.1, by = 0.01))


# plot p-val histogram
all_peptide_hist <- hist(lmer_anova_pval, breaks = 70, freq = F, 
                         xlab = "one-way ANOVA p-values", main = "p-values distribution for peptides")
polygon_ind <- which(all_peptide_hist$density >= lmer_qval_eta0)
for (i in polygon_ind){
  polygon( x = c(all_peptide_hist$breaks[i], all_peptide_hist$breaks[i+1], all_peptide_hist$breaks[i+1], all_peptide_hist$breaks[i]),
           y = c(all_peptide_hist$density[i], all_peptide_hist$density[i], lmer_qval_eta0, lmer_qval_eta0),
           col = "red")
}
text(x=0.65,y=4, labels = paste( "estimated proportion of \nnon-null peptides =",
                                 round( 100*(1 - lmer_qval_eta0),2 ),"%" ))



