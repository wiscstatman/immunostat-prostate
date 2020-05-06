library(tidyverse) # make sure you have the latest tidyverse !
library(lme4) # linear mixed effects model
library(fdrtool)
library(Rtsne)
library(heatmap3)
library(ggplot2)

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


# typical step in ANOVA
anova_func <- function(anova_pval){
  # control FDR
  anova_BH <- p.adjust(anova_pval, method = "BH")
  anova_qval <- fdrtool(anova_pval, statistic = "pvalue", verbose = F, plot  = F)$qval
  anova_qval_eta0 <- unname(fdrtool(anova_pval, statistic = "pvalue", verbose = F, plot  = F)$param[,"eta0"])
  
  # plot histogram of p-values
  all_peptide_hist <- hist(anova_pval, breaks = 70, freq = F, xlab = "one-way ANOVA p-values", 
                           main = paste0("p-values distribution for ", length(anova_pval), " peptides"))
  polygon_ind <- which(all_peptide_hist$density >= anova_qval_eta0)
  for (i in polygon_ind){
    polygon( x = c(all_peptide_hist$breaks[i], all_peptide_hist$breaks[i+1], all_peptide_hist$breaks[i+1], all_peptide_hist$breaks[i]),
             y = c(all_peptide_hist$density[i], all_peptide_hist$density[i], anova_qval_eta0, anova_qval_eta0),
             col = "red")
  }
  text(x=0.65,y=4, labels = paste( "estimated proportion of \nnon-null peptides =",
                                   round( 100*(1 - anova_qval_eta0),2 ),"%" ))
  
  return(list(anova_BH = anova_BH, anova_qval = anova_qval, anova_qval_eta0 = anova_qval_eta0))
}


load("04_aggregated_calls.RData")
raw_data_median <- read_csv("raw_data_median.csv")

####################################################################################### 
#                                      Data Processing                                #
####################################################################################### 

array_id_key = read_csv("sample_key_project1.csv") %>%
  janitor::clean_names() %>% 
  rename(stage = condition,
         array_id = "file_name") %>%
  mutate(id = str_replace_all(id, " ", ""),
         id = str_to_lower(id),
         id = str_replace_all(id, "/", ""),
         stage = as_factor(stage),
         stage = fct_recode(stage,
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
patient_key = array_id_key %>%
  group_by(id, stage) %>%
  summarize(n = n()) %>%
  ungroup() %>% 
  filter(n >= 2) %>%
  select(-n)
array_id_key = array_id_key %>%
  filter(id %in% patient_key$id)

# check patient counts (NOT distinct patients)
array_id_key %>%
  group_by(id, stage) %>%
  summarize() %>%
  group_by(stage) %>%
  tally()

# there are patients who were measured at two different stages
# To ensure unique patients, remove the following ids

ids_to_remove = c("adt181",
                  "adt223",
                  "pdv008",
                  "pap123",
                  "pap067",
                  "adt143")

# drop patients' earlier records
array_id_key = array_id_key %>%
  filter(!(id %in% ids_to_remove))
patient_key = array_id_key %>%
  group_by(id, stage) %>%
  tally() %>%
  select(-n) 
patient_key$stage <- relevel(patient_key$stage, ref = "normal")

# check patient counts (distinct patients)
array_id_key %>%
  group_by(id, stage) %>%
  summarize() %>%
  group_by(stage) %>%
  tally()

#-----------------------------------------------------------------------------------------------
# read calls data

# read aggregated calls data
calls = read_csv("aggregated_calls_full_nmcspc.csv") 
length(unique(calls$probe_sequence)) == nrow(calls) # probe_sequence NOT unique

# probe_sequence NOT unique
# need to generate unique PROBE_ID
# PROBE_ID in raw_data_complete.csv is paste0(SEQ_ID, ";", POSITION)
calls$PROBE_ID = paste0(calls$seq_id, ";", calls$position)
calls <- calls %>% select(container:position, PROBE_ID, everything()) # rearrange columns

# keep only patients that appear in patient_sample_key
calls <- calls %>% select(container:PROBE_ID, any_of(patient_key$id))

# remove peptides that have zero calls in ALL subjects
calls <- calls[ apply( calls %>% select(any_of(patient_key$id)) , 1, function(x){ !all(x==0) } ) , ]

# see row sums of calls
calls_rowsums <- calls %>% select(any_of(patient_key$id)) %>% rowSums()
calls_rowsums_df <- as.data.frame(table(calls_rowsums)) %>%
  mutate(Freq = ifelse(as.numeric(calls_rowsums) >=10, sum(Freq[as.numeric(calls_rowsums) >=10]), Freq))  %>%
  filter(as.numeric(calls_rowsums) <=10) %>%
  mutate(calls_rowsums = ifelse(as.numeric(calls_rowsums) == 10 , "10 and above", calls_rowsums),
         calls_rowsums = factor(calls_rowsums, levels = c(1:9,"10 and above")))
ggplot(data=calls_rowsums_df, aes(x=calls_rowsums, y=Freq)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=paste0(Freq)), vjust=-.3, size=3.5)+
  labs(x = "Sum of calls among all patients for a peptide",
       title = paste0("Sum of Calls Among All Patients for Each Peptide (Total ", nrow(calls), " Peptides)")) +
  theme(plot.title = element_text(hjust = 0.5))


# in the end, how many peptides with at least one call among all patients
nrow(calls)

# plots for calls 
# call_counts_by_peptide <- t(apply( calls %>% select(any_of(sample_key$id)), 1, function(x){
#   tapply(x, sample_key$stage, FUN=sum)
# }))
# colnames(call_counts_by_peptide) <- levels(sample_key$stage)
# call_counts_by_peptide <- as.data.frame(call_counts_by_peptide)
# hist(call_counts_by_peptide$normal)
# hist(call_counts_by_peptide$new_dx)
# hist(call_counts_by_peptide$nmCSPC)
# hist(call_counts_by_peptide$nmCRPC)
# hist(call_counts_by_peptide$mCRPC)
## a lot of zeroes anyway...

#-----------------------------------------------------------------------------------------------
# get raw_data_complete.csv and compute median    

raw_data = read_csv("raw_data_complete.csv")

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

# extract all sequence id
# may need to use this for gene set analysis
all_seq_id <-  unique(raw_data$SEQ_ID)


# take log2 transformation
raw_data <- raw_data %>%
  mutate_at(vars(matches("dat")), log2)


# make sure array_id_key align with raw_data_complete
array_iii <- match( colnames( raw_data %>% select(any_of(array_id_key$array_id)) ), array_id_key$array_id )
array_id_key <- array_id_key[array_iii , ]
sum(as.numeric( colnames( raw_data %>% select(any_of(array_id_key$array_id)) ) == array_id_key$array_id )) == nrow(array_id_key)


# compute median
raw_data_median <- t( apply( select(raw_data, contains("dat")), 1, function(x) {
  tapply(x, array_id_key$id, FUN=median)
} ) )
raw_data_median <- bind_cols( select(raw_data, -contains("dat")), as.data.frame(raw_data_median) ) 

colnames(raw_data_median)[1:15] # check

# write.table(raw_data_median, file = "raw_data_median.csv", sep = ",", row.names = F)


####################################################################################### 
#      Evaluate Reproducibility of Replicates via Linear Mixed Effects Model          #         
#######################################################################################

ncol_raw <- ncol(raw_data) 
nrep <- nrow(array_id_key)

# initiate
lmer_result <- matrix(NA, nrow = nrow(raw_data), ncol = 4)
colnames(lmer_result) <- c("variance_id", "variance_residual", "lrstat", "singularity")

# check array_id_key align with raw_data_complete
sum(as.numeric( colnames( raw_data %>% select(any_of(array_id_key$array_id)) ) == array_id_key$array_id )) == nrow(array_id_key)


for(i in 1:nrow(raw_data)){
  y <- as.numeric(raw_data[i, (ncol_raw - nrep + 1) : ncol_raw])
  fit1 <- lmer(y ~ stage + (1|id), data = array_id_key)
  fit2 <- lmer(y ~ 1 + (1|id), data = array_id_key)
  lmer_result[i,] <- c(
    as.data.frame(VarCorr(fit1))$'vcov',
    as.numeric(-2*(logLik(fit2, REML=T) - logLik(fit1, REML=T))),
    ( isSingular(fit1) | isSingular(fit2) )
  )
  print(i)
}

# check how many singular fits
sum(as.numeric(lmer_result[,'singularity']))
max(as.numeric(lmer_result[,'variance_id'][ lmer_result[,'singularity']==T ]))
min(as.numeric(lmer_result[,'variance_residual'] ))

# get estimated proportion of variances
lmer_var_ratio <- lmer_result[,'variance_id'] / ( lmer_result[,'variance_id'] + lmer_result[,'variance_residual']  )
hist(lmer_var_ratio, breaks = 100, xlab = "estimated proportion of variances",
     main = "Histogram of peptide-level proportion of random-effect variance to total variance")


####################################################################################### 
#                            TEST !!!! -- logistic regression                         #
####################################################################################### 

# make sure patients' stages align with calls
calls_iii <- match( colnames( calls %>% select(any_of(patient_key$id)) ), patient_key$id )
calls_stages <- (patient_key$stage)[calls_iii]
sum(as.numeric( colnames( calls %>% select(any_of(patient_key$id)) ) == (patient_key$id)[calls_iii] )) == nrow(patient_key)

# initiate
logreg_pval <- rep(NA, nrow(calls))
names(logreg_pval) <- calls$PROBE_ID
ncol_calls <- ncol(calls)
n <- nrow(patient_key)

# compute deviance test p-values
for(i in 1:nrow(calls)){
  y <- as.numeric( calls[i, (ncol_calls - n + 1): ncol_calls] )
  fit1 <- glm(y ~ calls_stages, family = binomial(link = "logit"))
  logreg_pval[i] <- 1 - pchisq( fit1$null.deviance - fit1$deviance, df = (fit1$df.null - fit1$df.residual) )
  print(i)
}

# control FDR
logreg_BH <- p.adjust(logreg_pval, method = "BH")
logreg_qval <- fdrtool(logreg_pval, statistic = "pvalue", verbose = F, plot  = F)$qval
logreg_qval_eta0 <- unname(fdrtool(logreg_pval, statistic = "pvalue", verbose = F, plot  = F)$param[,"eta0"])

# plot histogram of p-values
hist(logreg_pval, breaks = 50, freq = T, main = "Logistic Regression Deviance Test p-values", xlab = "p-values ")

# peptide counts at various FDR thresholds
count.func(logreg_BH, seq(0.01, 0.05, by = 0.01))
count.func(logreg_qval, seq(0.01, 0.05, by = 0.01))

# Calls are conservative
# Logistic regression based on number of calls yield (almost) no signal even before FDR control


####################################################################################### 
#                           Another test -- one-way ANOVA                             #
####################################################################################### 

# make sure patients' stages align with raw_data_median
median_iii <- match( colnames(select( raw_data_median, any_of(array_id_key$id) )), patient_key$id )
median_stage <- patient_key$stage[median_iii]
sum(as.numeric( colnames(select( raw_data_median, any_of(array_id_key$id) )) == patient_key$id[median_iii] )) == nrow(patient_key)

# initiate
all_anova_pval <- rep(NA, nrow(raw_data_median))
names(all_anova_pval) <- raw_data_median$PROBE_ID
ncol_median <- ncol(raw_data_median)
n <- nrow(patient_key)

# compute one-way anova p-values
for(i in 1:nrow(raw_data_median)){
  all_anova_pval[i] <- unname( unlist(summary( aov( 
    as.numeric(raw_data_median[i, (ncol_median - n + 1) : ncol_median]) ~ median_stage ) ))["Pr(>F)1"] )
  print(i)
}

# get ANOVA p-values for only peptides with at least one non-zero call among all patients
calls_anova_pval <- all_anova_pval[names(all_anova_pval) %in% calls$PROBE_ID]

# get p-values histogram and FDR 
all_anova <- anova_func(all_anova_pval)
calls_anova <- anova_func(calls_anova_pval)

# peptide counts at various FDR thresholds
count.func(all_anova$anova_BH, seq(0.01, 0.1, by = 0.01))
count.func(all_anova$anova_qval, seq(0.01, 0.1, by = 0.01))
count.func(calls_anova$anova_BH, seq(0.01, 0.1, by = 0.01))
count.func(calls_anova$anova_qval, seq(0.01, 0.1, by = 0.01))

# comparing with calls
names(calls_rowsums) <- calls$PROBE_ID
sum(as.numeric(names(calls_anova$anova_BH[calls_anova$anova_BH<=0.05]) %in% names(calls_rowsums[calls_rowsums >=3])))
sum(as.numeric(names(all_anova$anova_BH[all_anova$anova_BH<=0.05]) %in% names(calls_rowsums[calls_rowsums >=3])))

####################################################################################### 
#                                ANOVA Visualization                                  #
#######################################################################################

BH_FDR_cutoff <- 0.05

# all_anova
anova_dat <- raw_data_median %>% 
  select(PROBE_ID, any_of(patient_key$id)) %>% 
  mutate(anova_BH = all_anova$anova_BH[raw_data_median$PROBE_ID]) %>%
  filter(anova_BH <= BH_FDR_cutoff) %>%
  select(-anova_BH)

# calls_anova
anova_dat <- raw_data_median %>%
  filter((PROBE_ID %in% calls$PROBE_ID)) %>%
  select(PROBE_ID, any_of(patient_key$id)) 
anova_dat <- anova_dat %>%
  mutate(anova_BH = calls_anova$anova_BH[anova_dat$PROBE_ID]) %>%
  filter(anova_BH <= BH_FDR_cutoff) %>%
  select(-anova_BH)


anova_dat_demean <- sweep(as.matrix(anova_dat %>% select(-PROBE_ID)), 1, 
                          rowMeans(as.matrix(anova_dat %>% select(-PROBE_ID))), "-") # centering by row

# make sure stage aligns with anova_dat_demean
visual_iii <- match( colnames(anova_dat_demean) , patient_key$id )
visual_stage <- patient_key$stage[visual_iii]

# colors and shapes for the visualization techniques
cols = pal[ match(visual_stage, names(pal)) ]
shapes = shp[  match(visual_stage, names(shp)) ]

#---------------------------------------------------------------------------------------
# PCA

# svd
sv.dat <- sweep(t(anova_dat_demean), 2, colMeans(t(anova_dat_demean)), "-") # centering
sv <- svd(sv.dat)
V <- sv$v
D <- sv$d
U <- sv$u

# variance explained
pca.var <- D^2/sum(D^2) 
pca.cumvar <- cumsum(pca.var)

# plot PCA
par(mfrow = c(1,2), pty = "s", mar = c(2.2,2.3,1.5,0.45), mgp = c(1.4,0.4,0),
    cex.axis = 0.84, cex.lab = 0.84, cex.main = 0.84, tcl = -0.4)
PCload.func(cols, shapes, U, D, 1, 2, pca.var, title = "PC2 vs PC1") # PC loadings (PC2 vs PC1)
legend('topright', pch = shp, col = pal, cex = 0.5,
       c("normal",  "new_dx", "nmCSPC", "mCSPC", "nmCRPC", "mCRPC") )
PCload.func(cols, shapes, U, D, 3, 2, pca.var, title = "PC2 vs PC3") # PC loadings (PC2 vs PC3)

dev.off()

#---------------------------------------------------------------------------------------
# t-SNE

# how to specify parameter
# refer: https://lvdmaaten.github.io/tsne/User_guide.pdf
tsnedat <- unname(t(anova_dat_demean)) # dim(X) = N samples by D dimensions 
initdim <- 90 
perplex <- 30 

set.seed(10)
tsne_anova <- Rtsne(tsnedat, initial_dims = initdim, perplexity = perplex,
                    theta = 0, check_duplicates = F, max_iter = 50000L, eta = 50)

# t-SNE plot
plot(tsne_anova$Y, ylab = "", xlab = "", col = cols, pch = shapes, main = "t-SNE plot")
legend('topleft', pch = shp, col = pal,
       c("normal",  "new_dx", "nmCSPC", "mCSPC", "nmCRPC", "mCRPC") )

dev.off()

#---------------------------------------------------------------------------------------
# Heatmap

# heatmap color palette
cls <- colorRampPalette(c("navy", "honeydew", "firebrick3", "brown"))(n = 1024)

heatmap3(anova_dat_demean, 
         col = cls, # specify colors 
         ColSideColors = cols, # specify patient color code
         labCol = visual_stage, # specify patient
         ColSideLabs = "stages", 
         labRow = "",
         xlab = "Patients",
         # legendfun=function() showLegend(col = c("navy", "cornflowerblue", "turquoise1", "orchid1", 
         #                                         "darkorange1", "firebrick1"),
         #                                 legend = c("normal",  "new_dx", "nmCSPC", "mCSPC", "nmCRPC", "mCRPC"),
         #                                 cex = 1.2,
         #                                 lwd = 5  )
)


####################################################################################### 
#                                      Linear Contrast                                #
#######################################################################################

# let's compare healthy subjects vs cancer patients!

# again, make sure stage aligns with anova_dat
contr_iii <- match( colnames(anova_dat %>% select(-PROBE_ID)) , patient_key$id )
contr_stage <- patient_key$stage[contr_iii]

# get group means and MSE
BH_FDR_normal <- apply(as.matrix(anova_dat %>% select(-PROBE_ID)), 1, 
                       function(x){mean(x[contr_stage == "normal"])})
BH_FDR_cancer <- apply(as.matrix(anova_dat %>% select(-PROBE_ID)), 1, 
                       function(x){mean(x[contr_stage != "normal"])})
BH_FDR_MSE <- apply(as.matrix(anova_dat %>% select(-PROBE_ID)), 1, function(x){
  fit1 <- lm(x ~ contr_stage)
  return(deviance(fit1)/df.residual(fit1))
})

# get linear contrast p-values
normal_counts <- length(contr_stage[contr_stage == "normal"])
cancer_counts <- length(contr_stage[contr_stage != "normal"])
BH_FDR_tstat <- (BH_FDR_cancer - BH_FDR_normal) / sqrt(BH_FDR_MSE * (1/cancer_counts + 1/normal_counts))
BH_FDR_tstat_pval <- 1 - pt( BH_FDR_tstat, df = (length(unique(patient_key$stage)) - 1) )

# plot linear contrast p-values
hist(BH_FDR_tstat_pval, breaks= 50)

# how many peptides with cancer patients significantly higher than normal subjects
length(BH_FDR_tstat_pval[BH_FDR_tstat_pval <= 0.05])

# output these probes and their t-stat p-val
contr_df <- data.frame( peptide_id = anova_dat$PROBE_ID[BH_FDR_tstat_pval <= 0.05],
                        difference = round((BH_FDR_cancer - BH_FDR_normal)[BH_FDR_tstat_pval <= 0.05],4),
                        tstat_pvalues = round(BH_FDR_tstat_pval[BH_FDR_tstat_pval <= 0.05],4))
contr_df <- contr_df[base::order(contr_df$tstat_pvalues),]

# boxplot specs
stage2 <- as.character(contr_stage)
stage2[stage2 != "normal"] <- "cancer"
stage2 <- factor(stage2)
boxplot_col <- c("red", "blue")
boxplot_func <- function(mat, col = boxplot_col, draw){
  # par(mar = c(2, 4.5, 2.3, 1),cex = 0.84)
  graphics::boxplot( as.numeric(mat[draw,]) ~ stage2, 
                     col=col, horizontal=TRUE, las=1, xlab = "log2(fluorescence)", ylab = "groups",
                     main= paste(row.names(mat)[draw]) )
}

# get boxplot matrix ready
boxplot_mat <- anova_dat[BH_FDR_tstat_pval <= 0.05,]
boxplot_mat <- as.matrix(boxplot_mat %>% select(-PROBE_ID))
row.names(boxplot_mat) <- anova_dat$PROBE_ID[BH_FDR_tstat_pval <= 0.05]
boxplot_mat <- boxplot_mat[order(BH_FDR_tstat_pval[BH_FDR_tstat_pval <= 0.05]),]

# boxplot
boxplot_func(mat = boxplot_mat, draw = 1)
boxplot_func(mat = boxplot_mat, draw = 2)
boxplot_func(mat = boxplot_mat, draw = 3)
boxplot_func(mat = boxplot_mat, draw = 4)
boxplot_func(mat = boxplot_mat, draw = 5)


####################################################################################### 
#                                   Gene Set Analysis                                 #
#######################################################################################

library(AnnotationDbi) # intraIDMapper: map uniprot id to gene entrez id
library(allez)

# read uniprot csv
uniprot_data = read_csv("complete_uniprot_data.csv", 
                        col_types = cols( 
                          seq_id = col_character(),
                          uniprot_id = col_character(),
                          status = col_character(),
                          protein_names = col_character(),
                          gene_names = col_character(),
                          length = col_double(),
                          yourlist_m201912135c475328cef75220c360d524e9d456ce648e46h = col_character(),
                          function_cc = col_character(),
                          subcellular_location_cc = col_character(),
                          tissue_specificity = col_character(),
                          involvement_in_disease = col_character(),
                          doug_name = col_character(),
                          query = col_character())
)
uniprot_data  <- uniprot_data %>% select(seq_id : query)

# get 5% BH-FDR sequence
BH_seq_id <- unique(raw_data_median$SEQ_ID[anova_BH <= BH_FDR_cutoff])

# get zero-one binary vector indicating whether protein sequence makes the BH-FDR cutoff
seq_id_ok <- as.numeric(all_seq_id %in% BH_seq_id)
length(seq_id_ok) == length(all_seq_id) # check, nice

# get all uniprot id based on protein seq id
all_seq_uniprot <- uniprot_data$uniprot_id[ match(all_seq_id, uniprot_data$seq_id) ]
sum(as.numeric(is.na(all_seq_uniprot))) == 0 # check, has some NA's

# trim away those NA's
seq_id_ok <- seq_id_ok[!(is.na(all_seq_uniprot))]
all_seq_uniprot <- all_seq_uniprot[!(is.na(all_seq_uniprot))]
length(unique(all_seq_uniprot)) == length(all_seq_uniprot) # NOT nice, there are repeats

# collapse those repeats
seq_id_ok_df <- data.frame(
  all_seq_uniprot,
  seq_id_ok
) 
seq_id_ok_df <- seq_id_ok_df %>%
  group_by(all_seq_uniprot) %>%
  summarise(ok2 = min(1, sum(seq_id_ok))) 
nrow(seq_id_ok_df) == length(unique(all_seq_uniprot)) # check, ok
table(seq_id_ok_df$ok2) # check, ok
all_seq_uniprot <- as.character(seq_id_ok_df$all_seq_uniprot)
length(all_seq_uniprot) == length(unique(all_seq_uniprot)) # all uniprot ids now unique and NO NA's
seq_id_ok <- seq_id_ok_df$ok2

# use bioconductor to match uniprot id to gene entrez id
all_seq_entrez <- AnnotationDbi::intraIDMapper(ids = all_seq_uniprot, species = 'HOMSA')
sum(as.numeric(is.na(all_seq_entrez))) == 0 # check, nice no NA's

# mappings got problem with the following checks
length(all_seq_entrez) == length(all_seq_uniprot) # NOT nice, has multiple matches
length(unique(all_seq_entrez)) == length(all_seq_entrez) # this is ok
sum(as.numeric( names(all_seq_entrez) %in% all_seq_uniprot )) == length(all_seq_entrez) # at least this is ok
data.frame(uniprot = names(all_seq_entrez)) %>% group_by(uniprot) %>% tally() %>% filter(n > 1)
# some uniprot ids are matched to multiple entrez ids, some not matched at all

# update binary vector to match entrez id
entrez_df <- data.frame(
  entrez = unname(all_seq_entrez),
  all_seq_uniprot = names(all_seq_entrez)
) %>% inner_join(seq_id_ok_df)
table(entrez_df$ok2) # check, ok
seq_id_ok <- entrez_df$ok2
names(seq_id_ok) <- entrez_df$entrez

# gene-set analysis via allez!
allez.go <- allez(seq_id_ok, lib = "org.Hs.eg", idtype = "ENTREZID", sets = "GO")
allez.kegg <- allez(seq_id_ok, lib = "org.Hs.eg", idtype = "ENTREZID", sets = "KEGG")

# let's see results
nom.alpha <- 0.05
min.num.gene <- 2

# Extract a table of top-ranked functional sets from allez output
allezTable(allez.go, symbol = T, n.cell = min.num.gene, nominal.alpha = nom.alpha)[,1:5]
allezTable(allez.kegg, symbol = T, n.cell = min.num.gene, nominal.alpha = nom.alpha)[,1:4]

# Display an image of gene scores by functional sets
allezPlot(allez.go, n.cell = min.num.gene, nominal.alpha = nom.alpha)
allezPlot(allez.kegg, n.cell = min.num.gene, nominal.alpha = nom.alpha)




# save results
# save(logreg_pval,
#      lmer_result,
#      all_anova_pval,
#      file = "04_aggregated_calls.RData")

