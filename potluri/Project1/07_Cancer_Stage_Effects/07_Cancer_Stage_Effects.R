# library(AnnotationDbi) # intraIDMapper: map uniprot id to gene entrez id
library(allez)
library(lme4) # linear mixed effects model
library(fdrtool)
library(Rtsne)
library(heatmap3)
library(ggplot2)
library(matrixStats) #rowMedians
library(tidyverse) # make sure you have the latest tidyverse !
library(xlsx)

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


# commonly-used ggplot theme 
myTheme <- theme(panel.background = element_rect(fill = "grey90"),
                 panel.grid.major = element_line(color = "white"),
                 panel.grid.minor = element_line(color = "white"),
                 axis.text.x = element_text(angle = 90, hjust = 1),
                 legend.position = "bottom",
                 plot.title = element_text(hjust = 0.5))


load("07_Cancer_Stage_Effects.RData")
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


# check dimensions of calls match dimensions of raw_data_median
(nrow(calls) == nrow(raw_data_median)) & 
  (ncol(calls %>% select(any_of(patient_key$id))) == ncol(raw_data_median %>% select(any_of(patient_key$id))))
ncol(calls %>% select(any_of(patient_key$id))) == nrow(patient_key)

# check if PROBE_ID (& patient_id) in calls appear in PROBE_ID (& patient_id) in raw_data_median
sum(as.numeric( calls$PROBE_ID %in% raw_data_median$PROBE_ID )) == nrow(calls)
sum(as.numeric( colnames( calls %>% select(any_of(patient_key$id)) ) %in% 
                  colnames( raw_data_median %>% select(any_of(patient_key$id)) ) )) == nrow(patient_key)

# get median_long & calls_long
median_long2 <- raw_data_median %>%
  select(PROBE_ID, any_of(array_id_key$id)) 
calls_long <- calls %>%
  select(PROBE_ID, any_of(array_id_key$id)) 
dim(median_long2) == dim(calls_long) # check

# rearrange rows and columns to match calls_long & median_long
calls_long <- calls_long[, match(colnames(median_long2), colnames(calls_long))]
calls_long <- calls_long[ match(median_long2$PROBE_ID, calls_long$PROBE_ID) , ]
sum(as.numeric( colnames( calls_long ) == colnames( median_long2 ) )) == nrow(patient_key) + 1 # check
sum(as.numeric( calls_long$PROBE_ID ==  median_long2$PROBE_ID  )) == nrow(median_long2) # check

# pivot_longer
median_long2 <- median_long2 %>% select(-PROBE_ID) %>%
  pivot_longer(cols = everything(), names_to = "id", values_to = "fluorescence") 
calls_long <- calls_long %>% select(-PROBE_ID) %>%
  pivot_longer(cols = everything(), names_to = "id", values_to = "calls") 

# check median_long2 & calls_long
nrow(median_long2) == nrow(raw_data_median) * nrow(patient_key)
head(median_long2)
nrow(calls_long) == nrow(calls) * nrow(patient_key)
head(calls_long)
dim(median_long2) == dim(median_long2)
sum(as.numeric(median_long2$id == calls_long$id)) == nrow(raw_data_median) * nrow(patient_key)

# plot fluorescence of calls vs no-calls
calls_fl_df <- data.frame(
  calls = factor(calls_long$calls),
  fluorescence = median_long2$fluorescence
)
janitor::tabyl(calls_fl_df$calls) %>% janitor::adorn_pct_formatting() %>%
  rename(calls = "calls_fl_df$calls", patient_peptide_counts = n)
ggplot(calls_fl_df, aes(x = calls, y = fluorescence, fill = calls)) +
  geom_boxplot(outlier.shape = ".") +
  labs(title = "Boxplots of Fluorescence Levels per Peptide per Patient", 
       x = "Calls per peptide per patient", y = "Median (across replicates) Fluorescence Levels on log2 scale") +
  theme(panel.background = element_rect(fill = "grey90"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        plot.title = element_text(hjust = 0.5))

# check calls_long
nrow(calls_long) == nrow(calls) * nrow(patient_key)
head(calls_long)

# boxplots to check fluorescence of positive calls vs fluorescence of zero calls
calls_median_df <- inner_join(calls_long, median_long, by = "id")

# free up memory
rm(median_long2, calls_fl_df,  calls_long); gc()


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


####################################################################################### 
#                         Check normalization of fluorescence data                    #         
#######################################################################################

median_long <- raw_data_median %>%
  select(any_of(array_id_key$id)) %>%
  pivot_longer(cols = everything(), names_to = "id", values_to = "fluorescence") 

# check
nrow(median_long) == nrow(raw_data_median) * nrow(patient_key)
head(median_long)

# set fill color
median_long$stage <- patient_key$stage[ match(median_long$id, patient_key$id) ]
# median_long$cols <- pal[ match(median_long$stage, names(pal)) ]

# sort order of patients in boxplot
median_long$id <- factor(median_long$id, levels = c(
  patient_key$id[patient_key$stage == "normal"],
  patient_key$id[patient_key$stage == "new_dx"],
  patient_key$id[patient_key$stage == "nmCSPC"],
  patient_key$id[patient_key$stage == "nmCRPC"],
  patient_key$id[patient_key$stage == "mCRPC"]
))

ggplot(median_long, aes(x = id, y = fluorescence, fill = stage)) +
  geom_boxplot(outlier.shape = ".") +
  # scale_fill_discrete(name="Stage") +
  # guides(fill=guide_legend(title="Stage")) +
  scale_fill_manual(name = "Stage", values = pal) +
  labs(title = "Boxplots of Peptide Fluorescence Levels for All Patients", 
       x = "Patient ID", y = "Median Fluorescence Levels on log2 scale") +
  myTheme

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
#                 Another test -- one-way ANOVA & Kruskal-Wallis Test                 #
####################################################################################### 

# make sure patients' stages align with raw_data_median
median_iii <- match( colnames(select( raw_data_median, any_of(array_id_key$id) )), patient_key$id )
median_stage <- patient_key$stage[median_iii]
median_stage <- factor(median_stage)
sum(as.numeric( colnames(select( raw_data_median, any_of(array_id_key$id) )) == patient_key$id[median_iii] )) == nrow(patient_key)

ncol_median <- ncol(raw_data_median)
n <- nrow(patient_key)

#---------------------------------------------------------------------------------------

# initiate ANOVA
all_anova_pval <- rep(NA, nrow(raw_data_median))
names(all_anova_pval) <- raw_data_median$PROBE_ID
all_anova_mse <- rep(NA, nrow(raw_data_median))
names(all_anova_mse) <- raw_data_median$PROBE_ID

# compute one-way anova p-values
for(i in 1:nrow(raw_data_median)){
  fit1 <- lm( as.numeric(raw_data_median[i, (ncol_median - n + 1) : ncol_median]) ~ median_stage )
  all_anova_pval[i] <- unname( unlist(summary(aov(fit1)))["Pr(>F)1"] )
  all_anova_mse[i] <- deviance(fit1)/df.residual(fit1)
  print(i)
}

# get ANOVA p-values for only peptides with at least one non-zero call among all patients
# calls_anova_pval <- all_anova_pval[names(all_anova_pval) %in% calls$PROBE_ID] # filter then FDR control

# get p-values histogram and FDR 
all_anova <- anova_func(all_anova_pval)
# calls_anova <- anova_func(calls_anova_pval) # filter then FDR control

# peptide counts at various FDR thresholds
count.func(all_anova$anova_BH, seq(0.01, 0.1, by = 0.01))
count.func(all_anova$anova_qval, seq(0.01, 0.1, by = 0.01))
# count.func(calls_anova$anova_BH, seq(0.01, 0.1, by = 0.01))
# count.func(calls_anova$anova_qval, seq(0.01, 0.1, by = 0.01))

# comparing with calls
names(calls_rowsums) <- calls$PROBE_ID
sum(as.numeric(names(calls_anova$anova_BH[calls_anova$anova_BH<=0.05]) %in% names(calls_rowsums[calls_rowsums >=3])))
sum(as.numeric(names(all_anova$anova_BH[all_anova$anova_BH<=0.05]) %in% names(calls_rowsums[calls_rowsums >=3])))

#---------------------------------------------------------------------------------------
# now do kruskal-wallis tests

# initiate kruskal-wallis(KW)
all_kw_pval <- rep(NA, nrow(raw_data_median))
names(all_kw_pval) <- raw_data_median$PROBE_ID

for(i in 1:nrow(raw_data_median)){
  all_kw_pval[i] <- kruskal.test( as.numeric(raw_data_median[i, (ncol_median - n + 1) : ncol_median]) ~ median_stage )$'p.value' 
  print(i)
}

# get p-values histogram and FDR 
all_kw <- anova_func(all_kw_pval)

# peptide counts at various FDR thresholds
count.func(all_kw$anova_BH, seq(0.01, 0.1, by = 0.01))

# signif for both ANOVA & Kruskal-Wallis
length(which(all_kw$anova_BH <= .05 & all_anova$anova_BH <= .05))

#---------------------------------------------------------------------------------------
# compare kruskal-wallis pval vs ANOVA pval

plot(x = all_anova_pval, y = all_kw_pval, pch = ".", xlab = "ANOVA p-values", ylab = "Kruskal-Wallis p-values")
lines(x = all_anova_pval[all_anova$anova_BH <= .05], 
      y = all_kw_pval[all_anova$anova_BH <= .05], 
      type = "p", pch = 20, col = "red")
lines(x = all_anova_pval[all_kw$anova_BH <= .05], 
      y = all_kw_pval[all_kw$anova_BH <= .05], 
      type = "p", pch = 20, col = "blue")
lines(x = all_anova_pval[all_anova$anova_BH <= .05 & all_kw$anova_BH <= .05], 
      y = all_kw_pval[all_anova$anova_BH <= .05 & all_kw$anova_BH <= .05], 
      type = "p", pch = 20, col = "green")

# check boxplot of peptide with very signif ANOVA pval but not kruskal-wallis pval
ANOVA_but_not_KW <- raw_data_median %>% 
  filter(all_anova_pval <.001 & all_kw_pval >.2) %>%
  select(PROBE_ID, any_of(array_id_key$id)) %>%
  as.matrix()

par(mfrow=c(3,1))
par(mar=c(5.1, 6.1, 4.1, 2.1))
graphics::boxplot(as.numeric(ANOVA_but_not_KW[1,-1])~median_stage, varwidth = T, horizontal = T,las= 2,
                  col = c("navy", "cornflowerblue", "turquoise1","darkorange1", "firebrick1"),
                  main = ANOVA_but_not_KW[1,"PROBE_ID"], xlab = "log2(fluorescence)", ylab = "")
graphics::boxplot(as.numeric(ANOVA_but_not_KW[2,-1])~median_stage, varwidth = T, horizontal = T,las= 2,
                  col = c("navy", "cornflowerblue", "turquoise1","darkorange1", "firebrick1"),
                  main = ANOVA_but_not_KW[2,"PROBE_ID"], xlab = "log2(fluorescence)", ylab = "")
graphics::boxplot(as.numeric(ANOVA_but_not_KW[3,-1])~median_stage, varwidth = T, horizontal = T,las= 2,
                  col = c("navy", "cornflowerblue", "turquoise1","darkorange1", "firebrick1"),
                  main = ANOVA_but_not_KW[2,"PROBE_ID"], xlab = "log2(fluorescence)", ylab = "")


# check boxplot of peptide with very signif kruskal-wallis pval but not ANOVA pval
KW_but_not_ANOVA <- raw_data_median %>% 
  filter(all_anova_pval > .7 & all_kw_pval < .001) %>%
  select(PROBE_ID, any_of(array_id_key$id)) %>%
  as.matrix()

par(mfrow=c(3,1))
par(mar=c(5.1, 6.1, 4.1, 2.1))
graphics::boxplot(as.numeric(KW_but_not_ANOVA[1,-1])~median_stage, varwidth = T, horizontal = T,las= 2,
                  col = c("navy", "cornflowerblue", "turquoise1","darkorange1", "firebrick1"),
                  main = KW_but_not_ANOVA[1,"PROBE_ID"], xlab = "log2(fluorescence)", ylab = "")
graphics::boxplot(as.numeric(KW_but_not_ANOVA[2,-1])~median_stage, varwidth = T, horizontal = T,las= 2,
                  col = c("navy", "cornflowerblue", "turquoise1","darkorange1", "firebrick1"),
                  main = KW_but_not_ANOVA[2,"PROBE_ID"], xlab = "log2(fluorescence)", ylab = "")
graphics::boxplot(as.numeric(KW_but_not_ANOVA[3,-1])~median_stage, varwidth = T, horizontal = T,las= 2,
                  col = c("navy", "cornflowerblue", "turquoise1","darkorange1", "firebrick1"),
                  main = KW_but_not_ANOVA[2,"PROBE_ID"], xlab = "log2(fluorescence)", ylab = "")

####################################################################################### 
#                                ANOVA Visualization                                  #
#######################################################################################

BH_FDR_cutoff <- 0.05

# all_anova
anova_dat <- raw_data_median %>% 
  select(PROBE_ID, any_of(patient_key$id)) %>% 
  mutate(anova_BH = all_anova$anova_BH[raw_data_median$PROBE_ID]) %>%
  # mutate(anova_BH = all_kw$anova_BH[raw_data_median$PROBE_ID]) %>%
  filter(anova_BH <= BH_FDR_cutoff) %>%
  select(-anova_BH)

# calls_anova (filter then FDR control)
# anova_dat <- raw_data_median %>%
#   filter((PROBE_ID %in% calls$PROBE_ID)) %>%
#   select(PROBE_ID, any_of(patient_key$id)) 
# anova_dat <- anova_dat %>%
#   mutate(anova_BH = calls_anova$anova_BH[anova_dat$PROBE_ID]) %>%
#   filter(anova_BH <= BH_FDR_cutoff) %>%
#   select(-anova_BH)

# anova_dat <- raw_data_median %>% 
#   select(PROBE_ID, any_of(patient_key$id)) %>% 
#   mutate(anova_BH = all_anova$anova_BH[raw_data_median$PROBE_ID]) %>%
#   filter(anova_BH <= BH_FDR_cutoff) %>%
#   filter(PROBE_ID %in% calls$PROBE_ID) %>%
#   select(-anova_BH)


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

#---------------------------------------------------------------------------------------
# Column-Reordered Heatmap

get_column_order.func <- function(stages){
  heat_map <- heatmap3(anova_dat_demean[,visual_stage == stages],col = cls, labRow = "")
  return( colnames(anova_dat_demean[,visual_stage == stages])[heat_map$colInd] )
}

normal_id_order <- get_column_order.func("normal")
newdx_id_order <- get_column_order.func("new_dx")
nmCSPC_id_order <- get_column_order.func("nmCSPC")
nmCRPC_id_order <- get_column_order.func("nmCRPC")
mCRPC_id_order <- get_column_order.func("mCRPC")

id_order <- match( c(normal_id_order, newdx_id_order, nmCSPC_id_order, nmCRPC_id_order, mCRPC_id_order),
                   colnames(anova_dat_demean) )

heatmap3(anova_dat_demean[,id_order], 
         col = cls, # specify colors 
         ColSideColors = cols[id_order], # specify patient color code
         Colv = NA,
         labCol = visual_stage[id_order], # specify patient
         ColSideLabs = "stages", 
         labRow = "",
         xlab = "Patients",
         # legendfun=function() showLegend(col = c("navy", "cornflowerblue", "turquoise1", "darkorange1", "firebrick1"),
         #                                 legend = c("normal",  "new_dx", "nmCSPC", "nmCRPC", "mCRPC"),
         #                                 cex = 1.2,
         #                                 lwd = 5  )
)



####################################################################################### 
#                                   Gene Set Analysis                                 #
#######################################################################################

# read uniprot_gene csv
uniprot_gene <- read_csv("uniprot_data_entrez.csv", col_types = cols_only(
  seq_id = col_character(),
  uniprot_id = col_character(),
  gene_symbol = col_character(),
  entrez_gene_id_pete = col_double(),
  gene_names = col_character(),
  protein_names = col_character()
)) %>% filter(!(is.na(seq_id)) & !(is.na(entrez_gene_id_pete))) %>%
  select(seq_id, uniprot_id, gene_symbol, entrez_gene_id_pete, gene_names, protein_names)

# check if seq_id & entrez_id unique
length(unique(uniprot_gene$seq_id)) == length(uniprot_gene$seq_id) # yes! unique!
length(unique(uniprot_gene$entrez_gene_id_pete)) == length(uniprot_gene$entrez_gene_id_pete) # NOT unique

uniprot_gene[ uniprot_gene$entrez_gene_id_pete %in%
               (uniprot_gene %>%
                  group_by(entrez_gene_id_pete) %>% 
                  tally() %>%
                  filter(n > 1) %>%
                  pull(entrez_gene_id_pete)) ,]


## manually changing entrez_gene_id for PCA10 & PRO29
uniprot_gene$entrez_gene_id_pete[uniprot_gene$seq_id == "PCA10"] <- 28912
uniprot_gene$entrez_gene_id_pete[uniprot_gene$seq_id == "PRO29"] <- NA
uniprot_gene <- uniprot_gene[!(is.na(uniprot_gene$entrez_gene_id_pete)),]

# CAREFUL!
# further filter uniprot_gene to limit them to proteins associated with filtered peptides based on calls
# uniprot_gene  <- uniprot_gene[ uniprot_gene$seq_id %in% (unique(calls$seq_id)) ,]

# get unique seq_id that are associated with significant peptides
signif_seq_id <- unique( raw_data_median$SEQ_ID[raw_data_median$PROBE_ID %in% anova_dat$PROBE_ID] )
signif_seq_id <- signif_seq_id[!(is.na(signif_seq_id))] # just in case

# convert these into binary vector
seq_id_ok <- as.numeric(uniprot_gene$seq_id %in% signif_seq_id)
names(seq_id_ok) <- uniprot_gene$entrez_gene_id_pete

# gene-set analysis via allez!
allez.go <- allez(seq_id_ok, lib = "org.Hs.eg", idtype = "ENTREZID", sets = "GO")
allez.kegg <- allez(seq_id_ok, lib = "org.Hs.eg", idtype = "ENTREZID", sets = "KEGG")

# let's see results
nom.alpha <- 0.05
min.num.gene <- 2

# Extract a table of top-ranked functional sets from allez output
allezTable(allez.go, symbol = T, n.cell = min.num.gene, nominal.alpha = nom.alpha, in.set = T)[,c(1:5,7)]
allezTable(allez.kegg, symbol = T, n.cell = min.num.gene, nominal.alpha = nom.alpha)[,1:4]

# Display an image of gene scores by functional sets
allezPlot(allez.go, n.cell = min.num.gene, nominal.alpha = nom.alpha)
allezPlot(allez.kegg, n.cell = min.num.gene, nominal.alpha = nom.alpha)

# tabulate enriched gene-set with signif genes only
allez.tab <- allezTable(allez.go, symbol = T, n.cell = min.num.gene, nominal.alpha = nom.alpha, in.set = T)
allez.tab$set.size <- paste(allez.tab$in.set, allez.tab$set.size, sep = "/")
allez.tab <- allez.tab %>% dplyr::select(-c(in.set, genes)) %>%
  mutate(in.genes = str_replace_all(in.genes, ";", "; "))


####################################################################################### 
#                                     Post-Hoc Analysis                               #
#######################################################################################

# we want the following contrasts:
# consecutive-group comparison: mCRPC-nmCRPC, nmCRPC-nmCSPS,nmCSPC-new_dx, new_dx-normal
# normal vs canceer
# mCRPC vs the others

# first make sure stage aligns with raw_data_median
sum(as.numeric(colnames(raw_data_median %>% select(-(PROBE_DESIGN_ID:DESIGN_ID))) == patient_key$id)) == nrow(patient_key)

# next get group means
group_means.func <- function(group, BH_filter){
  raw_data_median %>% 
    select(any_of(array_id_key$id)) %>% 
    select(which(patient_key$stage %in% group)) %>%
    filter(all_anova$anova_BH <= BH_filter) %>%
    rowMeans() 
}
mCRPC_mean <- group_means.func("mCRPC", .05)
nmCRPC_mean <- group_means.func("nmCRPC", .05)
nmCSPC_mean <- group_means.func("nmCSPC", .05)
newdx_mean <- group_means.func("new_dx", .05)
normal_mean <- group_means.func("normal", .05)
cancer_mean <- group_means.func(c("new_dx", "nmCSPC", "nmCRPC", "mCRPC"), .05)
NOT_mCRPC_mean <- group_means.func(c("normal", "new_dx", "nmCSPC", "nmCRPC"), .05)

# get group medians as well
group_median.func <- function(group, BH_filter){
  raw_data_median %>% 
    select(any_of(array_id_key$id)) %>% 
    select(which(patient_key$stage %in% group)) %>%
    filter(all_anova$anova_BH <= BH_filter) %>%
    as.matrix() %>%
    matrixStats::rowMedians() 
}
mCRPC_median <- group_median.func("mCRPC", .05)
nmCRPC_median <- group_median.func("nmCRPC", .05)
nmCSPC_median <- group_median.func("nmCSPC", .05)
newdx_median <- group_median.func("new_dx", .05)
normal_median <- group_median.func("normal", .05)
cancer_median <- group_median.func(c("new_dx", "nmCSPC", "nmCRPC", "mCRPC"), .05)
NOT_mCRPC_median <- group_median.func(c("normal", "new_dx", "nmCSPC", "nmCRPC"), .05)


# post-hoc analysis function
# posthoc_func <- function(posthoc_pval, groupmean_diff, countfunc_seq=seq(.03,.1,by=.01), effectsize_thresh=1, BH_thresh=.05){
#   hist(posthoc_pval, breaks = 70)
#   BH <- p.adjust(posthoc_pval, method = "BH")
#   BH_threshold_counts <- count.func(BH, countfunc_seq)
#   signif_peptide_counts <- length(which(abs(groupmean_diff) > effectsize_thresh & BH <= BH_thresh))
# }

#----------------------------------------------------------------------------------------------
# wilcoxon rank-sum tests

median_subset <- raw_data_median %>%
  filter(all_anova$anova_BH <= .05) %>%
  select(any_of(array_id_key$id)) %>%
  as.matrix()

# initiate wilcox-pval
mCRPC_nmCRPC_wilcox_pval <- rep(NA, nrow(median_subset))
nmCRPC_nmCSPC_wilcox_pval <- rep(NA, nrow(median_subset))
nmCSPC_newdx_wilcox_pval <- rep(NA, nrow(median_subset))
newdx_normal_wilcox_pval <- rep(NA, nrow(median_subset))
cancer_normal_wilcox_pval <- rep(NA, nrow(median_subset))
mCRPC_others_wilcox_pval <- rep(NA, nrow(median_subset))

# mCRPC_normal_wilcox_pval <- rep(NA, nrow(median_subset))
# nmCRPC_normal_wilcox_pval <- rep(NA, nrow(median_subset))
# nmCSPC_normal_wilcox_pval <- rep(NA, nrow(median_subset))
# mCRPC_newdx_wilcox_pval <- rep(NA, nrow(median_subset))
# nmCRPC_newdx_wilcox_pval <- rep(NA, nrow(median_subset))
# mCRPC_nmCSPC_wilcox_pval <- rep(NA, nrow(median_subset))

# get wilcox pval (2-sided)
for(i in 1: nrow(median_subset)){
  mCRPC_nmCRPC_wilcox_pval[i] <- wilcox.test(
    x = as.numeric(median_subset[i, patient_key$stage == "mCRPC"]),
    y = as.numeric(median_subset[i, patient_key$stage == "nmCRPC"]),
    alternative = "two.sided", exact = F, conf.level = .95
  )$'p.value'
  nmCRPC_nmCSPC_wilcox_pval[i] <- wilcox.test(
    x = as.numeric(median_subset[i, patient_key$stage == "nmCRPC"]),
    y = as.numeric(median_subset[i, patient_key$stage == "nmCSPC"]),
    alternative = "two.sided", exact = F, conf.level = .95
  )$'p.value'
  nmCSPC_newdx_wilcox_pval[i] <- wilcox.test(
    x = as.numeric(median_subset[i, patient_key$stage == "nmCSPC"]),
    y = as.numeric(median_subset[i, patient_key$stage == "new_dx"]),
    alternative = "two.sided", exact = F, conf.level = .95
  )$'p.value'
  newdx_normal_wilcox_pval[i] <- wilcox.test(
    x = as.numeric(median_subset[i, patient_key$stage == "new_dx"]),
    y = as.numeric(median_subset[i, patient_key$stage == "normal"]),
    alternative = "two.sided", exact = F, conf.level = .95
  )$'p.value'
  cancer_normal_wilcox_pval[i] <-  wilcox.test(
    x = as.numeric(median_subset[i, patient_key$stage %in% c("new_dx", "nmCSPC", "nmCRPC", "mCRPC")]),
    y = as.numeric(median_subset[i, patient_key$stage == "normal"]),
    alternative = "two.sided", exact = F, conf.level = .95
  )$'p.value'
  mCRPC_others_wilcox_pval[i] <- wilcox.test(
    x = as.numeric(median_subset[i, patient_key$stage == "mCRPC"]),
    y = as.numeric(median_subset[i, patient_key$stage %in% c("normal", "new_dx", "nmCSPC", "nmCRPC")]),
    alternative = "two.sided", exact = F, conf.level = .95
  )$'p.value'
  
  # mCRPC_normal_wilcox_pval[i] <- wilcox.test(
  #   x = as.numeric(median_subset[i, patient_key$stage == "mCRPC"]),
  #   y = as.numeric(median_subset[i, patient_key$stage == "normal"]),
  #   alternative = "two.sided", exact = F, conf.level = .95
  # )$'p.value'
  # nmCRPC_normal_wilcox_pval[i] <- wilcox.test(
  #   x = as.numeric(median_subset[i, patient_key$stage == "nmCRPC"]),
  #   y = as.numeric(median_subset[i, patient_key$stage == "normal"]),
  #   alternative = "two.sided", exact = F, conf.level = .95
  # )$'p.value'
  # nmCSPC_normal_wilcox_pval[i] <- wilcox.test(
  #   x = as.numeric(median_subset[i, patient_key$stage == "nmCSPC"]),
  #   y = as.numeric(median_subset[i, patient_key$stage == "normal"]),
  #   alternative = "two.sided", exact = F, conf.level = .95
  # )$'p.value'
  # mCRPC_newdx_wilcox_pval[i] <- wilcox.test(
  #   x = as.numeric(median_subset[i, patient_key$stage == "mCRPC"]),
  #   y = as.numeric(median_subset[i, patient_key$stage == "new_dx"]),
  #   alternative = "two.sided", exact = F, conf.level = .95
  # )$'p.value'
  # nmCRPC_newdx_wilcox_pval[i] <- wilcox.test(
  #   x = as.numeric(median_subset[i, patient_key$stage == "nmCRPC"]),
  #   y = as.numeric(median_subset[i, patient_key$stage == "new_dx"]),
  #   alternative = "two.sided", exact = F, conf.level = .95
  # )$'p.value'
  # mCRPC_nmCSPC_wilcox_pval[i] <- wilcox.test(
  #   x = as.numeric(median_subset[i, patient_key$stage == "mCRPC"]),
  #   y = as.numeric(median_subset[i, patient_key$stage == "nmCSPC"]),
  #   alternative = "two.sided", exact = F, conf.level = .95
  # )$'p.value'
  
  print(i)
}

# check pval histograms (restricted to signif peptides from ANOVA)
hist(mCRPC_nmCRPC_wilcox_pval, breaks = 70)
hist(nmCRPC_nmCSPC_wilcox_pval, breaks = 70)
hist(nmCSPC_newdx_wilcox_pval, breaks = 70)
hist(newdx_normal_wilcox_pval, breaks = 70)
hist(cancer_normal_wilcox_pval, breaks = 70)
hist(mCRPC_others_wilcox_pval, breaks = 70)

# hist(mCRPC_normal_wilcox_pval, breaks = 70)
# hist(nmCRPC_normal_wilcox_pval, breaks = 70)
# hist(nmCSPC_normal_wilcox_pval, breaks = 70)
# hist(mCRPC_newdx_wilcox_pval, breaks = 70)
# hist(nmCRPC_newdx_wilcox_pval, breaks = 70)
# hist(mCRPC_nmCSPC_wilcox_pval, breaks = 70)

# get BH-corrected pval (restricted to signif peptides from ANOVA)
mCRPC_nmCRPC_wilcox_BH <- p.adjust(mCRPC_nmCRPC_wilcox_pval, method = "BH")
nmCRPC_nmCSPC_wilcox_BH <- p.adjust(nmCRPC_nmCSPC_wilcox_pval, method = "BH")
nmCSPC_newdx_wilcox_BH <- p.adjust(nmCSPC_newdx_wilcox_pval, method = "BH")
newdx_normal_wilcox_BH <- p.adjust(newdx_normal_wilcox_pval, method = "BH")
cancer_normal_wilcox_BH <- p.adjust(cancer_normal_wilcox_pval, method = "BH")
mCRPC_others_wilcox_BH <- p.adjust(mCRPC_others_wilcox_pval, method = "BH")

# mCRPC_normal_wilcox_BH <- p.adjust(mCRPC_normal_wilcox_pval, method = "BH")
# nmCRPC_normal_wilcox_BH <- p.adjust(nmCRPC_normal_wilcox_pval, method = "BH")
# nmCSPC_normal_wilcox_BH <- p.adjust(nmCSPC_normal_wilcox_pval, method = "BH")
# mCRPC_newdx_wilcox_BH <- p.adjust(mCRPC_newdx_wilcox_pval, method = "BH")
# nmCRPC_newdx_wilcox_BH <- p.adjust(nmCRPC_newdx_wilcox_pval, method = "BH")
# mCRPC_nmCSPC_wilcox_BH <- p.adjust(mCRPC_nmCSPC_wilcox_pval, method = "BH")

# check peptide counts at different BH FDR thresholds
count.func(mCRPC_nmCRPC_wilcox_BH, seq(0.03, 0.15, by = 0.01))
count.func(nmCRPC_nmCSPC_wilcox_BH, seq(0.03, 0.15, by = 0.01))
count.func(nmCSPC_newdx_wilcox_BH, seq(0.03, 0.15, by = 0.01))
count.func(newdx_normal_wilcox_BH, seq(0.03, 0.15, by = 0.01))
count.func(cancer_normal_wilcox_BH, seq(0.03, 0.15, by = 0.01))
count.func(mCRPC_others_wilcox_BH, seq(0.03, 0.15, by = 0.01))

# count.func(mCRPC_normal_wilcox_BH, seq(0.03, 0.15, by = 0.01)) 
# count.func(nmCRPC_normal_wilcox_BH , seq(0.03, 0.15, by = 0.01))
# count.func(nmCSPC_normal_wilcox_BH , seq(0.03, 0.15, by = 0.01))
# count.func(mCRPC_newdx_wilcox_BH , seq(0.03, 0.15, by = 0.01))
# count.func(nmCRPC_newdx_wilcox_BH , seq(0.03, 0.15, by = 0.01))
# count.func(mCRPC_nmCSPC_wilcox_BH , seq(0.03, 0.15, by = 0.01))

# check how many peptides meet BH-FDR & effect-size thresholds for particular contrast
length(which(abs(mCRPC_mean - nmCRPC_mean) > 1 & mCRPC_nmCRPC_wilcox_BH <= 0.05))
length(which(abs(nmCRPC_mean - nmCSPC_mean) > 1 & nmCRPC_nmCSPC_wilcox_BH <= 0.05))
length(which(abs(nmCSPC_mean - newdx_mean) > 1 & nmCSPC_newdx_wilcox_BH <= 0.05))
length(which(abs(newdx_mean - normal_mean) > 1 & newdx_normal_wilcox_BH<= 0.05))
length(which(abs(cancer_mean - normal_mean) > 1 & cancer_normal_wilcox_BH <= 0.05))
length(which(abs(mCRPC_mean - NOT_mCRPC_mean) > 1 & mCRPC_others_wilcox_BH <= 0.05))

# length(which(abs(mCRPC_mean - normal_mean) > 1 & mCRPC_normal_wilcox_BH <= 0.05))
# length(which(abs(nmCRPC_mean - normal_mean) > 1 & nmCRPC_normal_wilcox_BH <= 0.05))
# length(which(abs(nmCSPC_mean - normal_mean) > 1 & nmCSPC_normal_wilcox_BH <= 0.05))
# length(which(abs(mCRPC_mean - newdx_mean) > 1 & mCRPC_newdx_wilcox_BH <= 0.05))
# length(which(abs(nmCRPC_mean - newdx_mean) > 1 & nmCRPC_newdx_wilcox_BH <= 0.05))
# length(which(abs(mCRPC_mean - nmCSPC_mean) > 1 & mCRPC_nmCSPC_wilcox_BH <= 0.05))


#----------------------------------------------------------------------------------------------
# t-tests

# get group counts
mCRPC_counts <- nrow(patient_key[patient_key$stage=="mCRPC",])
nmCRPC_counts <- nrow(patient_key[patient_key$stage=="nmCRPC",])
nmCSPC_counts <- nrow(patient_key[patient_key$stage=="nmCSPC",])
newdx_counts <- nrow(patient_key[patient_key$stage=="new_dx",])
normal_counts <- nrow(patient_key[patient_key$stage=="normal",])
cancer_counts <- nrow(patient_key[patient_key$stage !="normal",])
NOT_mCRPC_counts <- nrow(patient_key[patient_key$stage !="mCRPC",])

# get t-statistics
BH_filter_mse <- all_anova_mse[all_anova$anova_BH <= .05]
mCRPC_nmCRPC_tstat <- (mCRPC_mean - nmCRPC_mean) / sqrt(BH_filter_mse * (1/mCRPC_counts + 1/nmCRPC_counts))
nmCRPC_nmCSPC_tstat <- (nmCRPC_mean - nmCSPC_mean) / sqrt(BH_filter_mse * (1/nmCRPC_counts + 1/nmCSPC_counts))
nmCSPC_newdx_tstat <- (nmCSPC_mean - newdx_mean) / sqrt(BH_filter_mse * (1/nmCSPC_counts + 1/newdx_counts))
newdx_normal_tstat <- (newdx_mean - normal_mean) / sqrt(BH_filter_mse * (1/newdx_counts + 1/normal_counts))
cancer_normal_tstat <- (cancer_mean - normal_mean) / sqrt(BH_filter_mse * (1/cancer_counts + 1/normal_counts))
mCRPC_others_tstat <- (mCRPC_mean - NOT_mCRPC_mean) / sqrt(BH_filter_mse * (1/mCRPC_counts + 1/NOT_mCRPC_counts))

# mCRPC_normal_tstat <- (mCRPC_mean - normal_mean) / sqrt(BH_filter_mse * (1/mCRPC_counts + 1/normal_counts))
# nmCRPC_normal_tstat <- (nmCRPC_mean - normal_mean) / sqrt(BH_filter_mse * (1/nmCRPC_counts + 1/normal_counts))
# nmCSPC_normal_tstat <- (nmCSPC_mean - normal_mean) / sqrt(BH_filter_mse * (1/nmCSPC_counts + 1/normal_counts))
# mCRPC_newdx_tstat <- (mCRPC_mean - newdx_mean) / sqrt(BH_filter_mse * (1/mCRPC_counts + 1/newdx_counts))
# nmCRPC_newdx_tstat <- (nmCRPC_mean - newdx_mean) / sqrt(BH_filter_mse * (1/nmCRPC_counts + 1/newdx_counts))
# mCRPC_nmCSPC_tstat <- (mCRPC_mean - nmCSPC_mean) / sqrt(BH_filter_mse * (1/mCRPC_counts + 1/nmCSPC_counts))

# get tstat-pval (two-sided)
mCRPC_nmCRPC_pval <- 2 * pt( -abs(mCRPC_nmCRPC_tstat), df = nrow(patient_key) - length(unique(patient_key$stage)) )
nmCRPC_nmCSPC_pval <- 2 * pt( -abs(nmCRPC_nmCSPC_tstat), df = nrow(patient_key) - length(unique(patient_key$stage)) )
nmCSPC_newdx_pval <- 2 * pt( -abs(nmCSPC_newdx_tstat), df = nrow(patient_key) - length(unique(patient_key$stage)) )
newdx_normal_pval <- 2 * pt( -abs(newdx_normal_tstat), df = nrow(patient_key) - length(unique(patient_key$stage)) )
cancer_normal_pval <- 2 * pt( -abs(cancer_normal_tstat), df = nrow(patient_key) - length(unique(patient_key$stage)) )
mCRPC_others_pval <- 2 * pt( -abs(mCRPC_others_tstat), df = nrow(patient_key) - length(unique(patient_key$stage)) )

# mCRPC_normal_pval <-2 * pt( -abs(mCRPC_normal_tstat), df = nrow(patient_key) - length(unique(patient_key$stage)) )
# nmCRPC_normal_pval <- 2 * pt( -abs(nmCRPC_normal_tstat), df = nrow(patient_key) - length(unique(patient_key$stage)) )
# nmCSPC_normal_pval <- 2 * pt( -abs(nmCSPC_normal_tstat), df = nrow(patient_key) - length(unique(patient_key$stage)) )
# mCRPC_newdx_pval <- 2 * pt( -abs(mCRPC_newdx_tstat), df = nrow(patient_key) - length(unique(patient_key$stage)) )
# nmCRPC_newdx_pval <- 2 * pt( -abs(nmCRPC_newdx_tstat), df = nrow(patient_key) - length(unique(patient_key$stage)) )
# mCRPC_nmCSPC_pval <- 2 * pt( -abs(mCRPC_nmCSPC_tstat), df = nrow(patient_key) - length(unique(patient_key$stage)) )

# check pval histograms (restricted to signif peptides from ANOVA)
hist(mCRPC_nmCRPC_pval, breaks = 70)
hist(nmCRPC_nmCSPC_pval, breaks = 70)
hist(nmCSPC_newdx_pval, breaks = 70)
hist(newdx_normal_pval, breaks = 70)
hist(cancer_normal_pval, breaks = 70)
hist(mCRPC_others_pval, breaks = 70)

hist(mCRPC_normal_pval, breaks = 70)
hist(nmCRPC_normal_pval, breaks = 70)
hist(nmCSPC_normal_pval, breaks = 70)
hist(mCRPC_newdx_pval, breaks = 70)
hist(nmCRPC_newdx_pval, breaks = 70)
hist(mCRPC_nmCSPC_pval, breaks = 70)

# get BH-corrected pval (restricted to signif peptides from ANOVA)
mCRPC_nmCRPC_BH <- p.adjust(mCRPC_nmCRPC_pval, method = "BH")
nmCRPC_nmCSPC_BH <- p.adjust(nmCRPC_nmCSPC_pval, method = "BH")
nmCSPC_newdx_BH <- p.adjust(nmCSPC_newdx_pval, method = "BH")
newdx_normal_BH <- p.adjust(newdx_normal_pval, method = "BH")
cancer_normal_BH <- p.adjust(cancer_normal_pval, method = "BH")
mCRPC_others_BH <- p.adjust(mCRPC_others_pval, method = "BH")

# mCRPC_normal_BH <- p.adjust(mCRPC_normal_pval, method = "BH")
# nmCRPC_normal_BH <- p.adjust(nmCRPC_normal_pval, method = "BH")
# nmCSPC_normal_BH <- p.adjust(nmCSPC_normal_pval, method = "BH")
# mCRPC_newdx_BH <- p.adjust(mCRPC_newdx_pval, method = "BH")
# nmCRPC_newdx_BH <- p.adjust(nmCRPC_newdx_pval, method = "BH")
# mCRPC_nmCSPC_BH <- p.adjust(mCRPC_nmCSPC_pval, method = "BH")

# check peptide counts at different BH FDR thresholds
count.func(mCRPC_nmCRPC_BH, seq(0.03, 0.15, by = 0.01))
count.func(nmCRPC_nmCSPC_BH, seq(0.03, 0.15, by = 0.01))
count.func(nmCSPC_newdx_BH, seq(0.03, 0.15, by = 0.01))
count.func(newdx_normal_BH, seq(0.03, 0.15, by = 0.01))
count.func(cancer_normal_BH, seq(0.03, 0.15, by = 0.01))
count.func(mCRPC_others_BH, seq(0.03, 0.15, by = 0.01))

# count.func(mCRPC_normal_BH, seq(0.03, 0.15, by = 0.01))
# count.func(nmCRPC_normal_BH, seq(0.03, 0.15, by = 0.01))
# count.func(nmCSPC_normal_BH, seq(0.03, 0.15, by = 0.01))
# count.func(mCRPC_newdx_BH , seq(0.03, 0.15, by = 0.01))
# count.func(nmCRPC_newdx_BH, seq(0.03, 0.15, by = 0.01))
# count.func(mCRPC_nmCSPC_BH, seq(0.03, 0.15, by = 0.01))

# check how many peptides meet BH-FDR & effect-size thresholds for particular contrast
length(which(abs(mCRPC_mean - nmCRPC_mean) > 1 & mCRPC_nmCRPC_BH <= 0.05))
length(which(abs(nmCRPC_mean - nmCSPC_mean) > 1 & nmCRPC_nmCSPC_BH <= 0.05))
length(which(abs(nmCSPC_mean - newdx_mean) > 1 & nmCSPC_newdx_BH <= 0.05))
length(which(abs(newdx_mean - normal_mean) > 1 & newdx_normal_BH<= 0.05))
length(which(abs(cancer_mean - normal_mean) > 1 & cancer_normal_BH <= 0.05))
length(which(abs(mCRPC_mean - NOT_mCRPC_mean) > 1 & mCRPC_others_BH <= 0.05))

# length(which(abs(mCRPC_mean - normal_mean) > 1 & mCRPC_normal_BH <= 0.05))
# length(which(abs(nmCRPC_mean - normal_mean) > 1 & nmCRPC_normal_BH <= 0.05))
# length(which(abs(nmCSPC_mean - normal_mean) > 1 & nmCSPC_normal_BH <= 0.05))
# length(which(abs(mCRPC_mean - newdx_mean) > 1 & mCRPC_newdx_BH <= 0.05))
# length(which(abs(nmCRPC_mean - newdx_mean) > 1 & nmCRPC_newdx_BH <= 0.05))
# length(which(abs(mCRPC_mean - nmCSPC_mean) > 1 & mCRPC_nmCSPC_BH <= 0.05))

#----------------------------------------------------------------------------------------------
# volcano plots of contrasts

# make sure contrast_diff and contrast_pval and contrast_BH of same length !!
contrast_volcano_plot.func <- function(contrast_diff, contrast_pval, contrast_BH, test_type, measure, contrast_title){
  plot(x = contrast_diff, y = -log10(contrast_pval), pch = 20,
       xlab = paste0("difference of ", measure, " log2(fluorescence)"), 
       ylab = paste0("-log10(contrast ",test_type, " p-values)"),
       main = paste0("Volcano plot of contrast: ", contrast_title))
  lines(x = contrast_diff[ contrast_BH <= .05 & abs(contrast_diff) >= 1 ], 
        y = -log10(contrast_pval[ contrast_BH <= .05 & abs(contrast_diff) >= 1]),
        type = "p", pch = 20, col = "blue")
}

contrast_volcano_plot.func(mCRPC_median - NOT_mCRPC_median, 
                           mCRPC_others_wilcox_pval, 
                           mCRPC_others_wilcox_BH, 
                           "Wilcoxon test",
                           "median",
                           "mCRPC vs others")

#----------------------------------------------------------------------------------------------
# Write to Excel

anova_df <- data.frame(
  PROBE_ID = raw_data_median$PROBE_ID[all_anova$anova_BH <= .05],
  SEQ_ID = raw_data_median$SEQ_ID[all_anova$anova_BH <= .05],
  ANOVA_BH_FDR = all_anova$anova_BH[all_anova$anova_BH <= .05]
)
anova_df <- anova_df %>% 
  remove_rownames() %>%
  mutate( ANOVA_BH_FDR = round(ANOVA_BH_FDR, 4) )

contrast_df.func <- function(contrast_BH, contrast_diff){
  df <- data.frame(
    PROBE_ID = raw_data_median$PROBE_ID[all_anova$anova_BH <= .05][contrast_BH <= .05],
    SEQ_ID = raw_data_median$SEQ_ID[all_anova$anova_BH <= .05][contrast_BH <= .05],
    Effect_size = contrast_diff[contrast_BH <= .05],
    ANOVA_BH_FDR = all_anova$anova_BH[all_anova$anova_BH <= .05][contrast_BH <= .05]
  )
  df <- df %>% 
    filter(abs(Effect_size)>1) %>%
    arrange(desc(Effect_size)) %>%
    remove_rownames() %>%
    mutate(Effect_size = round(Effect_size, 4),
           ANOVA_BH_FDR = round(ANOVA_BH_FDR, 4))
  return(df)
}

# based on t-test
# mCRPC_others_df <- contrast_df.func(mCRPC_others_BH, mCRPC_mean - NOT_mCRPC_mean)
# cancer_normal_df <- contrast_df.func(cancer_normal_BH, cancer_mean - normal_mean)
# mCRPC_nmCRPC_df <- contrast_df.func(mCRPC_nmCRPC_BH, mCRPC_mean - nmCRPC_mean) 
# nmCRPC_nmCSPC_df <- contrast_df.func(nmCRPC_nmCSPC_BH, nmCRPC_mean - nmCSPC_mean) 
# nmCSPC_newdx_df <- contrast_df.func(nmCSPC_newdx_BH, nmCSPC_mean - newdx_mean) 
# newdx_normal_df <- contrast_df.func(newdx_normal_BH, newdx_mean - normal_mean)

# based on wilcox test
mCRPC_others_df <- contrast_df.func(mCRPC_others_wilcox_BH, mCRPC_mean - NOT_mCRPC_mean)
cancer_normal_df <- contrast_df.func(cancer_normal_wilcox_BH, cancer_mean - normal_mean)
mCRPC_nmCRPC_df <- contrast_df.func(mCRPC_nmCRPC_wilcox_BH, mCRPC_mean - nmCRPC_mean) 
nmCRPC_nmCSPC_df <- contrast_df.func(nmCRPC_nmCSPC_wilcox_BH, nmCRPC_mean - nmCSPC_mean) 
nmCSPC_newdx_df <- contrast_df.func(nmCSPC_newdx_wilcox_BH, nmCSPC_mean - newdx_mean) 
newdx_normal_df <- contrast_df.func(newdx_normal_wilcox_BH, newdx_mean - normal_mean)

summary_page_df <- data.frame(
  contrasts = c("mCRPC vs others", "cancer vs normal", "mCRPC vs nmCRPC",
                "nmCRPC vs nmCSPC", "nmCSPC vs new_dx", "new_dx vs normal"),
  total_peptide_counts = c( nrow(mCRPC_others_df), nrow(cancer_normal_df), nrow(mCRPC_nmCRPC_df),
                            nrow(nmCRPC_nmCSPC_df), nrow(nmCSPC_newdx_df), nrow(newdx_normal_df) ),
  peptides_with_positive_effect_size = c(mCRPC_others_df %>% filter(Effect_size > 0) %>% nrow(),
                                         cancer_normal_df %>% filter(Effect_size > 0) %>% nrow(),
                                         mCRPC_nmCRPC_df %>% filter(Effect_size > 0) %>% nrow(),
                                         nmCRPC_nmCSPC_df %>% filter(Effect_size > 0) %>% nrow(),
                                         nmCSPC_newdx_df %>% filter(Effect_size > 0) %>% nrow(),
                                         newdx_normal_df %>% filter(Effect_size > 0) %>% nrow())
)
summary_page_df <- summary_page_df %>% 
  mutate(peptides_with_negative_effect_size = total_peptide_counts - peptides_with_positive_effect_size)

wb = createWorkbook()
sheet = createSheet(wb, "Contrast Summaries")
addDataFrame(summary_page_df, sheet=sheet, row.names=FALSE, startRow = 2, startColumn = 2)

sheet = createSheet(wb, "One-Way ANOVA")
addDataFrame(anova_df, sheet=sheet, row.names=FALSE)

sheet = createSheet(wb, "mCRPC vs others")
addDataFrame(mCRPC_others_df, sheet=sheet, row.names=FALSE)

sheet = createSheet(wb, "cancer vs normal")
addDataFrame(cancer_normal_df, sheet=sheet, row.names=FALSE)

sheet = createSheet(wb, "mCRPC vs nmCRPC")
addDataFrame(mCRPC_nmCRPC_df, sheet=sheet, row.names=FALSE)

sheet = createSheet(wb, "nmCRPC vs nmCSPC")
addDataFrame(nmCRPC_nmCSPC_df, sheet=sheet, row.names=FALSE)

sheet = createSheet(wb, "nmCSPC vs newdx")
addDataFrame(nmCSPC_newdx_df, sheet=sheet, row.names=FALSE)

sheet = createSheet(wb, "newdx vs normal")
addDataFrame(newdx_normal_df, sheet=sheet, row.names=FALSE)

saveWorkbook(wb, "07_CancerStages_SignifPeptides.xlsx")

#----------------------------------------------------------------------------------------------
# Boxplot of interesting peptides 

# first make sure stage aligns with raw_data_median
sum(as.numeric(colnames(raw_data_median %>% select(-(PROBE_DESIGN_ID:DESIGN_ID))) == patient_key$id)) == nrow(patient_key)

contrast_df.func <- function(contrast_BH, contrast_diff){
  df <- data.frame(
    PROBE_ID = raw_data_median$PROBE_ID[all_anova$anova_BH <= .05][contrast_BH <= .05],
    SEQ_ID = raw_data_median$SEQ_ID[all_anova$anova_BH <= .05][contrast_BH <= .05],
    Effect_size = contrast_diff[contrast_BH <= .05],
    ANOVA_BH_FDR = all_anova$anova_BH[all_anova$anova_BH <= .05][contrast_BH <= .05]
  )
  df <- df %>% 
    filter(abs(Effect_size)>1) %>%
    arrange(desc(Effect_size)) %>%
    remove_rownames() %>%
    mutate(Effect_size = round(Effect_size, 4),
           ANOVA_BH_FDR = round(ANOVA_BH_FDR, 4))
  return(df)
}

# based on wilcox test
mCRPC_others_df <- contrast_df.func(mCRPC_others_wilcox_BH, mCRPC_mean - NOT_mCRPC_mean)
cancer_normal_df <- contrast_df.func(cancer_normal_wilcox_BH, cancer_mean - normal_mean)
mCRPC_nmCRPC_df <- contrast_df.func(mCRPC_nmCRPC_wilcox_BH, mCRPC_mean - nmCRPC_mean) 
nmCRPC_nmCSPC_df <- contrast_df.func(nmCRPC_nmCSPC_wilcox_BH, nmCRPC_mean - nmCSPC_mean) 
nmCSPC_newdx_df <- contrast_df.func(nmCSPC_newdx_wilcox_BH, nmCSPC_mean - newdx_mean) 
newdx_normal_df <- contrast_df.func(newdx_normal_wilcox_BH, newdx_mean - normal_mean)

boxplot_func <- function(ref_df, grp1, grp2, draw, rename_grp1, rename_grp2, col = c("red", "blue")){
  mat <- raw_data_median %>%
    filter(PROBE_ID %in% ref_df$PROBE_ID) 
  mat <- mat[match(ref_df$PROBE_ID, mat$PROBE_ID), ] %>%
    select(-(PROBE_DESIGN_ID:DESIGN_ID)) %>%
    select(which(patient_key$stage %in% c(grp1, grp2))) %>%
    as.matrix()
  row.names(mat) <- ref_df$PROBE_ID
  mat_stage <- as.character(patient_key$stage[patient_key$stage %in% c(grp1, grp2)])
  if (length(grp1) > 1){
    mat_stage[mat_stage %in% grp1] = rename_grp1
  }
  if (length(grp2) > 1){
    mat_stage[mat_stage %in% grp2] = rename_grp2
  }
  mat_stage = factor(mat_stage)
  par(mar = c(4.1, 5.5, 2.8, 1),cex = 0.84)
  graphics::boxplot( as.numeric(mat[draw,]) ~ mat_stage, varwidth = T,
                       col=col, horizontal=TRUE, las=1, xlab = "log2(fluorescence)", ylab = "groups",
                       main= paste(row.names(mat)[draw]) )
    
}

png("07_1a_Cancer_higher_than_Normal.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=1)
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=2)
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=3)
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=4)
dev.off()

png("07_1b_Normal_higher_than_Cancer.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=nrow(cancer_normal_df))
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=nrow(cancer_normal_df)-1)
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=nrow(cancer_normal_df)-2)
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=nrow(cancer_normal_df)-3)
dev.off()

png("07_2a_mCRPC_higher_than_others.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=1)
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=2)
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=3)
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=4)
dev.off()

png("07_2b_others_higher_than_mCRPC.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=nrow(mCRPC_others_df))
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=nrow(mCRPC_others_df)-1)
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=nrow(mCRPC_others_df)-2)
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=nrow(mCRPC_others_df)-3)
dev.off()

png("07_3a_Newdx_higher_than_Normal.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = 1)
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = 2)
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = 3)
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = 4)
dev.off()

png("07_3b_Normal_higher_than_Newdx.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = nrow(newdx_normal_df))
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = nrow(newdx_normal_df)-1)
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = nrow(newdx_normal_df)-2)
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = nrow(newdx_normal_df)-3)
dev.off()

png("07_4a_nmCSPC_higher_than_Newdx.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = 1)
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = 2)
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = 3)
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = 4)
dev.off()

png("07_4b_Newdx_higher_than_nmCSPC.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = nrow(nmCSPC_newdx_df))
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = nrow(nmCSPC_newdx_df)-1)
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = nrow(nmCSPC_newdx_df)-2)
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = nrow(nmCSPC_newdx_df)-3)
dev.off()

png("07_5a_nmCRPC_higher_than_nmCSPC.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=1)
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=2)
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=3)
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=4)
dev.off()

png("07_5b_nmCSPC_higher_than_nmCRPC.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=nrow(nmCRPC_nmCSPC_df))
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=nrow(nmCRPC_nmCSPC_df)-1)
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=nrow(nmCRPC_nmCSPC_df)-2)
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=nrow(nmCRPC_nmCSPC_df)-3)
dev.off()

png("07_6a_mCRPC_higher_than_nmCRPC.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=1)
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=2)
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=3)
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=4)
dev.off()

png("07_6b_nmCRPC_higher_than_mCRPC.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=nrow(mCRPC_nmCRPC_df))
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=nrow(mCRPC_nmCRPC_df)-1)
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=nrow(mCRPC_nmCRPC_df)-2)
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=nrow(mCRPC_nmCRPC_df)-3)
dev.off()

#----------------------------------------------------------------------------------------------
# replot subset of ANOVA residual heatmap based on effect-size threshold from post-hoc analysis

posthoc_signif_crit <- ( abs(mCRPC_median - nmCRPC_median) > 1 & mCRPC_nmCRPC_wilcox_BH <= 0.05 ) |
( abs(nmCRPC_median - nmCSPC_median) > 1 & nmCRPC_nmCSPC_wilcox_BH <= 0.05 ) |
( abs(nmCSPC_median - newdx_median) > 1 & nmCSPC_newdx_wilcox_BH <= 0.05 ) |
( abs(newdx_median - normal_median) > 1 & newdx_normal_wilcox_BH<= 0.05 ) |
( abs(cancer_median - normal_median) > 1 & cancer_normal_wilcox_BH <= 0.05 ) |
( abs(mCRPC_median - NOT_mCRPC_median) > 1 & mCRPC_others_wilcox_BH <= 0.05 )

# check
length(posthoc_signif_crit)

anova_dat <- raw_data_median %>% 
  select(PROBE_ID, any_of(patient_key$id)) %>% 
  mutate(anova_BH = all_anova$anova_BH[raw_data_median$PROBE_ID]) %>%
  # mutate(anova_BH = all_kw$anova_BH[raw_data_median$PROBE_ID]) %>%
  filter(anova_BH <= .05) %>%
  select(-anova_BH) %>%
  filter(posthoc_signif_crit)

anova_dat_demean <- sweep(as.matrix(anova_dat %>% select(-PROBE_ID)), 1, 
                          rowMeans(as.matrix(anova_dat %>% select(-PROBE_ID))), "-") # centering by row

# make sure stage aligns with anova_dat_demean
visual_iii <- match( colnames(anova_dat_demean) , patient_key$id )
visual_stage <- patient_key$stage[visual_iii]

# colors and shapes for the visualization techniques
cols = pal[ match(visual_stage, names(pal)) ]
cls <- colorRampPalette(c("navy", "honeydew", "firebrick3", "brown"))(n = 1024)

get_column_order.func <- function(stages){
  heat_map <- heatmap3(anova_dat_demean[,visual_stage == stages],col = cls, labRow = "")
  return( colnames(anova_dat_demean[,visual_stage == stages])[heat_map$colInd] )
}

normal_id_order <- get_column_order.func("normal")
newdx_id_order <- get_column_order.func("new_dx")
nmCSPC_id_order <- get_column_order.func("nmCSPC")
nmCRPC_id_order <- get_column_order.func("nmCRPC")
mCRPC_id_order <- get_column_order.func("mCRPC")

id_order <- match( c(normal_id_order, newdx_id_order, nmCSPC_id_order, nmCRPC_id_order, mCRPC_id_order),
                   colnames(anova_dat_demean) )

heatmap3(anova_dat_demean[,id_order], 
         col = cls, # specify colors 
         ColSideColors = cols[id_order], # specify patient color code
         Colv = NA,
         labCol = visual_stage[id_order], # specify patient
         ColSideLabs = "stages", 
         labRow = "",
         xlab = "Patients",
         # legendfun=function() showLegend(col = c("navy", "cornflowerblue", "turquoise1", "darkorange1", "firebrick1"),
         #                                 legend = c("normal",  "new_dx", "nmCSPC", "nmCRPC", "mCRPC"),
         #                                 cex = 1.2,
         #                                 lwd = 5  )
)



#----------------------------------------------------------------------------------------------

# save results
# save(logreg_pval,
#      lmer_result,
#      all_anova_pval,
#      all_kw_pval,
#      all_anova_mse,
#      file = "07_Cancer_Stage_Effects.RData")

