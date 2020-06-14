library(xlsx)
library(allez)
library(gridExtra)
library(ggplot2)
library(matrixStats) # rowMedians
library(lme4) # linear mixed effects model
library(lmerTest)
library(fdrtool)
library(heatmap3)
library(tidyverse) # make sure you have the latest tidyverse !

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
  plot(Z[,x], Z[,y], col = cols, pch = shapes, main = title, las = 1,
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
anova_func <- function(anova_pval, xlab){
  # control FDR
  anova_BH <- p.adjust(anova_pval, method = "BH")
  anova_qval <- fdrtool(anova_pval, statistic = "pvalue", verbose = F, plot  = F)$qval
  anova_qval_eta0 <- unname(fdrtool(anova_pval, statistic = "pvalue", verbose = F, plot  = F)$param[,"eta0"])
  
  # plot histogram of p-values
  all_peptide_hist <- hist(anova_pval, breaks = 70, freq = F, xlab = xlab, las = 1,  
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

# post-process allez table for kable output
get_alleztable.func <- function(allez_go_input){
  allez.tab <- allezTable(allez_go_input, symbol = T, nominal.alpha = nom.alpha, n.upp = max_gene_in_set, in.set = T)
  allez.tab$set.size <- paste(allez.tab$in.set, allez.tab$set.size, sep = "/")
  allez.tab <- allez.tab %>% dplyr::select(-c(in.set, genes)) %>%
    mutate(in.genes = str_replace_all(in.genes, ";", "; "))
  return(allez.tab)
}


load("09_Cancer_Stage_Effects.RData")
load("09_LMER_results.RData")
raw_data_median <- read_csv("raw_data_median.csv")
raw_data_median_proj2 <- read_csv("raw_data_median_proj2.csv")

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

# raw_data = read_csv("raw_data_complete.csv")
# 
# sum(as.numeric(raw_data$X == raw_data$COL_NUM)) == nrow(raw_data) # X = COL_NUM
# sum(as.numeric(raw_data$Y == raw_data$ROW_NUM)) == nrow(raw_data) # Y = ROW_NUM
# sum(as.numeric(raw_data$MATCH_INDEX == raw_data$FEATURE_ID)) == nrow(raw_data) # MATCH_INDEX = FEATURE_ID
# unique(raw_data$DESIGN_NOTE) # only NA
# unique(raw_data$SELECTION_CRITERIA) # only NA
# unique(raw_data$MISMATCH) # only 0
# unique(raw_data$PROBE_CLASS) # only NA
# # we can drop X, Y, MATCH_INDEX, DESIGN_NOTE, SELECTION_CRITERIA, MISMATCH, PROBE_CLASS
# 
# raw_data = raw_data %>%
#   select(PROBE_DESIGN_ID:Y, any_of(array_id_key$array_id)) %>% # drop patients with records at different stages 
#   select( -c(X, Y, MATCH_INDEX, DESIGN_NOTE, SELECTION_CRITERIA, MISMATCH, PROBE_CLASS) ) # drop unhelpful columns
# 
# colnames(raw_data)[1:15] # check
# 
# # extract all sequence id
# # may need to use this for gene set analysis
# all_seq_id <-  unique(raw_data$SEQ_ID)
# 
# 
# # take log2 transformation
# raw_data <- raw_data %>%
#   mutate_at(vars(matches("dat")), log2)
# 
# 
# # make sure array_id_key align with raw_data_complete
# array_iii <- match( colnames( raw_data %>% select(any_of(array_id_key$array_id)) ), array_id_key$array_id )
# array_id_key <- array_id_key[array_iii , ]
# sum(as.numeric( colnames( raw_data %>% select(any_of(array_id_key$array_id)) ) == array_id_key$array_id )) == nrow(array_id_key)
# 
# 
# # compute median
# raw_data_median <- t( apply( select(raw_data, contains("dat")), 1, function(x) {
#   tapply(x, array_id_key$id, FUN=median)
# } ) )
# raw_data_median <- bind_cols( select(raw_data, -contains("dat")), as.data.frame(raw_data_median) ) 
# 
# colnames(raw_data_median)[1:15] # check
# 
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

# get calls_long for later
calls_long <- calls %>%
  select(PROBE_ID, any_of(array_id_key$id)) 

# remove peptides that have zero calls in ALL subjects
calls <- calls[ apply( calls %>% select(any_of(patient_key$id)) , 1, function(x){ !all(x==0) } ) , ]

# in the end, how many peptides with at least one call among all patients
nrow(calls)

#----------------------------------------------------------------------------------------------

# get median_long 
median_long2 <- raw_data_median %>%
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

# free up memory
rm(median_long2, calls_fl_df,  calls_long); gc()

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
  scale_fill_manual(name = "Stage", values = pal) +
  labs(title = "Boxplots of Peptide Fluorescence Levels for All Patients", 
       x = "Patient ID", y = "Median Fluorescence Levels on log2 scale") +
  theme(panel.background = element_rect(fill = "grey90"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 4.5),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))

####################################################################################### 
#      Evaluate Reproducibility of Replicates via Linear Mixed Effects Model          #         
#######################################################################################

# ncol_raw <- ncol(raw_data) 
# nrep <- nrow(array_id_key)
# 
# # initiate
# lmer_result <- matrix(NA, nrow = nrow(raw_data), ncol = 4)
# colnames(lmer_result) <- c("variance_id", "variance_residual", "lrstat", "singularity")
# 
# # check array_id_key align with raw_data_complete
# sum(as.numeric( colnames( raw_data %>% select(any_of(array_id_key$array_id)) ) == array_id_key$array_id )) == nrow(array_id_key)
# 
# 
# for(i in 1:nrow(raw_data)){
#   y <- as.numeric(raw_data[i, (ncol_raw - nrep + 1) : ncol_raw])
#   fit1 <- lmer(y ~ stage + (1|id), data = array_id_key)
#   fit2 <- lmer(y ~ 1 + (1|id), data = array_id_key)
#   lmer_result[i,] <- c(
#     as.data.frame(VarCorr(fit1))$'vcov',
#     as.numeric(-2*(logLik(fit2, REML=T) - logLik(fit1, REML=T))),
#     ( isSingular(fit1) | isSingular(fit2) )
#   )
#   print(i)
# }

# check how many singular fits
sum(as.numeric(lmer_result[,'singularity']))
max(as.numeric(lmer_result[,'variance_id'][ lmer_result[,'singularity']==T ]))
min(as.numeric(lmer_result[,'variance_residual'] ))

# get estimated proportion of variances
lmer_var_ratio <- lmer_result[,'variance_id'] / ( lmer_result[,'variance_id'] + lmer_result[,'variance_residual']  )
hist(lmer_var_ratio, breaks = 100, xlab = "estimated proportion of variances",
     main = "Histogram of peptide-level proportion of random-effect variance to total variance")


####################################################################################### 
#                               TEST  -- Logistic Regression                          #
####################################################################################### 

# # make sure patients' stages align with calls
# calls_iii <- match( colnames( calls %>% select(any_of(patient_key$id)) ), patient_key$id )
# calls_stages <- (patient_key$stage)[calls_iii]
# sum(as.numeric( colnames( calls %>% select(any_of(patient_key$id)) ) == (patient_key$id)[calls_iii] )) == nrow(patient_key)
# 
# # initiate
# logreg_pval <- rep(NA, nrow(calls))
# names(logreg_pval) <- calls$PROBE_ID
# ncol_calls <- ncol(calls)
# n <- nrow(patient_key)
# 
# # compute deviance test p-values
# for(i in 1:nrow(calls)){
#   y <- as.numeric( calls[i, (ncol_calls - n + 1): ncol_calls] )
#   fit1 <- glm(y ~ calls_stages, family = binomial(link = "logit"))
#   logreg_pval[i] <- 1 - pchisq( fit1$null.deviance - fit1$deviance, df = (fit1$df.null - fit1$df.residual) )
#   print(i)
# }

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
#                 Another Tests -- one-way ANOVA & Kruskal-Wallis Test                #
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
# all_anova_pval <- rep(NA, nrow(raw_data_median))
# names(all_anova_pval) <- raw_data_median$PROBE_ID
# all_anova_mse <- rep(NA, nrow(raw_data_median))
# names(all_anova_mse) <- raw_data_median$PROBE_ID

# compute one-way anova p-values
# for(i in 1:nrow(raw_data_median)){
#   fit1 <- lm( as.numeric(raw_data_median[i, (ncol_median - n + 1) : ncol_median]) ~ median_stage )
#   all_anova_pval[i] <- unname( unlist(summary(aov(fit1)))["Pr(>F)1"] )
#   all_anova_mse[i] <- deviance(fit1)/df.residual(fit1)
#   print(i)
# }

# get p-values histogram and FDR 
all_anova <- anova_func(all_anova_pval, "one-way ANOVA p-values")

# peptide counts at various FDR thresholds
count.func(all_anova$anova_BH, seq(0.01, 0.1, by = 0.01))
count.func(all_anova$anova_qval, seq(0.01, 0.1, by = 0.01))

#---------------------------------------------------------------------------------------
# some peptides may violate ANOVA assumption 

par(mar=c(3.1, 5.1, 2.1, 2.1), mgp=c(2, 0.8, 0))
graphics::boxplot(as.numeric( raw_data_median %>%
                                filter(PROBE_ID == "1324_KIAA1430_57587;185") %>%
                                select(any_of(array_id_key$id)) ) ~ median_stage, 
                  varwidth = T, horizontal = T, las = 1,
                  col = c("navy", "cornflowerblue", "turquoise1","darkorange1", "firebrick1"),
                  main = "Peptide ID: 1324_KIAA1430_57587;185", xlab = "log2(fluorescence)", ylab = "")

par(mar=c(3.1, 5.1, 2.1, 2.1), mgp=c(2, 0.8, 0))
graphics::boxplot(as.numeric( raw_data_median %>%
                                filter(PROBE_ID == "459_CLTC_1213;1421") %>%
                                select(any_of(array_id_key$id)) ) ~ median_stage, 
                  varwidth = T, horizontal = T, las= 1,
                  col = c("navy", "cornflowerblue", "turquoise1","darkorange1", "firebrick1"),
                  main = "Peptide ID: 459_CLTC_1213;1421", xlab = "log2(fluorescence)", ylab = "")

dev.off()

#---------------------------------------------------------------------------------------
# now do kruskal-wallis tests

# initiate kruskal-wallis(KW)
# all_kw_pval <- rep(NA, nrow(raw_data_median))
# names(all_kw_pval) <- raw_data_median$PROBE_ID
# 
# for(i in 1:nrow(raw_data_median)){
#   all_kw_pval[i] <- kruskal.test( as.numeric(raw_data_median[i, (ncol_median - n + 1) : ncol_median]) ~ median_stage )$'p.value' 
#   print(i)
# }

# get p-values histogram and FDR 
all_kw <- anova_func(all_kw_pval, "Kruskal-Wallis p-values")

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

#---------------------------------------------------------------------------------------
# PCA after Kruskal-Wallis

anova_dat <- raw_data_median %>% 
  select(PROBE_ID, any_of(patient_key$id)) %>% 
  mutate(anova_BH = all_kw$anova_BH[raw_data_median$PROBE_ID]) %>%
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
par(mfrow = c(1,2), pty = "s", mar = c(2.2,2.3,1.5,0.45), mgp = c(1.6,0.4,0),
    cex.axis = 0.84, cex.lab = 0.84, cex.main = 0.84, tcl = -0.4)
PCload.func(cols, shapes, U, D, 1, 2, pca.var, title = "PC2 vs PC1") # PC loadings (PC2 vs PC1)
legend('topright', pch = shp, col = pal, cex = 0.5,
       c("normal",  "new_dx", "nmCSPC", "mCSPC", "nmCRPC", "mCRPC") )
PCload.func(cols, shapes, U, D, 3, 2, pca.var, title = "PC2 vs PC3") # PC loadings (PC2 vs PC3)

dev.off()


####################################################################################### 
#                                     Pairwise Comparisons                            #
#######################################################################################

# we want the following contrasts:
# consecutive-group comparison: mCRPC-nmCRPC, nmCRPC-nmCSPS,nmCSPC-new_dx, new_dx-normal
# normal vs canceer
# mCRPC vs the others

# first make sure stage aligns with raw_data_median
sum(as.numeric(colnames(raw_data_median %>% select(-(PROBE_DESIGN_ID:DESIGN_ID))) == patient_key$id)) == 
  nrow(patient_key)

# set BH-FDR cutoff
BH_FDR_cutoff <- .05

# get group medians
group_median.func <- function(group, BH_filter){
  raw_data_median %>% 
    select(any_of(array_id_key$id)) %>% 
    select(which(patient_key$stage %in% group)) %>%
    # filter(all_anova$anova_BH <= BH_filter) %>%
    filter(all_kw$anova_BH <= BH_filter) %>%
    as.matrix() %>%
    matrixStats::rowMedians() 
}
mCRPC_median <- group_median.func("mCRPC", BH_FDR_cutoff)
nmCRPC_median <- group_median.func("nmCRPC", BH_FDR_cutoff)
nmCSPC_median <- group_median.func("nmCSPC", BH_FDR_cutoff)
newdx_median <- group_median.func("new_dx", BH_FDR_cutoff)
normal_median <- group_median.func("normal", BH_FDR_cutoff)
cancer_median <- group_median.func(c("new_dx", "nmCSPC", "nmCRPC", "mCRPC"), BH_FDR_cutoff)
NOT_mCRPC_median <- group_median.func(c("normal", "new_dx", "nmCSPC", "nmCRPC"), BH_FDR_cutoff)

# might need group means?
# group_means.func <- function(group, BH_filter){
#   raw_data_median %>% 
#     select(any_of(array_id_key$id)) %>% 
#     select(which(patient_key$stage %in% group)) %>%
#     # filter(all_anova$anova_BH <= BH_filter) %>%
#     filter(all_kw$anova_BH <= BH_filter) %>%
#     rowMeans() 
# }

#----------------------------------------------------------------------------------------------
# wilcoxon rank-sum tests

# median_subset <- raw_data_median %>%
#   filter(all_kw$anova_BH <= BH_FDR_cutoff) %>%  # change KW BH threshold here!
#   select(any_of(array_id_key$id)) %>%
#   as.matrix()
# 
# # initiate wilcox-pval
# mCRPC_nmCRPC_wilcox_pval <- rep(NA, nrow(median_subset))
# nmCRPC_nmCSPC_wilcox_pval <- rep(NA, nrow(median_subset))
# nmCSPC_newdx_wilcox_pval <- rep(NA, nrow(median_subset))
# newdx_normal_wilcox_pval <- rep(NA, nrow(median_subset))
# cancer_normal_wilcox_pval <- rep(NA, nrow(median_subset))
# mCRPC_others_wilcox_pval <- rep(NA, nrow(median_subset))
# 
# # get wilcox pval (2-sided)
# for(i in 1: nrow(median_subset)){
#   mCRPC_nmCRPC_wilcox_pval[i] <- wilcox.test(
#     x = as.numeric(median_subset[i, patient_key$stage == "mCRPC"]),
#     y = as.numeric(median_subset[i, patient_key$stage == "nmCRPC"]),
#     alternative = "two.sided", exact = T
#   )$'p.value'
#   nmCRPC_nmCSPC_wilcox_pval[i] <- wilcox.test(
#     x = as.numeric(median_subset[i, patient_key$stage == "nmCRPC"]),
#     y = as.numeric(median_subset[i, patient_key$stage == "nmCSPC"]),
#     alternative = "two.sided", exact = T
#   )$'p.value'
#   nmCSPC_newdx_wilcox_pval[i] <- wilcox.test(
#     x = as.numeric(median_subset[i, patient_key$stage == "nmCSPC"]),
#     y = as.numeric(median_subset[i, patient_key$stage == "new_dx"]),
#     alternative = "two.sided", exact = T
#   )$'p.value'
#   newdx_normal_wilcox_pval[i] <- wilcox.test(
#     x = as.numeric(median_subset[i, patient_key$stage == "new_dx"]),
#     y = as.numeric(median_subset[i, patient_key$stage == "normal"]),
#     alternative = "two.sided", exact = T
#   )$'p.value'
#   cancer_normal_wilcox_pval[i] <-  wilcox.test(
#     x = as.numeric(median_subset[i, patient_key$stage %in% c("new_dx", "nmCSPC", "nmCRPC", "mCRPC")]),
#     y = as.numeric(median_subset[i, patient_key$stage == "normal"]),
#     alternative = "two.sided", exact = T
#   )$'p.value'
#   mCRPC_others_wilcox_pval[i] <- wilcox.test(
#     x = as.numeric(median_subset[i, patient_key$stage == "mCRPC"]),
#     y = as.numeric(median_subset[i, patient_key$stage %in% c("normal", "new_dx", "nmCSPC", "nmCRPC")]),
#     alternative = "two.sided", exact = T
#   )$'p.value'
#   print(i)
# }

#----------------------------------------------------------------------------------------------
# get the kruskal-wallis 5% BH FDR peptide counts 
kw_FDR5prct <- length(all_kw$anova_BH[all_kw$anova_BH <= .05])

# pval histograms 
wilcox_pval_df <- data.frame(
  group_pair = c(
    rep("mCRPC vs nmCRPC", kw_FDR5prct),
    rep("nmCRPC vs nmCSPC", kw_FDR5prct),
    rep("nmCSPC vs new_dx", kw_FDR5prct),
    rep("new_dx vs normal", kw_FDR5prct),
    rep("cancer vs normal", kw_FDR5prct),
    rep("mCRPC vs others", kw_FDR5prct)
  ),
  p_values = c(
    mCRPC_nmCRPC_wilcox_pval,
    nmCRPC_nmCSPC_wilcox_pval,
    nmCSPC_newdx_wilcox_pval,
    newdx_normal_wilcox_pval,
    cancer_normal_wilcox_pval,
    mCRPC_others_wilcox_pval
  )
) 
wilcox_pval_df$group_pair <- factor(wilcox_pval_df$group_pair, levels = c(
  "cancer vs normal", "mCRPC vs others",
  "mCRPC vs nmCRPC", "nmCRPC vs nmCSPC",
  "nmCSPC vs new_dx", "new_dx vs normal"
))
ggplot(wilcox_pval_df, aes(x = p_values)) +
  geom_histogram(aes(y=..density..), bins = 50) +
  facet_wrap(. ~ group_pair, ncol=2) +
  labs(x = "Wilcoxon p-values", title = paste0("Density Histograms of the ", kw_FDR5prct, " Wilcoxon p-values")) + 
  theme(plot.title = element_text(hjust = 0.5))


#----------------------------------------------------------------------------------------------

# get BH-corrected pval (restricted to signif peptides from Kruskal-Wallis)
mCRPC_nmCRPC_wilcox_BH <- p.adjust(mCRPC_nmCRPC_wilcox_pval, method = "BH")
nmCRPC_nmCSPC_wilcox_BH <- p.adjust(nmCRPC_nmCSPC_wilcox_pval, method = "BH")
nmCSPC_newdx_wilcox_BH <- p.adjust(nmCSPC_newdx_wilcox_pval, method = "BH")
newdx_normal_wilcox_BH <- p.adjust(newdx_normal_wilcox_pval, method = "BH")
cancer_normal_wilcox_BH <- p.adjust(cancer_normal_wilcox_pval, method = "BH")
mCRPC_others_wilcox_BH <- p.adjust(mCRPC_others_wilcox_pval, method = "BH")

wilcox_peptide_counts_df <- data.frame(
  pairwise_comparison = c(
    "cancer vs normal",
    "mCRPC vs others",
    "mCRPC vs nmCRPC", 
    "nmCRPC vs nmCSPC",
    "nmCSPC vs new_dx", 
    "new_dx vs normal"
  ),
  peptide_counts = c(
    length(which(abs(cancer_median - normal_median) > 1 & cancer_normal_wilcox_BH <= BH_FDR_cutoff)),
    length(which(abs(mCRPC_median - NOT_mCRPC_median) > 1 & mCRPC_others_wilcox_BH <= BH_FDR_cutoff)),
    length(which(abs(mCRPC_median - nmCRPC_median) > 1 & mCRPC_nmCRPC_wilcox_BH <= BH_FDR_cutoff)),
    length(which(abs(nmCRPC_median - nmCSPC_median) > 1 & nmCRPC_nmCSPC_wilcox_BH <= BH_FDR_cutoff)),
    length(which(abs(nmCSPC_median - newdx_median) > 1 & nmCSPC_newdx_wilcox_BH <= BH_FDR_cutoff)),
    length(which(abs(newdx_median - normal_median) > 1 & newdx_normal_wilcox_BH<= BH_FDR_cutoff))
  )
)

#----------------------------------------------------------------------------------------------
# volcano plots of contrasts

# make sure contrast_diff and contrast_pval and contrast_BH of same length !!

contrast_volcano_plot.func <- function(contrast_diff, contrast_pval, contrast_BH, test_type, measure, contrast_title, xlim = c(-7,6), ylim = c(0,8)){
  signif_counts = length(which(abs(contrast_diff) > 1 & contrast_BH <= BH_FDR_cutoff))
  horizontal_pval = max(contrast_pval[contrast_BH<= BH_FDR_cutoff])
  plot(x = contrast_diff, y = -log10(contrast_pval), pch = 20, xlim = xlim, ylim = ylim, las = 1,
       xlab = paste0("difference of ", measure, " log2(fluorescence)"), 
       ylab = paste0("-log10(contrast ",test_type, " p-values)"),
       main = paste0("Volcano plot of contrast: ", contrast_title))
  text(x = 4, y = 7, paste0(signif_counts, " signif peptides"))
  lines(x = contrast_diff[ contrast_BH <= BH_FDR_cutoff & abs(contrast_diff) >= 1 ], 
        y = -log10(contrast_pval[ contrast_BH <= BH_FDR_cutoff & abs(contrast_diff) >= 1]),
        type = "p", pch = 20, col = "red")
  abline(v = 1, col = "blue", lty = 2, lwd = 2)
  abline(v = -1, col = "blue", lty = 2, lwd = 2)
  abline(h = -log10(horizontal_pval), col = "blue", lty = 2, lwd = 2)
}

# png("09_contrast_volcano_plots(KW_wilcox_diffmedian_cutoff).png", width = 1024, height = 1024)
par(mfrow = c(3,2))
# cancer vs normal
contrast_volcano_plot.func(cancer_median - normal_median, 
                           cancer_normal_wilcox_pval, 
                           cancer_normal_wilcox_BH, 
                           "Wilcoxon",
                           "median",
                           "cancer vs normal")

# mCRPC vs others
contrast_volcano_plot.func(mCRPC_median - NOT_mCRPC_median, 
                           mCRPC_others_wilcox_pval, 
                           mCRPC_others_wilcox_BH, 
                           "Wilcoxon",
                           "median",
                           "mCRPC vs others" )

# mCRPC vs nmCRPC
contrast_volcano_plot.func(mCRPC_median - nmCRPC_median, 
                           mCRPC_nmCRPC_wilcox_pval, 
                           mCRPC_nmCRPC_wilcox_BH, 
                           "Wilcoxon",
                           "median",
                           "mCRPC vs nmCRPC")

# nmCRPC vs nmCSPC
contrast_volcano_plot.func(nmCRPC_median - nmCSPC_median, 
                           nmCRPC_nmCSPC_wilcox_pval, 
                           nmCRPC_nmCSPC_wilcox_BH, 
                           "Wilcoxon",
                           "median",
                           "nmCRPC vs nmCSPC")

# nmCSPC vs new_dx
contrast_volcano_plot.func(nmCSPC_median - newdx_median, 
                           nmCSPC_newdx_wilcox_pval, 
                           nmCSPC_newdx_wilcox_BH, 
                           "Wilcoxon",
                           "median",
                           "nmCSPC vs new_dx")

# new_dx vs normal
contrast_volcano_plot.func(newdx_median - normal_median, 
                           newdx_normal_wilcox_pval, 
                           newdx_normal_wilcox_BH, 
                           "Wilcoxon",
                           "median",
                           "new_dx vs normal")

dev.off()


####################################################################################### 
#                                   Visualization                                     #
#######################################################################################

# replot subset of ANOVA residual heatmap based on effect-size threshold from post-hoc analysis

posthoc_signif_crit <- ( abs(mCRPC_median - nmCRPC_median) > 1 & mCRPC_nmCRPC_wilcox_BH <= BH_FDR_cutoff ) |
  ( abs(nmCRPC_median - nmCSPC_median) > 1 & nmCRPC_nmCSPC_wilcox_BH <= BH_FDR_cutoff ) |
  ( abs(nmCSPC_median - newdx_median) > 1 & nmCSPC_newdx_wilcox_BH <= BH_FDR_cutoff ) |
  ( abs(newdx_median - normal_median) > 1 & newdx_normal_wilcox_BH<= BH_FDR_cutoff ) |
  ( abs(cancer_median - normal_median) > 1 & cancer_normal_wilcox_BH <= BH_FDR_cutoff ) |
  ( abs(mCRPC_median - NOT_mCRPC_median) > 1 & mCRPC_others_wilcox_BH <= BH_FDR_cutoff )

# check
length(posthoc_signif_crit)
sum(as.numeric(posthoc_signif_crit))
secondary_cutoff_counts <- sum(as.numeric(posthoc_signif_crit))

anova_dat <- raw_data_median %>% 
  select(PROBE_ID, any_of(patient_key$id)) %>% 
  # mutate(anova_BH = all_anova$anova_BH[raw_data_median$PROBE_ID]) %>%
  mutate(anova_BH = all_kw$anova_BH[raw_data_median$PROBE_ID]) %>%
  filter(anova_BH <= BH_FDR_cutoff) %>%
  select(-anova_BH) %>%
  filter(posthoc_signif_crit)

anova_dat_demean <- sweep(as.matrix(anova_dat %>% select(-PROBE_ID)), 1, 
                          rowMeans(as.matrix(anova_dat %>% select(-PROBE_ID))), "-") # centering by row

# check
dim(anova_dat_demean)

# make sure stage aligns with anova_dat_demean
visual_iii <- match( colnames(anova_dat_demean) , patient_key$id )
visual_stage <- patient_key$stage[visual_iii]

#--------------------------------------------------------------------------------------------
# heatmap

# colors and shapes for the visualization techniques
cols = pal[ match(visual_stage, names(pal)) ]
# cls <- colorRampPalette(c("navy", "honeydew", "firebrick3", "brown"))(n = 1024)
# cls <- colorRampPalette(c("navy", "honeydew", "sienna1", "firebrick3", "brown"))(n = 1024)
cls <- colorRampPalette(c("navy", "honeydew", "firebrick1"))(n = 1024)

get_column_order.func <- function(stages){
  heat_map <- heatmap3(anova_dat_demean[,visual_stage == stages],col = cls, labRow = "", scale = "none", 
                       showColDendro = F, showRowDendro = F)
  return( colnames(anova_dat_demean[,visual_stage == stages])[heat_map$colInd] )
}

normal_id_order <- get_column_order.func("normal")
newdx_id_order <- get_column_order.func("new_dx")
nmCSPC_id_order <- get_column_order.func("nmCSPC")
nmCRPC_id_order <- get_column_order.func("nmCRPC")
mCRPC_id_order <- get_column_order.func("mCRPC")

id_order <- match( c(normal_id_order, newdx_id_order, nmCSPC_id_order, nmCRPC_id_order, mCRPC_id_order),
                   colnames(anova_dat_demean) )

# winsorize
quantile(as.numeric(anova_dat_demean), probs = c(.01, .05, .1, .15, .2, .25, .3, .7, .75, .8, .85, .9, .95)) # check
ecdf(as.numeric(anova_dat_demean))(c(-3.5, -2.6, -2, 2, 2.6, 3.5)) # check
anova_dat_demean_winsorize <- t( apply(anova_dat_demean, 1, function(x){
  x[x < -2] = -2
  x[x > 2] = 2
  return(x)
}) )

# png("09_residual_heatmap(KW_wilcox_diffmedian_cutoff).png", width = 1024, height = 1024)
heatmap3(anova_dat_demean_winsorize[,id_order], 
         col = cls, # specify colors 
         ColSideColors = cols[id_order], # specify patient color code
         Colv = NA,
         scale = "none", # no scaling by row
         labCol = visual_stage[id_order], # specify patient
         ColSideLabs = "stages", 
         labRow = "",
         xlab = "Patients",
         legendfun=function() showLegend(col = c("navy", "cornflowerblue", "turquoise1", "darkorange1", "firebrick1"),
                                         legend = c("normal",  "new_dx", "nmCSPC", "nmCRPC", "mCRPC"),
                                         cex = 1.2,
                                         lwd = 5  )
)
dev.off()

#--------------------------------------------------------------------------------------------
# PCA

# colors and shapes for the visualization techniques
cols = pal[ match(visual_stage, names(pal)) ]
shapes = shp[  match(visual_stage, names(shp)) ]

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
par(mfrow = c(1,2), pty = "s", mar = c(2.2,2.3,1.5,0.45), mgp = c(1.6,0.4,0),
    cex.axis = 0.84, cex.lab = 0.84, cex.main = 0.84, tcl = -0.4)
PCload.func(cols, shapes, U, D, 1, 2, pca.var, title = "PC2 vs PC1") # PC loadings (PC2 vs PC1)
legend('topright', pch = shp, col = pal, cex = 0.5,
       c("normal",  "new_dx", "nmCSPC", "mCSPC", "nmCRPC", "mCRPC") )
PCload.func(cols, shapes, U, D, 3, 2, pca.var, title = "PC2 vs PC3") # PC loadings (PC2 vs PC3)

dev.off()



####################################################################################### 
#              Gene Set Analyses After Kruskal-Wallis & Wilcoxon Tests                #
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
uniprot_gene <- uniprot_gene[!(is.na(uniprot_gene$entrez_gene_id_pete)),] # just in case
length(unique(uniprot_gene$seq_id)) == length(uniprot_gene$seq_id) # yes! unique!
length(unique(uniprot_gene$entrez_gene_id_pete)) == length(uniprot_gene$entrez_gene_id_pete) # NOT unique

# which seq_id has repeated gene_symbol
genesymb_repeat <- as.data.frame( uniprot_gene[ uniprot_gene$entrez_gene_id_pete %in%  
                                   (uniprot_gene %>% 
                                      group_by(entrez_gene_id_pete) %>% 
                                      tally() %>% 
                                      filter(n > 1) %>% 
                                      pull(entrez_gene_id_pete)) , ] )

entrez_id_repeat <- as.character( unique(genesymb_repeat$entrez_gene_id_pete) )

# # check if these seq_id make the Kruskal-Wallis 5% BH FDR cutoff ?
# genesymb_repeat$seq_id %in% raw_data_median$SEQ_ID[all_kw$anova_BH <= BH_FDR_cutoff] # yes...fine...
# # check if these seq_id make at least one of the secondary cutoffs ?
# genesymb_repeat$seq_id %in% raw_data_median$SEQ_ID[all_kw$anova_BH <= BH_FDR_cutoff][posthoc_signif_crit] # yes...fine
# # check each contrast one by one
# genesymb_repeat$seq_id %in% raw_data_median$SEQ_ID[all_kw$anova_BH <= BH_FDR_cutoff][( abs(mCRPC_median - nmCRPC_median) > 1 & mCRPC_nmCRPC_wilcox_BH <= BH_FDR_cutoff ) ]
# genesymb_repeat$seq_id %in% raw_data_median$SEQ_ID[all_kw$anova_BH <= BH_FDR_cutoff][( abs(nmCRPC_median - nmCSPC_median) > 1 & nmCRPC_nmCSPC_wilcox_BH <= BH_FDR_cutoff ) ]
# genesymb_repeat$seq_id %in% raw_data_median$SEQ_ID[all_kw$anova_BH <= BH_FDR_cutoff][( abs(nmCSPC_median - newdx_median) > 1 & nmCSPC_newdx_wilcox_BH <= BH_FDR_cutoff ) ]
# genesymb_repeat$seq_id %in% raw_data_median$SEQ_ID[all_kw$anova_BH <= BH_FDR_cutoff][( abs(newdx_median - normal_median) > 1 & newdx_normal_wilcox_BH<= BH_FDR_cutoff ) ]
# genesymb_repeat$seq_id %in% raw_data_median$SEQ_ID[all_kw$anova_BH <= BH_FDR_cutoff][( abs(cancer_median - normal_median) > 1 & cancer_normal_wilcox_BH <= BH_FDR_cutoff ) ]
# genesymb_repeat$seq_id %in% raw_data_median$SEQ_ID[all_kw$anova_BH <= BH_FDR_cutoff][( abs(mCRPC_median - NOT_mCRPC_median) > 1 & mCRPC_others_wilcox_BH <= BH_FDR_cutoff )]


# DECISION: treat them as repeats in the microarray
# protein deemed signif if either one of the repeated seq_id makes the cutoff

get_SeqID.func <- function(diff, BH){
  signif_seq_id <- unique( raw_data_median$SEQ_ID[all_kw$anova_BH <= BH_FDR_cutoff][abs(diff) > 1 & BH <= BH_FDR_cutoff] )
  seq_id_ok <- as.numeric(uniprot_gene$seq_id %in% signif_seq_id)
  names(seq_id_ok) <- uniprot_gene$entrez_gene_id_pete
  seq_id_ok2 <- seq_id_ok[!(names(seq_id_ok) %in% entrez_id_repeat)]
  seq_id_ok2 <- c(seq_id_ok2, sapply( entrez_id_repeat, function(x){max(seq_id_ok[which(names(seq_id_ok) == x)])} ) )
}

cancer_normal_SeqID <- get_SeqID.func(cancer_median - normal_median, cancer_normal_wilcox_BH)
mCRPC_others_SeqID <- get_SeqID.func(mCRPC_median - NOT_mCRPC_median, mCRPC_others_wilcox_BH)
mCRPC_nmCRPC_SeqID <- get_SeqID.func(mCRPC_median - nmCRPC_median, mCRPC_nmCRPC_wilcox_BH)
nmCRPC_nmCSPC_SeqID <- get_SeqID.func(nmCRPC_median - nmCSPC_median, nmCRPC_nmCSPC_wilcox_BH)
nmCSPC_newdx_SeqID <- get_SeqID.func(nmCSPC_median - newdx_median, nmCSPC_newdx_wilcox_BH)
newdx_normal_SeqID <- get_SeqID.func(newdx_median - normal_median, newdx_normal_wilcox_BH)

# gene-set analysis via allez!
cancer_normal_allez.go <- allez(cancer_normal_SeqID, lib = "org.Hs.eg", idtype = "ENTREZID", sets = "GO")
mCRPC_others_allez.go <- allez(mCRPC_others_SeqID, lib = "org.Hs.eg", idtype = "ENTREZID", sets = "GO")
mCRPC_nmCRPC_allez.go <- allez(mCRPC_nmCRPC_SeqID, lib = "org.Hs.eg", idtype = "ENTREZID", sets = "GO")
nmCRPC_nmCSPC_allez.go <- allez(nmCRPC_nmCSPC_SeqID, lib = "org.Hs.eg", idtype = "ENTREZID", sets = "GO")
nmCSPC_newdx_allez.go <- allez(nmCSPC_newdx_SeqID, lib = "org.Hs.eg", idtype = "ENTREZID", sets = "GO")
newdx_normal_allez.go <- allez(newdx_normal_SeqID, lib = "org.Hs.eg", idtype = "ENTREZID", sets = "GO")

#---------------------------------------------------------------------------------------------------
# get allez results

nom.alpha <- 0.05
max_gene_in_set <- 300

# Extract a table of top-ranked functional sets from allez output
# Display an image of gene scores by functional sets

# cancel vs normal
allezTable(cancer_normal_allez.go, symbol = T, nominal.alpha = nom.alpha, n.upp = max_gene_in_set, in.set = T)[,c(1:5,7)]
get_alleztable.func(cancer_normal_allez.go)
allezPlot(cancer_normal_allez.go, nominal.alpha = nom.alpha, n.upp = max_gene_in_set)


# mCRPC vs others
allezTable(mCRPC_others_allez.go, symbol = T, nominal.alpha = nom.alpha, n.upp = max_gene_in_set, in.set = T)[,c(1:5,7)]
get_alleztable.func(mCRPC_others_allez.go)
allezPlot(mCRPC_others_allez.go, nominal.alpha = nom.alpha, n.upp = max_gene_in_set)


# mCRPC vs nmCRPC
allezTable(mCRPC_nmCRPC_allez.go, symbol = T, nominal.alpha = nom.alpha, n.upp = max_gene_in_set, in.set = T)[,c(1:5,7)]
get_alleztable.func(mCRPC_nmCRPC_allez.go)
allezPlot(mCRPC_nmCRPC_allez.go, nominal.alpha = nom.alpha, n.upp = max_gene_in_set)


# nmCRPC vs nmCSPC
allezTable(nmCRPC_nmCSPC_allez.go, symbol = T, nominal.alpha = nom.alpha, n.upp = max_gene_in_set, in.set = T)[,c(1:5,7)]
get_alleztable.func(nmCRPC_nmCSPC_allez.go)
allezPlot(nmCRPC_nmCSPC_allez.go, nominal.alpha = nom.alpha, n.upp = max_gene_in_set)


# nmCSPC vs new_dx
allezTable(nmCSPC_newdx_allez.go, symbol = T, nominal.alpha = nom.alpha, n.upp = max_gene_in_set, in.set = T)[,c(1:5,7)]
get_alleztable.func(nmCSPC_newdx_allez.go)
allezPlot(nmCSPC_newdx_allez.go, nominal.alpha = nom.alpha, n.upp = max_gene_in_set)


# new_dx vs normal
allezTable(newdx_normal_allez.go, symbol = T, nominal.alpha = nom.alpha, n.upp = max_gene_in_set, in.set = T)[,c(1:5,7)]
get_alleztable.func(newdx_normal_allez.go)
allezPlot(newdx_normal_allez.go, nominal.alpha = nom.alpha, n.upp = max_gene_in_set)


####################################################################################### 
#                               Data Processing -- part II                            #
####################################################################################### 
pal_proj2 <- c("turquoise1", "cornflowerblue","navy", "orchid1", "darkorange1", "firebrick1")
names(pal_proj2) <- c("ADT_time:0", "ADT_time:3", "ADT_time:6", "PAP_time:0", "PAP_time:3", "PAP_time:6")

array_id_key_proj2 = read_tsv("sample_key_project2.txt") %>%
  janitor::clean_names() %>% 
  rename(treatment = condition,
         array_id = "file_name") %>%
  mutate(id = tolower(id))

sample_key_proj2 = array_id_key_proj2 %>%
  group_by(id, time) %>% 
  summarize(n = n()) %>%
  ungroup() %>%
  filter(n>=2) 

# check
table(sample_key_proj2$n)
## all patients associated with 3 replicates for each of the 3 time points

# remove patients with no replicates
array_id_key_proj2 = array_id_key_proj2 %>% 
  filter(id %in% sample_key_proj2$id)

# make sure sample_key_proj2 align with raw_data_median_proj2
sample_key_proj2 <- sample_key_proj2 %>%
  select(-n) %>%
  arrange(time, id)
sample_key_proj2 <- sample_key_proj2 %>% 
  left_join(array_id_key_proj2 %>% select(id, treatment) %>% distinct())

# check
sum(as.numeric( colnames(raw_data_median_proj2 %>% select(contains("time:"))) ==
                  paste( paste0("id:", sample_key_proj2$id), paste0("time:", sample_key_proj2$time), sep = "_" ) )) == 
  ncol(raw_data_median_proj2 %>% select(contains("time:")))

#----------------------------------------------------------------------------------------
# read fluorescence data

# raw_data = read_csv("raw_data_complete.csv")
# 
# raw_data = raw_data %>%
#   select(PROBE_DESIGN_ID:Y, any_of(array_id_key_proj2$array_id)) %>% # drop patients with records at different stages 
#   select( -c(X, Y, MATCH_INDEX, DESIGN_NOTE, SELECTION_CRITERIA, MISMATCH, PROBE_CLASS) ) # drop unhelpful columns
# 
# # check
# colnames(raw_data)[1:15]
# 
# # take log2 transformation
# raw_data <- raw_data %>%
#   mutate_at(vars(matches("dat")), log2)
# 
# # check any NA's
# raw_data %>% 
#   select_if(function(x) any(is.na(x))) 
# 
# # make sure array_id_key_proj2 align with raw_data_complete
# array_iii <- match( colnames( raw_data %>% select(any_of(array_id_key_proj2$array_id)) ), array_id_key_proj2$array_id )
# array_id_key_proj2 <- array_id_key_proj2[array_iii , ]
# sum(as.numeric( colnames( raw_data %>% select(any_of(array_id_key_proj2$array_id)) ) 
#                 == array_id_key_proj2$array_id )) == nrow(array_id_key_proj2)
# 
# 
# # compute median
# raw_data_median_proj2 <- t( apply( select(raw_data, contains("dat")), 1, function(x) {
#   as.vector( tapply( as.numeric(x), list(array_id_key_proj2$id, array_id_key_proj2$time), FUN=median) )
# } ) )
# raw_data_median_proj2 <- bind_cols( select(raw_data, -contains("dat")), as.data.frame(raw_data_median_proj2) ) 
# 
# colnames(raw_data_median_proj2)[1:15] # check
# 
# # want to rename column
# example_row <- tapply( as.numeric(select(raw_data, contains("dat"))[1,]), 
#                        list(array_id_key_proj2$id, array_id_key_proj2$time), FUN=median) %>% 
#   as.data.frame() %>% 
#   rownames_to_column("id") %>% 
#   pivot_longer(-id, names_to = "time", values_to = "fluorescence") %>% 
#   arrange(time, id)
# sum(as.numeric( example_row$fluorescence == select(raw_data_median_proj2, contains("V"))[1,] )) == nrow(example_row) 
# example_row$column_name <- paste( paste0("id:", example_row$id), paste0("time:", example_row$time), sep = "_" )
# raw_data_median_proj2 <- raw_data_median_proj2 %>% rename_at(vars(contains("V")), ~ example_row$column_name) 
# 
# # make sure sample_key_proj2 align with raw_data_median_proj2
# sample_key_proj2 <- sample_key_proj2 %>%
#   select(-n) %>%
#   arrange(time, id)
# sum(as.numeric(sample_key_proj2$id == example_row$id)) == nrow(sample_key_proj2) # check
# sum(as.numeric(sample_key_proj2$time == example_row$time)) == nrow(sample_key_proj2) # check
# sample_key_proj2 <- sample_key_proj2 %>% 
#   left_join(array_id_key_proj2 %>% select(id, treatment) %>% distinct())
# 
# write.table(raw_data_median_proj2, file = "raw_data_median_proj2.csv", sep = ",", row.names = F)

####################################################################################### 
#                           Check Fluorescence Normalization                          #
####################################################################################### 

median_long_proj2 <- raw_data_median_proj2 %>%
  select(contains("id:")) %>%
  pivot_longer(cols = everything(), names_to = "id_time", values_to = "fluorescence") 

# check
nrow(median_long_proj2) == nrow(raw_data_median_proj2) * nrow(sample_key_proj2)
head(median_long_proj2)

# set fill color
median_long_proj2$treat_time <- rep( paste(sample_key_proj2$treatment, sample_key_proj2$time, sep = "_"), 177604)
median_long_proj2  <- median_long_proj2 %>%
  mutate(treat_time = as_factor(treat_time),
         treat_time = fct_recode(treat_time,
                            "ADT_time:0" = "ADT_0",
                            "ADT_time:3" = "ADT_3",
                            "ADT_time:6" = "ADT_6",
                            "PAP_time:0" = "Vaccine_0",
                            "PAP_time:3" = "Vaccine_3",
                            "PAP_time:6" = "Vaccine_6"))

# sort order of patients in boxplot
median_long_proj2$id_time <- factor(median_long_proj2$id_time, levels = c(
  unique(median_long_proj2$id_time[median_long_proj2$treat_time == "ADT_time:0"]),
  unique(median_long_proj2$id_time[median_long_proj2$treat_time == "ADT_time:3"]),
  unique(median_long_proj2$id_time[median_long_proj2$treat_time == "ADT_time:6"]),
  unique(median_long_proj2$id_time[median_long_proj2$treat_time == "PAP_time:0"]),
  unique(median_long_proj2$id_time[median_long_proj2$treat_time == "PAP_time:3"]),
  unique(median_long_proj2$id_time[median_long_proj2$treat_time == "PAP_time:6"])
))

ggplot(median_long_proj2, aes(x = id_time, y = fluorescence, fill = treat_time)) +
  geom_boxplot(outlier.shape = ".") +
  scale_fill_manual(name = "treatment_time", values = pal_proj2[levels(median_long_proj2$treat_time)]) +
  labs(title = "Boxplots of Peptide Fluorescence Levels for Patients at 3 time points", 
       x = "Patient_Time", y = "Median log2 Fluorescence Levels") +
  theme(panel.background = element_rect(fill = "grey90"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 4.3),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))

####################################################################################### 
#                    Linear Mixed Model to Assess Time Effect (REML = TRUE)           #
#                       Separate Models for PAP and ADT Groups                        #
####################################################################################### 

# ncol_med_proj2 <- ncol(raw_data_median_proj2) 
# n_med_proj2 <- nrow(sample_key_proj2)
# 
# # initiate
# PAP_resid_fit0 <- matrix(NA, nrow = nrow(raw_data_median_proj2), 
#                          ncol = nrow(sample_key_proj2%>%filter(treatment == "Vaccine")))
# ADT_resid_fit0 <- matrix(NA, nrow = nrow(raw_data_median_proj2), 
#                          ncol = nrow(sample_key_proj2%>%filter(treatment == "ADT")))
# PAP_result <- matrix(NA, nrow = nrow(raw_data_median_proj2), ncol = 8)
# ADT_result <- matrix(NA, nrow = nrow(raw_data_median_proj2), ncol = 8)
# 
# colnames(PAP_resid_fit0) <- colnames(raw_data_median_proj2 %>% select(contains("id:pap")))
# colnames(ADT_resid_fit0) <- colnames(raw_data_median_proj2 %>% select(contains("id:adt")))
# colnames(PAP_result) <- paste0("PAP_", c(
#   "time_effect", 
#   "time_tstat",
#   "KR_df",
#   "KR_Ftest_pval",
#   "Satterthwaite_df",
#   "Satterthwaite_Ftest_pval",
#   "zval_1sided_KR",
#   "zval_1sided_Satterthwaite"
# ))
# colnames(ADT_result) <- paste0("ADT_", c(
#   "time_effect", 
#   "time_tstat",
#   "KR_df",
#   "KR_Ftest_pval",
#   "Satterthwaite_df",
#   "Satterthwaite_Ftest_pval",
#   "zval_1sided_KR",
#   "zval_1sided_Satterthwaite"
# ))
# 
# Test_Time.func <- function(y, treat_type){
#   resp <- y[sample_key_proj2$treatment == treat_type]
#   fit1 <- lmer(resp ~ time + (1 + time | id), REML = T, 
#                data = sample_key_proj2[sample_key_proj2$treatment == treat_type,])
#   fit0 <- lmer(resp ~ 1 + (1 + time | id), REML = T, 
#                data = sample_key_proj2[sample_key_proj2$treatment == treat_type,])
#   resid_fit0 <- unname(round(resid(fit0),4))
#   effect_tstat <- coef(summary(fit1))['time',c('Estimate', 't value')] 
#   KR_df_pval <- contest(fit1, c(0,1), ddf = "Kenward-Roger")[c('DenDF', 'Pr(>F)')] 
#   Satterthwaite_df_pval <- contest(fit1, c(0,1))[c('DenDF', 'Pr(>F)')] 
#   zval_1sided_KR <- qnorm(pt( as.numeric(effect_tstat['t value']), 
#                               df = as.numeric(KR_df_pval['DenDF']) ,  lower.tail = T ))
#   zval_1sided_Satterthwaite <- qnorm(pt( as.numeric(effect_tstat['t value']), 
#                                          df = as.numeric(Satterthwaite_df_pval['DenDF']) ,  lower.tail = T ))
#   result <- c(
#     as.numeric(effect_tstat), 
#     as.numeric(KR_df_pval),
#     as.numeric(Satterthwaite_df_pval),
#     zval_1sided_KR,
#     zval_1sided_Satterthwaite
#   )
#   return( list(
#     resid_fit0 = resid_fit0,
#     result = result
#   ) )
# }
# 
# 
# for(i in 1:nrow(raw_data_median_proj2)){
#   y <- as.numeric(raw_data_median_proj2[i, (ncol_med_proj2 - n_med_proj2 + 1) : ncol_med_proj2])
#   PAP_test <- Test_Time.func(y, "Vaccine")
#   ADT_test <- Test_Time.func(y, "ADT")
#   
#   PAP_resid_fit0[i,] <- PAP_test$resid_fit0
#   PAP_result[i,] <- PAP_test$result
#   
#   ADT_resid_fit0[i,] <- ADT_test$resid_fit0
#   ADT_result[i,] <- ADT_test$result
#   
#   if(i %% 100 == 0){
#     print(i)
#   }
# }
# 
# save(PAP_resid_fit0, ADT_resid_fit0, PAP_result, ADT_result,
#      file = "08_LMER_results.RData")

PAP_Satterth_Ftest_pval <- PAP_result[,"PAP_Satterthwaite_Ftest_pval"]
PAP_KR_Ftest_pval <- PAP_result[,"PAP_KR_Ftest_pval"]
ADT_Satterth_Ftest_pval <- ADT_result[,"ADT_Satterthwaite_Ftest_pval"]
ADT_KR_Ftest_pval <- ADT_result[,"ADT_KR_Ftest_pval"]

PAP_Ftest_KR_BH <- p.adjust(PAP_KR_Ftest_pval, method="BH")
PAP_Ftest_Satterthwaite_BH <- p.adjust(PAP_Satterth_Ftest_pval, method="BH")
ADT_Ftest_KR_BH <- p.adjust(ADT_KR_Ftest_pval,method="BH")
ADT_Ftest_Satterthwaite_BH <- p.adjust(ADT_Satterth_Ftest_pval,method="BH")

#---------------------------------------------------------------------------------------------

# F-test p-values based on KR adjustments more conservative than Satterthwaite

par(mfrow=c(1,2))
plot(PAP_Satterth_Ftest_pval[PAP_Satterth_Ftest_pval <= .2 & PAP_KR_Ftest_pval <= .2], 
     PAP_KR_Ftest_pval[PAP_Satterth_Ftest_pval <= .2 & PAP_KR_Ftest_pval <= .2], 
     pch = ".", xlim = c(0,.2), ylim = c(0,.2), las = 1, 
     main = "Time Fixed Effect p-values \nfor PAP patients",
     xlab = "Satterthwaite F-test p-values", ylab = "Kenward-Roger (KR) F-test p-values")
abline(a=0, b=1, col = "red", lty=2, lwd = 2)

plot(ADT_Satterth_Ftest_pval[ADT_Satterth_Ftest_pval <= .2 & ADT_KR_Ftest_pval <= .2], 
     ADT_KR_Ftest_pval[ADT_Satterth_Ftest_pval <= .2 & ADT_KR_Ftest_pval <= .2], 
     pch = ".", xlim = c(0,.2), ylim = c(0,.2), las = 1, 
     main = "Time Fixed Effect p-values \nfor ADT patients",
     xlab = "Satterthwaite F-test p-values", ylab = "Kenward-Roger (KR) F-test p-values")
abline(a=0, b=1, col = "red", lty=2, lwd = 2)

dev.off()


count.func(PAP_Ftest_KR_BH, seq(.01,.05,by=.01))
count.func(PAP_Ftest_Satterthwaite_BH, seq(.01,.05,by=.01))
count.func(ADT_Ftest_KR_BH, seq(.63,.7,by=.01))
count.func(ADT_Ftest_Satterthwaite_BH, seq(.63,.7,by=.01))

# check
sum(as.numeric( raw_data_median_proj2$PROBE_ID[PAP_Ftest_KR_BH <= .01] %in% 
                  raw_data_median_proj2$PROBE_ID[PAP_Ftest_Satterthwaite_BH <= .01] ))

# tabulate
Ftest_pval_counts <- data.frame(
  BH_FDR_thresholds = seq(.01, .05, by = .01),
  Peptide_counts_KR = count.func(PAP_Ftest_KR_BH, seq(.01,.05,by=.01))[2,2:6],
  Peptide_counts_Satterthwaite = count.func(PAP_Ftest_Satterthwaite_BH, seq(.01,.05,by=.01))[2,2:6]
)

#---------------------------------------------------------------------------------------------
#p-value density histograms

KR_Ftest_pval_df <- data.frame(
  treatment = c(
    rep("PAP", 177604*2),
    rep("ADT", 177604*2)
  ), 
  method = rep( rep( c("Satterthwaite", "Kenward-Roger"), each = 177604 ), 2) ,
  p_values = c(
    PAP_Satterth_Ftest_pval,
    PAP_KR_Ftest_pval,
    ADT_Satterth_Ftest_pval,
    ADT_KR_Ftest_pval
  )
)

ggplot(KR_Ftest_pval_df, aes(x = p_values, fill = method)) +
  geom_histogram(aes(y=..density..), bins = 100, position = "identity", alpha = .4) +
  facet_grid(. ~ treatment) +
  labs(x = "F-test p-values", title = paste0("Density Histograms of F-test p-values")) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

####################################################################################### 
#                                   Visualization                                     #
#######################################################################################

BH_FDR_cutoff_proj2 <- .05
signif_crit_proj2 <- (PAP_Ftest_KR_BH <= BH_FDR_cutoff_proj2) & 
  (PAP_Ftest_Satterthwaite_BH <= BH_FDR_cutoff_proj2) &
  (PAP_result[,"PAP_time_effect"] > .3333)

sum(as.numeric(signif_crit_proj2))
proj2_signif_count <- sum(as.numeric(signif_crit_proj2))

#---------------------------------------------------------------------------------------------
# volcano plots

proj2_volcano_plot.func <- function(time_effect, pval, BH, title){
  plot(x = time_effect, y = -log10(pval), pch = ".", las = 1,
       ylim = c(0,9), xlim = c(-.4, .9),
       xlab = "coefficient of time fixed effect", ylab = "-log10(KR F-test p-values)",
       main = title)
  lines(x = time_effect[ BH <= .01 & time_effect >= .3333 ], 
        y = -log10(pval[ BH <= .01 & time_effect >= .3333]),
        type = "p", pch = ".", col = "red")
}

par(mfrow=c(1,2))
proj2_volcano_plot.func(PAP_result[,"PAP_time_effect"], PAP_result[,"PAP_KR_Ftest_pval"], PAP_Ftest_KR_BH,
                  "PAP's volcano plot")
abline(v = .3333, lty = 2, lwd = 1.5, col = "blue")
proj2_volcano_plot.func(ADT_result[,"ADT_time_effect"], ADT_result[,"ADT_KR_Ftest_pval"], ADT_Ftest_KR_BH,
                  "ADT's volcano plot")

#---------------------------------------------------------------------------------------------
# heatmap

proj2_resid <- PAP_resid_fit0[signif_crit_proj2,]
dim(proj2_resid) # check
proj2_resid_time <- sample_key_proj2$time[sample_key_proj2$treatment=="Vaccine"]

# specify color scheme
cls <- colorRampPalette(c("navy", "honeydew", "firebrick1"))(n = 1024)
proj2_heatmap_pal <- c("lightgoldenrod1", "darkorange1", "brown")
names(proj2_heatmap_pal) <- c(0,3,6)
cols <- proj2_heatmap_pal[ match(proj2_resid_time, names(proj2_heatmap_pal)) ]

# winsorize
quantile(as.numeric(proj2_resid), probs = c(.01, .05, .1, .15, .2, .25, .3, .7, .75, .8, .85, .9, .95)) # check
ecdf(as.numeric(proj2_resid))(c(-1.7, -1.4, 1.4, 1.7)) # check
proj2_resid_winsorize <- t( apply(proj2_resid, 1, function(x){
  x[x < -1.7] = -1.7
  x[x > 1.7] = 1.7
  return(x)
}) )


# get column order of time 6
proj2_get_column_order.func <- function(resid_mat, Time){
  heat_map <- heatmap3(resid_mat[,proj2_resid_time == Time],col = cls, labRow = "", scale = "none", 
                       showColDendro = F, showRowDendro = F)
  return( colnames(resid_mat[,proj2_resid_time == Time])[heat_map$colInd] )
}
time6_order <- proj2_get_column_order.func(proj2_resid_winsorize,6)
time0_order <- gsub("time:6", "time:0", time6_order)
time3_order <- gsub("time:6", "time:3", time6_order)
time_order <- match( c(time0_order, time3_order, time6_order), colnames(proj2_resid) )

heatmap3(proj2_resid_winsorize[,time_order], 
         col = cls, # specify colors 
         ColSideColors = cols[time_order], # specify time color code
         Colv = NA,
         scale = "none", # no scaling by row
         labCol = colnames(proj2_resid)[time_order], # specify patient_time
         ColSideLabs = "Time", 
         labRow = "",
         xlab = "Patient_Time",
         legendfun=function() showLegend(col = proj2_heatmap_pal,
                                         legend = c("Time 0",  "Time 3", "Time 6"),
                                         cex = 1.2,
                                         lwd = 5  )
)

#--------------------------------------------------------------------------------------------
# longitudinal boxplots

# get pap_df from Write to Excel code chunks
proj2_signif_boxplot_df <- raw_data_median_proj2 %>% 
  select(PROBE_ID, contains("id:")) %>%
  filter( PROBE_ID %in% pap_df$PROBE_ID[1:6] )

proj2_signif_boxplot.func <- function(signif_mat, draw){
  signif_mat2 <- signif_mat[,-1]
  signif_df <- data.frame(
    treatment = factor( toupper( substr(colnames(signif_mat2), 4,6) ) ),
    time = factor( str_sub( colnames(signif_mat2), -1,-1 ) ), 
    fluorescence = as.numeric(signif_mat2[draw,])
  )
  ggplot(signif_df, aes(x = time, y = fluorescence, fill = treatment)) +
    geom_boxplot(width = 0.5, position=position_dodge2(width = 0.5)) +
    labs(title = paste0("Boxplots of Fluorescence Levels for \nPeptide: ", 
                        signif_mat$PROBE_ID[draw]), 
         x = "Time", y = "log2 Median Fluorescence") +
    scale_fill_manual(values=c("#F8766D", "#00BFC4")) +
    ylim(c(2,13.2)) +
    theme(panel.background = element_rect(fill = "grey90"),
          panel.grid.major = element_line(color = "white"),
          panel.grid.minor = element_line(color = "white"),
          # legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))
}

grid.arrange(
  proj2_signif_boxplot.func(proj2_signif_boxplot_df,1),
  proj2_signif_boxplot.func(proj2_signif_boxplot_df,2),
  proj2_signif_boxplot.func(proj2_signif_boxplot_df,3),
  proj2_signif_boxplot.func(proj2_signif_boxplot_df,4),
  proj2_signif_boxplot.func(proj2_signif_boxplot_df,5),
  proj2_signif_boxplot.func(proj2_signif_boxplot_df,6),
  ncol = 2
)

####################################################################################### 
#                               Gene Set Analyses After LMER                          #
#######################################################################################

proj2_get_SeqID.func <- function(coeff, BH){
  signif_seq_id <- unique( raw_data_median_proj2$SEQ_ID[coeff > .3333 & BH <= BH_FDR_cutoff_proj2] )
  seq_id_ok <- as.numeric(uniprot_gene$seq_id %in% signif_seq_id)
  names(seq_id_ok) <- uniprot_gene$entrez_gene_id_pete
  seq_id_ok2 <- seq_id_ok[!(names(seq_id_ok) %in% entrez_id_repeat)]
  seq_id_ok2 <- c(seq_id_ok2, sapply( entrez_id_repeat, function(x){max(seq_id_ok[which(names(seq_id_ok) == x)])} ) )
}

# deploy allez
PAP_SeqID <- proj2_get_SeqID.func(PAP_result[,"PAP_time_effect"], pmin(PAP_Ftest_KR_BH,PAP_Ftest_Satterthwaite_BH))
PAP_allez.go <- allez(PAP_SeqID, lib = "org.Hs.eg", idtype = "ENTREZID", sets = "GO")

# get allez result
nom.alpha <- 0.05
max_gene_in_set <- 300

allezTable(PAP_allez.go, symbol = T, nominal.alpha = nom.alpha, n.upp = max_gene_in_set, in.set = T)[,c(1:5,7)]
get_alleztable.func(PAP_allez.go)
allezPlot(PAP_allez.go, nominal.alpha = nom.alpha, n.upp = max_gene_in_set)


####################################################################################### 
#                               Boxplot of interesting peptides                       #
#######################################################################################

# first make sure stage aligns with raw_data_median
sum(as.numeric(colnames(raw_data_median %>% select(-(PROBE_DESIGN_ID:DESIGN_ID))) == patient_key$id)) == 
  nrow(patient_key)

contrast_df.func <- function(contrast_BH, contrast_diff){
  df <- data.frame(
    PROBE_ID = raw_data_median$PROBE_ID[all_kw$anova_BH <= BH_FDR_cutoff][contrast_BH <= BH_FDR_cutoff],
    SEQ_ID = raw_data_median$SEQ_ID[all_kw$anova_BH <= BH_FDR_cutoff][contrast_BH <= BH_FDR_cutoff],
    Effect_size = contrast_diff[contrast_BH <= BH_FDR_cutoff],
    KW_BH_FDR = all_kw$anova_BH[all_kw$anova_BH <= BH_FDR_cutoff][contrast_BH <= BH_FDR_cutoff]
  )
  df <- df %>% 
    filter(abs(Effect_size)>1) %>%
    arrange(desc(Effect_size)) %>%
    remove_rownames() %>%
    mutate(Effect_size = round(Effect_size, 4),
           KW_BH_FDR = round(KW_BH_FDR, 4))
  return(df)
}

# based on wilcox test
mCRPC_others_df <- contrast_df.func(mCRPC_others_wilcox_BH, mCRPC_median - NOT_mCRPC_median)
cancer_normal_df <- contrast_df.func(cancer_normal_wilcox_BH, cancer_median - normal_median)
mCRPC_nmCRPC_df <- contrast_df.func(mCRPC_nmCRPC_wilcox_BH, mCRPC_median - nmCRPC_median) 
nmCRPC_nmCSPC_df <- contrast_df.func(nmCRPC_nmCSPC_wilcox_BH, nmCRPC_median - nmCSPC_median) 
nmCSPC_newdx_df <- contrast_df.func(nmCSPC_newdx_wilcox_BH, nmCSPC_median - newdx_median) 
newdx_normal_df <- contrast_df.func(newdx_normal_wilcox_BH, newdx_median - normal_median)

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

png("09_1a_Cancer_higher_than_Normal.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=1)
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=2)
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=3)
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=4)
dev.off()

png("09_1b_Normal_higher_than_Cancer.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=nrow(cancer_normal_df))
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=nrow(cancer_normal_df)-1)
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=nrow(cancer_normal_df)-2)
boxplot_func(ref_df=cancer_normal_df,grp1=c("new_dx","nmCSPC","nmCRPC","mCRPC"),grp2="normal",rename_grp1="cancer",draw=nrow(cancer_normal_df)-3)
dev.off()

png("09_2a_mCRPC_higher_than_others.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=1)
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=2)
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=3)
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=4)
dev.off()

png("09_2b_others_higher_than_mCRPC.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=nrow(mCRPC_others_df))
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=nrow(mCRPC_others_df)-1)
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=nrow(mCRPC_others_df)-2)
boxplot_func(ref_df=mCRPC_others_df,grp1="mCRPC",grp2=c("new_dx","nmCSPC","nmCRPC","normal"),rename_grp2="others",draw=nrow(mCRPC_others_df)-3)
dev.off()

png("09_3a_Newdx_higher_than_Normal.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = 1)
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = 2)
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = 3)
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = 4)
dev.off()

png("09_3b_Normal_higher_than_Newdx.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = nrow(newdx_normal_df))
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = nrow(newdx_normal_df)-1)
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = nrow(newdx_normal_df)-2)
boxplot_func(ref_df=newdx_normal_df,grp1="new_dx",grp2="normal",draw = nrow(newdx_normal_df)-3)
dev.off()

png("09_4a_nmCSPC_higher_than_Newdx.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = 1)
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = 2)
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = 3)
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = 4)
dev.off()

png("09_4b_Newdx_higher_than_nmCSPC.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = nrow(nmCSPC_newdx_df))
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = nrow(nmCSPC_newdx_df)-1)
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = nrow(nmCSPC_newdx_df)-2)
boxplot_func(ref_df=nmCSPC_newdx_df,grp1="nmCSPC",grp2="new_dx",draw = nrow(nmCSPC_newdx_df)-3)
dev.off()

png("09_5a_nmCRPC_higher_than_nmCSPC.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=1)
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=2)
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=3)
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=4)
dev.off()

png("09_5b_nmCSPC_higher_than_nmCRPC.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=nrow(nmCRPC_nmCSPC_df))
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=nrow(nmCRPC_nmCSPC_df)-1)
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=nrow(nmCRPC_nmCSPC_df)-2)
boxplot_func(ref_df=nmCRPC_nmCSPC_df,grp1="nmCRPC",grp2="nmCSPC",draw=nrow(nmCRPC_nmCSPC_df)-3)
dev.off()

png("09_6a_mCRPC_higher_than_nmCRPC.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=1)
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=2)
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=3)
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=4)
dev.off()

png("09_6b_nmCRPC_higher_than_mCRPC.png", height = 1024, width = 512)
par(mfrow=c(4,1))
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=nrow(mCRPC_nmCRPC_df))
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=nrow(mCRPC_nmCRPC_df)-1)
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=nrow(mCRPC_nmCRPC_df)-2)
boxplot_func(ref_df=mCRPC_nmCRPC_df,grp1="mCRPC",grp2="nmCRPC",draw=nrow(mCRPC_nmCRPC_df)-3)
dev.off()


####################################################################################### 
#                             Write to Excel and Save Results                         #
#######################################################################################

anova_df <- data.frame(
  PROBE_ID = raw_data_median$PROBE_ID[all_kw$anova_BH <= BH_FDR_cutoff],
  SEQ_ID = raw_data_median$SEQ_ID[all_kw$anova_BH <= BH_FDR_cutoff],
  KW_BH_FDR = all_kw$anova_BH[all_kw$anova_BH <= BH_FDR_cutoff]
)
anova_df <- anova_df %>% 
  remove_rownames() %>%
  mutate( KW_BH_FDR = round(KW_BH_FDR, 4) )

contrast_df.func <- function(contrast_BH, contrast_diff){
  df <- data.frame(
    PROBE_ID = raw_data_median$PROBE_ID[all_kw$anova_BH <= BH_FDR_cutoff][contrast_BH <= BH_FDR_cutoff],
    SEQ_ID = raw_data_median$SEQ_ID[all_kw$anova_BH <= BH_FDR_cutoff][contrast_BH <= BH_FDR_cutoff],
    Effect_size = contrast_diff[contrast_BH <= BH_FDR_cutoff],
    KW_BH_FDR = all_kw$anova_BH[all_kw$anova_BH <= BH_FDR_cutoff][contrast_BH <= BH_FDR_cutoff],
    contrast_BH_FDR = contrast_BH[contrast_BH <= BH_FDR_cutoff]
  )
  df <- df %>% 
    filter(abs(Effect_size)>1) %>%
    arrange(desc(Effect_size)) %>%
    remove_rownames() %>%
    mutate(Effect_size = round(Effect_size, 4),
           KW_BH_FDR = round(KW_BH_FDR, 4),
           contrast_BH_FDR = round(contrast_BH_FDR, 4)
    )
  return(df)
}

# based on wilcox test
cancer_normal_df <- contrast_df.func(cancer_normal_wilcox_BH, cancer_median - normal_median)
mCRPC_others_df <- contrast_df.func(mCRPC_others_wilcox_BH, mCRPC_median - NOT_mCRPC_median)
mCRPC_nmCRPC_df <- contrast_df.func(mCRPC_nmCRPC_wilcox_BH, mCRPC_median - nmCRPC_median) 
nmCRPC_nmCSPC_df <- contrast_df.func(nmCRPC_nmCSPC_wilcox_BH, nmCRPC_median - nmCSPC_median) 
nmCSPC_newdx_df <- contrast_df.func(nmCSPC_newdx_wilcox_BH, nmCSPC_median - newdx_median) 
newdx_normal_df <- contrast_df.func(newdx_normal_wilcox_BH, newdx_median - normal_median)

summary_page_df <- data.frame(
  contrasts = c("cancer vs normal", "mCRPC vs others", "mCRPC vs nmCRPC",
                "nmCRPC vs nmCSPC", "nmCSPC vs new_dx", "new_dx vs normal"),
  total_peptide_counts = c( nrow(cancer_normal_df), nrow(mCRPC_others_df), nrow(mCRPC_nmCRPC_df),
                            nrow(nmCRPC_nmCSPC_df), nrow(nmCSPC_newdx_df), nrow(newdx_normal_df) ),
  peptides_with_positive_effect_size = c(cancer_normal_df %>% filter(Effect_size > 0) %>% nrow(), 
                                         mCRPC_others_df %>% filter(Effect_size > 0) %>% nrow(),
                                         mCRPC_nmCRPC_df %>% filter(Effect_size > 0) %>% nrow(),
                                         nmCRPC_nmCSPC_df %>% filter(Effect_size > 0) %>% nrow(),
                                         nmCSPC_newdx_df %>% filter(Effect_size > 0) %>% nrow(),
                                         newdx_normal_df %>% filter(Effect_size > 0) %>% nrow())
)
summary_page_df <- summary_page_df %>% 
  mutate(peptides_with_negative_effect_size = total_peptide_counts - peptides_with_positive_effect_size)

# longitudinal
pap_df <- data.frame(
  PROBE_ID = raw_data_median_proj2$PROBE_ID[signif_crit_proj2],
  SEQ_ID = raw_data_median_proj2$SEQ_ID[signif_crit_proj2],
  Time_Effect = PAP_result[,"PAP_time_effect"][signif_crit_proj2],
  KR_BH = PAP_Ftest_KR_BH[signif_crit_proj2],
  Satterth_BH = PAP_Ftest_Satterthwaite_BH[signif_crit_proj2]
) %>% 
  mutate(BH_FDR = pmin(KR_BH, Satterth_BH),
         BH_FDR = round(BH_FDR, 8), 
         Time_Effect = round(Time_Effect, 4)) %>%
  select(PROBE_ID, SEQ_ID, BH_FDR, Time_Effect) %>%
  arrange(desc(Time_Effect))


wb = createWorkbook()
sheet = createSheet(wb, "Contrast Summaries")
addDataFrame(summary_page_df, sheet=sheet, row.names=FALSE, startRow = 2, startColumn = 2)

sheet = createSheet(wb, "Kruskal-Wallis")
addDataFrame(anova_df, sheet=sheet, row.names=FALSE)

sheet = createSheet(wb, "cancer vs normal")
addDataFrame(cancer_normal_df, sheet=sheet, row.names=FALSE)

sheet = createSheet(wb, "mCRPC vs others")
addDataFrame(mCRPC_others_df, sheet=sheet, row.names=FALSE)

sheet = createSheet(wb, "mCRPC vs nmCRPC")
addDataFrame(mCRPC_nmCRPC_df, sheet=sheet, row.names=FALSE)

sheet = createSheet(wb, "nmCRPC vs nmCSPC")
addDataFrame(nmCRPC_nmCSPC_df, sheet=sheet, row.names=FALSE)

sheet = createSheet(wb, "nmCSPC vs newdx")
addDataFrame(nmCSPC_newdx_df, sheet=sheet, row.names=FALSE)

sheet = createSheet(wb, "newdx vs normal")
addDataFrame(newdx_normal_df, sheet=sheet, row.names=FALSE)

sheet = createSheet(wb, "PAP_Longitudinal")
addDataFrame(pap_df, sheet=sheet, row.names=FALSE)

saveWorkbook(wb, "09_Significant_Peptides.xlsx")


#----------------------------------------------------------------------------------------------

# save results
# save(logreg_pval,
#      lmer_result,
#      all_anova_pval,
#      all_kw_pval,
#      all_anova_mse,
#      file = "09_Cancer_Stage_Effects.RData")

# save(cancer_normal_allez.go,
#      mCRPC_others_allez.go,
#      mCRPC_nmCRPC_allez.go,
#      nmCRPC_nmCSPC_allez.go,
#      nmCSPC_newdx_allez.go,
#      newdx_normal_allez.go,
#      PAP_allez.go,
#      file = "09_allez_results.RData"
# )

