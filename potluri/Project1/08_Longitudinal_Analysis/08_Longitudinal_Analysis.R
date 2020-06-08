library(allez)
# library(pbkrtest) # Kenward-Roger approx
library(lme4) # linear mixed effects model
library(lmerTest) # Satterthwaite approx
library(fdrtool)
library(locfdr)
library(Rtsne)
library(heatmap3)
library(ggplot2)
library(gridExtra)
library(tidyverse) # make sure you have the latest tidyverse ! and install the package "janitor"
library(xlsx)

####################################################################################### 
#                           Some Common Variables and Functions                       #
####################################################################################### 

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
  countab <- rbind(c("locFDR threshold", thresh.vec),
                   c("Peptide counts", counter))
  return(countab)
}


# if you dont have the following data files, download it from shared Box folder:
# https://uwmadison.app.box.com/folder/110549175126 

raw_data_median_proj2 <- read_csv("raw_data_median_proj2.csv")
load("08_LMER_results.RData")
load("08_RandomEffects_LRstat.RData")

####################################################################################### 
#                                      Data Processing                                #
####################################################################################### 

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

raw_data = read_csv("raw_data_complete.csv")

raw_data = raw_data %>%
  select(PROBE_DESIGN_ID:Y, any_of(array_id_key_proj2$array_id)) %>% # drop patients with records at different stages 
  select( -c(X, Y, MATCH_INDEX, DESIGN_NOTE, SELECTION_CRITERIA, MISMATCH, PROBE_CLASS) ) # drop unhelpful columns

# check
colnames(raw_data)[1:15]

# take log2 transformation
raw_data <- raw_data %>%
  mutate_at(vars(matches("dat")), log2)

# check any NA's
raw_data %>% 
  select_if(function(x) any(is.na(x))) 

# make sure array_id_key_proj2 align with raw_data_complete
array_iii <- match( colnames( raw_data %>% select(any_of(array_id_key_proj2$array_id)) ), array_id_key_proj2$array_id )
array_id_key_proj2 <- array_id_key_proj2[array_iii , ]
sum(as.numeric( colnames( raw_data %>% select(any_of(array_id_key_proj2$array_id)) ) 
                == array_id_key_proj2$array_id )) == nrow(array_id_key_proj2)


# compute median
raw_data_median_proj2 <- t( apply( select(raw_data, contains("dat")), 1, function(x) {
  as.vector( tapply( as.numeric(x), list(array_id_key_proj2$id, array_id_key_proj2$time), FUN=median) )
} ) )
raw_data_median_proj2 <- bind_cols( select(raw_data, -contains("dat")), as.data.frame(raw_data_median_proj2) ) 

colnames(raw_data_median_proj2)[1:15] # check

# want to rename column
example_row <- tapply( as.numeric(select(raw_data, contains("dat"))[1,]), 
                       list(array_id_key_proj2$id, array_id_key_proj2$time), FUN=median) %>% 
  as.data.frame() %>% 
  rownames_to_column("id") %>% 
  pivot_longer(-id, names_to = "time", values_to = "fluorescence") %>% 
  arrange(time, id)
sum(as.numeric( example_row$fluorescence == select(raw_data_median_proj2, contains("V"))[1,] )) == nrow(example_row) 
example_row$column_name <- paste( paste0("id:", example_row$id), paste0("time:", example_row$time), sep = "_" )
raw_data_median_proj2 <- raw_data_median_proj2 %>% rename_at(vars(contains("V")), ~ example_row$column_name) 

# make sure sample_key_proj2 align with raw_data_median_proj2
sample_key_proj2 <- sample_key_proj2 %>%
  select(-n) %>%
  arrange(time, id)
sum(as.numeric(sample_key_proj2$id == example_row$id)) == nrow(sample_key_proj2) # check
sum(as.numeric(sample_key_proj2$time == example_row$time)) == nrow(sample_key_proj2) # check
sample_key_proj2 <- sample_key_proj2 %>% 
  left_join(array_id_key_proj2 %>% select(id, treatment) %>% distinct())

write.table(raw_data_median_proj2, file = "raw_data_median_proj2.csv", sep = ",", row.names = F)


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

# sort order of patients in boxplot
median_long_proj2$id_time <- factor(median_long_proj2$id_time, levels = c(
  unique(median_long_proj2$id_time[median_long_proj2$treat_time == "ADT_0"]),
  unique(median_long_proj2$id_time[median_long_proj2$treat_time == "ADT_3"]),
  unique(median_long_proj2$id_time[median_long_proj2$treat_time == "ADT_6"]),
  unique(median_long_proj2$id_time[median_long_proj2$treat_time == "Vaccine_0"]),
  unique(median_long_proj2$id_time[median_long_proj2$treat_time == "Vaccine_3"]),
  unique(median_long_proj2$id_time[median_long_proj2$treat_time == "Vaccine_6"])
))

ggplot(median_long_proj2, aes(x = id_time, y = fluorescence, fill = treat_time)) +
  geom_boxplot(outlier.shape = ".") +
  # scale_fill_discrete(name="Stage") +
  # guides(fill=guide_legend(title="Stage")) +
  scale_fill_manual(name = "treatment_time", values = pal_proj2) +
  labs(title = "Boxplots of Peptide Fluorescence Levels for Patients at 3 time points", 
       x = "Patient_Time", y = "Median Fluorescence Levels on log2 scale") +
  theme(panel.background = element_rect(fill = "grey90"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))


####################################################################################### 
#                                 Test Time Random Effects                            #
####################################################################################### 

# want to know if random time slope is needed in the random effects?

ncol_med_proj2 <- ncol(raw_data_median_proj2) 
n_med_proj2 <- nrow(sample_key_proj2)

# sample_key_proj2$treatment <- as.factor(sample_key_proj2$treatment)
# sample_key_proj2$treatment <- relevel(sample_key_proj2$treatment, ref = "ADT")
# lrstat_bigmodel <- rep(NA, nrow(raw_data_median_proj2))

lrstat_smallmodel_pap <- rep(NA, nrow(raw_data_median_proj2))
lrstat_smallmodel_adt <- rep(NA, nrow(raw_data_median_proj2))

smallmodel_testrandom.func <- function(y, treat_type){
  resp <- y[sample_key_proj2$treatment == treat_type]
  fit1 <- lmer(resp ~ time + (1 + time | id), REML = T, 
               data = sample_key_proj2[sample_key_proj2$treatment == treat_type,])
  fit0 <- lmer(resp ~ time + (1 | id), REML = T, 
               data = sample_key_proj2[sample_key_proj2$treatment == treat_type,])
  return( -2*( as.numeric(logLik(fit0)) - as.numeric(logLik(fit1)) ) )
}

for(i in 1:nrow(raw_data_median_proj2)){
  y <- as.numeric(raw_data_median_proj2[i, (ncol_med_proj2 - n_med_proj2 + 1) : ncol_med_proj2])
  # fit1_bigmodel <- lmer(y ~ treatment*time + (1 + time | id), data = sample_key_proj2, REML = T)
  # fit0_bigmodel <- lmer(y ~ treatment*time + (1 | id), data = sample_key_proj2, REML = T)
  # lrstat_bigmodel[i] <- -2*( as.numeric(logLik(fit0_bigmodel)) - as.numeric(logLik(fit1_bigmodel)) )
  lrstat_smallmodel_pap[i] <- smallmodel_testrandom.func(y = y, treat_type = "Vaccine")
  lrstat_smallmodel_adt[i] <- smallmodel_testrandom.func(y = y, treat_type = "ADT") 
  print(i)
}

save(lrstat_bigmodel, lrstat_smallmodel_pap, lrstat_smallmodel_adt, 
     file = "08_RandomEffects_LRstat.RData")

# chisq2_pval_bigmodel <- pchisq(lrstat_bigmodel, df = 2, lower.tail = F)
chisq2_pval_smallmodel_pap <- pchisq(lrstat_smallmodel_pap, df = 2, lower.tail = F)
chisq2_pval_smallmodel_adt <- pchisq(lrstat_smallmodel_adt, df = 2, lower.tail = F)

# chisq2_BH_bigmodel <- p.adjust(chisq2_pval_bigmodel, method = "BH")
chisq2_BH_smallmodel_pap <- p.adjust(chisq2_pval_smallmodel_pap, method = "BH")
chisq2_BH_smallmodel_adt <- p.adjust(chisq2_pval_smallmodel_adt, method = "BH")

# hist(chisq2_pval_bigmodel,breaks= 100)
hist(chisq2_pval_smallmodel_pap,breaks= 100)
hist(chisq2_pval_smallmodel_adt,breaks= 100)

# count.func(chisq2_BH_bigmodel, seq(.05,.1,by=.01))
count.func(chisq2_BH_smallmodel_pap, seq(.05,.1,by=.01))
count.func(chisq2_BH_smallmodel_adt, seq(.05,.1,by=.01))

## Even under very conservative tests, there are many peptides that have signif random time slope
## Hemanth's analysis also suggests changes in individual antibody response over time vary across different patients
## With only 3 time points including baseline, we refrain from fitting more complicated time random effects


####################################################################################### 
#             Linear Mixed Model to Assess Treatment Effect (REML = TRUE)             #
#                       Separate Models for PAP and ADT Groups                        #
####################################################################################### 

ncol_med_proj2 <- ncol(raw_data_median_proj2) 
n_med_proj2 <- nrow(sample_key_proj2)

# initiate
PAP_resid_fit0 <- matrix(NA, nrow = nrow(raw_data_median_proj2), 
                         ncol = nrow(sample_key_proj2%>%filter(treatment == "Vaccine")))
ADT_resid_fit0 <- matrix(NA, nrow = nrow(raw_data_median_proj2), 
                         ncol = nrow(sample_key_proj2%>%filter(treatment == "ADT")))
PAP_result <- matrix(NA, nrow = nrow(raw_data_median_proj2), ncol = 8)
ADT_result <- matrix(NA, nrow = nrow(raw_data_median_proj2), ncol = 8)

colnames(PAP_resid_fit0) <- colnames(raw_data_median_proj2 %>% select(contains("id:pap")))
colnames(ADT_resid_fit0) <- colnames(raw_data_median_proj2 %>% select(contains("id:adt")))
colnames(PAP_result) <- paste0("PAP_", c(
  "time_effect", 
  "time_tstat",
  "KR_df",
  "KR_Ftest_pval",
  "Satterthwaite_df",
  "Satterthwaite_Ftest_pval",
  "zval_1sided_KR",
  "zval_1sided_Satterthwaite"
))
colnames(ADT_result) <- paste0("ADT_", c(
  "time_effect", 
  "time_tstat",
  "KR_df",
  "KR_Ftest_pval",
  "Satterthwaite_df",
  "Satterthwaite_Ftest_pval",
  "zval_1sided_KR",
  "zval_1sided_Satterthwaite"
))

Test_Time.func <- function(y, treat_type){
  resp <- y[sample_key_proj2$treatment == treat_type]
  fit1 <- lmer(resp ~ time + (1 + time | id), REML = T, 
               data = sample_key_proj2[sample_key_proj2$treatment == treat_type,])
  fit0 <- lmer(resp ~ 1 + (1 + time | id), REML = T, 
               data = sample_key_proj2[sample_key_proj2$treatment == treat_type,])
  resid_fit0 <- unname(round(resid(fit0),4))
  effect_tstat <- coef(summary(fit1))['time',c('Estimate', 't value')] 
  KR_df_pval <- contest(fit1, c(0,1), ddf = "Kenward-Roger")[c('DenDF', 'Pr(>F)')] 
  Satterthwaite_df_pval <- contest(fit1, c(0,1))[c('DenDF', 'Pr(>F)')] 
  zval_1sided_KR <- qnorm(pt( as.numeric(effect_tstat['t value']), 
                             df = as.numeric(KR_df_pval['DenDF']) ,  lower.tail = T ))
  zval_1sided_Satterthwaite <- qnorm(pt( as.numeric(effect_tstat['t value']), 
                                        df = as.numeric(Satterthwaite_df_pval['DenDF']) ,  lower.tail = T ))
  result <- c(
    as.numeric(effect_tstat), 
    as.numeric(KR_df_pval),
    as.numeric(Satterthwaite_df_pval),
    zval_1sided_KR,
    zval_1sided_Satterthwaite
  )
  return( list(
    resid_fit0 = resid_fit0,
    result = result
  ) )
}


for(i in 1:nrow(raw_data_median_proj2)){
  y <- as.numeric(raw_data_median_proj2[i, (ncol_med_proj2 - n_med_proj2 + 1) : ncol_med_proj2])
  PAP_test <- Test_Time.func(y, "Vaccine")
  ADT_test <- Test_Time.func(y, "ADT")
  
  PAP_resid_fit0[i,] <- PAP_test$resid_fit0
  PAP_result[i,] <- PAP_test$result
  
  ADT_resid_fit0[i,] <- ADT_test$resid_fit0
  ADT_result[i,] <- ADT_test$result

  if(i %% 100 == 0){
    print(i)
  }
}

save(PAP_resid_fit0, ADT_resid_fit0, PAP_result, ADT_result,
     file = "08_LMER_results.RData")


####################################################################################### 
#                                           Local FDR                                 #
#######################################################################################

# F-test p-values based on KR adjustments more conservative than Satterthwaite

plot(PAP_result[,"PAP_Satterthwaite_Ftest_pval"], PAP_result[,"PAP_KR_Ftest_pval"], pch = ".", xlim = c(0,.2), ylim = c(0,.2))
abline(a=0, b=1, col = "red", lty=2, lwd = 2)

plot(ADT_result[,"ADT_Satterthwaite_Ftest_pval"], ADT_result[,"ADT_KR_Ftest_pval"], pch = ".", xlim = c(0,.2), ylim = c(0,.2))
abline(a=0, b=1, col = "red", lty=2, lwd = 2)

PAP_Ftest_KR_BH <- p.adjust(PAP_result[,"PAP_KR_Ftest_pval"],method="BH")
PAP_Ftest_Satterthwaite_BH <- p.adjust(PAP_result[,"PAP_Satterthwaite_Ftest_pval"],method="BH")
ADT_Ftest_KR_BH <- p.adjust(ADT_result[,"ADT_KR_Ftest_pval"],method="BH")
ADT_Ftest_Satterthwaite_BH <- p.adjust(ADT_result[,"ADT_Satterthwaite_Ftest_pval"],method="BH")

## previous analysis concurs that Likelihood ratio test p-values (fitted with REML = F) are very liberal

png("08_pvalues_density_histograms.png", width = 1024, height = 1024)
par(mfrow=c(2,2))
hist(PAP_result[,"PAP_KR_Ftest_pval"], breaks=100, freq = F, ylim = c(0,32),
     main = "PAP's (KR F-test) p-values density histogram")
text(x=.7, y=20, paste0(sum(as.numeric(PAP_Ftest_KR_BH <= .01))," peptides at 1% BH FDR"), cex = 1)

hist(PAP_result[,"PAP_Satterthwaite_Ftest_pval"], breaks=100, freq = F,ylim = c(0,32),
     main = "PAP's (Satterthwaite F-test) p-values density histogram ")
text(x=.7, y=20, paste0(sum(as.numeric(PAP_Ftest_Satterthwaite_BH <= .01))," peptides at 1% BH FDR"), cex = 1)

hist(ADT_result[,"ADT_KR_Ftest_pval"], breaks=100, freq = F,ylim = c(0,32),
     main = "ADT's (KR F-test) p-values density histogram ")
text(x=.7, y=20, "No signif peptides at 20% BH FDR", cex = 1)

hist(ADT_result[,"ADT_Satterthwaite_Ftest_pval"], breaks=100, freq = F,ylim = c(0,32),
     main = "ADT's (Satterthwaite F-test) p-values density histogram ")
text(x=.7, y=20, "No signif peptides at 20% BH FDR", cex = 1)
dev.off()

# BH threshold peptide counts
count.func(PAP_Ftest_KR_BH, seq(.01,.1,by=.01))
count.func(PAP_Ftest_Satterthwaite_BH, seq(.01,.1,by=.01))
count.func(ADT_Ftest_KR_BH, seq(.63,.7,by=.01))
count.func(ADT_Ftest_Satterthwaite_BH, seq(.63,.7,by=.01))

# volcano plot (might be useful in setting initial values to estimate f0 in locfdr)
volcano_plot.func <- function(time_effect, pval, BH, title){
  plot(x = time_effect, y = -log10(pval), pch = ".", ylim = c(0,12.5), xlim = c(-.4, .85),
       xlab = "coefficient of time fixed effect", ylab = "-log10(F-test p-values)",
       main = title)
  lines(x = time_effect[ BH <= .01 & time_effect >= .3333 ], 
        y = -log10(pval[ BH <= .01 & time_effect >= .3333]),
        type = "p", pch = ".", col = "blue")
}

png("08_volcano_plots.png", width = 1024, height = 1024)
par(mfrow=c(2,2))
volcano_plot.func(PAP_result[,"PAP_time_effect"], PAP_result[,"PAP_KR_Ftest_pval"], PAP_Ftest_KR_BH,
                  "PAP's volcano plot (KR F-test)")
volcano_plot.func(PAP_result[,"PAP_time_effect"], PAP_result[,"PAP_Satterthwaite_Ftest_pval"], PAP_Ftest_Satterthwaite_BH,
                  "PAP's volcano plot (Satterthwaite F-test)")
volcano_plot.func(ADT_result[,"ADT_time_effect"], ADT_result[,"ADT_KR_Ftest_pval"], ADT_Ftest_KR_BH,
                  "ADT's volcano plot (KR F-test)")
volcano_plot.func(ADT_result[,"ADT_time_effect"], ADT_result[,"ADT_Satterthwaite_Ftest_pval"], ADT_Ftest_Satterthwaite_BH,
                  "ADT's volcano plot (Satterthwaite F-test)")
dev.off()

## Note that no time effect of ADT exceeds .3333, ie. one fold-change after 3 months

#--------------------------------------------------------------------------------------
# check PAP's locFDR based on z-scores of LMER tstat with KR-adjusted df 

PAP_locfdr_KR <- locfdr(PAP_result[,"PAP_zval_1sided_KR"], 
                        df = 17, # to fit estimated f(z) ,
                        mlests = c(-2, 1.8),
                        main = "PAP's locFDR based on LMER t-stat 1-sided pval with KR-adjusted df")

## green solid curve: natural spline estimate of mixture density f(z)
## blue dashed curve: ML estimated null subdensity p0f0
## colored hist bars: estimated non-null counts

PAP_locfdr_KR$fp0['mlest', 'p0'] # ML estimate of proportion of null peptides, p0 

PAP_locfdr_KR$fp0['mlest', c('delta', 'sigma')] # ML estimate of null mean and sd
## mean z-scores of null peptides so negative -- lower antibody response over time on null peptides?

PAP_locfdr_KR$Efdr["Efdr"] # expected locFDR for non-null peptides : 

locfdr(PAP_result[,"PAP_zval_1sided_KR"],df=17,plot=3,main="PAP's locFDR based on LMER t-stat 1-sided pval with KR-adjusted df")$call
## proportion of non-null peptides less than a given level of locFDR

# let's get what we want
locFDR_PAP_KR <- PAP_locfdr_KR$fdr
count.func(locFDR_PAP_KR, seq(.01,.1,by=.01))

length( which(locFDR_PAP_KR <= .01 & PAP_result[,"PAP_time_effect"] >= .3333) )

#--------------------------------------------------------------------------------------
# check PAP's locFDR based on z-scores of LMER tstat with Satterthwaite-adjusted df 

PAP_locfdr_Satterthwaite <- locfdr(PAP_result[,"PAP_zval_1sided_Satterthwaite"], 
                                   df = 14,  # to fit estimated f(z) 
                                   # mlests = c(-1.627, 1.378), # initial value for mean & sd of estimating f(0)
                                   mlests = c(-2, 1.4), # initial value for mean & sd of estimating f(0)
                                   main = "PAP's locFDR based on LMER t-stat 1-sided pval with Satterthwaite-adjusted df")

PAP_locfdr_Satterthwaite$fp0['mlest', 'p0'] # ML estimate of proportion of null peptides, p0 

PAP_locfdr_Satterthwaite$fp0['mlest', c('delta', 'sigma')] # ML estimate of null mean and sd
## mean z-scores of null peptides so negative -- lower antibody response over time on null peptides?

PAP_locfdr_Satterthwaite$Efdr["Efdr"] # expected locFDR for non-null peptides 

# let's get what we want
locFDR_PAP_Satterthwaite <- PAP_locfdr_Satterthwaite$fdr
count.func(locFDR_PAP_Satterthwaite, seq(.01,.1,by=.01))

length( which(locFDR_PAP_Satterthwaite <= .01 & PAP_result[,"PAP_time_effect"] >= .3333) )

#--------------------------------------------------------------------------------------
# check ADT's locFDR based on z-scores of LMER tstat with KR-adjusted df 

ADT_locfdr_KR <- locfdr(ADT_result[,"ADT_zval_1sided_KR"], 
                        df = 10, # to fit estimated f(z) 
                        mlests = c(-1.5,1), # initial value for mean & sd of estimating f(0)
                        # mlests = c(.05,2.3), # initial value for mean & sd of estimating f(0)
                        main = "ADT's locFDR based on LMER t-stat 1-sided pval with KR-adjusted df")

# let's get what we want
locFDR_ADT_KR <- ADT_locfdr_KR$fdr
count.func(locFDR_ADT_KR, seq(.01,.1,by=.01))

#--------------------------------------------------------------------------------------
# check ADT's locFDR based on z-scores of LMER tstat with Satterthwaite-adjusted df 

ADT_locfdr_Satterthwaite <- locfdr(ADT_result[,"ADT_zval_1sided_Satterthwaite"], 
                                   df = 10, # to fit estimated f(z) 
                                   mlests = c(-1.5,1), # initial value for mean & sd of estimating f(0)
                                   # mlests = c(.05,2.3), # initial value for mean & sd of estimating f(0)
                                   main = "ADT's locFDR based on LMER t-stat 1-sided pval with Satterthwaite-adjusted df")

# let's get what we want
locFDR_ADT_Satterthwaite <- ADT_locfdr_Satterthwaite$fdr
count.func(locFDR_ADT_Satterthwaite, seq(.01,.1,by=.01))

#--------------------------------------------------------------------------------------
# let's see if we could plot all locfdr graphs together

png("08_locFDR_plots.png", width = 1024, height = 1024)
par(mfrow=c(2,2))
locfdr(PAP_result[,"PAP_zval_1sided_KR"], 
       df = 17, # to fit estimated f(z) ,
       mlests = c(-2, 1.8),
       main = "PAP's locFDR based on LMER t-stat 1-sided pval \nwith KR-adjusted df")
text(x = 3, y = 3000, paste0(sum(as.numeric(locFDR_PAP_KR <= .01)), " peptides at 1% locFDR"), cex = 1)

locfdr(PAP_result[,"PAP_zval_1sided_Satterthwaite"], 
       df = 14,  # to fit estimated f(z) 
       # mlests = c(-1.627, 1.378), # initial value for mean & sd of estimating f(0)
       mlests = c(-2, 1.4), # initial value for mean & sd of estimating f(0)
       main = "PAP's locFDR based on LMER t-stat 1-sided pval \nwith Satterthwaite-adjusted df")
text(x = 2.5, y = 3000, paste0(sum(as.numeric(locFDR_PAP_Satterthwaite <= .01)), " peptides at 1% locFDR"), cex = 1)

locfdr(ADT_result[,"ADT_zval_1sided_KR"], 
       df = 10, # to fit estimated f(z) 
       mlests = c(-1.5,1), # initial value for mean & sd of estimating f(0)
       # mlests = c(.05,2.3), # initial value for mean & sd of estimating f(0)
       main = "ADT's locFDR based on LMER t-stat 1-sided pval \nwith KR-adjusted df")
text(x = 3, y = 2500, paste0(sum(as.numeric(locFDR_ADT_KR <= .01)), " peptides at 1% locFDR"), cex = 1)

locfdr(ADT_result[,"ADT_zval_1sided_Satterthwaite"], 
       df = 10, # to fit estimated f(z) 
       mlests = c(-1.5,1), # initial value for mean & sd of estimating f(0)
       # mlests = c(.05,2.3), # initial value for mean & sd of estimating f(0)
       main = "ADT's locFDR based on LMER t-stat 1-sided pval \nwith Satterthwaite-adjusted df")
text(x = 3, y = 3000, paste0(sum(as.numeric(locFDR_ADT_Satterthwaite <= .01)), " peptides at 1% locFDR"), cex = 1)
dev.off()