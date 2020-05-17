library(allez)
library(pbkrtest) # Kenward-Roger approx
library(lme4) # linear mixed effects model
library(lmerTest) # Satterthwaite approx
library(fdrtool)
library(Rtsne)
library(heatmap3)
library(ggplot2)
library(gridExtra)
library(tidyverse) # make sure you have the latest tidyverse !


####################################################################################### 
#                           Some Common Variables and Functions                       #
####################################################################################### 

# specified color and shape scheme 
pal_proj2 <- c("turquoise1", "cornflowerblue","navy", "orchid1", "darkorange1", "firebrick1")
names(pal_proj2) <- c("ADT_0", "ADT_3", "ADT_6", "PAP_0", "PAP_3", "PAP_6")
shp_proj2 <- c(8, 15, 16, 3, 17, 18)
shp_proj2 <- c(8, 16, 18, 8, 16, 18)
names(shp_proj2) <- names(pal_proj2)

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

# typical step in FDR adjustment for LMM pval
LMM_func_proj2 <- function(LMM_pval){
  # control FDR
  LMM_BH <- p.adjust(LMM_pval, method = "BH")
  LMM_qval <- fdrtool(LMM_pval, statistic = "pvalue", verbose = F, plot  = F)$qval
  LMM_qval_eta0 <- unname(fdrtool(LMM_pval, statistic = "pvalue", verbose = F, plot  = F)$param[,"eta0"])
  
  # plot histogram of p-values
  all_peptide_hist <- hist(LMM_pval, breaks = 70, freq = F, xlab = "Linear Mixed Models (LMM) p-values", 
                           main = paste0("p-values distribution for ", length(LMM_pval), " peptides"))
  polygon_ind <- which(all_peptide_hist$density >= LMM_qval_eta0)
  for (i in polygon_ind){
    polygon( x = c(all_peptide_hist$breaks[i], all_peptide_hist$breaks[i+1], all_peptide_hist$breaks[i+1], all_peptide_hist$breaks[i]),
             y = c(all_peptide_hist$density[i], all_peptide_hist$density[i], LMM_qval_eta0, LMM_qval_eta0),
             col = "red")
  }
  text(x=0.65,y=4, labels = paste( "estimated proportion of \nnon-null peptides =",
                                   round( 100*(1 - LMM_qval_eta0),2 ),"%" ))
  
  return(list(LMM_BH = LMM_BH, LMM_qval = LMM_qval, LMM_qval_eta0 = LMM_qval_eta0))
}


load("KR_lmer_proj2_REML.RData")
load("LMM_proj2_NoREML.RData")
raw_data_median_proj2 <- read_csv("raw_data_median_proj2.csv")

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
sum(as.numeric( example_row$fluorescence == select(raw_data_median_proj2, contains("V"))[1,] )) == nrow(example_row) # check
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
#             Linear Mixed Model to Assess Treatment Effect (REML = TRUE)             #
####################################################################################### 

ncol_med_proj2 <- ncol(raw_data_median_proj2) 
n_med_proj2 <- nrow(sample_key_proj2)

sample_key_proj2$treatment <- as.factor(sample_key_proj2$treatment)
sample_key_proj2$treatment <- relevel(sample_key_proj2$treatment, ref = "ADT")
  
# initiate
pval_proj2 <- rep(NA, nrow(raw_data_median_proj2))
check_singular <- rep(NA, nrow(raw_data_median_proj2))
pval2_proj2 <- rep(NA, nrow(raw_data_median_proj2))
check_singular2 <- rep(NA, nrow(raw_data_median_proj2))
lrstat <- rep(NA, nrow(raw_data_median_proj2))

for(i in 1:nrow(raw_data_median_proj2)){
  y <- as.numeric(raw_data_median_proj2[i, (ncol_med_proj2 - n_med_proj2 + 1) : ncol_med_proj2])
  
  fit0 <- lmer(y ~ time + (1 + time | id), data = sample_key_proj2, REML = T)
  fit1 <- lmer(y ~ treatment*time + (1 + time | id), data = sample_key_proj2, REML = T)
  lrstat[i] <- -2*( logLik(fit0) - logLik(fit1) )
  
  # next time also try these
  # lmer(y ~ 1 + (1|id))
  # lmer(y ~ treatment + (1|id))
  
  fit2 <- lmer(y ~ time + (1  | id), data = sample_key_proj2, REML = T)
  fit3 <- lmer(y ~ treatment*time + (1  | id), data = sample_key_proj2, REML = T)
  pval_proj2[i] <- unlist(KRmodcomp(fit1, fit0))$stats.p.value
  pval2_proj2[i] <- unlist(KRmodcomp(fit3, fit2))$stats.p.value
  check_singular[i] <- ( isSingular(fit0) | isSingular(fit1) )
  check_singular2[i] <- ( isSingular(fit2) | isSingular(fit3) )
  print(i)
}


KR_lmer_proj2 <- data.frame(
  PROBE_ID = raw_data_median_proj2$PROBE_ID,
  pval_proj2 = pval_proj2,
  check_singular = check_singular,
  pval2_proj2 = pval2_proj2,
  check_singular2 = check_singular2
)

save(KR_lmer_proj2, file = "KR_lmer_proj2_REML.RData")


# add satterthwaite approx

# initiate 
satterthwaite_pval <- rep(NA, nrow(raw_data_median_proj2))

for(i in 1:nrow(raw_data_median_proj2)){
  y <- as.numeric(raw_data_median_proj2[i, (ncol_med_proj2 - n_med_proj2 + 1) : ncol_med_proj2])
  fit1 <- lmer(y ~ treatment*time + (1 + time | id), data = sample_key_proj2, REML = T)
  satterthwaite_pval[i] <- contest(fit1, L = diag(4)[c(2,4),] , ddf = "Satterthwaite")$'Pr(>F)'
  # print(i)
}

KR_lmer_proj2 <- cbind( KR_lmer_proj2, data.frame(satterthwaite_pval) )
save(KR_lmer_proj2, file = "KR_lmer_proj2_REML.RData")


# FDR adjustment of p-values
KR_LMM <- LMM_func_proj2(KR_lmer_proj2$pval_proj2)
KR_LMM2 <- LMM_func_proj2(KR_lmer_proj2$pval2_proj2)
Satterthwaite_LMM <- LMM_func_proj2(KR_lmer_proj2$satterthwaite_pval)

# peptide counts at various FDR thresholds
count.func(KR_LMM$LMM_BH, seq(0.01, 0.1, by = 0.01))
count.func(KR_LMM$LMM_qval, seq(0.01, 0.1, by = 0.01))
count.func(KR_LMM2$LMM_BH, seq(0.01, 0.1, by = 0.01))
count.func(KR_LMM2$LMM_qval, seq(0.01, 0.1, by = 0.01))
count.func(Satterthwaite_LMM$LMM_BH, seq(0.01, 0.1, by = 0.01))
count.func(Satterthwaite_LMM$LMM_qval, seq(0.01, 0.1, by = 0.01))

# compare KR p-val from 2 models
KR_rand_df <- data.frame(
  model = rep(c("simpler", "more"), each = 177604),
  p_values = c(KR_lmer_proj2$pval2_proj2, KR_lmer_proj2$pval_proj2)
)
ggplot(KR_rand_df, aes(x=p_values, fill=model, color=model)) +
  geom_histogram( position="identity", alpha=0.2, bins = 70) +
  labs(title = "density histogram of p-values for peptides \nwith simpler or more random effects") +
  theme(plot.title = element_text(hjust = 0.5))

# how many singular solutions
sum(as.numeric(KR_lmer_proj2$check_singular))
sum(as.numeric(KR_lmer_proj2$check_singular2))

# careful check if distribution of p-val different for singular solutions
pval_singular.func <- function(pval_vec, singular_vec){
  pval_singular_df <- data.frame(
    singularity = as.factor( singular_vec ),
    p_values = pval_vec
  )
  pval_singular_p1 <- ggplot(pval_singular_df, aes(x=p_values, fill=singularity, color=singularity)) +
    geom_histogram(aes(y=..density..), position="identity", alpha=0.2, bins = 70) +
    labs(title = "density histogram of p-values for peptides \nwith singular or non-singular LMM solutions") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
  pval_singular_p2 <- ggplot(pval_singular_df, aes(x=p_values, fill=singularity, color=singularity)) +
    geom_histogram(position="identity", alpha=0.2, bins = 70)  +
    labs(title = "frequency histogram of p-values for peptides \nwith singular or non-singular LMM solutions") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
  grid.arrange(pval_singular_p1, pval_singular_p2, ncol = 2)
  
}
pval_singular.func(KR_lmer_proj2$pval_proj2, KR_lmer_proj2$check_singular)
pval_singular.func(KR_lmer_proj2$pval2_proj2, KR_lmer_proj2$check_singular2)


# compare Satterthwaite approx p-val vs KR approx F-test pval
plot_pval_thresh <- 1
plot(x = KR_lmer_proj2$satterthwaite_pval[KR_lmer_proj2$satterthwaite_pval <= plot_pval_thresh & 
                                            KR_lmer_proj2$pval_proj2 <= plot_pval_thresh], 
     y = KR_lmer_proj2$pval_proj2[KR_lmer_proj2$satterthwaite_pval <= plot_pval_thresh & 
                                    KR_lmer_proj2$pval_proj2 <= plot_pval_thresh], 
     pch = ".", xlab = "Satterthwaite approx pval (REML=TRUE)", ylab = "KR F-test pval (REML=TRUE)")
abline(a= 0, b = 1, col = "blue", lty = 2)

# check if KR approx F-test signif peptides nested within Satterthwaite signif peptides
sum(as.numeric( KR_lmer_proj2$PROBE_ID[KR_LMM$LMM_BH <= 0.05] %in% 
                  KR_lmer_proj2$PROBE_ID[Satterthwaite_LMM$LMM_BH <= 0.05] )) == 
  length( KR_lmer_proj2$PROBE_ID[KR_LMM$LMM_BH <= 0.05] )


# to get residuals with REML = TRUE fit

resid_fit0 <- matrix(NA, nrow = nrow(raw_data_median_proj2), ncol = nrow(sample_key_proj2))

for(i in 1:nrow(raw_data_median_proj2)){
  y <- as.numeric(raw_data_median_proj2[i, (ncol_med_proj2 - n_med_proj2 + 1) : ncol_med_proj2])
  fit0 <- lmer(y ~ time + (1 + time | id), data = sample_key_proj2, REML = T)
  resid_fit0[i,] <- unname(round(resid(fit0),4))
  # print(i)
}

colnames(resid_fit0) <- colnames(raw_data_median_proj2 %>% select(contains("id:")))
KR_lmer_proj2 <- cbind( KR_lmer_proj2, as.data.frame(resid_fit0) )
save(KR_lmer_proj2, file = "KR_lmer_proj2_REML.RData")

####################################################################################### 
#             Linear Mixed Model to Assess Treatment Effect (REML = FALSE)            #
#######################################################################################

ncol_med_proj2 <- ncol(raw_data_median_proj2) 
n_med_proj2 <- nrow(sample_key_proj2)

sample_key_proj2$treatment <- as.factor(sample_key_proj2$treatment)
sample_key_proj2$treatment <- relevel(sample_key_proj2$treatment, ref = "ADT")

# initiate
loglik_fit0 <- rep(NA, nrow(raw_data_median_proj2))
loglik_fit1 <- rep(NA, nrow(raw_data_median_proj2))
treat_fit1 <- matrix(NA, nrow = nrow(raw_data_median_proj2), ncol = 2)
LMM_KR_pval <- rep(NA, nrow(raw_data_median_proj2))
LMM_singular <- rep(NA, nrow(raw_data_median_proj2))
resid_fit0 <- matrix(NA, nrow = nrow(raw_data_median_proj2), ncol = nrow(sample_key_proj2))


for(i in 1:nrow(raw_data_median_proj2)){
  y <- as.numeric(raw_data_median_proj2[i, (ncol_med_proj2 - n_med_proj2 + 1) : ncol_med_proj2])
  fit0 <- lmer(y ~ time + (1 + time | id), data = sample_key_proj2, REML = F)
  fit1 <- lmer(y ~ treatment*time + (1 + time | id), data = sample_key_proj2, REML = F)
  loglik_fit0[i] <- as.numeric( logLik(fit0) )
  loglik_fit1[i] <- as.numeric( logLik(fit1) )
  treat_fit1[i,] <- unname(round( coef(summary(fit1))[,'Estimate'][c('treatmentVaccine', 'treatmentVaccine:time')], 4 ))
  resid_fit0[i,] <- unname(round(resid(fit0),4))
  LMM_KR_pval[i] <- unlist(KRmodcomp(fit1, fit0))$stats.p.value
  LMM_singular[i] <- ( isSingular(fit0) | isSingular(fit1) )
  # print(i)
}

colnames(treat_fit1) <- c("Treat_Main", "Treat_Time") 
colnames(resid_fit0) <- colnames(raw_data_median_proj2 %>% select(contains("id:")))

LMM_proj2_NoREML <- data.frame(
  PROBE_ID = raw_data_median_proj2$PROBE_ID,
  loglik_fit0 = loglik_fit0,
  loglik_fit1 = loglik_fit1,
  LMM_KR_pval = LMM_KR_pval,
  LMM_singular = LMM_singular
)

LMM_proj2_NoREML <- cbind(
  LMM_proj2_NoREML,
  as.data.frame(treat_fit1),
  as.data.frame(resid_fit0)
)

save(LMM_proj2_NoREML, file = "LMM_proj2_NoREML.RData")


# FDR adjustment of p-values
LR_NoREML_pval <- pchisq( -2 * ( LMM_proj2_NoREML$loglik_fit0 - LMM_proj2_NoREML$loglik_fit1 ), df = 2, lower.tail = F ) 
LR_NoREML <- LMM_func_proj2( LR_NoREML_pval )

# peptide counts at various FDR thresholds
count.func(LR_NoREML$LMM_BH, seq(0.01, 0.1, by = 0.01))
count.func(LR_NoREML$LMM_qval, seq(0.01, 0.1, by = 0.01))

# compare LRT chisq p-val vs KR approx F-test pval
plot_pval_thresh <- 1
plot(x = LR_NoREML_pval[LR_NoREML_pval <= plot_pval_thresh & KR_lmer_proj2$pval_proj2 <= plot_pval_thresh], 
     y = KR_lmer_proj2$pval_proj2[LR_NoREML_pval <= plot_pval_thresh & KR_lmer_proj2$pval_proj2 <= plot_pval_thresh], 
     pch = ".", xlab = "LRT chisq pval (REML=FALSE)", ylab = "KR F-test pval (REML=TRUE)")
abline(a= 0, b = 1, col = "blue", lty = 2)

# check if KR approx F-test signif peptides nested within LRT signif peptides
sum(as.numeric( KR_lmer_proj2$PROBE_ID[KR_LMM$LMM_BH <= 0.05] %in% 
                  LMM_proj2_NoREML$PROBE_ID[LR_NoREML$LMM_BH <= 0.05] )) == 
  length( KR_lmer_proj2$PROBE_ID[KR_LMM$LMM_BH <= 0.05] )

# compare LRT vs Satterthwaite
plot(x = LR_NoREML_pval[KR_lmer_proj2$satterthwaite_pval <= plot_pval_thresh & 
                          LR_NoREML_pval <= plot_pval_thresh], 
     y = KR_lmer_proj2$satterthwaite_pval[KR_lmer_proj2$satterthwaite_pval <= plot_pval_thresh & 
                                            LR_NoREML_pval <= plot_pval_thresh], 
     pch = ".", xlab = "LRT chisq pval (REML=FALSE)", ylab = "Satterthwaite approx pval (REML=TRUE)")
abline(a= 0, b = 1, col = "blue", lty = 2)


####################################################################################### 
#                         ANOVA Visualization (Overall)                               #
#######################################################################################

BH_FDR_cutoff <- 0.05
# Eff_time0 <- 1
Eff_time3 <- 1
Eff_time6 <- 1
signif_crit <- (KR_LMM$LMM_BH <= BH_FDR_cutoff) & 
  (Satterthwaite_LMM$LMM_BH <= BH_FDR_cutoff) &
  # ( LMM_proj2_NoREML$Treat_Main  >= Eff_time0 ) &
  ( (LMM_proj2_NoREML$Treat_Main + 3 * LMM_proj2_NoREML$Treat_Time) >= Eff_time3 ) &
  ( (LMM_proj2_NoREML$Treat_Main + 6 * LMM_proj2_NoREML$Treat_Time) >= Eff_time6 )
sum(as.numeric(signif_crit))  # check

anova_dat_demean <- as.matrix( KR_lmer_proj2 %>% select(contains("id:")) )[ signif_crit , ]
treat_time <- paste(
  toupper( substr(colnames(anova_dat_demean), 4,6) ),
  str_sub( colnames(anova_dat_demean), -1,-1 ), sep = "_"
) 
anova_dat_demean <- sweep(anova_dat_demean, 1,  rowMeans(anova_dat_demean), "-") # centering by row


# colors and shapes for the visualization techniques
cols = pal_proj2[ match(treat_time, names(pal_proj2)) ]
shapes = shp_proj2[  match(treat_time, names(shp_proj2)) ]

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
legend('topright', pch = shp_proj2, col = pal_proj2, cex = 0.5, names(pal_proj2) )
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
legend('topright', pch = shp_proj2, col = pal_proj2, names(pal_proj2))

dev.off()

#---------------------------------------------------------------------------------------
# Heatmap

# heatmap color palette
cls <- colorRampPalette(c("navy", "blue1", "white", "firebrick1", "brown4"))(n = 1024)

heatmap3(anova_dat_demean, 
         col = cls, # specify colors 
         ColSideColors = cols, # specify patient color code
         labCol = colnames(anova_dat_demean), # specify patient
         ColSideLabs = "treatment_time", 
         labRow = "",
         xlab = "Patient_time",
         # legendfun=function() showLegend(col = pal_proj2,
         #                                 legend = names(pal_proj2),
         #                                 cex = 1.2,
         #                                 lwd = 5  )
)


####################################################################################### 
#                         ANOVA Visualization (Time 3)                                #
#######################################################################################

BH_FDR_cutoff <- 0.01
# Eff_time0 <- 1
Eff_time3 <- 1
Eff_time6 <- 1
signif_crit <- (KR_LMM$LMM_BH <= BH_FDR_cutoff) & 
  (Satterthwaite_LMM$LMM_BH <= BH_FDR_cutoff) &
  # ( LMM_proj2_NoREML$Treat_Main  >= Eff_time0 ) &
  ( (LMM_proj2_NoREML$Treat_Main + 3 * LMM_proj2_NoREML$Treat_Time) >= Eff_time3 ) &
  ( (LMM_proj2_NoREML$Treat_Main + 6 * LMM_proj2_NoREML$Treat_Time) >= Eff_time6 )
sum(as.numeric(signif_crit))  # check

anova_dat_demean <- as.matrix( KR_lmer_proj2 %>% select(contains("id:")) )[ signif_crit , ]
treat_time <- paste(
  toupper( substr(colnames(anova_dat_demean), 4,6) ),
  str_sub( colnames(anova_dat_demean), -1,-1 ), sep = "_"
) 
anova_dat_demean3 <- anova_dat_demean[ , (treat_time == "ADT_3") | (treat_time == "PAP_3") ]
anova_dat_demean3 <- sweep(anova_dat_demean3, 1,  rowMeans(anova_dat_demean3), "-") # centering by row

# colors and shapes for the visualization techniques
cols_time3 = pal_proj2[ match(treat_time[(treat_time == "ADT_3") | (treat_time == "PAP_3")], names(pal_proj2)) ]
shapes_time3 = shp_proj2[  match(treat_time[(treat_time == "ADT_3") | (treat_time == "PAP_3")], names(shp_proj2)) ]
pal_time3 <- c("blue", "red")
names(pal_time3) <- c("ADT_3", "PAP_3")
cols_time3 = pal_time3[ match(treat_time[(treat_time == "ADT_3") | (treat_time == "PAP_3")], names(pal_time3)) ]

#---------------------------------------------------------------------------------------
# PCA -- focus on heatmap at time 3?

# svd
sv.dat <- sweep(t(anova_dat_demean3), 2, colMeans(t(anova_dat_demean3)), "-") # centering
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
PCload.func(cols_time3, shapes_time3, U, D, 1, 2, pca.var, title = "PC2 vs PC1") # PC loadings (PC2 vs PC1)
legend('topright', cex = 0.5, pch = shp_proj2[c('ADT_3','PAP_3')], col = pal_proj2[c('ADT_3','PAP_3')], c('ADT_3','PAP_3') )
PCload.func(cols_time3, shapes_time3, U, D, 3, 2, pca.var, title = "PC2 vs PC3") # PC loadings (PC2 vs PC3)

dev.off()

#---------------------------------------------------------------------------------------
# t-SNE -- focus on heatmap at time 3?

# how to specify parameter
# refer: https://lvdmaaten.github.io/tsne/User_guide.pdf

tsnedat <- unname(t(anova_dat_demean3)) # dim(X) = N samples by D dimensions 
initdim <- 40 
perplex <- 13 

set.seed(10)
tsne_anova <- Rtsne(tsnedat, initial_dims = initdim, perplexity = perplex,
                    theta = 0, check_duplicates = F, max_iter = 50000L, eta = 50)

# t-SNE plot
plot(tsne_anova$Y, ylab = "", xlab = "", col = cols_time3, pch = shapes_time3, main = "t-SNE plot")
# legend('topleft', pch = shp_proj2[c('ADT_3','PAP_3')], col = pal_proj2[c('ADT_3','PAP_3')], c('ADT_3','PAP_3'))
legend('topleft', pch = shp_proj2[c('ADT_3','PAP_3')], col = pal_time3, c('ADT_3','PAP_3'))

dev.off()

#---------------------------------------------------------------------------------------
# Heatmap -- focus on heatmap at time 3?


# heatmap color palette
cls <- colorRampPalette(c("navy", "royalblue1", "white", "firebrick1", "brown4"))(n = 1024)
cls <- colorRampPalette(c("navy", "white", "firebrick3"))(n = 1024)

heatmap3(anova_dat_demean3, 
         col = cls, # specify colors 
         ColSideColors = cols_time3, # specify patient color code
         labCol = colnames(anova_dat_demean3), # specify patient
         ColSideLabs = "treatment_time", 
         labRow = "",
         xlab = "Patient_time",
         # legendfun=function() showLegend(col = pal_proj2[c('ADT_3','PAP_3')],
         #                                 legend = c('ADT_3','PAP_3'),
         #                                 cex = 1.2,
         #                                 lwd = 5  )
)

####################################################################################### 
#                         ANOVA Visualization (Time 6)                                #
#######################################################################################

BH_FDR_cutoff <- 0.01
# Eff_time0 <- 1
Eff_time3 <- 1
Eff_time6 <- 1
signif_crit <- (KR_LMM$LMM_BH <= BH_FDR_cutoff) & 
  (Satterthwaite_LMM$LMM_BH <= BH_FDR_cutoff) &
  # ( LMM_proj2_NoREML$Treat_Main  >= Eff_time0 ) &
  ( (LMM_proj2_NoREML$Treat_Main + 3 * LMM_proj2_NoREML$Treat_Time) >= Eff_time3 ) &
  ( (LMM_proj2_NoREML$Treat_Main + 6 * LMM_proj2_NoREML$Treat_Time) >= Eff_time6 )
sum(as.numeric(signif_crit))  # check

anova_dat_demean <- as.matrix( KR_lmer_proj2 %>% select(contains("id:")) )[ signif_crit , ]
treat_time <- paste(
  toupper( substr(colnames(anova_dat_demean), 4,6) ),
  str_sub( colnames(anova_dat_demean), -1,-1 ), sep = "_"
) 
anova_dat_demean6 <- anova_dat_demean[ , (treat_time == "ADT_6") | (treat_time == "PAP_6") ]
anova_dat_demean6 <- sweep(anova_dat_demean6, 1,  rowMeans(anova_dat_demean6), "-") # centering by row

# colors and shapes for the visualization techniques
cols_time6 = pal_proj2[ match(treat_time[(treat_time == "ADT_6") | (treat_time == "PAP_6")], names(pal_proj2)) ]
shapes_time6 = shp_proj2[  match(treat_time[(treat_time == "ADT_6") | (treat_time == "PAP_6")], names(shp_proj2)) ]
pal_time6 <- c("blue", "red")
names(pal_time6) <- c("ADT_6", "PAP_6")
cols_time6 = pal_time6[ match(treat_time[(treat_time == "ADT_6") | (treat_time == "PAP_6")], names(pal_time6)) ]

#---------------------------------------------------------------------------------------
# PCA -- focus on heatmap at time 6?

# svd
sv.dat <- sweep(t(anova_dat_demean6), 2, colMeans(t(anova_dat_demean6)), "-") # centering
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
PCload.func(cols_time6, shapes_time6, U, D, 1, 2, pca.var, title = "PC2 vs PC1") # PC loadings (PC2 vs PC1)
legend('topright', cex = 0.5, pch = shp_proj2[c('ADT_6','PAP_6')], col = pal_proj2[c('ADT_6','PAP_6')], c('ADT_6','PAP_6') )
PCload.func(cols_time6, shapes_time6, U, D, 3, 2, pca.var, title = "PC2 vs PC3") # PC loadings (PC2 vs PC3)

dev.off()

#---------------------------------------------------------------------------------------
# t-SNE -- focus on heatmap at time 6?

# how to specify parameter
# refer: https://lvdmaaten.github.io/tsne/User_guide.pdf

tsnedat <- unname(t(anova_dat_demean6)) # dim(X) = N samples by D dimensions 
initdim <- 40 
perplex <- 13 

set.seed(10)
tsne_anova <- Rtsne(tsnedat, initial_dims = initdim, perplexity = perplex,
                    theta = 0, check_duplicates = F, max_iter = 50000L, eta = 50)

# t-SNE plot
plot(tsne_anova$Y, ylab = "", xlab = "", col = cols_time6, pch = shapes_time6, main = "t-SNE plot")
# legend('topleft', pch = shp_proj2[c('ADT_6','PAP_6')], col = pal_proj2[c('ADT_6','PAP_6')], c('ADT_6','PAP_6'))
legend('topleft', pch = shp_proj2[c('ADT_6','PAP_6')], col = pal_time6, c('ADT_6','PAP_6'))

dev.off()

#---------------------------------------------------------------------------------------
# Heatmap -- focus on heatmap at time 6?


# heatmap color palette
cls <- colorRampPalette(c("navy", "royalblue1", "white", "firebrick1", "brown4"))(n = 1024)

heatmap3(anova_dat_demean6, 
         col = cls, # specify colors 
         ColSideColors = cols_time6, # specify patient color code
         labCol = colnames(anova_dat_demean6), # specify patient
         ColSideLabs = "treatment_time", 
         labRow = "",
         xlab = "Patient_time",
         # legendfun=function() showLegend(col = pal_proj2[c('ADT_6','PAP_6')],
         #                                 legend = c('ADT_6','PAP_6'),
         #                                 cex = 1.2,
         #                                 lwd = 5  )
)

dev.off()


####################################################################################### 
#                          Boxplots of signif peptides                                #
#######################################################################################

BH_FDR_cutoff <- 0.0001
# Eff_time0 <- 1
Eff_time3 <- 2
Eff_time6 <- 2
signif_crit <- (KR_LMM$LMM_BH <= BH_FDR_cutoff) & 
  (Satterthwaite_LMM$LMM_BH <= BH_FDR_cutoff) &
  # ( LMM_proj2_NoREML$Treat_Main  >= Eff_time0 ) &
  ( (LMM_proj2_NoREML$Treat_Main + 3 * LMM_proj2_NoREML$Treat_Time) >= Eff_time3 ) &
  ( (LMM_proj2_NoREML$Treat_Main + 6 * LMM_proj2_NoREML$Treat_Time) >= Eff_time6 )
sum(as.numeric(signif_crit))  # check

signif_median <- raw_data_median_proj2 %>% 
  select(PROBE_ID, contains("id:")) %>%
  filter( signif_crit )

signif_total_eff <- 2*(LMM_proj2_NoREML$Treat_Main)[signif_crit] + 9*(LMM_proj2_NoREML$Treat_Time)[signif_crit]
signif_median <- signif_median[ order(signif_total_eff, decreasing = T),]
head(signif_median) # check

signif_boxplot.func <- function(signif_mat, draw){
  signif_mat2 <- signif_mat[,-1]
  signif_df <- data.frame(
    treatment = factor( toupper( substr(colnames(signif_mat2), 4,6) ) ),
    time = factor( str_sub( colnames(signif_mat2), -1,-1 ) ), 
    fluorescence = as.numeric(signif_mat2[draw,])
  )
  ggplot(signif_df, aes(x = time, y = fluorescence, fill = treatment)) +
    geom_boxplot(width = 0.5, position=position_dodge2(width = 0.5)) +
    labs(title = paste0("Boxplots of Fluorescence Levels for \nSignificant Peptide: ", 
                        signif_mat$PROBE_ID[draw]), 
         x = "Time", y = "log2 Median Fluorescence") +
    scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
    # ylim(c(2,16)) +
    theme(panel.background = element_rect(fill = "grey90"),
          panel.grid.major = element_line(color = "white"),
          panel.grid.minor = element_line(color = "white"),
          # legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))
}

grid.arrange(
  signif_boxplot.func(signif_median,1),
  signif_boxplot.func(signif_median,2),
  signif_boxplot.func(signif_median,3),
  signif_boxplot.func(signif_median,4),
  signif_boxplot.func(signif_median,5),
  signif_boxplot.func(signif_median,6),
  signif_boxplot.func(signif_median,7),
  signif_boxplot.func(signif_median,8),
  signif_boxplot.func(signif_median,9),
  ncol = 3
)

####################################################################################### 
#                                   Gene Set Analysis                                 #
#######################################################################################

BH_FDR_cutoff <- 0.01
# Eff_time0 <- 1
Eff_time3 <- 1
Eff_time6 <- 1
signif_crit <- (KR_LMM$LMM_BH <= BH_FDR_cutoff) & 
  (Satterthwaite_LMM$LMM_BH <= BH_FDR_cutoff) &
  # ( LMM_proj2_NoREML$Treat_Main  >= Eff_time0 ) &
  ( (LMM_proj2_NoREML$Treat_Main + 3 * LMM_proj2_NoREML$Treat_Time) >= Eff_time3 ) &
  ( (LMM_proj2_NoREML$Treat_Main + 6 * LMM_proj2_NoREML$Treat_Time) >= Eff_time6 )
sum(as.numeric(signif_crit))  # check

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

# what are the repeated genes?
repeat_gene <- uniprot_gene[ uniprot_gene$entrez_gene_id_pete %in%
                (uniprot_gene %>%
                   group_by(entrez_gene_id_pete) %>% 
                   tally() %>%
                   filter(n > 1) %>%
                   pull(entrez_gene_id_pete)) ,] %>% select(-protein_names)

# are these repeated genes among the ones that are signif?
repeat_gene$seq_id %in% unique( raw_data_median_proj2$SEQ_ID[signif_crit] )

# how many signif proteins?
length(unique( raw_data_median_proj2$SEQ_ID[signif_crit] ))

## manually changing entrez_gene_id for PCA10 & PRO29
uniprot_gene$entrez_gene_id_pete[uniprot_gene$seq_id == "PCA10"] <- 28912
uniprot_gene$entrez_gene_id_pete[uniprot_gene$seq_id == "PRO29"] <- NA
uniprot_gene <- uniprot_gene[!(is.na(uniprot_gene$entrez_gene_id_pete)),]

# get unique seq_id that are associated with significant peptides
signif_seq_id <- unique( raw_data_median_proj2$SEQ_ID[signif_crit] )
signif_seq_id <- signif_seq_id[!(is.na(signif_seq_id))] # just in case

# convert these into binary vector
seq_id_ok <- as.numeric(uniprot_gene$seq_id %in% signif_seq_id)
names(seq_id_ok) <- uniprot_gene$entrez_gene_id_pete

# check 
length(seq_id_ok)
sum(seq_id_ok)

# gene-set analysis via allez!
allez.go <- allez(seq_id_ok, lib = "org.Hs.eg", idtype = "ENTREZID", sets = "GO")

# let's see results
nom.alpha <- 0.05
min.num.gene <- 2

# Extract a table of top-ranked functional sets from allez output
allezTable(allez.go, symbol = T, n.cell = min.num.gene, nominal.alpha = nom.alpha, in.set = T)[,c(1:5,7)]

# Display an image of gene scores by functional sets
allezPlot(allez.go, n.cell = min.num.gene, nominal.alpha = nom.alpha)

# tabulate enriched gene-set with signif genes only
allez.tab <- allezTable(allez.go, symbol = T, n.cell = min.num.gene, nominal.alpha = nom.alpha, in.set = T)
allez.tab$set.size <- paste(allez.tab$in.set, allez.tab$set.size, sep = "/")
allez.tab <- allez.tab %>% dplyr::select(-c(in.set, genes)) %>%
  mutate(in.genes = str_replace_all(in.genes, ";", "; "))
