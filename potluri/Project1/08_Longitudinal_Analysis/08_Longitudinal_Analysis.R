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
  countab <- rbind(c("FDR threshold", thresh.vec),
                   c("Peptide counts", counter))
  return(countab)
}

raw_data_median_proj2 <- read_csv("raw_data_median_proj2.csv")
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
#                                 Test Time Random Effects                            #
####################################################################################### 

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

