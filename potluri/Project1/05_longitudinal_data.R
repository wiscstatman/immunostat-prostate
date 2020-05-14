library(lme4) # linear mixed effects model
library(pbkrtest) # LMM LRT
library(fdrtool)
library(Rtsne)
library(heatmap3)
library(ggplot2)
library(gridExtra)
library(tidyverse) # make sure you have the latest tidyverse !


####################################################################################### 
#                           Some Common Variables and Functions                       #
####################################################################################### 

pal_proj2 <- c("turquoise1", "cornflowerblue","navy", "orchid1", "darkorange1", "firebrick1")

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


load("KR_lmer_proj2.RData")
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
#                     Linear Mixed Model to Assess Treatment Effect                   #
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

for(i in 1:nrow(raw_data_median_proj2)){
  y <- as.numeric(raw_data_median_proj2[i, (ncol_med_proj2 - n_med_proj2 + 1) : ncol_med_proj2])
  fit0 <- lmer(y ~ time + (1 + time | id), data = sample_key_proj2)
  fit1 <- lmer(y ~ treatment*time + (1 + time | id), data = sample_key_proj2)
  fit2 <- lmer(y ~ time + (1  | id), data = sample_key_proj2)
  fit3 <- lmer(y ~ treatment*time + (1  | id), data = sample_key_proj2)
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

save(KR_lmer_proj2, file = "KR_lmer_proj2.RData")

# FDR adjustment of p-values
KR_LMM <- LMM_func_proj2(KR_lmer_proj2$pval_proj2)

# peptide counts at various FDR thresholds
count.func(KR_LMM$LMM_BH, seq(0.01, 0.1, by = 0.01))
count.func(KR_LMM$LMM_qval, seq(0.01, 0.1, by = 0.01))

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
