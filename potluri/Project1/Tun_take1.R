library(tidyverse) # make sure you have the latest tidyverse !
library(fdrtool)

load("median_data.RData")
load("anova_pval.RData")

####################################################################################### 
#                                      Data Processing                                #
####################################################################################### 

# Hemanth excluded patients without replicates

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

# check patient counts (NOT distinct patients)
array_id_key %>%
  group_by(id, stage) %>%
  summarize() %>%
  group_by(stage) %>%
  tally()

sample_key = array_id_key %>%
  group_by(id, stage) %>%
  summarize(n = n()) %>%  # assign (id,stage) group counts to the variable "n"
  ungroup() %>% 
  filter(n >= 2) %>%
  select(-n) # drop the column "n"

array_id_key = array_id_key %>%
  filter(id %in% sample_key$id)

# after dropping patients without replicates, check patient counts (NOT distinct patients)
array_id_key %>%
  group_by(id, stage) %>%
  summarize() %>%
  group_by(stage) %>%
  tally()

# drop unnecessary columns from array_id_key
length(array_id_key$rep==1) == length(array_id_key$rep==2) # replicates column no longer helpful
table(array_id_key$notes); subset(array_id_key, is.na(notes)==F) # need to keep notes
array_id_key = array_id_key %>%
  select( -c(index, rep) ) # Hemanth says index column not helpful
  
#-----------------------------------------------------------------------------------------------
# Ensuring only distinct patients

# 11 patients with replicates in 2 stages
repeated_patient_id_pair <- t( matrix( c(
    "mcv005",      "pdv004",      
    "mcv119",      "arv003",      
    "mcv083",      "pap092",      
    "mcv003",      "pap032",      
    "adt034",      "pap079",      
    "adt181",      "pap001",      
    "adt223",      "pap054",      
    "vpd037",      "pdv008",      
    "pap123",      "vpd041",      
    "pap067",      "vpd016",      
    "adt143",      "pdv104"), nrow=2) ) 

cbind( match( repeated_patient_id_pair[,1], array_id_key$id ),
       match( repeated_patient_id_pair[,2], array_id_key$id ) )

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

array_id_key = array_id_key %>%
  filter(!(id %in% ids_to_remove))

cbind( match( repeated_patient_id_pair[,1], array_id_key$id ),
       match( repeated_patient_id_pair[,2], array_id_key$id ) )

# check patient counts (finally distinct patients)
distinct_patients_with_reps <- array_id_key %>%
  group_by(id, stage) %>%
  summarize() %>%
  group_by(stage) %>%
  tally()

#-----------------------------------------------------------------------------------------------
# Now read in fluorescence intensity

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
  select(PROBE_DESIGN_ID:Y, any_of(array_id_key$array_id)) %>% # drop patients without replicates 
  select( -c(X, Y, MATCH_INDEX, DESIGN_NOTE, SELECTION_CRITERIA, MISMATCH, PROBE_CLASS) ) # drop unhelpful columns


#-----------------------------------------------------------------------------------------------
# find median of replicates

iii <- match( colnames( select(raw_data, contains("dat")) ), array_id_key$array_id )
array_id_key <- array_id_key[iii,]
sum(as.numeric( colnames(select(raw_data, contains("dat"))) == array_id_key$array_id )) == nrow(array_id_key)
# now array_id_key and raw_data align ! 

array_id_key$id <- as.factor(array_id_key$id)
raw_data_median <- t( apply( select(raw_data, contains("dat")), 1, function(x) {
  tapply(x, array_id_key$id, FUN=median)
} ) )
raw_data_median <- log2(raw_data_median)
raw_data_median <- bind_cols( select(raw_data, -contains("dat")), as.data.frame(raw_data_median) ) 
rm(iii); gc()

####################################################################################### 
#                                       One-Way ANOVA                                 #
#######################################################################################

# how many distinct patients with replicates
n <- sum(distinct_patients_with_reps$n); n

# making sure patient ID's are unique in raw_data_median
length(unique( tail(colnames(raw_data_median), n) )) == length( tail(colnames(raw_data_median), n) )

# making sure patient ID's are unique in sample_key
length(unique( sample_key$id )) == nrow(sample_key)

# get matching stage for patients
stage_iii <- match( tail(colnames(raw_data_median), n), sample_key$id )
stage <- sample_key$stage[stage_iii]
stage <- factor(stage) # remove unused levels
sum(as.numeric( tail(colnames(raw_data_median), n) == sample_key$id[stage_iii] )) == n
rm(stage_iii); gc()

# initiate
one_way_anova_pval <- rep(NA, nrow(raw_data_median))

# punch it!

median_ncol <- ncol(raw_data_median)
for(i in 1:nrow(raw_data_median)){
  one_way_anova_pval[i] <- unname( unlist(summary( aov( 
    as.numeric(raw_data_median[i, (median_ncol-n+1) : median_ncol]) ~ stage ) ))["Pr(>F)1"] )
  print(i)
}

# save(one_way_anova_pval, file = "anova_pval.RData")

# plot p-values histogram
hist(one_way_anova_pval, breaks = 70, freq = F, xlab = "one-way ANOVA p-values", main = "p-values distribution \nfor peptides")

one_way_anova_BH <- p.adjust(one_way_anova_pval, method = "BH")
one_way_anova_qval <- fdrtool(one_way_anova_pval, statistic = "pvalue", verbose = F, plot  = F)$qval
qval_eta0 <- unname(fdrtool(one_way_anova_pval, statistic = "pvalue", verbose = F, plot  = F)$param[,"eta0"])

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

# peptide counts at various FDR thresholds
count.func(one_way_anova_BH, seq(0.01, 0.1, by = 0.01))
count.func(one_way_anova_qval, seq(0.01, 0.1, by = 0.01))

# significant peptides
raw_data_median[ which(one_way_anova_qval <= 0.05) , 1: (median_ncol-n)]

# plot histogram with proportion of non-null peptides
all_peptide_hist <- hist(one_way_anova_pval, breaks = 70, freq = F, 
                         xlab = "one-way ANOVA p-values", main = "p-values distribution for peptides")
polygon_ind <- which(all_peptide_hist$density >= qval_eta0)
for (i in polygon_ind){
  polygon( x = c(all_peptide_hist$breaks[i], all_peptide_hist$breaks[i+1], all_peptide_hist$breaks[i+1], all_peptide_hist$breaks[i]),
           y = c(all_peptide_hist$density[i], all_peptide_hist$density[i], qval_eta0, qval_eta0),
           col = "red")
}

####################################################################################### 
#                   Marginal Variance Filtering within each protein                   #
#######################################################################################

raw_data_median$marginal_var <- apply( raw_data_median[, (median_ncol-n+1) : median_ncol], 1, var )

raw_data_median <- raw_data_median %>%
  group_by(SEQ_ID) %>%
  mutate(rank_variance = n() - rank(marginal_var) + 1) %>%
  ungroup()

# number of proteins
length(unique(raw_data_median$SEQ_ID)) # 1611 protein

# check ranking counts
table(raw_data_median$rank_variance[raw_data_median$rank_variance<=7])
# because default ties.method = average

# remove those .5 ranking
raw_data_median$rank_variance = floor(raw_data_median$rank_variance)
table(raw_data_median$rank_variance[raw_data_median$rank_variance<=7])

# save(raw_data_median, file = "median_data.RData")


# marginal variance filtering function
marginal_var_filter.func <- function(pval_vector, topnum){
  pval_filter = pval_vector[raw_data_median$rank_variance <= topnum]
  BH_filter = p.adjust(pval_filter, method = "BH")
  qval_filter = fdrtool(pval_filter, statistic = "pvalue", verbose = F, plot  = F)$qval
  qval_eta0_filter = unname(fdrtool(pval_filter, statistic = "pvalue", verbose = F, plot  = F)$param[,"eta0"])
  hist_plot = hist(pval_filter, breaks = 70, freq = F, xlab = "one-way ANOVA p-values", 
                   main = paste0("p-values distribution for top ", topnum, " peptides with largest variance within protein") )
  return(list(
    pval_filter = pval_filter,
    BH_filter = BH_filter,
    qval_filter = qval_filter,
    qval_eta0_filter = qval_eta0_filter,
    hist = hist_plot
  ))
}

top <- marginal_var_filter.func(one_way_anova_pval, 5)
top$qval_eta0_filter
count.func(top$BH_filter, seq(0.01, 0.1, by = 0.01))
count.func(top$qval_filter, seq(0.01, 0.1, by = 0.01))
