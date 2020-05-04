library(tidyverse) # make sure you have the latest tidyverse !
library(fdrtool)
library(Rtsne)
library(heatmap3)
library(xlsx)

# median_double_log2 <- read.table("ANOVA_median_log2log2.csv", header = T, row.names = NULL, sep = ",")
load("remove_rep1.RData")
raw_data_median <- read.table("ANOVA_median_remove_rep1.csv", header = T, row.names = NULL, sep = ",")

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

# We keep patients with rep = 1

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

# check patient counts (NOT distinct patients)
array_id_key %>%
  group_by(id, stage) %>%
  summarize() %>%
  group_by(stage) %>%
  tally()

# if you want to remove patients whose rep == 1
sample_key = array_id_key %>%
  group_by(id, stage) %>%
  summarize(n = n()) %>%  # assign (id,stage) group counts to the variable "n"
  ungroup() %>% 
  filter(n >= 2) %>%
  select(-n) # drop the column "n"
array_id_key = array_id_key %>%
  filter(id %in% sample_key$id)

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

# check all these pairs are present in the array_id_key dataframe
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

# drop patients' earlier records
array_id_key = array_id_key %>%
  filter(!(id %in% ids_to_remove))

# now each of these patients has 1 set of records
cbind( match( repeated_patient_id_pair[,1], array_id_key$id ),
       match( repeated_patient_id_pair[,2], array_id_key$id ) )

# check patient counts (finally distinct patients)
distinct_patients_with_reps <- array_id_key %>%
  group_by(id, stage) %>%
  summarize() %>%
  group_by(stage) %>%
  tally()
distinct_patients_with_reps

# how many distinct patients with replicates
n <- sum(distinct_patients_with_reps$n); n

# group unique patients by their id
sample_key = array_id_key %>%
  group_by(id, stage) %>%
  summarize(rep = n()) %>%  # assign (id,stage) group counts to the variable "n"
  ungroup()

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
  select(PROBE_DESIGN_ID:Y, any_of(array_id_key$array_id)) %>% # drop patients with records at different stages 
  select( -c(X, Y, MATCH_INDEX, DESIGN_NOTE, SELECTION_CRITERIA, MISMATCH, PROBE_CLASS) ) # drop unhelpful columns

colnames(raw_data)[1:15] # check

#-----------------------------------------------------------------------------------------------
# take log2 transformation
# then compute median of replicates

# verify array_id's are unique
length(unique(colnames( select(raw_data, contains("dat")) ))) == length(colnames( select(raw_data, contains("dat")) ))
length(unique(array_id_key$array_id )) == length(array_id_key$array_id )
length(colnames( select(raw_data, contains("dat")) )) == length(array_id_key$array_id )

# align array_id in array_id_key with the ones in raw_data
iii <- match( colnames( select(raw_data, contains("dat")) ), array_id_key$array_id )
array_id_key <- array_id_key[iii,]
sum(as.numeric( colnames(select(raw_data, contains("dat"))) == array_id_key$array_id )) == nrow(array_id_key)
# now array_id_key and raw_data align ! 

# now take log2 transformation
raw_data <- raw_data %>%
  mutate_at(vars(matches("dat")), log2)

# also try double_log2 transformation
double_log2 <-  raw_data %>%
  mutate_at(vars(matches("dat")), log2)

# now compute median
array_id_key$id <- as.factor(array_id_key$id)
raw_data_median <- t( apply( select(raw_data, contains("dat")), 1, function(x) {
  tapply(x, array_id_key$id, FUN=median)
} ) )
raw_data_median <- bind_cols( select(raw_data, -contains("dat")), as.data.frame(raw_data_median) ) 

# also compute median for double_log2 
median_double_log2 <- t( apply( select(double_log2, contains("dat")), 1, function(x) {
  tapply(x, array_id_key$id, FUN=median)
} ) )
median_double_log2 <- bind_cols( select(double_log2, -contains("dat")), as.data.frame(median_double_log2) ) 

####################################################################################### 
#                                       One-Way ANOVA                                 #
#######################################################################################

# making sure patient ID's are unique in raw_data_median
length(unique( tail(colnames(raw_data_median), n) )) == length( tail(colnames(raw_data_median), n) )
length(unique( tail(colnames(median_double_log2), n) )) == length( tail(colnames(median_double_log2), n) )

# making sure patient ID's are unique in sample_key
length(unique( sample_key$id )) == nrow(sample_key)

# get matching stage for patients
stage_iii <- match( tail(colnames(raw_data_median), n), sample_key$id )
sum(as.numeric( sample_key$id[stage_iii] == tail(colnames(raw_data_median), n) )) == n # check
stage <- sample_key$stage[stage_iii]
stage <- factor(stage, levels = c("normal","new_dx","nmCSPC","nmCRPC","mCRPC")) # order levels of stages, matter for Tukey HSD 

stage_iii_double_log2 <- match( tail(colnames(median_double_log2), n), sample_key$id )
sum(as.numeric( sample_key$id[stage_iii_double_log2] == tail(colnames(median_double_log2), n) )) == n # check
stage_double_log2 <- sample_key$stage[stage_iii_double_log2]
stage_double_log2 <- factor(stage_double_log2) # remove unused levels
sum(as.numeric(stage_double_log2 == stage)) == nrow(sample_key)

#----------------------------------------------------------------------------------------
# initiate
one_way_anova_pval <- rep(NA, nrow(raw_data_median))
anova_pval_double_log2 <- rep(NA, nrow(median_double_log2))

# punch it!

median_ncol <- ncol(raw_data_median)
for(i in 1:nrow(raw_data_median)){
  one_way_anova_pval[i] <- unname( unlist(summary( aov( 
    as.numeric(raw_data_median[i, (median_ncol-n+1) : median_ncol]) ~ stage ) ))["Pr(>F)1"] )
  print(i)
}
names(one_way_anova_pval) <- raw_data_median$PROBE_ID

for(i in 1:nrow(median_double_log2)){
  anova_pval_double_log2[i] <- unname( unlist(summary( aov( 
    as.numeric(median_double_log2[i, (median_ncol-n+1) : median_ncol]) ~ stage_double_log2 ) ))["Pr(>F)1"] )
  print(i)
}
names(anova_pval_double_log2) <- raw_data_median$PROBE_ID

# plot p-values histogram
hist(one_way_anova_pval, breaks = 70, freq = F, xlab = "one-way ANOVA p-values", main = "p-values distribution \nfor peptides")
hist(anova_pval_double_log2, breaks = 70, freq = F, xlab = "one-way ANOVA p-values", main = "p-values distribution \nfor peptides")

#----------------------------------------------------------------------------------------
one_way_anova_BH <- p.adjust(one_way_anova_pval, method = "BH")
one_way_anova_qval <- fdrtool(one_way_anova_pval, statistic = "pvalue", verbose = F, plot  = F)$qval
qval_eta0 <- unname(fdrtool(one_way_anova_pval, statistic = "pvalue", verbose = F, plot  = F)$param[,"eta0"])

anova_BH_double_log2 <- p.adjust(anova_pval_double_log2, method = "BH")
anova_qval_double_log2 <- fdrtool(anova_pval_double_log2, statistic = "pvalue", verbose = F, plot  = F)$qval
qval_eta0_double_log2 <- unname(fdrtool(anova_pval_double_log2, statistic = "pvalue", verbose = F, plot  = F)$param[,"eta0"])

# peptide counts at various FDR thresholds
count.func(one_way_anova_BH, seq(0.01, 0.1, by = 0.01))
count.func(one_way_anova_qval, seq(0.01, 0.1, by = 0.01))
count.func(anova_BH_double_log2, seq(0.01, 0.1, by = 0.01))
count.func(anova_qval_double_log2, seq(0.01, 0.1, by = 0.01))

#----------------------------------------------------------------------------------------
# plot histogram with proportion of non-null peptides
all_peptide_hist <- hist(one_way_anova_pval, breaks = 70, freq = F, 
                         xlab = "one-way ANOVA p-values", main = "p-values distribution for peptides")
polygon_ind <- which(all_peptide_hist$density >= qval_eta0)
for (i in polygon_ind){
  polygon( x = c(all_peptide_hist$breaks[i], all_peptide_hist$breaks[i+1], all_peptide_hist$breaks[i+1], all_peptide_hist$breaks[i]),
           y = c(all_peptide_hist$density[i], all_peptide_hist$density[i], qval_eta0, qval_eta0),
           col = "red")
}
text(x=0.65,y=4, labels = paste( "estimated proportion of \nnon-null peptides =",
                                 round( 100*(1 - qval_eta0),2 ),"%" ))

# get corresponding protein id
signif_protein <- unique( raw_data_median$SEQ_ID[one_way_anova_BH <= 0.05]  )
signif_protein <- unique( raw_data_median$SEQ_ID[one_way_anova_qval <= 0.05]  )

####################################################################################### 
#                   Marginal Variance Filtering within each protein                   #
#######################################################################################

# compute marginal variance
raw_data_median$marginal_var <- apply( raw_data_median[, (median_ncol-n+1) : median_ncol], 1, var )
median_double_log2$marginal_var <- apply( median_double_log2[, (median_ncol-n+1) : median_ncol], 1, var )

# check number of proteins
length(unique(raw_data_median$SEQ_ID)) # 1611 protein

# rank peptides by marginal variance within each protein
raw_data_median <- raw_data_median %>%
  group_by(SEQ_ID) %>%
  mutate(rank_variance = n() - rank(marginal_var) + 1) %>%
  ungroup()

median_double_log2 <- median_double_log2 %>%
  group_by(SEQ_ID) %>%
  mutate(rank_variance = n() - rank(marginal_var) + 1) %>%
  ungroup()

# check ranking counts
table(raw_data_median$rank_variance[raw_data_median$rank_variance<=7])
table(median_double_log2$rank_variance[median_double_log2$rank_variance<=7])
# because default ties.method = average

# remove those .5 ranking
raw_data_median$rank_variance = floor(raw_data_median$rank_variance)
median_double_log2$rank_variance = floor(median_double_log2$rank_variance)
table(raw_data_median$rank_variance[raw_data_median$rank_variance<=7])
table(median_double_log2$rank_variance[median_double_log2$rank_variance<=7])

#----------------------------------------------------------------------------------------
# rearrange columns
raw_data_median <- raw_data_median[, c(1:10, (ncol(raw_data_median)-1):ncol(raw_data_median), 
                                       11:(ncol(raw_data_median)-2))]
colnames(raw_data_median)[1:15];  tail(colnames(raw_data_median)) # check
median_double_log2 <- median_double_log2[, c(1:10, (ncol(raw_data_median)-1):ncol(raw_data_median), 
                                             11:(ncol(raw_data_median)-2))]
colnames(median_double_log2)[1:15];  tail(colnames(median_double_log2)) # check

# write table
raw_data_median <- raw_data_median %>% 
  add_column(p_values = one_way_anova_pval,
             BH_FDR = one_way_anova_BH,
             Storey_Qvalues = one_way_anova_qval)
raw_data_median <- raw_data_median[, c(1:12, (ncol(raw_data_median)-2):ncol(raw_data_median), 
                                       13:(ncol(raw_data_median)-3))]
colnames(raw_data_median)[1:20];  tail(colnames(raw_data_median));  raw_data_median[1:7,1:5] # check
write.table(raw_data_median, file = "ANOVA_median.csv", sep = ",", row.names = F, col.names = T)

median_double_log2 <- median_double_log2 %>% 
  add_column(p_values = anova_pval_double_log2,
             BH_FDR = anova_BH_double_log2,
             Storey_Qvalues = anova_qval_double_log2)
median_double_log2 <- median_double_log2[, c(1:12, (ncol(raw_data_median)-2):ncol(raw_data_median), 
                                             13:(ncol(raw_data_median)-3))]
colnames(median_double_log2)[1:20];  tail(colnames(median_double_log2));  median_double_log2[1:7,1:5] # check
write.table(median_double_log2, file = "ANOVA_median_log2log2.csv", sep = ",", row.names = F, col.names = T)

#----------------------------------------------------------------------------------------
# marginal variance filtering function
marginal_var_filter.func <- function(pval_vector, topnum){
  # pval_filter = pval_vector[raw_data_median$rank_variance <= topnum]
  pval_filter = pval_vector[rank_variance <= topnum]
  BH_filter = p.adjust(pval_filter, method = "BH")
  qval_filter = fdrtool(pval_filter, statistic = "pvalue", verbose = F, plot  = F)$qval
  qval_eta0_filter = unname(fdrtool(pval_filter, statistic = "pvalue", verbose = F, plot  = F)$param[,"eta0"])
  return(list(
    pval_filter = pval_filter,
    BH_filter = BH_filter,
    qval_filter = qval_filter,
    qval_eta0_filter = qval_eta0_filter
  ))
}

# results
topnum <- 5

top <- marginal_var_filter.func(raw_data_median$p_values, topnum)
1 - top$qval_eta0_filter
count.func(top$BH_filter, seq(0.01, 0.1, by = 0.01))
count.func(top$qval_filter, seq(0.01, 0.1, by = 0.01))

# histogram 
top_plot <- hist(top$pval_filter, breaks = 70, freq = F, xlab = "one-way ANOVA p-values", 
                 main = paste("p-values distribution for top ", topnum, 
                              "peptide(s) with largest variance within protein") )
polygon_ind <- which(top_plot$density >= top$qval_eta0_filter)
for (i in polygon_ind){
  polygon( x = c(top_plot$breaks[i], top_plot$breaks[i+1], top_plot$breaks[i+1], top_plot$breaks[i]),
           y = c(top_plot$density[i], top_plot$density[i], top$qval_eta0_filter, top$qval_eta0_filter),
           col = "red")
}
text(x=0.65,y=4, labels = paste( "estimated proportion of \nnon-null peptides =",
                                 round( 100*(1 - top$qval_eta0_filter),2 ),"%" ))

top_protein <- unique( raw_data_median$SEQ_ID[top$BH_filter <= 0.05]  )



####################################################################################### 
#                                        ANOVA Visualization                          #
#######################################################################################

cols = pal[ match(stage, names(pal)) ]
shapes = shp[  match(stage, names(shp)) ]

# anova_dat <- raw_data_median[ top$BH_filter <= 0.05 ,  16:ncol(raw_data_median)]
anova_dat <- as.matrix(raw_data_median[ (raw_data_median$BH_FDR <= 0.05) , 16:ncol(raw_data_median)])
row.names(anova_dat) <- raw_data_median$PROBE_ID[ (raw_data_median$BH_FDR <= 0.05) ]

anova_dat <- sweep(anova_dat, 1, rowMeans(anova_dat), "-") # centering by row

cols = pal[ match(stage, names(pal)) ]
shapes = shp[  match(stage, names(shp)) ]

# in case you want to remove a certain stage from visualization
# anova_dat <- anova_dat[, (stage != "normal")]
# cols = pal[ match(stage[stage != "normal"], names(pal)) ]
# shapes <- shp[ match(stage[stage != "normal"], names(shp)) ]

#---------------------------------------------------------------------------------------
# PCA

# svd
sv.dat <- sweep(t(anova_dat), 2, colMeans(t(anova_dat)), "-") # centering
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
PCload.func(cols, shapes, U, D, 3, 2, pca.var, title = "PC2 vs PC3") # PC loadings (PC2 vs PC3)
legend('bottomleft', pch = shp, col = pal,
      c("normal",  "new_dx", "nmCSPC", "mCSPC", "nmCRPC", "mCRPC") )
dev.off()

#---------------------------------------------------------------------------------------
# t-SNE

# how to specify parameter
# refer: https://lvdmaaten.github.io/tsne/User_guide.pdf
tsnedat <- unname(t(anova_dat)) # dim(X) = N samples by D dimensions 
initdim <- 90 
perplex <- 30 

set.seed(10)
tsne_anova <- Rtsne(tsnedat, initial_dims = initdim, perplexity = perplex,
                  theta = 0, check_duplicates = F, max_iter = 50000L, eta = 50)

# t-SNE plot
plot(tsne_anova$Y, ylab = "", xlab = "", col = cols, pch = shapes, main = "t-SNE plot")
legend('topright', pch = shp, col = pal,
       c("normal",  "new_dx", "nmCSPC", "mCSPC", "nmCRPC", "mCRPC") )
dev.off()

#---------------------------------------------------------------------------------------
# Heatmap

# heatmap color scheme
# cls <- colorRampPalette(c("navy", "white","firebrick3"))(n = 2048)
cls <- colorRampPalette(c("navy", "white", "firebrick1", "firebrick3"))(n = 1024)

# heatmap_dat <- t(apply(anova_dat, 1, function(x){(x-mean(x))/sd(x)}))

heatmap3(anova_dat, 
         col = cls, # specify colors 
         ColSideColors = cols, # specify patient color code
         labCol = stage, # specify patient
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
# should we also exclude some rows and/or columns with huge outliers?


anova_dat <- raw_data_median[ (raw_data_median$BH_FDR <= 0.05) , 16:ncol(raw_data_median)]

# get 5-num summary
fivenum_summary <- t( apply( anova_dat, 1, function(x){unname(summary(as.numeric(x)))} ) )
colnames(fivenum_summary) <- c("min", "first_q", "median", "mean", "third_q", "max")
row.names(fivenum_summary) <- raw_data_median$PROBE_ID[(raw_data_median$BH_FDR <= 0.05)]

# get 95th percentile and 00th percentile
fivenum_summary <- as.data.frame(fivenum_summary)
fivenum_summary$percentile95 <- apply(anova_dat, 1, function(x){quantile(x, probs = .95, names = F)})
fivenum_summary$percentile99 <- apply(anova_dat, 1, function(x){quantile(x, probs = .99, names = F)})

hist(fivenum_summary$max, main = "Histogram of maximum log2(fluorescence) for peptides at BH FDR 5%")
length(fivenum_summary$max[fivenum_summary$max >= 13])
hist(fivenum_summary$percentile95, main = "Histogram of 95th percentile log2(fluorescence) for peptides at BH FDR 5%")
hist(fivenum_summary$percentile99, main = "Histogram of 99th percentile log2(fluorescence) for peptides at BH FDR 5%")

# hist(fivenum_summary$third_q, main = "Histogram of 3rd quartile log2(fluorescence) for peptides at BH FDR 5%")
# hist(fivenum_summary$min, main = "Histogram of minimum log2(fluorescence) for peptides at BH FDR 5%")
# hist(fivenum_summary$first_q, main = "Histogram of 1st quartile log2(fluorescence) for peptides at BH FDR 5%")

fivenum_summary %>% 
  filter((max - percentile99) <=  (third_q - first_q)) %>%
  filter(max <= 14) %>%
  tally()

# now plot heatmap
heatmap_dat <- anova_dat[
  (fivenum_summary$max <= 13)
  & ((fivenum_summary$max - fivenum_summary$percentile99) <=
     (fivenum_summary$third_q - fivenum_summary$first_q))
  , 
  # (stage != "normal") & (stage != "nmCRPC")
] 
heatmap_dat <- sweep(heatmap_dat, 1, rowMeans(heatmap_dat), "-") # centering by row

cls <- colorRampPalette(c("navy", "blue", "white", "orchid", "firebrick3"))(n=1024)

heatmap3(heatmap_dat, 
         col = cls, # specify colors 
         cexCol = 1.1, # font size of column labels,
         cexRow = 1.1, # font size of row labels,
         ColSideColors = cols,
         ColSideLabs = "stages", 
         labRow = "",
         xlab = "Patients",
         labCol = stage,
         # legendfun=function() showLegend(col = c("navy", "cornflowerblue", "turquoise1", "orchid1", 
         #                                         "darkorange1", "firebrick1"),
         #                                 legend = c("normal",  "new_dx", "nmCSPC", "mCSPC", "nmCRPC", "mCRPC"),
         #                                 cex = 1.2,
         #                                 lwd = 5  )
)

####################################################################################### 
#                                 Tukey HSD Pairwise Comparison                       #
#######################################################################################

median_subset <- raw_data_median[ (raw_data_median$BH_FDR <= 0.05) , ]

tukey_adjusted_pval <- matrix(NA, nrow = nrow(median_subset), ncol = choose(nlevels(stage),2))
tukey_effect_estimate <- matrix(NA, nrow = nrow(median_subset), ncol = choose(nlevels(stage),2))
tukey_CI_lower <- matrix(NA, nrow = nrow(median_subset), ncol = choose(nlevels(stage),2))
tukey_CI_upper <- matrix(NA, nrow = nrow(median_subset), ncol = choose(nlevels(stage),2))
anova_mse <- rep( NA, nrow(median_subset) ) # just in case we need it 

for (i in 1: nrow(median_subset)){
  fit1 <- aov( as.numeric(median_subset[i, (ncol(median_subset)-n+1) : ncol(median_subset)]) ~ stage )
  anova_mse[i] <- unname( unlist(summary(fit1))['Mean Sq2'] )
  fit1_TukeyHSD <- TukeyHSD(fit1)$stage
  tukey_adjusted_pval[i,] <- unname(fit1_TukeyHSD[,'p adj'])
  tukey_effect_estimate[i,] <- unname(fit1_TukeyHSD[,'diff'])
  tukey_CI_lower[i,] <- unname(fit1_TukeyHSD[,'lwr'])
  tukey_CI_upper[i,] <- unname(fit1_TukeyHSD[,'upr'])
  print(i)
}
colnames(tukey_adjusted_pval) = colnames(tukey_effect_estimate) = colnames(tukey_CI_lower) = colnames(tukey_CI_upper) <- names(fit1_TukeyHSD[,'p adj'])
row.names(tukey_adjusted_pval) = row.names(tukey_effect_estimate) = row.names(tukey_CI_lower) = row.names(tukey_CI_upper) <- median_subset$PROBE_ID
names(anova_mse) <- median_subset$PROBE_ID

tukey_cutoff <- 0.05

# column (pairwise) summary
apply(tukey_adjusted_pval,2, function(x){length(x[x <= tukey_cutoff])})

# pairwise comparison pattern
# 1: at most cutoff, 0: more than cutoff
tukey_pairwise_pattern <- as.data.frame(sign( -1 * sign(tukey_adjusted_pval - tukey_cutoff) + 1 )) 

# 1: positive difference; 0: exactly same; -1: negative difference
tukey_pairwise_diff <- as.data.frame(sign( tukey_effect_estimate ) ) 

# element-wise multiplication: 1: both significant and positive
significant_positive <- tukey_pairwise_pattern * tukey_pairwise_diff

#----------------------------------------------------------------------------------------
# extract some funny pattern

which(
  tukey_pairwise_pattern$"new_dx-normal" == 1 &
    tukey_pairwise_pattern$"nmCSPC-normal" == 1 &
    tukey_pairwise_pattern$"nmCRPC-normal" == 1 &
    tukey_pairwise_pattern$"mCRPC-normal"  == 1 &
    tukey_pairwise_pattern$"nmCSPC-new_dx" == 0 &
    tukey_pairwise_pattern$"nmCRPC-new_dx" == 0 &
    tukey_pairwise_pattern$"mCRPC-new_dx"  == 0 &
    tukey_pairwise_pattern$"nmCRPC-nmCSPC" == 0 &
    tukey_pairwise_pattern$"mCRPC-nmCSPC"  == 0 &
    tukey_pairwise_pattern$"mCRPC-nmCRPC"  == 0 &
    median_subset$rank_variance <= 3
  )

which(
  significant_positive$"new_dx-normal" == 1 &
    significant_positive$"nmCSPC-normal" == 1 &
    significant_positive$"nmCRPC-normal" == 1 &
    significant_positive$"mCRPC-normal"  == 1
)

which(
  ( significant_positive$"new_dx-normal" == 1 &
    significant_positive$"nmCSPC-normal" == 1 &
    significant_positive$"nmCRPC-normal" == 1 &
    significant_positive$"mCRPC-normal"  == 1 ) &
   ( significant_positive$"nmCSPC-new_dx" != 1 &
      significant_positive$"nmCRPC-new_dx" != 1 &
      significant_positive$"mCRPC-new_dx"  != 1 &
      significant_positive$"nmCRPC-nmCSPC" != 1 &
      significant_positive$"mCRPC-nmCSPC"  != 1 &
      significant_positive$"mCRPC-nmCRPC"  != 1 ) 
)

which(
  significant_positive$`mCRPC-normal` == 1 &
    significant_positive$`mCRPC-new_dx`== 1 &
    significant_positive$`mCRPC-nmCSPC` == 1 &
    significant_positive$`mCRPC-nmCRPC` == 1
)


#----------------------------------------------------------------------------------------
# write signif peptides and signif proteins into multiple worksheets of an Excel file

wb = createWorkbook()

# All_Peptides sheet
sheet = createSheet(wb, "All_Peptides")

# all peptides
all_peptides_df <- as.data.frame( bind_cols(
  select(raw_data, -contains("dat"))[one_way_anova_BH <= 0.05,],
  rank_variance = rank_variance[one_way_anova_BH <= 0.05],
  p_values = one_way_anova_pval[one_way_anova_BH <= 0.05],
  BH_FDR = one_way_anova_BH[one_way_anova_BH <= 0.05],
  q_values = one_way_anova_qval[one_way_anova_BH <= 0.05]
) )
head(all_peptides_df) # check
addDataFrame(all_peptides_df, sheet=sheet, row.names=FALSE)

all_peptides_unique_protein_df <- data.frame(
  unique_protein = unique(all_peptides_df$SEQ_ID)
)
addDataFrame(all_peptides_unique_protein_df, sheet=sheet, startColumn = ncol(all_peptides_df) + 2, row.names=FALSE)


# Tukey HSD sheet

sheet = createSheet(wb, "Tukey_HSD")
addDataFrame(as.data.frame(tukey_adjusted_pval), sheet=sheet)


# Top_5 Sheet
sheet = createSheet(wb, "Top_5")

# marginal filtering top 5 
topnum <- 5
top <- marginal_var_filter.func(one_way_anova_pval, topnum)
top5_df <- as.data.frame( bind_cols(
  select(raw_data, -contains("dat"))[rank_variance <= topnum,][ top$BH_filter <= 0.05 ,],
  rank_variance = rank_variance[rank_variance <= topnum][ top$BH_filter <= 0.05],
  p_values = top$pval_filter[ top$BH_filter <= 0.05],
  BH_FDR = top$BH_filter[ top$BH_filter <= 0.05],
  q_values = top$qval_filter[ top$BH_filter <= 0.05]
) )
head(top5_df) # check
addDataFrame(top5_df, sheet=sheet, row.names=FALSE)

top5_unique_protein_df <- data.frame(
  unique_protein = unique(top5_df$SEQ_ID)
)
addDataFrame(top5_unique_protein_df, sheet=sheet, startColumn = ncol(top5_df) + 2, row.names=FALSE)


# Top_3 Sheet
sheet = createSheet(wb, "Top_3")

# marginal filtering top 3 
topnum <- 3
top <- marginal_var_filter.func(one_way_anova_pval, topnum)
top3_df <- as.data.frame( bind_cols(
  select(raw_data, -contains("dat"))[rank_variance <= topnum,][ top$BH_filter <= 0.05 ,],
  rank_variance = rank_variance[rank_variance <= topnum][ top$BH_filter <= 0.05],
  p_values = top$pval_filter[ top$BH_filter <= 0.05],
  BH_FDR = top$BH_filter[ top$BH_filter <= 0.05],
  q_values = top$qval_filter[ top$BH_filter <= 0.05]
) )
head(top3_df) # check
addDataFrame(top3_df, sheet=sheet, row.names=FALSE)

top3_unique_protein_df <- data.frame(
  unique_protein = unique(top3_df$SEQ_ID)
)
addDataFrame(top3_unique_protein_df, sheet=sheet, startColumn = ncol(top3_df) + 2, row.names=FALSE)


# Top_1 Sheet
sheet = createSheet(wb, "Top_1")

# marginal filtering top 3 
topnum <- 1
top <- marginal_var_filter.func(one_way_anova_pval, topnum)
top1_df <- as.data.frame( bind_cols(
  select(raw_data, -contains("dat"))[rank_variance <= topnum,][ top$BH_filter <= 0.05 ,],
  rank_variance = rank_variance[rank_variance <= topnum][ top$BH_filter <= 0.05],
  p_values = top$pval_filter[ top$BH_filter <= 0.05],
  BH_FDR = top$BH_filter[ top$BH_filter <= 0.05],
  q_values = top$qval_filter[ top$BH_filter <= 0.05]
) )
head(top1_df) # check
addDataFrame(top1_df, sheet=sheet, row.names=FALSE)

top1_unique_protein_df <- data.frame(
  unique_protein = unique(top1_df$SEQ_ID)
)
addDataFrame(top1_unique_protein_df, sheet=sheet, startColumn = ncol(top1_df) + 2, row.names=FALSE)


saveWorkbook(wb, "BH_FDR_5prct.xlsx")

#---------------------------------------------------------------------------------------
# Boxplots of a few peptides

boxplot_col <- unname( pal[ match(levels(stage), names(pal)) ] )

boxplot_func <- function(mat, col = boxplot_col, draw){
  par(mar = c(2, 4.5, 2.3, 1),cex = 0.84)
  boxplot( as.numeric(mat[draw,]) ~ stage, 
           col=col, horizontal=TRUE, las=1, xlab = NA, ylab = NA,
           main= paste(boxplot_peptides[draw]) )
  abline( h=2.5, col="grey" )
  abline( h=4.5, col="grey" )
  abline( h=6.5, col="grey" )
}

boxplot_peptides <- c(
  "1152_HIST1H4L_8368;9",
  "1241_C21orf33_8209;125",
  "50_HERPUD1_9709;189",
  "48_RPS6_6194;133",
  "439_TAX1BP1_8887;493",
  "1220_PKP3_11187;577"
)

boxplot_mat <- as.matrix(raw_data_median[match(boxplot_peptides, raw_data_median$PROBE_ID),16:ncol(raw_data_median)])  
boxplot_func(mat = boxplot_mat, draw = 1)
boxplot_func(mat = boxplot_mat, draw = 2)
boxplot_func(mat = boxplot_mat, draw = 3)
boxplot_func(mat = boxplot_mat, draw = 4)
boxplot_func(mat = boxplot_mat, draw = 5)
boxplot_func(mat = boxplot_mat, draw = 6)

# save results
# save(one_way_anova_pval,
#      rank_variance,
#      anova_dat,
#      anova_mse,
#      tukey_adjusted_pval,
#      tukey_effect_estimate,
#      tukey_CI_lower,
#      tukey_CI_upper,
#      file = "remove_rep1.RData")

