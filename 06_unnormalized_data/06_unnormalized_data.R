library(tidyverse) # make sure you have the latest tidyverse !
library(ggplot2)

#-------------------------------------------------------------------------------------------
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

unnormalized_raw_median <- read_csv("unnormalized_raw_median.csv")

#-------------------------------------------------------------------------------------------
# get unnormalized data

unnormalized_raw <- read_csv("non_normalized_data.csv")

# check
sum(as.numeric(unnormalized_raw$X == unnormalized_raw$COL_NUM)) == nrow(unnormalized_raw) # X = COL_NUM
sum(as.numeric(unnormalized_raw$Y == unnormalized_raw$ROW_NUM)) == nrow(unnormalized_raw) # Y = ROW_NUM
sum(as.numeric(unnormalized_raw$MATCH_INDEX == unnormalized_raw$FEATURE_ID)) == nrow(unnormalized_raw) # MATCH_INDEX = FEATURE_ID
unique(unnormalized_raw$DESIGN_NOTE) # only NA
unique(unnormalized_raw$SELECTION_CRITERIA) # only NA
unique(unnormalized_raw$MISMATCH) # only 0
unique(unnormalized_raw$PROBE_CLASS) # only NA
# we can drop X, Y, MATCH_INDEX, DESIGN_NOTE, SELECTION_CRITERIA, MISMATCH, PROBE_CLASS

unnormalized_raw = unnormalized_raw %>%
  select(PROBE_DESIGN_ID:Y, any_of(array_id_key$array_id)) %>% # drop patients with records at different stages 
  select( -c(X, Y, MATCH_INDEX, DESIGN_NOTE, SELECTION_CRITERIA, MISMATCH, PROBE_CLASS) ) # drop unhelpful columns

colnames(unnormalized_raw)[1:15] # check

# take log2 transformation
unnormalized_raw <- unnormalized_raw %>%
  mutate_at(vars(matches("dat")), log2)

# make sure array_id_key align with raw_data_complete
array_iii <- match( colnames( unnormalized_raw %>% select(any_of(array_id_key$array_id)) ), array_id_key$array_id )
array_id_key <- array_id_key[array_iii , ]
sum(as.numeric( colnames( unnormalized_raw %>% select(any_of(array_id_key$array_id)) ) == 
                  array_id_key$array_id )) == nrow(array_id_key)

# compute median
unnormalized_raw_median <- t( apply( select(unnormalized_raw, contains("dat")), 1, function(x) {
  tapply(x, array_id_key$id, FUN=median)
} ) )
unnormalized_raw_median <- bind_cols( select(unnormalized_raw, -contains("dat")), 
                                      as.data.frame(unnormalized_raw_median) ) 

colnames(unnormalized_raw_median)[1:15] # check

write.table(unnormalized_raw_median, file = "unnormalized_raw_median.csv", sep = ",", row.names = F)

#-------------------------------------------------------------------------------------------

unnormalized_median_long <- unnormalized_raw_median %>%
  select(any_of(array_id_key$id)) %>%
  pivot_longer(cols = everything(), names_to = "id", values_to = "fluorescence") 

# check
nrow(unnormalized_median_long) == nrow(unnormalized_raw_median) * nrow(patient_key)
head(unnormalized_median_long)

# set fill color
unnormalized_median_long$stage <- patient_key$stage[ match(unnormalized_median_long$id, patient_key$id) ]

# sort according to group and median
unnormalized_sort <- unnormalized_median_long %>%
  group_by(stage,id) %>%
  mutate(median = median(fluorescence)) %>%
  select(-fluorescence) %>%
  distinct() %>%
  arrange(stage, median) # arrange according to stage and median

# sort order of patients in boxplot
# unnormalized_sort <- unnormalized_sort %>% arrange(median) # arrange only according to median
unnormalized_median_long$id <- factor(unnormalized_median_long$id, levels = unnormalized_sort$id)

pal <- c("navy", "cornflowerblue", "turquoise1", "orchid1", "darkorange1", "firebrick1")
names(pal) <- c("normal",  "new_dx", "nmCSPC", "mCSPC", "nmCRPC", "mCRPC")

ggplot(unnormalized_median_long, aes(x = id, y = fluorescence, fill = stage)) +
  geom_boxplot(outlier.shape = ".") +
  # scale_fill_discrete(name="Stage") +
  # guides(fill=guide_legend(title="Stage")) +
  scale_fill_manual(name = "Stage", values = pal) +
  labs(title = "Boxplots of Peptide (Pre-normalized) Fluorescence Levels for All Patients", 
       x = "Patient ID", y = "Median Fluorescence Levels on log2 scale") +
  theme(panel.background = element_rect(fill = "grey90"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 4.5),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))
