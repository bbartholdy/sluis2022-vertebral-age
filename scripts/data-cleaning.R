# Stats for Iris

library(tidyverse)
library(here)

# Upload data -------------------------------------------------------------

raw_data <- readr::read_csv(here("data/raw_data/raw-data.csv"))

# Data cleaning -----------------------------------------------------------

raw_data_clean <- raw_data
raw_data_clean[,5:ncol(raw_data_clean)] <- apply(raw_data_clean[,5:ncol(raw_data_clean)], 
                                                 2, function(x) ifelse(x == "M" | x == "9", NA, x))

raw_data_clean <- raw_data_clean %>%
  mutate(across(C2_1:L5_3, as.numeric)) %>%
  rename(sex = "Sex",
         known_age = "Archival age",
         vertebral_position = "Side vertebral body"
         )

write_csv(raw_data_clean, here("data/derived_data/raw-data_clean.csv"))

# convert data to long form (easier to calculate method- and region-specific means, etc.)

raw_data_long <- raw_data_clean %>%
  pivot_longer(cols = C2_1:L5_3, values_to = "score") %>%
  mutate(method = case_when(grepl("_1", name) ~ "snodgrass", # create separate column for methods
                            grepl("_2", name) ~ "watanabe",
                            grepl("_3", name) ~ "praneatpolgrang")) %>%
  mutate(region = case_when(grepl("C", name) ~ "cervical", # create separate column for vertebral region
                            grepl("T", name) ~ "thoracic",
                            grepl("L", name) ~ "lumbar")) %>%
  mutate(name = str_replace(name, "_.", "")) %>% # convert to contain element name
  rename(element = name,
         position = vertebral_position)

write_csv(raw_data_long, here("data/derived_data/raw-data_long.csv"))
