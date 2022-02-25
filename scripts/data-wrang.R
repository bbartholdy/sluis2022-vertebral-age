library(tidyverse)
library(here)
library(corrplot)
library(broom)
#library(infer)
#library(rstatix)

# Upload data -------------------------------------------------------------

raw_data <- readr::read_csv(here("data/derived_data/raw-data_clean.csv"))
raw_data_long <- readr::read_csv(here("data/derived_data/raw-data_long.csv"))

raw_data <- raw_data %>%
  mutate(sex = as_factor(sex),
         vertebral_position = as_factor(vertebral_position)) %>%
  filter(id != "MB14") # individual identified as an extreme outlier

raw_data_long <- raw_data_long %>%
  mutate(across(c(sex, position, element, method, region), as_factor)) %>%
  filter(id != "MB14") # individual identified as an extreme outlier


# Data wrangling ----------------------------------------------------------

# combine superior and inferior margin scores using mean
element_data_long <- raw_data_long %>% 
  group_by(id, method, element) %>%
  mutate(score = mean(score, na.rm = T)) %>%
  filter(position == "Superior") %>%
  dplyr::select(!position)

# EDA ---------------------------------------------------------------------

region_means <- raw_data_long %>%
  group_by(id, known_age, sex, region, method) %>%
  summarise(region_mean = mean(score, na.rm = T)) 

region_means2 <- raw_data_long %>%
  group_by(id, known_age, sex, method) %>%
  summarise(region_mean = mean(score, na.rm = T)) %>%
  mutate(region = "all")

region_means_all <- full_join(region_means, region_means2)

region_means_all <- region_means_all %>%
  ungroup() %>%
  mutate(age_group = case_when(known_age >= 0 & known_age <= 19 ~ 1,
            known_age >= 20 & known_age <= 29 ~ 2,
            known_age >= 30 & known_age <= 39 ~ 3,
            known_age >= 40 & known_age <= 49 ~ 4,
            known_age >= 50 & known_age <= 59 ~ 5,
            known_age >= 60 & known_age <= 69 ~ 6,
            known_age >= 70 & known_age <= 79 ~ 7,
            known_age >= 80 & known_age <= 89 ~ 8)) %>%
  mutate(age_group = as_factor(age_group))

# correlation between osteophyte score and age

  # scatter plots grouped by region means
region_means %>%
  #group_by(feat, known_age, sex, region, method) %>%
  #summarise(region_mean = mean(score, na.rm = T)) %>%
  ggplot(aes(x = known_age, region_mean, col = sex)) +
  geom_point() +
  scale_color_viridis_d() +
  facet_wrap(~ region + method) +
  theme_bw()

  # by element
cor_matrix <- raw_data %>%
  filter(vertebral_position == "Superior") %>%
  dplyr::select(!c(id, sex, vertebral_position)) %>%
  cor(use = "pairwise.complete.obs")

