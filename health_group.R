#################################
# install.packages("eurostat")
library(eurostat)
library(dplyr)
library(tidyr)

# set the countries codes
geo_dic <- get_eurostat_dic("geo")
country_codes <- geo_dic$code_name[
  grepl("^[A-Z]{2}$", geo_dic$code_name)
]

# load the datasets from eurostat
health_data <- get_eurostat(
  id = "hlth_cd_asdr2",
  filters = list(
    geo = country_codes,
    time = c(2011, 2019),
    sex = "T",
    age = "TOTAL",
    unit = "RT" # stanadrized death
  ),
  time_format = "num"
)

# filter countries
health_data <- health_data %>%
  filter(grepl("^[A-Z]{2}$", geo))

# lot of NA --> skip
health_data <- health_data %>%
  filter(!geo %in% c("UK", "MD", "GE"))

# search for level 1 death causes
sort(unique(health_data$icd10))

# choose the categories
level1_causes <- c(
  "A_B", "C00-D48", "D50-D89", "E", "F",
  "G_H", "I", "J", "K", "L", "M", "N", "V01-Y89")

# filter for the categories
health_data_level1 <- health_data %>%
  filter(icd10 %in% level1_causes)

# TOTAL and ARVY problem: TOTAL(2011); ARVY(2019) --> eurostat changed the code
all_causes <- health_data %>%
  filter(icd10 %in% c("TOTAL", "A-R_V-Y")) %>%
  filter(!is.na(values)) %>%
  group_by(geo, time) %>%
  slice(1) %>% ungroup() %>%
  mutate(icd10 = "ALL_CAUSES") %>%
  select(geo, time, icd10, values)

# same order of the columns to merge the datasets
all_causes <- all_causes %>%
  select(geo, time, icd10, values)

health_data_level1 <- health_data_level1 %>%
  select(geo, time, icd10, values)

# merge datasets: 
health_final <- rbind(health_data_level1, all_causes)

########------------------------------
comparison_data <- health_final %>%
  filter(
    geo %in% c("HU", "DK"))

wide_data <- comparison_data %>%
  pivot_wider(
    names_from = geo,
    values_from = values
  )

# V
ratio_data <- wide_data %>%
  mutate(
    V = HU / DK
  ) %>%
  select(time, icd10, HU, DK, V)

# k = V1 - V0
convergence <- ratio_data %>%
  select(time, icd10, V) %>%
  pivot_wider(
    names_from = time,
    values_from = V,
    names_prefix = "V_"
  ) %>%
  mutate(
    k = V_2019 - V_2011
  )

# cause (NOT TOTAL)
cause_data <- wide_data %>%
  filter(icd10 != "ALL_CAUSES")

# kitagawa weights
kitagawa_data <- cause_data %>%
  group_by(time) %>%
  mutate(
    comparator_total = sum(DK, na.rm = TRUE),
    weight = DK / comparator_total
  ) %>%
  ungroup() %>%
  select(time, icd10, HU, DK, weight) %>%
  pivot_wider(
    names_from = time,
    values_from = c(HU, DK, weight),
    names_sep = "_"
  )

# Kitagawa results
kitagawa_results <- kitagawa_data %>%
  mutate(
    V_2011 = HU_2011 / DK_2011,
    V_2019 = HU_2019 / DK_2019,
    
    k = V_2019 - V_2011,
    
    Kp_2011 = weight_2011 * k,
    Kp_2019 = weight_2019 * k,
    Kp = (Kp_2011 + Kp_2019) / 2,
    
    Kpp_2011 =
      (weight_2019 * V_2011) -
      (weight_2011 * V_2011),
    
    Kpp_2019 =
      (weight_2019 * V_2019) -
      (weight_2011 * V_2019),
    
    Kpp = (Kpp_2011 + Kpp_2019) / 2
  )

# SUMMARY
kitagawa_summary <- kitagawa_results %>%
  summarise(
    sum_Kp = sum(Kp, na.rm = TRUE),
    composition_effect = sum(Kpp, na.rm = TRUE)
  )

total_k <- convergence %>% filter(icd10 == "ALL_CAUSES") %>% pull(k)
kitagawa_summary$sum_Kp + kitagawa_summary$composition_effect 
