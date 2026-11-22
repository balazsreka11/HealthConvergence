#############################################################
# Health convergence: Hungary vs. comparator country
# Kitagawa decomposition + parametric bootstrap (Poisson)
#############################################################

# install.packages(c("eurostat", "dplyr", "tidyr"))
library(eurostat)
library(dplyr)
library(tidyr)

#------------------------------------------------------------
# 0) SETTINGS
#------------------------------------------------------------
investigated <- "HU"   # investigated country
comparator   <- "DK"   # comparator country (change it freely, e.g. "AT")
base_year    <- 2011
ref_year     <- 2019

B     <- 10000         # number of bootstrap replications
alpha <- 0.05          # 95% confidence level
set.seed(2026)         # reproducibility


#------------------------------------------------------------
# 1) DATA DOWNLOAD AND CLEANING (same as before)
#------------------------------------------------------------
geo_dic <- get_eurostat_dic("geo")
country_codes <- geo_dic$code_name[
  grepl("^[A-Z]{2}$", geo_dic$code_name)
]

health_data <- get_eurostat(
  id = "hlth_cd_asdr2",
  filters = list(
    geo = country_codes,
    time = c(base_year, ref_year),
    sex = "T",
    age = "TOTAL",
    unit = "RT"          # age-standardised death rate
  ),
  time_format = "num"
)

# keep only countries, drop those with lots of NA
health_data <- health_data %>%
  filter(grepl("^[A-Z]{2}$", geo)) %>%
  filter(!geo %in% c("UK", "MD", "GE"))

# check the available ICD-10 codes
sort(unique(health_data$icd10))

# the 13 Level 1 causes of death
level1_causes <- c(
  "A_B", "C00-D48", "D50-D89", "E", "F",
  "G_H", "I", "J", "K", "L", "M", "N", "V01-Y89")

health_data_level1 <- health_data %>%
  filter(icd10 %in% level1_causes) %>%
  select(geo, time, icd10, values)

# all-cause rate: Eurostat changed the code (TOTAL in 2011, A-R_V-Y in 2019)
all_causes <- health_data %>%
  filter(icd10 %in% c("TOTAL", "A-R_V-Y")) %>%
  filter(!is.na(values)) %>%
  group_by(geo, time) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(icd10 = "ALL_CAUSES") %>%
  select(geo, time, icd10, values)

health_final <- rbind(health_data_level1, all_causes)

#------------------------------------------------------------
# 2) WIDE TABLE: one row per cause, 4 rates
#    inv_2011, inv_2019, comp_2011, comp_2019
#------------------------------------------------------------
rates <- health_final %>%
  filter(geo %in% c(investigated, comparator),
         time %in% c(base_year, ref_year)) %>%
  mutate(group = ifelse(geo == investigated, "inv", "comp"),
         col   = paste0(group, "_", time)) %>%
  select(icd10, col, values) %>%
  pivot_wider(names_from = col, values_from = values) %>%
  select(icd10,
         inv_0  = all_of(paste0("inv_",  base_year)),
         inv_1  = all_of(paste0("inv_",  ref_year)),
         comp_0 = all_of(paste0("comp_", base_year)),
         comp_1 = all_of(paste0("comp_", ref_year)))


#------------------------------------------------------------
# 3) CONVERGENCE MEASURE AND KITAGAWA DECOMPOSITION
#    (point estimates)
#------------------------------------------------------------
# V and k for all rows (including ALL_CAUSES)
convergence <- rates %>%
  mutate(V_0 = inv_0 / comp_0,
         V_1 = inv_1 / comp_1,
         k   = V_1 - V_0)

# Kitagawa only on the 13 causes (without the total)
kitagawa_results <- rates %>%
  filter(icd10 != "ALL_CAUSES") %>%
  mutate(
    w_0 = comp_0 / sum(comp_0),      # comparator cause shares, base year
    w_1 = comp_1 / sum(comp_1),      # comparator cause shares, reference year
    V_0 = inv_0 / comp_0,
    V_1 = inv_1 / comp_1,
    k   = V_1 - V_0,
    Kp_0  = w_0 * k,
    Kp_1  = w_1 * k,
    Kp    = (Kp_0 + Kp_1) / 2,
    Kpp_0 = w_1 * V_0 - w_0 * V_0,
    Kpp_1 = w_1 * V_1 - w_0 * V_1,
    Kpp   = (Kpp_0 + Kpp_1) / 2
  )

kitagawa_summary <- kitagawa_results %>%
  summarise(sum_Kp = sum(Kp),
            composition_effect = sum(Kpp))

# check: sum(Kp) + sum(Kpp) should be close to the total change in V
total_k <- convergence %>% filter(icd10 == "ALL_CAUSES") %>% pull(k)
kitagawa_summary$sum_Kp + kitagawa_summary$composition_effect
total_k

#------------------------------------------------------------
# 4) PARAMETRIC BOOTSTRAP (logic of compstat_09.R)
#
# Assumption: the number of deaths follows a Poisson distribution
# whose parameter is the observed rate:
#   Investigated_0 ~ Poi(inv_0),  Investigated_1 ~ Poi(inv_1)
#   Comparator_0   ~ Poi(comp_0), Comparator_1   ~ Poi(comp_1)
# We simulate B new "samples" for each cause and each of the 4 rates,
# then take the 2.5% and 97.5% quantiles (percentile CI).
#------------------------------------------------------------
rate_cols <- c("inv_0", "inv_1", "comp_0", "comp_1")

# list: one B x 4 matrix of simulated rates for each cause
bs_sim <- lapply(seq_len(nrow(rates)), function(i) {
  sapply(rate_cols, function(col) {
    rpois(n = B, lambda = rates[[col]][i])
  })
})
names(bs_sim) <- rates$icd10

# 4a) 95% percentile CIs of the mortality rates
ci_rates <- do.call(rbind, lapply(rates$icd10, function(cause) {
  sims <- bs_sim[[cause]]
  data.frame(
    icd10 = cause,
    series = rate_cols,
    observed = unlist(rates[rates$icd10 == cause, rate_cols]),
    lower = apply(sims, 2, quantile, p = alpha / 2),
    upper = apply(sims, 2, quantile, p = 1 - alpha / 2),
    row.names = NULL
  )
}))

# wide format: one row per cause, lower/upper for all 4 series
ci_wide <- ci_rates %>%
  pivot_wider(names_from = series,
              values_from = c(observed, lower, upper),
              names_glue = "{series}_{.value}")

#------------------------------------------------------------
# 5) SIGNIFICANT DIFFERENCES (investigated vs comparator)
#    A difference is significant if the two 95% CIs do not overlap
#    in the given year.
#------------------------------------------------------------
significance <- ci_wide %>%
  mutate(
    sig_base = inv_0_upper < comp_0_lower | inv_0_lower > comp_0_upper,
    sig_ref  = inv_1_upper < comp_1_lower | inv_1_lower > comp_1_upper,
    dir_base = case_when(!sig_base ~ "n.s.",
                         inv_0_observed > comp_0_observed ~ "investigated higher",
                         TRUE ~ "investigated lower"),
    dir_ref  = case_when(!sig_ref ~ "n.s.",
                         inv_1_observed > comp_1_observed ~ "investigated higher",
                         TRUE ~ "investigated lower")
  ) %>%
  select(icd10, sig_base, dir_base, sig_ref, dir_ref)

print(significance)

#------------------------------------------------------------
# 5b) Multiple comparisons (compstat_09.R, Bonferroni idea):
#     14 causes are tested in each year, so the family-wise error
#     rate grows. Bonferroni-adjusted level: alpha / g.
#------------------------------------------------------------
g <- nrow(rates)
alpha_adj <- alpha / g

ci_adj <- do.call(rbind, lapply(rates$icd10, function(cause) {
  sims <- bs_sim[[cause]]
  data.frame(
    icd10 = cause,
    series = rate_cols,
    lower_adj = apply(sims, 2, quantile, p = alpha_adj / 2),
    upper_adj = apply(sims, 2, quantile, p = 1 - alpha_adj / 2),
    row.names = NULL
  )
})) %>%
  pivot_wider(names_from = series,
              values_from = c(lower_adj, upper_adj),
              names_glue = "{series}_{.value}")

significance_bonf <- ci_adj %>%
  mutate(
    sig_base_bonf = inv_0_upper_adj < comp_0_lower_adj | inv_0_lower_adj > comp_0_upper_adj,
    sig_ref_bonf  = inv_1_upper_adj < comp_1_lower_adj | inv_1_lower_adj > comp_1_upper_adj
  ) %>%
  select(icd10, sig_base_bonf, sig_ref_bonf)

significance_all <- significance %>% left_join(significance_bonf, by = "icd10")
print(significance_all)

#------------------------------------------------------------
# 6) BOOTSTRAP CIs OF V, k, Kp, Kpp
#    The same simulated rates are pushed through the formulas, so
#    each replication gives one full decomposition.
#------------------------------------------------------------
causes13 <- rates$icd10[rates$icd10 != "ALL_CAUSES"]

# B x 13 matrices of simulated rates
get_mat <- function(col) sapply(causes13, function(cs) bs_sim[[cs]][, col])
inv0 <- get_mat("inv_0");  inv1 <- get_mat("inv_1")
comp0 <- get_mat("comp_0"); comp1 <- get_mat("comp_1")

# weights (row-wise shares of the comparator's total)
w0 <- comp0 / rowSums(comp0)
w1 <- comp1 / rowSums(comp1)

V0 <- inv0 / comp0
V1 <- inv1 / comp1
k_sim <- V1 - V0

Kp_sim  <- (w0 * k_sim + w1 * k_sim) / 2
Kpp_sim <- ((w1 * V0 - w0 * V0) + (w1 * V1 - w0 * V1)) / 2

# replications with 0 deaths give Inf/NaN -> they are ignored (na.rm)
Kp_sim[!is.finite(Kp_sim)]   <- NA
Kpp_sim[!is.finite(Kpp_sim)] <- NA
k_sim[!is.finite(k_sim)]     <- NA

q_ci <- function(m) {
  t(apply(m, 2, quantile, p = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE))
}

ci_kitagawa <- data.frame(
  icd10 = causes13,
  k_obs = kitagawa_results$k,   k_low = q_ci(k_sim)[, 1],   k_up = q_ci(k_sim)[, 2],
  Kp_obs = kitagawa_results$Kp, Kp_low = q_ci(Kp_sim)[, 1], Kp_up = q_ci(Kp_sim)[, 2],
  Kpp_obs = kitagawa_results$Kpp, Kpp_low = q_ci(Kpp_sim)[, 1], Kpp_up = q_ci(Kpp_sim)[, 2]
)

# CIs of the totals
sum_Kp_sim  <- rowSums(Kp_sim)
sum_Kpp_sim <- rowSums(Kpp_sim)
ci_totals <- data.frame(
  what = c("sum_Kp", "composition_effect"),
  observed = c(kitagawa_summary$sum_Kp, kitagawa_summary$composition_effect),
  lower = c(quantile(sum_Kp_sim,  alpha / 2, na.rm = TRUE),
            quantile(sum_Kpp_sim, alpha / 2, na.rm = TRUE)),
  upper = c(quantile(sum_Kp_sim,  1 - alpha / 2, na.rm = TRUE),
            quantile(sum_Kpp_sim, 1 - alpha / 2, na.rm = TRUE))
)

print(ci_kitagawa)
print(ci_totals)

#------------------------------------------------------------
# 7) QUICK VISUAL CHECK of one simulated distribution
#------------------------------------------------------------
hist(bs_sim[["I"]][, "inv_1"],
     main = paste("Simulated rate, cause I,", investigated, ref_year),
     xlab = "rate per 100 000")
abline(v = quantile(bs_sim[["I"]][, "inv_1"], c(alpha / 2, 1 - alpha / 2)),
       col = "red")
