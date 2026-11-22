#############################################################
# Health convergence of Hungary towards a comparator country
# R Shiny app: Kitagawa decomposition, waterfall chart,
# parametric (Poisson) bootstrap CIs, significance tests
#
# Run with:  shiny::runApp("app.R")
# First start downloads the Eurostat data (table hlth_cd_asdr2)
# and caches it in "health_cache.rds" next to this file.
#############################################################

# install.packages(c("shiny", "ggplot2", "dplyr", "tidyr", "eurostat"))
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(eurostat)

#------------------------------------------------------------
# 0) GLOBAL SETTINGS
#------------------------------------------------------------
investigated <- "HU"      # investigated country (fixed)
base_year    <- 2011      # base year   (index 0)
ref_year     <- 2019      # reference year (index 1)
alpha        <- 0.05      # 95% confidence level

# the 13 Level 1 causes of death (Eurostat ICD-10 codes)
level1_causes <- c(
  "A_B", "C00-D48", "D50-D89", "E", "F",
  "G_H", "I", "J", "K", "L", "M", "N", "V01-Y89")

# readable labels (also defines the display order)
cause_labels <- c(
  "ALL_CAUSES" = "All causes",
  "A_B"        = "Infectious & parasitic diseases",
  "C00-D48"    = "Neoplasms",
  "D50-D89"    = "Blood & immune disorders",
  "E"          = "Endocrine & metabolic diseases",
  "F"          = "Mental & behavioural disorders",
  "G_H"        = "Nervous system & sense organs",
  "I"          = "Circulatory system diseases",
  "J"          = "Respiratory system diseases",
  "K"          = "Digestive system diseases",
  "L"          = "Skin diseases",
  "M"          = "Musculoskeletal diseases",
  "N"          = "Genitourinary diseases",
  "V01-Y89"    = "External causes (accidents etc.)")

#------------------------------------------------------------
# 1) DATA: download from Eurostat (cached) and clean
#------------------------------------------------------------
load_health_data <- function(cache = "health_cache.rds") {
  if (file.exists(cache)) return(readRDS(cache))

  # country codes (two-letter codes only, no EU aggregates)
  geo_dic <- get_eurostat_dic("geo")
  country_codes <- geo_dic$code_name[grepl("^[A-Z]{2}$", geo_dic$code_name)]

  raw <- get_eurostat(
    id = "hlth_cd_asdr2",
    filters = list(
      geo  = country_codes,
      time = c(base_year, ref_year),
      sex  = "T",
      age  = "TOTAL",
      unit = "RT"                 # age-standardised death rate / 100 000
    ),
    time_format = "num"
  )

  raw <- raw %>%
    filter(grepl("^[A-Z]{2}$", geo)) %>%
    filter(!geo %in% c("UK", "MD", "GE"))     # too many NAs

  # Level 1 causes
  level1 <- raw %>%
    filter(icd10 %in% level1_causes) %>%
    select(geo, time, icd10, values)

  # all-cause rate: Eurostat changed the code (TOTAL in 2011, A-R_V-Y in 2019)
  all_causes <- raw %>%
    filter(icd10 %in% c("TOTAL", "A-R_V-Y"), !is.na(values)) %>%
    group_by(geo, time) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(icd10 = "ALL_CAUSES") %>%
    select(geo, time, icd10, values)

  health_final <- rbind(level1, all_causes)

  # country names for the drop-down menu
  geo_names <- setNames(geo_dic$full_name, geo_dic$code_name)

  out <- list(health_final = health_final, geo_names = geo_names)
  saveRDS(out, cache)
  out
}

dat <- load_health_data()
health_final <- dat$health_final
geo_names    <- dat$geo_names

# countries with complete data: 14 causes x 2 years, no NA
complete_geo <- health_final %>%
  filter(!is.na(values)) %>%
  group_by(geo) %>%
  summarise(n = n(), n_cause = n_distinct(icd10), .groups = "drop") %>%
  filter(n == 28, n_cause == 14) %>%
  pull(geo)

comparator_choices <- setdiff(complete_geo, investigated)
comparator_choices <- setNames(
  comparator_choices,
  paste0(geo_names[comparator_choices], " (", comparator_choices, ")"))
comparator_choices <- comparator_choices[order(names(comparator_choices))]

#------------------------------------------------------------
# 2) CALCULATION FUNCTIONS
#------------------------------------------------------------

# 2a) one row per cause with the 4 rates:
#     inv_0, inv_1 (investigated, base/ref year), comp_0, comp_1 (comparator)
prep_rates <- function(hf, inv, comp) {
  hf %>%
    filter(geo %in% c(inv, comp)) %>%
    mutate(series = paste0(ifelse(geo == inv, "inv", "comp"), "_",
                           ifelse(time == base_year, 0, 1))) %>%
    select(icd10, series, values) %>%
    pivot_wider(names_from = series, values_from = values) %>%
    select(icd10, inv_0, inv_1, comp_0, comp_1) %>%
    arrange(match(icd10, names(cause_labels)))
}

# 2b) convergence measure V0, V1, k (all causes + 13 causes)
convergence_table <- function(rates) {
  rates %>%
    mutate(V_0 = inv_0 / comp_0,
           V_1 = inv_1 / comp_1,
           k   = V_1 - V_0)
}

# 2c) Kitagawa decomposition on the 13 causes
kitagawa <- function(rates) {
  detail <- rates %>%
    filter(icd10 != "ALL_CAUSES") %>%
    mutate(
      w_0 = comp_0 / sum(comp_0),        # comparator cause shares, base year
      w_1 = comp_1 / sum(comp_1),        # comparator cause shares, ref year
      V_0 = inv_0 / comp_0,
      V_1 = inv_1 / comp_1,
      k   = V_1 - V_0,
      Kp_0  = w_0 * k,
      Kp_1  = w_1 * k,
      Kp    = (Kp_0 + Kp_1) / 2,
      Kpp_0 = w_1 * V_0 - w_0 * V_0,
      Kpp_1 = w_1 * V_1 - w_0 * V_1,
      Kpp   = (Kpp_0 + Kpp_1) / 2)
  list(detail = detail,
       sum_Kp = sum(detail$Kp),
       comp_effect = sum(detail$Kpp))
}

# 2d) data for the waterfall chart:
#     ratio 2011 -> 13 cause contributions (Kp) -> composition effect
#     -> residual -> ratio 2019
make_waterfall <- function(conv, kit) {
  V0 <- conv$V_0[conv$icd10 == "ALL_CAUSES"]
  V1 <- conv$V_1[conv$icd10 == "ALL_CAUSES"]

  steps <- kit$detail %>% arrange(Kp)
  # the total all-cause ratio is not exactly the sum of the 13 causes,
  # so a small residual closes the chart honestly
  resid <- (V1 - V0) - kit$sum_Kp - kit$comp_effect

  vals   <- c(steps$Kp, kit$comp_effect, resid)
  labels <- c(unname(cause_labels[steps$icd10]),
              "Composition effect", "Residual (rounding)")
  ends   <- V0 + cumsum(vals)
  starts <- c(V0, head(ends, -1))

  df <- tibble(
    label = c(paste0("Ratio ", base_year), labels, paste0("Ratio ", ref_year)),
    start = c(0, starts, 0),
    end   = c(V0, ends, V1),
    type  = c("Total", ifelse(vals >= 0, "Divergence (worse)", "Convergence (better)"),
              "Total"),
    value = c(V0, vals, V1)
  )
  df$label <- factor(df$label, levels = df$label)
  df$id <- seq_len(nrow(df))
  df
}

#------------------------------------------------------------
# 3) PARAMETRIC BOOTSTRAP (Poisson)
#
# Death counts are assumed to be Poisson with the observed rate as
# parameter (as in the course file compstat_09.R, parametric bootstrap):
#   Inv_0 ~ Poi(inv_0), Inv_1 ~ Poi(inv_1),
#   Comp_0 ~ Poi(comp_0), Comp_1 ~ Poi(comp_1)
# For each cause B new values are simulated; the 2.5% and 97.5%
# quantiles give the 95% percentile CI.
#------------------------------------------------------------
bootstrap_all <- function(rates, B, alpha) {
  rate_cols <- c("inv_0", "inv_1", "comp_0", "comp_1")

  # list: one B x 4 matrix of simulated rates for each cause
  bs_sim <- lapply(seq_len(nrow(rates)), function(i) {
    sapply(rate_cols, function(col) rpois(n = B, lambda = rates[[col]][i]))
  })
  names(bs_sim) <- rates$icd10

  # 3a) CIs of the 4 rates for each cause (also Bonferroni-adjusted)
  g <- nrow(rates)                 # number of causes tested per year
  alpha_adj <- alpha / g
  ci_long <- do.call(rbind, lapply(rates$icd10, function(cause) {
    sims <- bs_sim[[cause]]
    data.frame(
      icd10 = cause,
      series = rate_cols,
      observed  = unlist(rates[rates$icd10 == cause, rate_cols]),
      lower     = apply(sims, 2, quantile, p = alpha / 2),
      upper     = apply(sims, 2, quantile, p = 1 - alpha / 2),
      lower_adj = apply(sims, 2, quantile, p = alpha_adj / 2),
      upper_adj = apply(sims, 2, quantile, p = 1 - alpha_adj / 2),
      row.names = NULL)
  }))

  # 3b) significance: CIs of the two countries do not overlap
  ci_w <- ci_long %>%
    pivot_wider(names_from = series,
                values_from = c(observed, lower, upper, lower_adj, upper_adj),
                names_glue = "{series}_{.value}")

  sig_one <- function(y, adj = FALSE) {
    sfx <- if (adj) "_adj" else ""
    iu <- ci_w[[paste0("inv_", y, "_upper", sfx)]]
    il <- ci_w[[paste0("inv_", y, "_lower", sfx)]]
    cu <- ci_w[[paste0("comp_", y, "_upper", sfx)]]
    cl <- ci_w[[paste0("comp_", y, "_lower", sfx)]]
    iu < cl | il > cu
  }
  dir_one <- function(y, sig) {
    case_when(!sig ~ "no",
              ci_w[[paste0("inv_", y, "_observed")]] >
                ci_w[[paste0("comp_", y, "_observed")]] ~ "yes (HU higher)",
              TRUE ~ "yes (HU lower)")
  }
  sig <- data.frame(icd10 = ci_w$icd10)
  sig$base      <- dir_one(0, sig_one(0))
  sig$ref       <- dir_one(1, sig_one(1))
  sig$base_bonf <- dir_one(0, sig_one(0, TRUE))
  sig$ref_bonf  <- dir_one(1, sig_one(1, TRUE))

  # 3c) bootstrap CIs of k, Kp, Kpp (each replication = full decomposition)
  causes13 <- rates$icd10[rates$icd10 != "ALL_CAUSES"]
  get_mat <- function(col) sapply(causes13, function(cs) bs_sim[[cs]][, col])
  inv0 <- get_mat("inv_0");   inv1 <- get_mat("inv_1")
  comp0 <- get_mat("comp_0"); comp1 <- get_mat("comp_1")

  w0 <- comp0 / rowSums(comp0)
  w1 <- comp1 / rowSums(comp1)
  V0 <- inv0 / comp0
  V1 <- inv1 / comp1
  k_sim   <- V1 - V0
  Kp_sim  <- (w0 * k_sim + w1 * k_sim) / 2
  Kpp_sim <- ((w1 * V0 - w0 * V0) + (w1 * V1 - w0 * V1)) / 2

  # replications with 0 simulated deaths give Inf/NaN -> ignored
  k_sim[!is.finite(k_sim)]     <- NA
  Kp_sim[!is.finite(Kp_sim)]   <- NA
  Kpp_sim[!is.finite(Kpp_sim)] <- NA

  q_ci <- function(m) {
    t(apply(m, 2, quantile, p = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE))
  }
  ci_kit <- data.frame(
    icd10 = causes13,
    k_low = q_ci(k_sim)[, 1],     k_up = q_ci(k_sim)[, 2],
    Kp_low = q_ci(Kp_sim)[, 1],   Kp_up = q_ci(Kp_sim)[, 2],
    Kpp_low = q_ci(Kpp_sim)[, 1], Kpp_up = q_ci(Kpp_sim)[, 2],
    row.names = NULL)

  sum_Kp_sim  <- rowSums(Kp_sim)
  sum_Kpp_sim <- rowSums(Kpp_sim)
  ci_tot <- list(
    sum_Kp = quantile(sum_Kp_sim, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE),
    comp   = quantile(sum_Kpp_sim, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE),
    n_dropped = sum(is.na(sum_Kp_sim)))

  list(ci_long = ci_long, sig = sig, ci_kit = ci_kit, ci_tot = ci_tot)
}

#------------------------------------------------------------
# 4) USER INTERFACE
#------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Health convergence of Hungary: Kitagawa decomposition"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("comparator", "Comparator country",
                  choices = comparator_choices,
                  selected = if ("DK" %in% comparator_choices) "DK"
                             else comparator_choices[1]),
      numericInput("B", "Bootstrap replications", value = 5000,
                   min = 500, max = 50000, step = 500),
      numericInput("seed", "Random seed", value = 2026),
      actionButton("run", "Re-run bootstrap"),
      hr(),
      helpText("Convergence measure: ratio of age-standardised death rates",
               "(Hungary / comparator). k = ratio in 2019 minus ratio in 2011;",
               "negative k = Hungary is catching up."),
      helpText("Bootstrap: deaths ~ Poisson(observed rate), 95% percentile CIs.")
    ),
    mainPanel(
      width = 9,
      htmlOutput("summary_text"),
      tabsetPanel(
        tabPanel("Waterfall chart", plotOutput("waterfall", height = "550px")),
        tabPanel("Decomposition table", tableOutput("kit_table")),
        tabPanel("Rates with 95% CIs",
                 plotOutput("rate_plot", height = "600px"),
                 tableOutput("rate_table")),
        tabPanel("Significant differences",
                 helpText("A difference is significant if the 95% bootstrap CIs",
                          "of the two countries do not overlap in that year.",
                          "The 'Bonferroni' columns use alpha / 14 per cause."),
                 tableOutput("sig_table"))
      )
    )
  )
)

#------------------------------------------------------------
# 5) SERVER
#------------------------------------------------------------
server <- function(input, output, session) {

  comp_name <- reactive(geo_names[[input$comparator]])

  # point estimates
  rates_r <- reactive({
    req(input$comparator)
    prep_rates(health_final, investigated, input$comparator)
  })
  conv_r <- reactive(convergence_table(rates_r()))
  kit_r  <- reactive(kitagawa(rates_r()))

  # bootstrap: re-run at start, when the button is clicked or the country changes
  boot_r <- eventReactive(list(input$run, input$comparator), {
    req(input$B >= 100)
    set.seed(input$seed)
    withProgress(message = "Running parametric bootstrap...", value = 0.5, {
      bootstrap_all(rates_r(), input$B, alpha)
    })
  }, ignoreNULL = FALSE)

  # --- summary box ---
  output$summary_text <- renderUI({
    conv <- conv_r(); kit <- kit_r(); bt <- boot_r()
    tot <- conv[conv$icd10 == "ALL_CAUSES", ]
    f <- function(x) sprintf("%.3f", x)
    ci <- function(q) sprintf("[%.3f; %.3f]", q[1], q[2])
    HTML(paste0(
      "<p><b>Hungary vs ", comp_name(), "</b><br>",
      "All-cause ratio: ", f(tot$V_0), " (", base_year, ") &rarr; ",
      f(tot$V_1), " (", ref_year, "), k = ", f(tot$k), "<br>",
      "Sum of cause contributions (&Sigma;Kp): ", f(kit$sum_Kp), " ",
      ci(bt$ci_tot$sum_Kp), "<br>",
      "Composition effect (&Sigma;Kpp): ", f(kit$comp_effect), " ",
      ci(bt$ci_tot$comp), "</p>"))
  })

  # --- waterfall chart ---
  output$waterfall <- renderPlot({
    wf <- make_waterfall(conv_r(), kit_r())
    pad <- 0.05 * diff(range(c(wf$start[wf$type != "Total"],
                               wf$end)))
    ylo <- min(wf$start[wf$type != "Total"], wf$end[wf$type != "Total"]) - 2 * pad
    yhi <- max(wf$start[wf$type != "Total"], wf$end[wf$type != "Total"]) + 2 * pad

    ggplot(wf) +
      geom_rect(aes(xmin = id - 0.4, xmax = id + 0.4,
                    ymin = pmin(start, end), ymax = pmax(start, end),
                    fill = type)) +
      geom_text(aes(x = id, y = pmax(start, end),
                    label = ifelse(type == "Total",
                                   sprintf("%.3f", value),
                                   sprintf("%+.3f", value))),
                vjust = -0.5, size = 3.3) +
      scale_x_continuous(breaks = wf$id, labels = levels(wf$label)) +
      scale_fill_manual(values = c("Total" = "grey55",
                                   "Divergence (worse)" = "#c0392b",
                                   "Convergence (better)" = "#27ae60")) +
      coord_cartesian(ylim = c(ylo, yhi)) +
      labs(title = paste0("Hungary / ", comp_name(),
                          ": decomposition of the change in the all-cause mortality ratio"),
           x = NULL, y = "Mortality ratio (Hungary / comparator)", fill = NULL) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1),
            legend.position = "top")
  })

  # --- decomposition table ---
  output$kit_table <- renderTable({
    kit <- kit_r()$detail; bt <- boot_r()
    out <- kit %>%
      left_join(bt$ci_kit, by = "icd10") %>%
      transmute(
        Cause = cause_labels[icd10],
        `Ratio base year` = V_0,
        `Ratio ref. year` = V_1,
        `k` = k,
        `Kp` = Kp,
        `Kp 95% CI` = sprintf("[%.3f; %.3f]", Kp_low, Kp_up),
        `Kpp` = Kpp,
        `Kpp 95% CI` = sprintf("[%.3f; %.3f]", Kpp_low, Kpp_up))
    out
  }, digits = 3)

  # --- rate plot and table ---
  ci_long_r <- reactive({
    boot_r()$ci_long %>%
      mutate(
        Country = ifelse(startsWith(series, "inv"), "Hungary", comp_name()),
        Year    = ifelse(endsWith(series, "0"), base_year, ref_year),
        Cause   = cause_labels[icd10])
  })

  output$rate_plot <- renderPlot({
    d <- ci_long_r() %>% filter(icd10 != "ALL_CAUSES")
    ggplot(d, aes(x = observed, y = reorder(Cause, observed), colour = Country)) +
      geom_pointrange(aes(xmin = lower, xmax = upper),
                      position = position_dodge(width = 0.6)) +
      facet_wrap(~ Year) +
      labs(x = "Age-standardised death rate per 100 000 (95% bootstrap CI)",
           y = NULL, colour = NULL) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "top")
  })

  output$rate_table <- renderTable({
    ci_long_r() %>%
      transmute(Cause, Country, Year,
                `Rate [95% CI]` = sprintf("%.1f [%.1f; %.1f]", observed, lower, upper)) %>%
      pivot_wider(names_from = c(Country, Year), values_from = `Rate [95% CI]`) %>%
      arrange(match(Cause, cause_labels))
  })

  # --- significance table ---
  output$sig_table <- renderTable({
    boot_r()$sig %>%
      transmute(Cause = cause_labels[icd10],
                !!paste0("Signif. ", base_year) := base,
                !!paste0("Signif. ", ref_year)  := ref,
                !!paste0("Bonferroni ", base_year) := base_bonf,
                !!paste0("Bonferroni ", ref_year)  := ref_bonf)
  })
}

shinyApp(ui, server)
