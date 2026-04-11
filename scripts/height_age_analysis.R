# ============================================================
# Height–Age analysis (expects a clean plantation-only input dataset)
#
# Core analysis (PAPER):
#   - Modelling window: age ≤ 40 years
#   - Two-component strategy:
#       1) Baseline trajectory: Chapman–Richards fitted to non-fast-juvenile records (≤40)
#       2) Fast-juvenile pathway: power curve fitted over ages 1–10 (subset within ≤40)
#   - Weighting compared (≤40 only): unweighted (paper), capped-weighted (p99), weighted
#   - "All-records" curve shown only as within-window sensitivity (≤40), not >40
#
# Repository paths:
#   - Input data:  data/heightage.xlsx
#   - Expected worksheet: data
#   - Figures:     outputs/PDF/H/
#   - Tables:      outputs/xlsx/resultspaper.xlsx
#
# Last updated: 2026-01-24
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(ggplot2)
  library(minpack.lm)
  library(purrr)
  library(tibble)
  library(openxlsx)
  library(gridExtra)
  library(here)
})

options(na.print = "NA")

# ---------------------------
# PATHS (repository-relative)
# ---------------------------
data_path    <- file.path("data", "heightage.xlsx")
results_path <- file.path("outputs", "xlsx", "resultspaper.xlsx")
fig_dir      <- file.path("outputs", "PDF", "H")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(results_path), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_path)) {
  stop("Input file not found: ", data_path)
}

if (!"data" %in% readxl::excel_sheets(data_path)) {
  stop("Worksheet 'data' not found in: ", data_path)
}

# ---------------------------
# SETTINGS (match manuscript)
# ---------------------------
age_max_main     <- 40      # everything in the paper is within ≤40
A_fast_cluster   <- 15      # classify taxa using mean height at ages ≤15
A_fast_fit       <- 10      # fit power model on ages 1–10
age_min_fastfit  <- 1       # avoid age = 0 in power model
early_plot_max   <- 25      # x-limit for the early window figure

k_classes        <- 2
min_total_obs    <- 5
min_early_obs    <- 5

do_PI_supp       <- FALSE   # PI not reported in paper

# ---------------------------
# 1) Read data
# ---------------------------
raw <- as.data.frame(read_excel(data_path, sheet = "data", na = c("NA", "")))

needed <- c("ageexa","ageave","heightexa","heightave","taxon")
missing_needed <- setdiff(needed, names(raw))
if (length(missing_needed) > 0) {
  stop("Missing expected columns in height dataset: ", paste(missing_needed, collapse = ", "))
}

num_cols <- intersect(
  c("ageexa","ageave","heightexa","heightave","heightsd","heightmin","heightmax","notrees"),
  names(raw)
)
raw[num_cols] <- lapply(raw[num_cols], function(x) suppressWarnings(as.numeric(x)))

fac_cols <- intersect(c("taxon","level"), names(raw))
raw[fac_cols] <- lapply(raw[fac_cols], as.factor)

# ---------------------------
# 2) Harmonise variables + weights (ALL comparisons within ≤40)
# ---------------------------
dat <- raw %>%
  mutate(
    agecom = coalesce(ageexa, ageave),
    hcom   = coalesce(heightexa, heightave),
    taxon  = as.character(taxon)
  ) %>%
  dplyr::filter(!is.na(agecom), !is.na(hcom))

df_main <- dat %>% dplyr::filter(agecom <= age_max_main)

df_main <- df_main %>%
  mutate(
    weight = case_when(
      ("level" %in% names(.)) & !is.na(level) & level == "T" ~ 1,
      ("heightsd" %in% names(.)) & ("notrees" %in% names(.)) &
        !is.na(heightsd) & !is.na(notrees) & heightsd > 0 ~ notrees / (heightsd^2),
      ("heightmin" %in% names(.)) & ("heightmax" %in% names(.)) & ("notrees" %in% names(.)) &
        !is.na(heightmin) & !is.na(heightmax) & !is.na(notrees) &
        (heightmax - heightmin) > 0 ~ 16 * notrees / ((heightmax - heightmin)^2),
      ("notrees" %in% names(.)) & !is.na(notrees) ~ notrees,
      TRUE ~ 1
    )
  )

w_cap <- quantile(df_main$weight, 0.99, na.rm = TRUE)
df_main <- df_main %>% mutate(weight_fit = pmin(weight, w_cap))

# Sanity checks (these should match your manuscript if you’re running the same sheet)
cat("\nHeight dataset (all ages): n =", nrow(dat), " max age =", max(dat$agecom, na.rm=TRUE), "\n")
cat("Height modelling subset (age ≤ 40): n =", nrow(df_main), " max age =", max(df_main$agecom, na.rm=TRUE), "\n\n")

# ---------------------------
# 3) Identify fast-juvenile taxa (k-means on early mean height; within ≤40 subset)
# ---------------------------
tax_summ <- df_main %>%
  dplyr::filter(!is.na(taxon), taxon != "") %>%
  group_by(taxon) %>%
  summarise(
    n_tot = n(),
    n_early = sum(agecom <= A_fast_cluster, na.rm = TRUE),
    mean_h_early = mean(hcom[agecom <= A_fast_cluster], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(eligible = (n_tot >= min_total_obs & n_early >= min_early_obs))

early_summary <- tax_summ %>% dplyr::filter(eligible) %>%
  dplyr::select(taxon, n_tot, n_early, mean_h_early)

if (nrow(early_summary) < 2) {
  stop("Not enough eligible taxa for k-means. Lower thresholds or check data.")
}

set.seed(42)
km <- kmeans(scale(early_summary$mean_h_early), centers = k_classes, nstart = 25)
early_summary$cluster <- km$cluster

cluster_rank <- early_summary %>%
  group_by(cluster) %>%
  summarise(mu = mean(mean_h_early), .groups = "drop") %>%
  arrange(mu) %>%
  mutate(group = c("baseline", "fast_juvenile"))

early_summary <- early_summary %>% left_join(cluster_rank, by = "cluster")

df_main <- df_main %>%
  left_join(early_summary %>% dplyr::select(taxon, group), by = "taxon") %>%
  mutate(
    group = case_when(
      is.na(taxon) | taxon == "" ~ "unclassified",
      is.na(group) ~ "baseline",
      TRUE ~ as.character(group)
    ),
    group = factor(group, levels = c("baseline", "fast_juvenile", "unclassified")),
    # ---- UPDATED point legend (clear wording) ----
    point_group = ifelse(group == "fast_juvenile", "fast-juvenile", "non fast-juvenile")
  )

df_baseline <- df_main %>% dplyr::filter(group != "fast_juvenile")
df_fast_early <- df_main %>%
  dplyr::filter(group == "fast_juvenile", agecom >= age_min_fastfit, agecom <= A_fast_fit)

cat("Group counts (age ≤ 40):\n")
print(df_main %>% count(group) %>% mutate(pct = n / sum(n)))

fast_taxa <- early_summary %>%
  dplyr::filter(group == "fast_juvenile") %>%
  arrange(desc(mean_h_early))

# ---------------------------
# 4) Model helpers
# ---------------------------
fit_cr <- function(df, weights_col = NULL,
                   start = list(a = 30, b = 0.05, c = 1.1),
                   maxiter = 2000) {
  tryCatch({
    if (is.null(weights_col)) {
      nlsLM(
        hcom ~ a * (1 - exp(-b * agecom))^c,
        data = df,
        start = start,
        control = nls.lm.control(maxiter = maxiter)
      )
    } else {
      nlsLM(
        hcom ~ a * (1 - exp(-b * agecom))^c,
        data = df,
        start = start,
        weights = df[[weights_col]],
        control = nls.lm.control(maxiter = maxiter)
      )
    }
  }, error = function(e) NULL)
}

fit_power <- function(df, start = list(alpha = 6, beta = 0.5), maxiter = 2000) {
  tryCatch({
    nlsLM(
      hcom ~ alpha * (agecom^beta),
      data = df,
      start = start,
      control = nls.lm.control(maxiter = maxiter)
    )
  }, error = function(e) NULL)
}

extract_fit <- function(m, model_name = "") {
  if (is.null(m)) {
    return(tibble(Model = model_name, AIC = NA_real_, RSE = NA_real_, converged = FALSE))
  }
  tibble(
    Model = model_name,
    AIC = AIC(m),
    RSE = summary(m)$sigma,
    converged = TRUE
  ) %>% bind_cols(as_tibble(as.list(coef(m))))
}

predict_ci_cr <- function(model, age_seq, level = 0.95) {
  coefs <- coef(model); vc <- vcov(model)
  a <- coefs["a"]; b <- coefs["b"]; c <- coefs["c"]
  
  exp_term  <- exp(-b * age_seq)
  one_minus <- 1 - exp_term
  log_term  <- log1p(-exp_term)
  
  pred <- a * one_minus^c
  
  da <- one_minus^c
  db <- a * c * exp_term * age_seq * one_minus^(c - 1)
  dc <- a * one_minus^c * log_term
  
  X <- cbind(da, db, dc)
  se_mean <- sqrt(rowSums((X %*% vc) * X))
  
  z <- qnorm(1 - (1 - level) / 2)
  
  tibble(agecom = age_seq, fit = pred, lo = pred - z * se_mean, hi = pred + z * se_mean)
}

predict_ci_power <- function(model, age_seq, level = 0.95) {
  coefs <- coef(model); vc <- vcov(model)
  alpha <- coefs["alpha"]; beta <- coefs["beta"]
  
  pred <- alpha * (age_seq^beta)
  
  d_alpha <- age_seq^beta
  d_beta  <- alpha * (age_seq^beta) * log(age_seq)
  
  X <- cbind(d_alpha, d_beta)
  se_mean <- sqrt(rowSums((X %*% vc) * X))
  
  z <- qnorm(1 - (1 - level) / 2)
  
  tibble(agecom = age_seq, fit = pred, lo = pred - z * se_mean, hi = pred + z * se_mean)
}

save_diag_pdf <- function(df, model, out_pdf, title_prefix = "") {
  if (is.null(model) || nrow(df) < 10) return(invisible(FALSE))
  
  d <- df %>%
    mutate(fitted = as.numeric(predict(model)),
           resid  = hcom - fitted)
  
  p1 <- ggplot(d, aes(x = fitted, y = resid)) +
    geom_point(alpha = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = paste0(title_prefix, "Residuals vs fitted"), x = "Fitted", y = "Residuals") +
    theme_minimal()
  
  p2 <- ggplot(d, aes(sample = resid)) +
    stat_qq(alpha = 0.25) +
    stat_qq_line() +
    labs(title = paste0(title_prefix, "Normal Q–Q"), x = "Theoretical", y = "Sample") +
    theme_minimal()
  
  p3 <- ggplot(d, aes(x = fitted, y = sqrt(abs(resid)))) +
    geom_point(alpha = 0.25) +
    labs(title = paste0(title_prefix, "Scale–location"), x = "Fitted", y = "Sqrt(|resid|)") +
    theme_minimal()
  
  p4 <- ggplot(d, aes(x = agecom, y = resid)) +
    geom_point(alpha = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = paste0(title_prefix, "Residuals vs age"), x = "Age", y = "Residuals") +
    theme_minimal()
  
  grDevices::pdf(out_pdf, width = 10, height = 7)
  gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  grDevices::dev.off()
  invisible(TRUE)
}

# ---------------------------
# 5) Fits (ALL within age ≤ 40)
# ---------------------------
m_base_unw <- fit_cr(df_baseline, weights_col = NULL)
m_base_cap <- fit_cr(df_baseline, weights_col = "weight_fit")
m_base_wgt <- fit_cr(df_baseline, weights_col = "weight")

# within-window sensitivity (≤40): baseline vs including fast-juvenile
m_all_unw  <- fit_cr(df_main, weights_col = NULL)

m_fast_pow <- fit_power(df_fast_early)

tab_baseline <- bind_rows(
  extract_fit(m_base_unw, "Baseline CR – Unweighted (paper)"),
  extract_fit(m_base_cap, "Baseline CR – Capped-weighted (p99)"),
  extract_fit(m_base_wgt, "Baseline CR – Weighted")
)

tab_fast <- extract_fit(m_fast_pow, paste0("Fast-juvenile power – ages ", age_min_fastfit, "–", A_fast_fit))

tab_within40_sens <- bind_rows(
  extract_fit(m_base_unw, "Baseline CR (paper; ≤40; fast-juvenile excluded)"),
  extract_fit(m_all_unw,  "CR including fast-juvenile (≤40)")
)

print(tab_baseline)
print(tab_fast)
print(tab_within40_sens)

# ---------------------------
# 6) Figures (CI only; ≤40)
# ---------------------------
age_seq_fastcurve <- seq(age_min_fastfit, A_fast_fit, length.out = 220)
age_seq_cr_early  <- seq(0.5, early_plot_max, length.out = 220)
age_seq_full      <- seq(0.5, age_max_main, length.out = 260)

plot_early_CI <- function() {
  base <- ggplot(df_main %>% dplyr::filter(agecom <= early_plot_max),
                 aes(x = agecom, y = hcom, colour = point_group)) +
    geom_point(alpha = 0.35, size = 1.0) +
    theme_minimal() +
    labs(
      title = paste0("Height–age (age ≤ ", age_max_main, "): early window"),
      subtitle = paste0("Fast-juvenile taxa from mean height at ages ≤ ", A_fast_cluster,
                        "; power fit uses ages ", age_min_fastfit, "–", A_fast_fit),
      x = "Age (years)", y = "Height (m)", colour = "Points"
    ) +
    coord_cartesian(xlim = c(0, early_plot_max)) +
    # ---- UPDATED legend labels for points ----
  scale_colour_discrete(labels = c(
    "fast-juvenile" = "fast-juvenile",
    "non fast-juvenile" = "non fast-juvenile"
  ))
  
  p <- base
  
  if (!is.null(m_base_unw)) {
    bb <- predict_ci_cr(m_base_unw, age_seq_cr_early)
    p <- p +
      geom_ribbon(data = bb, aes(x = agecom, ymin = lo, ymax = hi),
                  inherit.aes = FALSE, alpha = 0.12, colour = NA) +
      geom_line(data = bb, aes(x = agecom, y = fit),
                inherit.aes = FALSE, linewidth = 1.1)
  }
  
  if (!is.null(m_fast_pow) && nrow(df_fast_early) >= 3) {
    bf <- predict_ci_power(m_fast_pow, age_seq_fastcurve)
    p <- p +
      geom_ribbon(data = bf, aes(x = agecom, ymin = lo, ymax = hi),
                  inherit.aes = FALSE, alpha = 0.12, colour = NA) +
      geom_line(data = bf, aes(x = agecom, y = fit),
                inherit.aes = FALSE, linewidth = 1.1)
  }
  
  p
}

plot_full_CI <- function() {
  base <- ggplot(df_main, aes(x = agecom, y = hcom, colour = point_group)) +
    geom_point(alpha = 0.25, size = 0.9) +
    theme_minimal() +
    labs(
      title = paste0("Height–age (age ≤ ", age_max_main, "): full window"),
      subtitle = "Baseline curve (paper) and within-window curve including fast-juvenile (≤40)",
      x = "Age (years)", y = "Height (m)", colour = "Points"
    ) +
    coord_cartesian(xlim = c(0, age_max_main)) +
    # ---- UPDATED legend labels for points ----
  scale_colour_discrete(labels = c(
    "fast-juvenile" = "fast-juvenile",
    "non fast-juvenile" = "non fast-juvenile"
  ))
  
  dd <- bind_rows(
    if (!is.null(m_base_unw)) predict_ci_cr(m_base_unw, age_seq_full) %>%
      mutate(Fit = "Baseline (fast-juvenile excluded)") else NULL,
    if (!is.null(m_all_unw))  predict_ci_cr(m_all_unw,  age_seq_full) %>%
      mutate(Fit = "Including fast-juvenile") else NULL
  )
  
  base +
    geom_ribbon(
      data = dd,
      aes(x = agecom, ymin = lo, ymax = hi, fill = Fit),
      inherit.aes = FALSE,
      alpha = 0.10,
      colour = NA
    ) +
    geom_line(
      data = dd,
      aes(x = agecom, y = fit, linetype = Fit),
      inherit.aes = FALSE,
      linewidth = 1.1
    ) +
    # ---- UPDATED legend labels for fits ----
  scale_linetype_manual(values = c(
    "Baseline (fast-juvenile excluded)" = "dashed",
    "Including fast-juvenile"          = "solid"
  )) +
    guides(fill = "none", linetype = guide_legend(title = "Fit"))
}

p_early_CI <- plot_early_CI()
p_full_CI  <- plot_full_CI()

print(p_early_CI)
print(p_full_CI)

ggsave(file.path(fig_dir, paste0("H_Fig_Early_CI_Aclust", A_fast_cluster, "_Afit", A_fast_fit, "_Amax", age_max_main, ".pdf")),
       p_early_CI, width = 7.0, height = 5.0)

ggsave(file.path(fig_dir, paste0("H_Fig_Full_CI_Aclust", A_fast_cluster, "_Afit", A_fast_fit, "_Amax", age_max_main, ".pdf")),
       p_full_CI, width = 7.0, height = 5.0)

message("Figures saved to: ", fig_dir)

# ---------------------------
# 7) Diagnostics
# ---------------------------
save_diag_pdf(df_baseline, m_base_unw,
              out_pdf = file.path(fig_dir, paste0("H_Diag_Baseline_CR_Amax", age_max_main, ".pdf")),
              title_prefix = "Baseline CR (paper; ≤40): ")

save_diag_pdf(df_main, m_all_unw,
              out_pdf = file.path(fig_dir, paste0("H_Diag_IncludingFast_CR_Amax", age_max_main, ".pdf")),
              title_prefix = "CR including fast-juvenile (≤40): ")

save_diag_pdf(df_fast_early, m_fast_pow,
              out_pdf = file.path(fig_dir, paste0("H_Diag_FastJuvenile_Power_Afit", A_fast_fit, ".pdf")),
              title_prefix = paste0("Fast-juvenile power (ages ", age_min_fastfit, "–", A_fast_fit, "; ≤40): "))

message("Diagnostics saved to: ", fig_dir)

# ---------------------------
# 8) Write outputs to Excel
# ---------------------------
write_sheet <- function(wb, name, df) {
  if (name %in% names(wb)) removeWorksheet(wb, name)
  addWorksheet(wb, name)
  writeData(wb, name, df, withFilter = TRUE)
}

wb <- if (file.exists(results_path)) loadWorkbook(results_path) else createWorkbook()

write_sheet(wb, paste0("H_Baseline_CR_A", age_max_main), tab_baseline)
write_sheet(wb, paste0("H_FastJuvenile_Power_A", age_min_fastfit, "_", A_fast_fit), tab_fast)
write_sheet(wb, paste0("H_Within40_Comparison_A", age_max_main), tab_within40_sens)

write_sheet(wb, paste0("H_Group_Counts_A", age_max_main),
            df_main %>% count(group) %>% mutate(pct = n / sum(n)))

write_sheet(wb, paste0("H_FastJuvenile_Taxa_A", age_max_main), fast_taxa)

write_sheet(wb, paste0("H_FastJuvenile_TaxonCounts_A", age_max_main),
            df_main %>% dplyr::filter(group == "fast_juvenile") %>% count(taxon, sort = TRUE))

saveWorkbook(wb, results_path, overwrite = TRUE)
message("Tables written to: ", results_path)
