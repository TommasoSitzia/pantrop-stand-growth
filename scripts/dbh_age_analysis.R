# ============================================================
# DBH–Age analysis (expects a clean plantation-only input dataset)
#
# Core analysis (PAPER):
#   - Modelling window: age ≤ 40 years
#   - Models: Chapman–Richards for pooled/slow/medium; power model for fast DBH (95% CI for mean curve)
#   - Fits compared (≤40 only): unweighted (paper), capped-weighted (p99), weighted
#   - Growth classes (≤40 only): k-means on taxon mean DBH at ages ≤ 15
#
# Repository paths:
#   - Input data:  data/dbhage.xlsx
#   - Expected worksheet: data
#   - Output root: outputs/
#   - Figures:     outputs/PDF/DBH/
#   - Tables:      outputs/xlsx/dbh_results.xlsx
#
# Last updated: 2026-07-22
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
})

options(na.print = "NA")

# ----------------------------- #
# PATHS (repository-relative)
# ----------------------------- #
data_path    <- file.path("data", "dbhage.xlsx")
output_root  <- "outputs"
results_path <- file.path(output_root, "xlsx", "dbh_results.xlsx")
fig_dir      <- file.path(output_root, "PDF", "DBH")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(results_path), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_path)) {
  stop("Input file not found: ", data_path)
}

if (!"data" %in% readxl::excel_sheets(data_path)) {
  stop("Worksheet 'data' not found in: ", data_path)
}

# ----------------------------- #
# SETTINGS
# ----------------------------- #
age_max_main   <- 40

# growth-class settings (taxon-based)
min_total_obs  <- 10
min_early_obs  <- 5
early_age_max  <- 15
k_classes      <- 3
kmeans_nstart  <- 25

# intervals
do_PI_supp     <- FALSE  # PI not reported in paper

# focused reviewer analyses
run_loso_cv       <- TRUE
run_class_sens    <- TRUE
run_fast_alt      <- TRUE
run_level_sens    <- TRUE
cv_age_breaks     <- c(-Inf, 5, 10, 20, 30, Inf)
cv_age_labels     <- c("<=5", "6-10", "11-20", "21-30", ">30")

# ----------------------------- #
# 1) Read data
# ----------------------------- #
raw <- as.data.frame(read_excel(data_path, sheet = "data", na = c("NA", "")))

needed <- c("ageexa","ageave","dbhexa","dbhave","taxon")
missing_needed <- setdiff(needed, names(raw))
if (length(missing_needed) > 0) {
  stop("Missing expected columns in DBH dataset: ", paste(missing_needed, collapse = ", "))
}

num_cols <- intersect(
  c("idarticle",
    "ageexa","ageave","agesd","agemin","agemax",
    "dbhexa","dbhave","dbhsd","dbhmin","dbhmax",
    "rain","temp","treeha","treespha","notrees","stand"),
  names(raw)
)
raw[num_cols] <- lapply(raw[num_cols], function(x) suppressWarnings(as.numeric(x)))

fac_cols <- intersect(c("taxon","country","location","level"), names(raw))
raw[fac_cols] <- lapply(raw[fac_cols], as.factor)

# ----------------------------- #
# 2) Harmonise variables + weights
# ----------------------------- #
dat <- raw %>%
  mutate(
    agecom = coalesce(ageexa, ageave),
    dbhcom = coalesce(dbhexa, dbhave),
    taxon  = as.character(taxon)
  ) %>%
  dplyr::filter(!is.na(agecom), !is.na(dbhcom)) %>%
  mutate(
    weight = case_when(
      ("level" %in% names(.)) & !is.na(level) & level == "T" ~ 1,
      ("dbhsd" %in% names(.)) & ("notrees" %in% names(.)) &
        !is.na(dbhsd) & !is.na(notrees) & dbhsd > 0 ~ notrees / (dbhsd^2),
      ("dbhmin" %in% names(.)) & ("dbhmax" %in% names(.)) & ("notrees" %in% names(.)) &
        !is.na(dbhmin) & !is.na(dbhmax) & !is.na(notrees) &
        (dbhmax - dbhmin) > 0 ~ 16 * notrees / ((dbhmax - dbhmin)^2),
      ("notrees" %in% names(.)) & !is.na(notrees) ~ notrees,
      TRUE ~ 1
    )
  )

# IMPORTANT: all fits/weights comparison are done on the modelling window (≤ 40)
df_main <- dat %>% dplyr::filter(agecom <= age_max_main)

w_cap <- quantile(df_main$weight, 0.99, na.rm = TRUE)
df_main <- df_main %>%
  mutate(weight_fit = pmin(weight, w_cap))

# Study-balanced weights: each source study has the same total weight.
# Rescaling to mean 1 does not change the NLS optimum, but keeps magnitudes readable.
if (!("idarticle" %in% names(df_main))) stop("Column 'idarticle' is required for study-aware analyses.")
df_main <- df_main %>%
  group_by(idarticle) %>%
  mutate(study_weight_raw = 1 / n()) %>%
  ungroup() %>%
  mutate(study_weight = study_weight_raw / mean(study_weight_raw, na.rm = TRUE))

# Sanity checks (these should match your paper counts if your sheet is the same one used for the manuscript)
cat("\nDBH dataset (all ages): n =", nrow(dat), " max age =", max(dat$agecom, na.rm=TRUE), "\n")
cat("DBH modelling subset (age ≤ 40): n =", nrow(df_main), " max age =", max(df_main$agecom, na.rm=TRUE), "\n\n")

age_seq_main <- seq(0.5, age_max_main, length.out = 220)

# ----------------------------- #
# 3) Chapman–Richards helpers
# ----------------------------- #
fit_cr <- function(df, weights_col = NULL,
                   start = list(a = 40, b = 0.05, c = 1.2),
                   maxiter = 1024) {
  tryCatch({
    if (is.null(weights_col)) {
      nlsLM(
        dbhcom ~ a * (1 - exp(-b * agecom))^c,
        data = df,
        start = start,
        control = nls.lm.control(maxiter = maxiter)
      )
    } else {
      nlsLM(
        dbhcom ~ a * (1 - exp(-b * agecom))^c,
        data = df,
        start = start,
        weights = df[[weights_col]],
        control = nls.lm.control(maxiter = maxiter)
      )
    }
  }, error = function(e) NULL)
}

extract_fit <- function(m) {
  if (is.null(m)) {
    return(tibble(AIC = NA_real_, RSE = NA_real_, a = NA_real_, b = NA_real_, c = NA_real_, converged = FALSE))
  }
  tibble(
    AIC = AIC(m),
    RSE = summary(m)$sigma,
    a = unname(coef(m)["a"]),
    b = unname(coef(m)["b"]),
    c = unname(coef(m)["c"]),
    converged = TRUE
  )
}

predict_ci_cr <- function(model, age_seq, level = 0.95) {
  coefs <- coef(model)
  vc <- vcov(model)
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

# Delta-method confidence interval for the fast-class power model:
# DBH(t) = alpha * t^beta
predict_ci_power <- function(model, age_seq, level = 0.95) {
  coefs <- coef(model)
  vc <- vcov(model)
  alpha <- unname(coefs["alpha"])
  beta  <- unname(coefs["beta"])

  pred <- alpha * age_seq^beta

  # Gradient with respect to alpha and beta
  d_alpha <- age_seq^beta
  d_beta  <- alpha * age_seq^beta * log(age_seq)
  X <- cbind(d_alpha, d_beta)

  se_mean <- sqrt(rowSums((X %*% vc) * X))
  z <- qnorm(1 - (1 - level) / 2)

  tibble(
    agecom = age_seq,
    fit = pred,
    lo = pred - z * se_mean,
    hi = pred + z * se_mean
  )
}

save_diag_pdf <- function(df, model, out_pdf, title_prefix = "") {
  if (is.null(model) || nrow(df) < 10) return(invisible(FALSE))
  
  d <- df %>% mutate(
    fitted = as.numeric(predict(model)),
    resid  = dbhcom - fitted
  )
  
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
  
  grDevices::pdf(out_pdf, width = 8.5, height = 6.5, onefile = TRUE)
  gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  grDevices::dev.off()
  
  invisible(TRUE)
}

# ----------------------------- #
# 3b) Focused reviewer-analysis helpers
# ----------------------------- #
extract_param_ci <- function(m, fit_label, level = 0.95) {
  if (is.null(m)) return(tibble(Fit = fit_label, parameter = NA_character_, estimate = NA_real_, SE = NA_real_, lower = NA_real_, upper = NA_real_))
  cc <- coef(m)
  se <- sqrt(diag(vcov(m)))
  z <- qnorm(1 - (1 - level) / 2)
  tibble(Fit = fit_label, parameter = names(cc), estimate = unname(cc), SE = unname(se),
         lower = estimate - z * SE, upper = estimate + z * SE)
}

prediction_metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  if (!length(obs)) return(tibble(n = 0L, RMSE = NA_real_, MAE = NA_real_, Bias = NA_real_))
  err <- pred - obs
  tibble(n = length(obs), RMSE = sqrt(mean(err^2)), MAE = mean(abs(err)), Bias = mean(err))
}

loso_cr <- function(df, study_col, response_col, fit_fun, age_breaks, age_labels) {
  studies <- sort(unique(df[[study_col]]))
  pred_rows <- purrr::map_dfr(studies, function(st) {
    train <- df[df[[study_col]] != st, , drop = FALSE]
    test  <- df[df[[study_col]] == st, , drop = FALSE]
    m <- fit_fun(train)
    if (is.null(m)) {
      return(tibble(study = st, agecom = test$agecom, observed = test[[response_col]], predicted = NA_real_, converged = FALSE))
    }
    tibble(study = st, agecom = test$agecom, observed = test[[response_col]],
           predicted = as.numeric(predict(m, newdata = test)), converged = TRUE)
  })
  by_study <- pred_rows %>%
    group_by(study) %>%
    group_modify(~ prediction_metrics(.x$observed, .x$predicted)) %>%
    ungroup()
  overall <- prediction_metrics(pred_rows$observed, pred_rows$predicted) %>% mutate(scope = "All ages")
  by_age <- pred_rows %>%
    mutate(age_band = cut(agecom, breaks = age_breaks, labels = age_labels, right = TRUE)) %>%
    group_by(age_band) %>%
    group_modify(~ prediction_metrics(.x$observed, .x$predicted)) %>%
    ungroup() %>% mutate(scope = as.character(age_band)) %>% select(-age_band)
  list(predictions = pred_rows, by_study = by_study, summary = bind_rows(overall, by_age))
}

adjusted_rand_index <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 2) return(NA_real_)
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  sum_nij <- sum(choose2(tab))
  sum_ai <- sum(choose2(rowSums(tab)))
  sum_bj <- sum(choose2(colSums(tab)))
  total <- choose2(n)
  expected <- sum_ai * sum_bj / total
  max_index <- 0.5 * (sum_ai + sum_bj)
  if (max_index == expected) return(1)
  (sum_nij - expected) / (max_index - expected)
}

classify_dbh_taxa <- function(df, min_total, min_early, early_max, k = 3) {
  eligible <- df %>% count(taxon, name = "n_tot") %>% filter(n_tot >= min_total)
  es <- df %>%
    filter(taxon %in% eligible$taxon, agecom <= early_max) %>%
    group_by(taxon) %>%
    summarise(mean_early = mean(dbhcom, na.rm = TRUE), n_early = n(), .groups = "drop") %>%
    filter(n_early >= min_early)
  if (nrow(es) < k) return(tibble(taxon = character(), class = character(), mean_early = numeric(), n_early = integer()))
  set.seed(42)
  kk <- kmeans(scale(es$mean_early), centers = k, nstart = kmeans_nstart)
  es$cluster <- kk$cluster
  lab <- es %>% group_by(cluster) %>% summarise(mu = mean(mean_early), .groups = "drop") %>%
    arrange(mu) %>% mutate(class = c("slow", "medium", "fast"))
  es %>% left_join(lab, by = "cluster") %>% select(taxon, class, mean_early, n_early)
}



# Alternative models for the fast DBH class
fit_power_dbh <- function(df, start = list(alpha = 5, beta = 0.7), maxiter = 1024) {
  tryCatch({
    nlsLM(
      dbhcom ~ alpha * (agecom^beta),
      data = df,
      start = start,
      control = nls.lm.control(maxiter = maxiter)
    )
  }, error = function(e) NULL)
}

fit_log_dbh <- function(df) {
  tryCatch(
    lm(dbhcom ~ log1p(agecom), data = df),
    error = function(e) NULL
  )
}

extract_fit_generic <- function(m, model_name = "") {
  if (is.null(m)) {
    return(tibble(Model = model_name, AIC = NA_real_, RSE = NA_real_, converged = FALSE))
  }
  out <- tibble(
    Model = model_name,
    AIC = AIC(m),
    RSE = if (inherits(m, "lm")) summary(m)$sigma else summary(m)$sigma,
    converged = TRUE
  )
  bind_cols(out, as_tibble(as.list(coef(m))))
}

within_sample_metrics <- function(df, model, response_col) {
  if (is.null(model)) return(tibble(n = nrow(df), RMSE = NA_real_, MAE = NA_real_, Bias = NA_real_))
  pred <- tryCatch(as.numeric(predict(model, newdata = df)), error = function(e) rep(NA_real_, nrow(df)))
  prediction_metrics(df[[response_col]], pred)
}

# Reporting-level sensitivity helper. Each fit changes only the included reporting levels.
reporting_level_sensitivity_dbh <- function(df) {
  if (!("level" %in% names(df))) {
    return(list(summary = tibble(), parameters = tibble(), predictions = tibble()))
  }

  d <- df %>% filter(!is.na(level)) %>% mutate(level_chr = as.character(level))
  levs <- sort(unique(d$level_chr))

  scenarios <- list(`All levels` = d)
  for (lv in levs) scenarios[[paste0("Only level ", lv)]] <- d %>% filter(level_chr == lv)
  for (lv in levs) scenarios[[paste0("Exclude level ", lv)]] <- d %>% filter(level_chr != lv)

  fits <- purrr::imap(scenarios, function(dd, nm) {
    if (nrow(dd) < 10 || dplyr::n_distinct(dd$agecom) < 4) return(NULL)
    fit_cr(dd, weights_col = NULL)
  })

  summary <- purrr::imap_dfr(fits, function(m, nm) {
    dd <- scenarios[[nm]]
    bind_cols(
      tibble(Scenario = nm, n = nrow(dd), n_studies = n_distinct(dd$idarticle), levels_included = paste(sort(unique(dd$level_chr)), collapse = ", ")),
      extract_fit(m)
    )
  })

  parameters <- purrr::imap_dfr(fits, function(m, nm) extract_param_ci(m, nm))

  age_seq <- seq(0.5, age_max_main, length.out = 220)
  predictions <- purrr::imap_dfr(fits, function(m, nm) {
    if (is.null(m)) return(NULL)
    tibble(Scenario = nm, agecom = age_seq, predicted = as.numeric(predict(m, newdata = data.frame(agecom = age_seq))))
  })

  list(summary = summary, parameters = parameters, predictions = predictions)
}

# ----------------------------- #
# 4) Pooled models (age ≤ 40 ONLY)
# ----------------------------- #
m_unw <- fit_cr(df_main, weights_col = NULL)
m_cap <- fit_cr(df_main, weights_col = "weight_fit")
m_wgt <- fit_cr(df_main, weights_col = "weight")
m_study_bal <- fit_cr(df_main, weights_col = "study_weight")

m_paper <- m_unw
paper_label <- "Unweighted"

tab_pooled <- bind_rows(
  extract_fit(m_unw) %>% mutate(Fit = "Unweighted"),
  extract_fit(m_cap) %>% mutate(Fit = "Capped-weighted (p99)"),
  extract_fit(m_wgt) %>% mutate(Fit = "Weighted"),
  extract_fit(m_study_bal) %>% mutate(Fit = "Study-balanced")
) %>%
  mutate(
    Dataset = paste0("Age ≤ ", age_max_main),
    Model = "Chapman–Richards: DBH(t)=a*(1-exp(-b*t))^c"
  ) %>%
  dplyr::select(Dataset, Model, Fit, AIC, RSE, a, b, c, converged)

print(tab_pooled)

# ----------------------------- #
# 5) Growth classes (taxon-based; age ≤ 40 ONLY)
# ----------------------------- #
# For classification we still use the DBH data available (not restricted to ≤40 is fine),
# but in practice your DBH dataset is dominated by young ages; to keep everything aligned
# with the manuscript framing, we classify using the same modelling subset df_main.
taxa_ok <- df_main %>%
  group_by(taxon) %>%
  summarise(n_tot = n(), .groups = "drop") %>%
  dplyr::filter(n_tot >= min_total_obs) %>%
  pull(taxon)

early_summary <- df_main %>%
  dplyr::filter(taxon %in% taxa_ok, agecom <= early_age_max) %>%
  group_by(taxon) %>%
  summarise(
    mean_dbh_early = mean(dbhcom, na.rm = TRUE),
    n_early = n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_early >= min_early_obs)

if (nrow(early_summary) < k_classes) {
  stop("Not enough eligible taxa for k-means. Lower thresholds or check data.")
}

set.seed(42)
km <- kmeans(scale(early_summary$mean_dbh_early), centers = k_classes, nstart = kmeans_nstart)
early_summary$cluster <- km$cluster

cluster_map <- early_summary %>%
  group_by(cluster) %>%
  summarise(mu = mean(mean_dbh_early), .groups = "drop") %>%
  arrange(mu) %>%
  mutate(growth_class = c("slow", "medium", "fast"))

early_summary <- early_summary %>% left_join(cluster_map, by = "cluster")

df_main_classed <- df_main %>%
  left_join(early_summary %>% dplyr::select(taxon, growth_class), by = "taxon") %>%
  mutate(
    growth_class = ifelse(is.na(growth_class), "unclassified", growth_class),
    growth_class = factor(growth_class, levels = c("slow","medium","fast","unclassified"))
  )

df_fitclass <- df_main_classed %>%
  dplyr::filter(growth_class %in% c("slow","medium","fast")) %>%
  droplevels()

counts_class <- df_fitclass %>%
  count(growth_class) %>%
  mutate(pct = n / sum(n))

print(counts_class)

models_class <- df_fitclass %>%
  dplyr::group_split(growth_class) %>%
  stats::setNames(levels(df_fitclass$growth_class)) %>%
  purrr::map(~ fit_cr(.x, weights_col = NULL, start = list(a = 40, b = 0.05, c = 1.0)))

# Preliminary Chapman–Richards class table (slow/medium plus fast CR alternative).
tab_class_cr <- purrr::imap_dfr(models_class, function(m, cls) {
  extract_fit(m) %>% dplyr::mutate(growth_class = cls)
}) %>%
  dplyr::mutate(
    Dataset = paste0("Age ≤ ", age_max_main),
    Fit = "Unweighted",
    Model = "Chapman–Richards by growth class"
  ) %>%
  dplyr::select(Dataset, Model, growth_class, Fit, AIC, RSE, a, b, c, converged)

print(tab_class_cr)

# ----------------------------- #
# 5b) Focused reviewer analyses
# ----------------------------- #
param_ci_pooled <- bind_rows(
  extract_param_ci(m_unw, "Unweighted"),
  extract_param_ci(m_cap, "Capped-weighted (p99)"),
  extract_param_ci(m_wgt, "Weighted"),
  extract_param_ci(m_study_bal, "Study-balanced")
)

weight_diagnostics <- df_main %>%
  summarise(
    n = n(), n_studies = n_distinct(idarticle),
    min = min(weight, na.rm = TRUE), p50 = quantile(weight, 0.50, na.rm = TRUE),
    p90 = quantile(weight, 0.90, na.rm = TRUE), p95 = quantile(weight, 0.95, na.rm = TRUE),
    p99 = quantile(weight, 0.99, na.rm = TRUE), max = max(weight, na.rm = TRUE),
    n_capped = sum(weight > w_cap, na.rm = TRUE)
  )

study_dominance <- df_main %>% count(idarticle, sort = TRUE, name = "n_records") %>%
  mutate(share = n_records / sum(n_records), cumulative_share = cumsum(share), rank = row_number())

if (run_loso_cv) {
  loso <- loso_cr(df_main, "idarticle", "dbhcom", function(d) fit_cr(d, weights_col = NULL),
                  cv_age_breaks, cv_age_labels)
  loso_summary <- loso$summary
  loso_by_study <- loso$by_study
  loso_predictions <- loso$predictions
} else {
  loso_summary <- loso_by_study <- loso_predictions <- tibble()
}

if (run_class_sens) {
  ref_class <- classify_dbh_taxa(df_main, min_total_obs, min_early_obs, early_age_max, k_classes)
  sens_grid <- tidyr::crossing(
    min_total = c(5, 10, 15),
    min_early = c(3, 5, 8),
    early_max = c(10, 15, 20)
  )
  class_sensitivity <- purrr::pmap_dfr(sens_grid, function(min_total, min_early, early_max) {
    z <- classify_dbh_taxa(df_main, min_total, min_early, early_max, k_classes)
    common <- inner_join(ref_class %>% select(taxon, ref = class), z %>% select(taxon, alt = class), by = "taxon")
    tibble(min_total = min_total, min_early = min_early, early_max = early_max,
           n_eligible = nrow(z), n_common = nrow(common),
           agreement = ifelse(nrow(common), mean(common$ref == common$alt), NA_real_),
           ARI = ifelse(nrow(common) > 1, adjusted_rand_index(common$ref, common$alt), NA_real_),
           n_slow = sum(z$class == "slow"), n_medium = sum(z$class == "medium"), n_fast = sum(z$class == "fast"))
  })

} else class_sensitivity <- tibble()

# Fast-class alternative-model comparison
if (run_fast_alt) {
  df_fast_dbh <- df_fitclass %>% filter(growth_class == "fast")
  m_fast_cr_alt <- models_class[["fast"]]
  m_fast_power <- fit_power_dbh(df_fast_dbh)
  m_fast_log <- fit_log_dbh(df_fast_dbh)

  fast_alt_models <- bind_rows(
    bind_cols(extract_fit_generic(m_fast_cr_alt, "Chapman-Richards"), within_sample_metrics(df_fast_dbh, m_fast_cr_alt, "dbhcom")),
    bind_cols(extract_fit_generic(m_fast_power, "Power"), within_sample_metrics(df_fast_dbh, m_fast_power, "dbhcom")),
    bind_cols(extract_fit_generic(m_fast_log, "Log-linear"), within_sample_metrics(df_fast_dbh, m_fast_log, "dbhcom"))
  ) %>% mutate(n = nrow(df_fast_dbh), age_min = min(df_fast_dbh$agecom), age_max = max(df_fast_dbh$agecom))

  fast_alt_predictions <- bind_rows(
    if (!is.null(m_fast_cr_alt)) tibble(Model = "Chapman-Richards", agecom = age_seq_main,
      predicted = as.numeric(predict(m_fast_cr_alt, newdata = data.frame(agecom = age_seq_main)))) else NULL,
    if (!is.null(m_fast_power)) tibble(Model = "Power", agecom = age_seq_main,
      predicted = as.numeric(predict(m_fast_power, newdata = data.frame(agecom = age_seq_main)))) else NULL,
    if (!is.null(m_fast_log)) tibble(Model = "Log-linear", agecom = age_seq_main,
      predicted = as.numeric(predict(m_fast_log, newdata = data.frame(agecom = age_seq_main)))) else NULL
  )
} else {
  df_fast_dbh <- tibble()
  fast_alt_models <- fast_alt_predictions <- tibble()
}

# Manuscript-aligned class-model table:
# slow and medium use Chapman–Richards; fast uses the selected power model.
tab_class <- bind_rows(
  extract_fit(models_class[["slow"]]) %>%
    transmute(
      Dataset = paste0("Age ≤ ", age_max_main),
      Model = "Chapman–Richards",
      growth_class = "slow",
      Fit = "Unweighted",
      AIC, RSE, a, b, c,
      alpha = NA_real_, beta = NA_real_,
      converged
    ),
  extract_fit(models_class[["medium"]]) %>%
    transmute(
      Dataset = paste0("Age ≤ ", age_max_main),
      Model = "Chapman–Richards",
      growth_class = "medium",
      Fit = "Unweighted",
      AIC, RSE, a, b, c,
      alpha = NA_real_, beta = NA_real_,
      converged
    ),
  {
    if (run_fast_alt && !is.null(m_fast_power)) {
      ff <- extract_fit_generic(m_fast_power, "Power")
      tibble(
        Dataset = paste0("Age ≤ ", age_max_main),
        Model = "Power",
        growth_class = "fast",
        Fit = "Unweighted",
        AIC = ff$AIC,
        RSE = ff$RSE,
        a = NA_real_, b = NA_real_, c = NA_real_,
        alpha = unname(coef(m_fast_power)["alpha"]),
        beta = unname(coef(m_fast_power)["beta"]),
        converged = TRUE
      )
    } else {
      tibble(
        Dataset = paste0("Age ≤ ", age_max_main),
        Model = "Power",
        growth_class = "fast",
        Fit = "Unweighted",
        AIC = NA_real_, RSE = NA_real_,
        a = NA_real_, b = NA_real_, c = NA_real_,
        alpha = NA_real_, beta = NA_real_,
        converged = FALSE
      )
    }
  }
)

print(tab_class)

# Reporting-level sensitivity
if (run_level_sens) {
  level_sens_dbh <- reporting_level_sensitivity_dbh(df_main)
} else {
  level_sens_dbh <- list(summary = tibble(), parameters = tibble(), predictions = tibble())
}

# ----------------------------- #
# 6) Figures (CI only; age ≤ 40)
# ----------------------------- #
plot_pooled_CI <- function(show_all_fits = TRUE) {
  base <- ggplot(df_main, aes(x = agecom, y = dbhcom)) +
    geom_point(alpha = 0.20, size = 0.8) +
    theme_minimal() +
    labs(
      title = paste0("DBH–age in plantations (age ≤ ", age_max_main, ")"),
      x = "Age (years)", y = "DBH (cm)"
    ) +
    coord_cartesian(xlim = c(0, age_max_main))
  
  dd <- bind_rows(
    if (!is.null(m_unw)) predict_ci_cr(m_unw, age_seq_main) %>% mutate(Fit = "Unweighted") else NULL,
    if (!is.null(m_cap)) predict_ci_cr(m_cap, age_seq_main) %>% mutate(Fit = "Capped-weighted (p99)") else NULL,
    if (!is.null(m_wgt)) predict_ci_cr(m_wgt, age_seq_main) %>% mutate(Fit = "Weighted") else NULL
  )
  
  if (!show_all_fits) {
    dd <- dd %>% dplyr::filter(Fit == paper_label)
  }
  
  base +
    geom_ribbon(
      data = dd %>% dplyr::filter(Fit == paper_label),
      aes(x = agecom, ymin = lo, ymax = hi),
      inherit.aes = FALSE,
      alpha = 0.12,
      colour = NA
    ) +
    geom_line(
      data = dd,
      aes(x = agecom, y = fit, linetype = Fit),
      inherit.aes = FALSE,
      linewidth = 1.1
    ) +
    guides(linetype = guide_legend(title = "Fit"))
}

plot_classes_CI <- function() {
  base <- ggplot(
    df_fitclass,
    aes(x = agecom, y = dbhcom, colour = growth_class)
  ) +
    geom_point(alpha = 0.20, size = 0.8) +
    theme_minimal() +
    labs(
      title = paste0(
        "DBH–age by growth class (age ≤ ",
        age_max_main,
        ")"
      ),
      x = "Age (years)",
      y = "DBH (cm)",
      colour = "Growth class"
    ) +
    coord_cartesian(xlim = c(0, age_max_main))

  # Slow and medium classes retain Chapman–Richards models.
  pred_cr <- purrr::map_dfr(c("slow", "medium"), function(cls) {
    m <- models_class[[cls]]
    if (is.null(m)) return(NULL)

    predict_ci_cr(m, age_seq_main) %>%
      dplyr::mutate(
        growth_class = cls,
        model = "Chapman–Richards"
      )
  })

  # The fast class uses the selected non-asymptotic power model and is
  # drawn only over the observed fast-class age range.
  pred_fast <- NULL
  if (run_fast_alt && !is.null(m_fast_power) && nrow(df_fast_dbh) > 0) {
    fast_age_seq <- seq(
      max(0.5, min(df_fast_dbh$agecom, na.rm = TRUE)),
      min(age_max_main, max(df_fast_dbh$agecom, na.rm = TRUE)),
      length.out = 220
    )

    pred_fast <- predict_ci_power(m_fast_power, fast_age_seq) %>%
      dplyr::mutate(
        growth_class = "fast",
        model = "Power"
      )
  }

  pred_df <- bind_rows(pred_cr, pred_fast) %>%
    mutate(growth_class = factor(growth_class, levels = c("slow", "medium", "fast")))

  class_colours <- c(
    slow   = "#F8766D",
    medium = "#00BA38",
    fast   = "#619CFF"
  )

  base +
    geom_ribbon(
      data = pred_df,
      aes(
        x = agecom,
        ymin = lo,
        ymax = hi,
        fill = growth_class
      ),
      inherit.aes = FALSE,
      alpha = 0.12,
      colour = NA,
      show.legend = FALSE
    ) +
    geom_line(
      data = pred_df,
      aes(
        x = agecom,
        y = fit,
        colour = growth_class
      ),
      inherit.aes = FALSE,
      linewidth = 1.1
    ) +
    scale_colour_manual(values = class_colours) +
    scale_fill_manual(values = class_colours)
}


p_main_pooled  <- plot_pooled_CI(show_all_fits = TRUE)
p_main_classes <- plot_classes_CI()

print(p_main_pooled)
print(p_main_classes)

ggsave(file.path(fig_dir, paste0("DBH_Fig_Main_Pooled_CI_A", age_max_main, ".pdf")),
       p_main_pooled, width = 7.0, height = 5.0)

ggsave(file.path(fig_dir, paste0("DBH_Fig_Main_Classes_CI_A", age_max_main, ".pdf")),
       p_main_classes, width = 7.0, height = 5.0)



if (run_fast_alt && nrow(fast_alt_predictions) > 0) {
  p_fast_alt <- ggplot(df_fast_dbh, aes(x = agecom, y = dbhcom)) +
    geom_point(alpha = 0.25, size = 0.9) +
    geom_line(data = fast_alt_predictions,
              aes(x = agecom, y = predicted, linetype = Model),
              inherit.aes = FALSE, linewidth = 1.1) +
    coord_cartesian(xlim = c(0, age_max_main)) +
    theme_minimal() +
    labs(title = "Fast DBH class: alternative within-range models",
         subtitle = "Curves are descriptive within the observed age range; asymptotic extrapolation is not implied",
         x = "Age (years)", y = "DBH (cm)", linetype = "Model")
  print(p_fast_alt)
  ggsave(file.path(fig_dir, "DBH_FastClass_AlternativeModels.pdf"), p_fast_alt, width = 7, height = 5)
}

if (run_level_sens && nrow(level_sens_dbh$predictions) > 0) {
  p_level_dbh <- ggplot(level_sens_dbh$predictions,
                        aes(x = agecom, y = predicted, linetype = Scenario)) +
    geom_line(linewidth = 1) +
    coord_cartesian(xlim = c(0, age_max_main)) +
    theme_minimal() +
    labs(title = "DBH reporting-level sensitivity",
         x = "Age (years)", y = "Predicted DBH (cm)", linetype = "Scenario")
  print(p_level_dbh)
  ggsave(file.path(fig_dir, "DBH_ReportingLevel_Sensitivity.pdf"), p_level_dbh, width = 8, height = 5.5)
}

# ----------------------------- #
# 7) Diagnostics PDFs
# ----------------------------- #
save_diag_pdf(df_main, m_paper,
              out_pdf = file.path(fig_dir, paste0("DBH_Diag_Pooled_", paper_label, "_Amax", age_max_main, ".pdf")),
              title_prefix = paste0("Pooled ", paper_label, " (Amax ", age_max_main, "): "))

# Diagnostics for the selected class models.
invisible(purrr::imap(models_class[c("slow", "medium")], function(m, cls) {
  d_cls <- df_fitclass %>% dplyr::filter(growth_class == cls)
  save_diag_pdf(
    d_cls, m,
    out_pdf = file.path(fig_dir, paste0("DBH_Diag_Class_", cls, "_Amax", age_max_main, ".pdf")),
    title_prefix = paste0("Class ", cls, " (Amax ", age_max_main, "): ")
  )
}))

# Keep the fast Chapman–Richards diagnostic only as an explicitly named alternative-model check.
if (run_fast_alt && !is.null(m_fast_cr_alt) && nrow(df_fast_dbh) >= 10) {
  save_diag_pdf(
    df_fast_dbh,
    m_fast_cr_alt,
    out_pdf = file.path(fig_dir, paste0("DBH_Diag_Class_fast_CR_alternative_Amax", age_max_main, ".pdf")),
    title_prefix = paste0("Fast DBH class, Chapman–Richards alternative (Amax ", age_max_main, "): ")
  )
}

# Diagnostics for the selected fast-class power model used in Figure 2.
if (run_fast_alt && !is.null(m_fast_power) && nrow(df_fast_dbh) >= 10) {
  save_diag_pdf(
    df_fast_dbh,
    m_fast_power,
    out_pdf = file.path(
      fig_dir,
      paste0("DBH_Diag_Class_fast_power_Amax", age_max_main, ".pdf")
    ),
    title_prefix = paste0(
      "Fast DBH class, power model (Amax ", age_max_main, "): "
    )
  )
}

message("Figures + diagnostics saved to: ", fig_dir)

# ----------------------------- #
# 8) Write outputs to Excel
# ----------------------------- #
write_sheet <- function(wb, name, df) {
  if (name %in% names(wb)) removeWorksheet(wb, name)
  addWorksheet(wb, name)
  writeData(wb, name, df, withFilter = TRUE)
}

wb <- createWorkbook()

write_sheet(wb, paste0("DBH_Pooled_A", age_max_main), tab_pooled)
write_sheet(wb, paste0("DBH_GrowthClass_A", age_max_main), tab_class)
write_sheet(wb, paste0("DBH_GrowthClass_Counts_A", age_max_main), counts_class)

write_sheet(wb, "DBH_GrowthClass_Taxa",
            early_summary %>% arrange(growth_class, desc(mean_dbh_early)))

write_sheet(wb, "DBH_Param_CI", param_ci_pooled)
write_sheet(wb, "DBH_Weight_Diag", weight_diagnostics)
write_sheet(wb, "DBH_Study_Dominance", study_dominance)
write_sheet(wb, "DBH_LOSO_Summary", loso_summary)
write_sheet(wb, "DBH_LOSO_ByStudy", loso_by_study)
write_sheet(wb, "DBH_LOSO_Pred", loso_predictions)
write_sheet(wb, "DBH_Class_Sensitivity", class_sensitivity)
write_sheet(wb, "DBH_Fast_AltModels", fast_alt_models)
write_sheet(wb, "DBH_Fast_AltPred", fast_alt_predictions)
write_sheet(wb, "DBH_Level_Sens", level_sens_dbh$summary)
write_sheet(wb, "DBH_Level_Sens_Param", level_sens_dbh$parameters)
write_sheet(wb, "DBH_Level_Sens_Pred", level_sens_dbh$predictions)

saveWorkbook(wb, results_path, overwrite = TRUE)
message("Tables written to: ", results_path)
