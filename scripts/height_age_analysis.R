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
# Paths:
#   - Input data:  data/heightage.xlsx
#   - Expected worksheet: data
#   - Output root: outputs/
#   - Figures:     outputs/PDF/H/
#   - Tables:      outputs/xlsx/height_results.xlsx
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

# ---------------------------
# PATHS
# ---------------------------
data_path    <- file.path("data", "heightage.xlsx")
output_root  <- "outputs"
results_path <- file.path(output_root, "xlsx", "height_results.xlsx")
fig_dir      <- file.path(output_root, "PDF", "H")

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
kmeans_nstart   <- 25

do_PI_supp       <- FALSE   # PI not reported in paper

# focused reviewer analyses
run_loso_cv       <- TRUE
run_class_sens    <- TRUE
run_level_sens    <- TRUE
cv_age_breaks     <- c(-Inf, 5, 10, 20, 30, Inf)
cv_age_labels     <- c("<=5", "6-10", "11-20", "21-30", ">30")

# ---------------------------
# 1) Read data
# ---------------------------
raw <- as.data.frame(read_excel(data_path, sheet = "data", na = c("NA", "")))

needed <- c("ageexa","ageave","heightexa","heightave","taxon")
missing_needed <- setdiff(needed, names(raw))
if (length(missing_needed) > 0) {
  stop("Missing expected columns in height dataset: ", paste(missing_needed, collapse = ", "))
}

# Keep the study identifier as text: values such as "39H" must never be
# coerced to numeric, otherwise they become NA and study-aware analyses fail.
num_cols <- intersect(
  c(
    "ageexa", "ageave",
    "heightexa", "heightave",
    "heightsd", "heightmin", "heightmax",
    "notrees"
  ),
  names(raw)
)
raw[num_cols] <- lapply(raw[num_cols], function(x) suppressWarnings(as.numeric(x)))

fac_cols <- intersect(
  c("taxon", "country", "location", "level"),
  names(raw)
)
raw[fac_cols] <- lapply(raw[fac_cols], as.factor)

if (!("idarticleh" %in% names(raw))) {
  stop("Column 'idarticleh' is required for study-aware analyses.")
}
raw$idarticleh <- trimws(as.character(raw$idarticleh))
raw$idarticleh[raw$idarticleh == ""] <- NA_character_

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
df_main <- df_main %>%
  mutate(weight_fit = pmin(weight, w_cap))

if (anyNA(df_main$idarticleh)) {
  warning(sum(is.na(df_main$idarticleh)),
          " height records have no study identifier and will be excluded from study-aware analyses.")
}

n_studies_main <- dplyr::n_distinct(df_main$idarticleh, na.rm = TRUE)
if (n_studies_main < 2) {
  stop(
    "Only ", n_studies_main,
    " non-missing study identifier(s) found in the age <= 40 dataset. ",
    "Check that 'idarticleh' has been preserved as character data."
  )
}

# Each study receives the same aggregate weight. Records without a study ID
# receive NA and are omitted only from the study-balanced fit.
df_main <- df_main %>%
  group_by(idarticleh) %>%
  mutate(
    study_n = ifelse(is.na(idarticleh), NA_integer_, dplyr::n()),
    study_weight_raw = ifelse(is.na(idarticleh), NA_real_, 1 / study_n)
  ) %>%
  ungroup() %>%
  mutate(
    study_weight = study_weight_raw / mean(study_weight_raw, na.rm = TRUE)
  )

cat("Height studies in modelling subset:", n_studies_main, "\n")
cat("Study-weight values:", dplyr::n_distinct(df_main$study_weight, na.rm = TRUE), "\n\n")

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
km <- kmeans(scale(early_summary$mean_h_early), centers = k_classes, nstart = kmeans_nstart)
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
                   maxiter = 1024) {
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

fit_power <- function(df, start = list(alpha = 6, beta = 0.5), maxiter = 1024) {
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
# 4b) Focused reviewer-analysis helpers
# ---------------------------
extract_param_ci <- function(m, fit_label, level = 0.95) {
  if (is.null(m)) return(tibble(Fit = fit_label, parameter = NA_character_, estimate = NA_real_, SE = NA_real_, lower = NA_real_, upper = NA_real_))
  cc <- coef(m); se <- sqrt(diag(vcov(m))); z <- qnorm(1 - (1 - level) / 2)
  tibble(Fit = fit_label, parameter = names(cc), estimate = unname(cc), SE = unname(se),
         lower = estimate - z * SE, upper = estimate + z * SE)
}

prediction_metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred); obs <- obs[ok]; pred <- pred[ok]
  if (!length(obs)) return(tibble(n = 0L, RMSE = NA_real_, MAE = NA_real_, Bias = NA_real_))
  err <- pred - obs
  tibble(n = length(obs), RMSE = sqrt(mean(err^2)), MAE = mean(abs(err)), Bias = mean(err))
}

loso_model <- function(df, study_col, response_col, fit_fun, age_breaks, age_labels) {
  df_loso <- df %>%
    dplyr::filter(
      !is.na(.data[[study_col]]),
      !is.na(.data[[response_col]]),
      is.finite(agecom)
    ) %>%
    dplyr::mutate(
      .study = as.character(.data[[study_col]]),
      .observed = as.numeric(.data[[response_col]])
    )

  studies <- sort(unique(df_loso$.study))
  if (length(studies) < 2) {
    stop("LOSO requires at least two non-missing studies; found ", length(studies), ".")
  }

  pred_rows <- purrr::map_dfr(studies, function(st) {
    train <- df_loso %>% dplyr::filter(.study != st)
    test  <- df_loso %>% dplyr::filter(.study == st)

    model <- fit_fun(train)

    if (is.null(model)) {
      return(tibble::tibble(
        study = st,
        agecom = test$agecom,
        observed = test$.observed,
        predicted = NA_real_,
        converged = FALSE
      ))
    }

    predictions <- tryCatch(
      as.numeric(stats::predict(model, newdata = test)),
      error = function(e) rep(NA_real_, nrow(test))
    )

    tibble::tibble(
      study = st,
      agecom = test$agecom,
      observed = test$.observed,
      predicted = predictions,
      converged = all(is.finite(predictions))
    )
  })

  if (nrow(pred_rows) == 0) {
    return(list(
      predictions = tibble::tibble(),
      by_study = tibble::tibble(),
      summary = tibble::tibble()
    ))
  }

  by_study <- pred_rows %>%
    dplyr::group_by(study) %>%
    dplyr::summarise(
      n = sum(is.finite(observed) & is.finite(predicted)),
      RMSE = {
        ok <- is.finite(observed) & is.finite(predicted)
        if (any(ok)) sqrt(mean((predicted[ok] - observed[ok])^2)) else NA_real_
      },
      MAE = {
        ok <- is.finite(observed) & is.finite(predicted)
        if (any(ok)) mean(abs(predicted[ok] - observed[ok])) else NA_real_
      },
      Bias = {
        ok <- is.finite(observed) & is.finite(predicted)
        if (any(ok)) mean(predicted[ok] - observed[ok]) else NA_real_
      },
      converged = all(converged),
      .groups = "drop"
    )

  overall <- prediction_metrics(pred_rows$observed, pred_rows$predicted) %>%
    dplyr::mutate(scope = "All ages")

  by_age <- pred_rows %>%
    dplyr::mutate(
      age_band = cut(
        agecom,
        breaks = age_breaks,
        labels = age_labels,
        right = TRUE
      )
    ) %>%
    dplyr::filter(!is.na(age_band)) %>%
    dplyr::group_by(age_band) %>%
    dplyr::summarise(
      n = sum(is.finite(observed) & is.finite(predicted)),
      RMSE = {
        ok <- is.finite(observed) & is.finite(predicted)
        if (any(ok)) sqrt(mean((predicted[ok] - observed[ok])^2)) else NA_real_
      },
      MAE = {
        ok <- is.finite(observed) & is.finite(predicted)
        if (any(ok)) mean(abs(predicted[ok] - observed[ok])) else NA_real_
      },
      Bias = {
        ok <- is.finite(observed) & is.finite(predicted)
        if (any(ok)) mean(predicted[ok] - observed[ok]) else NA_real_
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(scope = as.character(age_band)) %>%
    dplyr::select(-age_band)

  list(
    predictions = pred_rows,
    by_study = by_study,
    summary = dplyr::bind_rows(overall, by_age)
  )
}

adjusted_rand_index <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 2) return(NA_real_)
  tab <- table(x, y); choose2 <- function(z) z * (z - 1) / 2
  sum_nij <- sum(choose2(tab)); sum_ai <- sum(choose2(rowSums(tab))); sum_bj <- sum(choose2(colSums(tab))); total <- choose2(n)
  expected <- sum_ai * sum_bj / total; max_index <- 0.5 * (sum_ai + sum_bj)
  if (max_index == expected) return(1)
  (sum_nij - expected) / (max_index - expected)
}

classify_height_taxa <- function(df, min_total, min_early, early_max) {
  es <- df %>% filter(!is.na(taxon), taxon != "") %>% group_by(taxon) %>%
    summarise(n_tot = n(), n_early = sum(agecom <= early_max),
              mean_early = mean(hcom[agecom <= early_max], na.rm = TRUE), .groups = "drop") %>%
    filter(n_tot >= min_total, n_early >= min_early, is.finite(mean_early))
  if (nrow(es) < 2) return(tibble(taxon = character(), class = character(), mean_early = numeric(), n_tot = integer(), n_early = integer()))
  set.seed(42); kk <- kmeans(scale(es$mean_early), centers = 2, nstart = kmeans_nstart); es$cluster <- kk$cluster
  lab <- es %>% group_by(cluster) %>% summarise(mu = mean(mean_early), .groups = "drop") %>%
    arrange(mu) %>% mutate(class = c("baseline", "fast_juvenile"))
  es %>% left_join(lab, by = "cluster") %>% select(taxon, class, mean_early, n_tot, n_early)
}



# Reporting-level sensitivity helper for the baseline height model.
reporting_level_sensitivity_height <- function(df) {
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
      tibble(Scenario = nm, n = nrow(dd), n_studies = n_distinct(dd$idarticleh), levels_included = paste(sort(unique(dd$level_chr)), collapse = ", ")),
      extract_fit(m, nm) %>% select(-Model)
    )
  })

  parameters <- purrr::imap_dfr(fits, function(m, nm) extract_param_ci(m, nm))

  age_seq <- seq(0.5, age_max_main, length.out = 260)
  predictions <- purrr::imap_dfr(fits, function(m, nm) {
    if (is.null(m)) return(NULL)
    tibble(Scenario = nm, agecom = age_seq, predicted = as.numeric(predict(m, newdata = data.frame(agecom = age_seq))))
  })

  list(summary = summary, parameters = parameters, predictions = predictions)
}

# ---------------------------
# 5) Fits (ALL within age ≤ 40)
# ---------------------------
m_base_unw <- fit_cr(df_baseline, weights_col = NULL)
m_base_cap <- fit_cr(df_baseline, weights_col = "weight_fit")
m_base_wgt <- fit_cr(df_baseline, weights_col = "weight")
m_base_study <- fit_cr(df_baseline %>% dplyr::filter(is.finite(study_weight)), weights_col = "study_weight")

# within-window sensitivity (≤40): baseline vs including fast-juvenile
m_all_unw  <- fit_cr(df_main, weights_col = NULL)

m_fast_pow <- fit_power(df_fast_early)

tab_baseline <- bind_rows(
  extract_fit(m_base_unw, "Baseline CR – Unweighted (paper)"),
  extract_fit(m_base_cap, "Baseline CR – Capped-weighted (p99)"),
  extract_fit(m_base_wgt, "Baseline CR – Weighted"),
  extract_fit(m_base_study, "Baseline CR – Study-balanced")
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
# 5b) Focused reviewer analyses
# ---------------------------
param_ci_height <- bind_rows(
  extract_param_ci(m_base_unw, "Baseline unweighted"),
  extract_param_ci(m_base_cap, "Baseline capped-weighted (p99)"),
  extract_param_ci(m_base_wgt, "Baseline weighted"),
  extract_param_ci(m_base_study, "Baseline study-balanced"),
  extract_param_ci(m_all_unw, "All records unweighted"),
  extract_param_ci(m_fast_pow, "Fast-juvenile power")
)

weight_diagnostics <- df_main %>% summarise(
  n = n(),
  n_studies = n_distinct(idarticleh, na.rm = TRUE),
  n_missing_study = sum(is.na(idarticleh)),
  min = min(weight, na.rm = TRUE),
  p50 = quantile(weight, .50, na.rm = TRUE), p90 = quantile(weight, .90, na.rm = TRUE),
  p95 = quantile(weight, .95, na.rm = TRUE), p99 = quantile(weight, .99, na.rm = TRUE),
  max = max(weight, na.rm = TRUE), n_capped = sum(weight > w_cap, na.rm = TRUE)
)

study_dominance <- df_main %>%
  dplyr::filter(!is.na(idarticleh)) %>%
  count(idarticleh, sort = TRUE, name = "n_records") %>%
  mutate(
    share = n_records / sum(n_records),
    cumulative_share = cumsum(share),
    rank = row_number()
  )

# Conditional LOSO: uses the manuscript's fixed fast-juvenile labels, then leaves out entire studies.
# This evaluates the fitted trajectories without record-level leakage; classification sensitivity is assessed separately below.
if (run_loso_cv) {
  loso_base <- loso_model(df_baseline, "idarticleh", "hcom", function(d) fit_cr(d, weights_col = NULL), cv_age_breaks, cv_age_labels)
  loso_fast <- loso_model(df_fast_early, "idarticleh", "hcom", function(d) fit_power(d), c(-Inf, 3, 5, 8, Inf), c("1-3", "4-5", "6-8", "9-10"))
} else {
  loso_base <- loso_fast <- list(summary = tibble(), by_study = tibble(), predictions = tibble())
}

# Baseline composition sensitivity: only explicitly classified baseline taxa versus current operational baseline.
explicit_baseline_taxa <- early_summary %>% filter(group == "baseline") %>% pull(taxon)
df_baseline_explicit <- df_main %>% filter(taxon %in% explicit_baseline_taxa)
m_base_explicit <- fit_cr(df_baseline_explicit, weights_col = NULL)
baseline_composition <- bind_rows(
  extract_fit(m_base_unw, "Operational baseline (includes ineligible taxa)"),
  extract_fit(m_base_explicit, "Explicitly classified baseline taxa only")
) %>% mutate(n = c(nrow(df_baseline), nrow(df_baseline_explicit)))

if (run_class_sens) {
  ref_class <- classify_height_taxa(df_main, min_total_obs, min_early_obs, A_fast_cluster)
  sens_grid <- tidyr::crossing(min_total = c(3, 5, 10), min_early = c(3, 5, 8), early_max = c(10, 15, 20))
  class_sensitivity <- purrr::pmap_dfr(sens_grid, function(min_total, min_early, early_max) {
    z <- classify_height_taxa(df_main, min_total, min_early, early_max)
    common <- inner_join(ref_class %>% select(taxon, ref = class), z %>% select(taxon, alt = class), by = "taxon")
    tibble(min_total = min_total, min_early = min_early, early_max = early_max,
           n_eligible = nrow(z), n_common = nrow(common),
           agreement = ifelse(nrow(common), mean(common$ref == common$alt), NA_real_),
           ARI = ifelse(nrow(common) > 1, adjusted_rand_index(common$ref, common$alt), NA_real_),
           n_baseline = sum(z$class == "baseline"), n_fast = sum(z$class == "fast_juvenile"))
  })

  fast_fit_sensitivity <- purrr::map_dfr(c(8, 10, 12), function(max_age) {
    dd <- df_main %>% filter(group == "fast_juvenile", agecom >= age_min_fastfit, agecom <= max_age)
    mm <- fit_power(dd)
    extract_fit(mm, paste0("Fast power ages 1-", max_age)) %>% mutate(age_max = max_age, n = nrow(dd))
  })

} else {
  class_sensitivity <- fast_fit_sensitivity <- tibble()
}

# Reporting-level sensitivity is evaluated for the operational baseline,
# because this is the principal height curve used in the paper.
if (run_level_sens) {
  level_sens_height <- reporting_level_sensitivity_height(df_baseline)
} else {
  level_sens_height <- list(summary = tibble(), parameters = tibble(), predictions = tibble())
}

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



if (run_level_sens && nrow(level_sens_height$predictions) > 0) {
  p_level_height <- ggplot(level_sens_height$predictions,
                           aes(x = agecom, y = predicted, linetype = Scenario)) +
    geom_line(linewidth = 1) +
    coord_cartesian(xlim = c(0, age_max_main)) +
    theme_minimal() +
    labs(title = "Height reporting-level sensitivity",
         subtitle = "Operational baseline model fitted under alternative reporting-level subsets",
         x = "Age (years)", y = "Predicted height (m)", linetype = "Scenario")
  print(p_level_height)
  ggsave(file.path(fig_dir, "H_ReportingLevel_Sensitivity.pdf"), p_level_height, width = 8, height = 5.5)
}

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

wb <- createWorkbook()

write_sheet(wb, paste0("H_Baseline_CR_A", age_max_main), tab_baseline)
write_sheet(wb, paste0("H_FastJuvenile_Power_A", age_min_fastfit, "_", A_fast_fit), tab_fast)
write_sheet(wb, paste0("H_Within40_Comparison_A", age_max_main), tab_within40_sens)

write_sheet(wb, paste0("H_Group_Counts_A", age_max_main),
            df_main %>% count(group) %>% mutate(pct = n / sum(n)))

write_sheet(wb, paste0("H_FastJuvenile_Taxa_A", age_max_main), fast_taxa)

write_sheet(wb, paste0("H_FastJuvenile_TaxonCounts_A", age_max_main),
            df_main %>% dplyr::filter(group == "fast_juvenile") %>% count(taxon, sort = TRUE))

write_sheet(wb, "H_Param_CI", param_ci_height)
write_sheet(wb, "H_Weight_Diag", weight_diagnostics)
write_sheet(wb, "H_Study_Dominance", study_dominance)
write_sheet(wb, "H_LOSO_Base_Summary", loso_base$summary)
write_sheet(wb, "H_LOSO_Base_ByStudy", loso_base$by_study)
write_sheet(wb, "H_LOSO_Base_Pred", loso_base$predictions)
write_sheet(wb, "H_LOSO_Fast_Summary", loso_fast$summary)
write_sheet(wb, "H_LOSO_Fast_ByStudy", loso_fast$by_study)
write_sheet(wb, "H_LOSO_Fast_Pred", loso_fast$predictions)
write_sheet(wb, "H_Baseline_Composition", baseline_composition)
write_sheet(wb, "H_Class_Sensitivity", class_sensitivity)
write_sheet(wb, "H_FastFit_Sensitivity", fast_fit_sensitivity)
write_sheet(wb, "H_Level_Sens", level_sens_height$summary)
write_sheet(wb, "H_Level_Sens_Param", level_sens_height$parameters)
write_sheet(wb, "H_Level_Sens_Pred", level_sens_height$predictions)

saveWorkbook(wb, results_path, overwrite = TRUE)
message("Tables written to: ", results_path)
