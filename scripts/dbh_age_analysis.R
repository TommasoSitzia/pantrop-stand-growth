# ============================================================
# DBH–Age analysis (expects a clean plantation-only input dataset)
#
# Core analysis (PAPER):
#   - Modelling window: age ≤ 40 years
#   - Models: Chapman–Richards (NLS; 95% CI for mean curve)
#   - Fits compared (≤40 only): unweighted (paper), capped-weighted (p99), weighted
#   - Growth classes (≤40 only): k-means on taxon mean DBH at ages ≤ 15
#
# Repository paths:
#   - Input data:  data/dbhage.xlsx
#   - Expected worksheet: data
#   - Figures:     outputs/PDF/DBH/
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

# ----------------------------- #
# PATHS (repository-relative)
# ----------------------------- #
data_path    <- here("data", "dbhage.xlsx")
results_path <- here("outputs", "xlsx", "resultspaper.xlsx")
fig_dir      <- here("outputs", "PDF", "DBH")

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

# intervals
do_PI_supp     <- FALSE  # PI not reported in paper

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
  filter(!is.na(agecom), !is.na(dbhcom)) %>%
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
df_main <- dat %>% filter(agecom <= age_max_main)

w_cap <- quantile(df_main$weight, 0.99, na.rm = TRUE)
df_main <- df_main %>% mutate(weight_fit = pmin(weight, w_cap))

# Sanity checks (these should match your paper counts if your sheet is the same one used for the manuscript)
cat("\nDBH dataset (all ages): n =", nrow(dat), " max age =", max(dat$agecom, na.rm=TRUE), "\n")
cat("DBH modelling subset (age ≤ 40): n =", nrow(df_main), " max age =", max(df_main$agecom, na.rm=TRUE), "\n\n")

age_seq_main <- seq(0.5, age_max_main, length.out = 220)

# ----------------------------- #
# 3) Chapman–Richards helpers
# ----------------------------- #
fit_cr <- function(df, weights_col = NULL,
                   start = list(a = 40, b = 0.05, c = 1.2),
                   maxiter = 2000) {
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
# 4) Pooled models (age ≤ 40 ONLY)
# ----------------------------- #
m_unw <- fit_cr(df_main, weights_col = NULL)
m_cap <- fit_cr(df_main, weights_col = "weight_fit")
m_wgt <- fit_cr(df_main, weights_col = "weight")

m_paper <- m_unw
paper_label <- "Unweighted"

tab_pooled <- bind_rows(
  extract_fit(m_unw) %>% mutate(Fit = "Unweighted"),
  extract_fit(m_cap) %>% mutate(Fit = "Capped-weighted (p99)"),
  extract_fit(m_wgt) %>% mutate(Fit = "Weighted")
) %>%
  mutate(
    Dataset = paste0("Age ≤ ", age_max_main),
    Model = "Chapman–Richards: DBH(t)=a*(1-exp(-b*t))^c"
  ) %>%
  select(Dataset, Model, Fit, AIC, RSE, a, b, c, converged)

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
  filter(n_tot >= min_total_obs) %>%
  pull(taxon)

early_summary <- df_main %>%
  filter(taxon %in% taxa_ok, agecom <= early_age_max) %>%
  group_by(taxon) %>%
  summarise(
    mean_dbh_early = mean(dbhcom, na.rm = TRUE),
    n_early = n(),
    .groups = "drop"
  ) %>%
  filter(n_early >= min_early_obs)

if (nrow(early_summary) < k_classes) {
  stop("Not enough eligible taxa for k-means. Lower thresholds or check data.")
}

set.seed(42)
km <- kmeans(scale(early_summary$mean_dbh_early), centers = k_classes, nstart = 25)
early_summary$cluster <- km$cluster

cluster_map <- early_summary %>%
  group_by(cluster) %>%
  summarise(mu = mean(mean_dbh_early), .groups = "drop") %>%
  arrange(mu) %>%
  mutate(growth_class = c("slow", "medium", "fast"))

early_summary <- early_summary %>% left_join(cluster_map, by = "cluster")

df_main_classed <- df_main %>%
  left_join(early_summary %>% select(taxon, growth_class), by = "taxon") %>%
  mutate(
    growth_class = ifelse(is.na(growth_class), "unclassified", growth_class),
    growth_class = factor(growth_class, levels = c("slow","medium","fast","unclassified"))
  )

df_fitclass <- df_main_classed %>%
  filter(growth_class %in% c("slow","medium","fast")) %>%
  droplevels()

counts_class <- df_fitclass %>%
  count(growth_class) %>%
  mutate(pct = n / sum(n))

print(counts_class)

models_class <- df_fitclass %>%
  group_split(growth_class) %>%
  setNames(levels(df_fitclass$growth_class)) %>%
  map(~ fit_cr(.x, weights_col = NULL, start = list(a = 40, b = 0.05, c = 1.0)))

tab_class <- imap_dfr(models_class, function(m, cls) {
  extract_fit(m) %>% mutate(growth_class = cls)
}) %>%
  mutate(
    Dataset = paste0("Age ≤ ", age_max_main),
    Fit = "Unweighted",
    Model = "Chapman–Richards by growth class"
  ) %>%
  select(Dataset, Model, growth_class, Fit, AIC, RSE, a, b, c, converged)

print(tab_class)

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
    dd <- dd %>% filter(Fit == paper_label)
  }
  
  base +
    geom_ribbon(
      data = dd %>% filter(Fit == paper_label),
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
  base <- ggplot(df_fitclass, aes(x = agecom, y = dbhcom, colour = growth_class)) +
    geom_point(alpha = 0.20, size = 0.8) +
    theme_minimal() +
    labs(
      title = paste0("DBH–age by growth class (age ≤ ", age_max_main, ")"),
      x = "Age (years)", y = "DBH (cm)", colour = "Growth class"
    ) +
    coord_cartesian(xlim = c(0, age_max_main))
  
  pred_df <- imap_dfr(models_class, function(m, cls) {
    if (is.null(m)) return(NULL)
    predict_ci_cr(m, age_seq_main) %>% mutate(growth_class = cls)
  })
  
  base +
    geom_ribbon(
      data = pred_df,
      aes(x = agecom, ymin = lo, ymax = hi, fill = growth_class),
      inherit.aes = FALSE,
      alpha = 0.12,
      colour = NA,
      show.legend = FALSE
    ) +
    geom_line(
      data = pred_df,
      aes(x = agecom, y = fit, colour = growth_class),
      inherit.aes = FALSE,
      linewidth = 1.1
    )
}

p_main_pooled  <- plot_pooled_CI(show_all_fits = TRUE)
p_main_classes <- plot_classes_CI()

print(p_main_pooled)
print(p_main_classes)

ggsave(file.path(fig_dir, paste0("DBH_Fig_Main_Pooled_CI_A", age_max_main, ".pdf")),
       p_main_pooled, width = 7.0, height = 5.0)

ggsave(file.path(fig_dir, paste0("DBH_Fig_Main_Classes_CI_A", age_max_main, ".pdf")),
       p_main_classes, width = 7.0, height = 5.0)

# ----------------------------- #
# 7) Diagnostics PDFs
# ----------------------------- #
save_diag_pdf(df_main, m_paper,
              out_pdf = file.path(fig_dir, paste0("DBH_Diag_Pooled_", paper_label, "_Amax", age_max_main, ".pdf")),
              title_prefix = paste0("Pooled ", paper_label, " (Amax ", age_max_main, "): "))

invisible(imap(models_class, function(m, cls) {
  d_cls <- df_fitclass %>% filter(growth_class == cls)
  save_diag_pdf(d_cls, m,
                out_pdf = file.path(fig_dir, paste0("DBH_Diag_Class_", cls, "_Amax", age_max_main, ".pdf")),
                title_prefix = paste0("Class ", cls, " (Amax ", age_max_main, "): "))
}))

message("Figures + diagnostics saved to: ", fig_dir)

# ----------------------------- #
# 8) Write outputs to Excel
# ----------------------------- #
write_sheet <- function(wb, name, df) {
  if (name %in% names(wb)) removeWorksheet(wb, name)
  addWorksheet(wb, name)
  writeData(wb, name, df, withFilter = TRUE)
}

wb <- if (file.exists(results_path)) loadWorkbook(results_path) else createWorkbook()

write_sheet(wb, paste0("DBH_Pooled_A", age_max_main), tab_pooled)
write_sheet(wb, paste0("DBH_GrowthClass_A", age_max_main), tab_class)
write_sheet(wb, paste0("DBH_GrowthClass_Counts_A", age_max_main), counts_class)

write_sheet(wb, "DBH_GrowthClass_Taxa",
            early_summary %>% arrange(growth_class, desc(mean_dbh_early)))

saveWorkbook(wb, results_path, overwrite = TRUE)
message("Tables written to: ", results_path)
