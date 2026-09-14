# =============================================================================
# make_plots.R
# =============================================================================
# Redesigned plotting script — cleaner, easier to read.
#
# KEY IMPROVEMENTS over the original make_plots.R:
#   1. Causal and non-causal SNPs always plotted SEPARATELY — no zero spike
#   2. Bias plot split by causal/non-causal, with clearer x-axis limits
#   3. Overlapping density ("hist_all") replaced with a ridge plot so each
#      estimator has its own lane
#   4. Coverage and MSE annotated directly on the estimate density panels
#   5. Forest plot shows causal and non-causal side by side
#   6. CI width distribution replaces the confusing CI bounds overlay
#
# HOW TO RUN (from the AnchorGene project root):
#
#   Rscript make_plots.R
#
#
# Optional: run aggregate_results.R first if the CSVs are stale.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggridges)   
})

TREATMENT_SNP <- Sys.getenv("TREATMENT_SNP", "rs0254")
OUT_DIR <- file.path(Sys.getenv("RES_DIR", "results"), TREATMENT_SNP)

RES_DIR  <- Sys.getenv("RES_DIR", "results")
PLOT_DIR <- file.path(Sys.getenv("PLOT_DIR", "plots"), TREATMENT_SNP)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Consistent colour palette — one colour per estimator
ALG_COLOURS <- c(
  "TMLE"            = "#B279A2",
  "C-TMLE0"         = "#E45756",
  "C-TMLE1"         = "#54A24B",
  "Scalable-CTMLE"  = "#4C78A8",
  "Wang-TMLE"       = "#534949",
  "Wang-CTMLE"      = "#72B7B2",
  "C-TMLE1-OAL"     = "#7c2ba1"
)

# Estimator order: group the well-performing ones together
ALG_ORDER <- c("TMLE","C-TMLE0","C-TMLE1","Scalable-CTMLE", "Wang-TMLE", "Wang-CTMLE", "C-TMLE1-OAL")

ALLOWED_ALGORITHMS <- ALG_ORDER

base_theme <- function() {
  theme_minimal(base_size = 12, base_family = "sans") +
    theme(
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.background   = element_rect(fill = "white", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = "#E5E7EB", linewidth = 0.35),
      strip.background   = element_rect(fill = "#F3F4F6", colour = NA),
      strip.text         = element_text(face = "bold", size = 11, colour = "#111827"),
      plot.title         = element_text(face = "bold", size = 14, hjust = 0.5, colour = "#111827"),
      plot.subtitle      = element_text(size = 9.2, hjust = 0.5, colour = "#4B5563", lineheight = 0.95),
      axis.title         = element_text(face = "bold", colour = "#111827"),
      axis.text          = element_text(colour = "#374151"),
      legend.position    = "none"
    )
}

wrap_subtitle <- function(text, width = 115) {
  paste(strwrap(text, width = width), collapse = "\n")
}

save_plot <- function(p, path, w = 14, h = 9) {
  tryCatch(
    ggsave(path, p, width = w, height = h, dpi = 150, bg = "white"),
    error = function(e) message("WARNING: could not save ", path, " — ", e$message)
  )
  message("  → ", path)
}

# =============================================================================
# Helper: compute per-estimator metrics for annotation.
# Groups by IS_CAUSAL only when both TRUE and FALSE rows exist (sim1 style).
# In sim2 every row has IS_CAUSAL = TRUE, so grouping is by Algorithm alone.
# =============================================================================
compute_metrics <- function(df) {
  has_noncausal <- any(df$IS_CAUSAL == FALSE, na.rm = TRUE)
  grp_vars <- if (has_noncausal) c("Algorithm", "IS_CAUSAL") else "Algorithm"

  df %>%
    mutate(
      covered  = (true_ATE >= CI_Lower) & (true_ATE <= CI_Upper),
      sq_err   = (Estimate - true_ATE)^2,
      CI_width = CI_Upper - CI_Lower,
      # Ensure IS_CAUSAL exists as a column even if not used for grouping
      IS_CAUSAL = IS_CAUSAL
    ) %>%
    group_by(across(all_of(grp_vars))) %>%
    summarise(
      Coverage = mean(covered,   na.rm = TRUE),
      MSE      = mean(sq_err,    na.rm = TRUE),
      MeanBias = mean(Bias,      na.rm = TRUE),
      MeanCIw  = mean(CI_width,  na.rm = TRUE),
      N        = sum(!is.na(Estimate)),
      IS_CAUSAL = IS_CAUSAL[1],   # carry through for downstream filters
      .groups  = "drop"
    ) %>%
    mutate(
      label_cov  = sprintf("Cov: %.0f%%", Coverage * 100),
      label_mse  = sprintf("MSE: %.3f",   MSE),
      label_bias = sprintf("Bias: %+.3f", MeanBias)
    )
}

# Helper: extract replicate and SNP counts for subtitle labels
# =============================================================================
count_label <- function(df, causal_only = FALSE) {
  sub    <- if (causal_only) df[df$IS_CAUSAL == TRUE, ] else df
  n_reps <- length(unique(sub$bootstrap_id))
  n_snps <- length(unique(sub$SNP))
  if (causal_only) {
    sprintf("mean over %d replicates x %d causal SNPs", n_reps, n_snps)
  } else {
    sprintf("mean over %d replicates x %d SNPs", n_reps, n_snps)
  }
}

# Helper: build a scenario label from metadata columns
# =============================================================================
scenario_label <- function(df) {
  pick_one <- function(colname) {
    if (!(colname %in% colnames(df))) return(NA_character_)
    vals <- unique(df[[colname]])
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NA_character_)
    as.character(vals[1])
  }

  beta          <- pick_one("beta")
  fvu           <- pick_one("fvu")
  treatment_snp <- pick_one("treatment_snp")
  run_tag       <- pick_one("run_tag")

  parts <- c()
  if (!is.na(treatment_snp)) parts <- c(parts, paste0("treatment_snp=", treatment_snp))
  if (!is.na(beta))          parts <- c(parts, paste0("beta=",          beta))
  if (!is.na(fvu))           parts <- c(parts, paste0("fvu=",           fvu))
  if (!is.na(run_tag) && nzchar(run_tag)) parts <- c(parts, paste0("run_tag=", run_tag))

  if (length(parts) == 0) return("Scenario: metadata unavailable")
  paste("Scenario:", paste(parts, collapse = " | "))
}

# Helper: drop run_tag from a scenario label string (for compact subtitles)
# =============================================================================
scenario_brief <- function(scen_label) {
  pieces <- unlist(strsplit(scen_label, "\\s*\\|\\s*"))
  pieces <- pieces[!grepl("^run_tag=", pieces)]
  paste(pieces, collapse = " | ")
}

# =============================================================================
# PLOT 1: Estimate density — causal SNPs only, one panel per estimator.
#         Annotated with Coverage and MSE.
#         Vertical dashed line at 0 (no effect) AND the mean true ATE.
# =============================================================================
plot_estimates_causal <- function(df, metrics, sim_label, prefix, scen_label) {
  # Clip per estimator to 2nd–98th percentile so outliers don't flatten panels
  res <- df %>%
    filter(IS_CAUSAL == TRUE, !is.na(Estimate)) %>%
    group_by(Algorithm) %>%
    mutate(
      lo       = quantile(Estimate, 0.02, na.rm = TRUE),
      hi       = quantile(Estimate, 0.98, na.rm = TRUE),
      Estimate = ifelse(Estimate < lo | Estimate > hi, NA_real_, Estimate)
    ) %>%
    ungroup() %>%
    filter(!is.na(Estimate)) %>%
    mutate(Algorithm = factor(Algorithm, levels = ALG_ORDER))

  met <- metrics %>%
    filter(IS_CAUSAL == TRUE) %>%
    mutate(Algorithm = factor(Algorithm, levels = ALG_ORDER))

  # Mean true ATE across causal SNPs (non-zero rows only — non-causal rows
  # would drag this toward 0)
  mean_true <- df %>%
    filter(IS_CAUSAL == TRUE, true_ATE != 0) %>%
    group_by(Algorithm) %>%
    summarise(mt = mean(true_ATE, na.rm = TRUE), .groups = "drop") %>%
    mutate(Algorithm = factor(Algorithm, levels = ALG_ORDER))

  p <- ggplot(res, aes(x = Estimate, fill = Algorithm)) +
    geom_density(alpha = 0.6, colour = "white", na.rm = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "black", linewidth = 0.5) +
    geom_vline(data = mean_true, aes(xintercept = mt),
               linetype = "solid", colour = "black", linewidth = 0.7) +
    geom_text(data = met,
              aes(x = Inf, y = Inf,
                  label = paste0(label_cov, "\n", label_mse)),
              hjust = 1.05, vjust = 1.3, size = 3.2, colour = "grey20",
              inherit.aes = FALSE) +
    scale_fill_manual(values = ALG_COLOURS) +
    # free_x: each panel uses its own x range — essential when estimators
    # have very different spread
    facet_wrap(~ Algorithm, nrow = 2, scales = "free") +
    labs(
      title    = paste("ATE Estimates at CAUSAL SNPs —", sim_label),
      subtitle = wrap_subtitle(paste0(
        scenario_brief(scen_label), "  |  ",
        count_label(df, causal_only = TRUE), "  |  ",
        "Dashed = 0  |  Solid = mean true ATE  |  ",
        "Coverage/MSE annotated  |  Clipped to 2–98th pct"
      )),
      x = "Estimated ATE", y = "Density"
    ) +
    base_theme()

  save_plot(p, file.path(PLOT_DIR, paste0(prefix, "_est_causal.png")), w = 16, h = 8)
}

# =============================================================================
# PLOT 2: Estimate density — non-causal SNPs only.
#         True ATE = 0 for all non-causal SNPs, so tight concentration around
#         0 indicates a well-calibrated estimator.
# =============================================================================
plot_estimates_noncausal <- function(df, metrics, sim_label, prefix, scen_label) {
  res <- df %>%
    filter(IS_CAUSAL == FALSE, !is.na(Estimate)) %>%
    group_by(Algorithm) %>%
    mutate(
      lo       = quantile(Estimate, 0.02, na.rm = TRUE),
      hi       = quantile(Estimate, 0.98, na.rm = TRUE),
      Estimate = ifelse(Estimate < lo | Estimate > hi, NA_real_, Estimate)
    ) %>%
    ungroup() %>%
    filter(!is.na(Estimate)) %>%
    mutate(Algorithm = factor(Algorithm, levels = ALG_ORDER))

  met <- metrics %>%
    filter(IS_CAUSAL == FALSE) %>%
    mutate(Algorithm = factor(Algorithm, levels = ALG_ORDER))

  p <- ggplot(res, aes(x = Estimate, fill = Algorithm)) +
    geom_density(alpha = 0.6, colour = "white", na.rm = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "black", linewidth = 0.6) +
    geom_text(data = met,
              aes(x = Inf, y = Inf,
                  label = paste0(label_cov, "\n", label_mse)),
              hjust = 1.05, vjust = 1.3, size = 3.2, colour = "grey20",
              inherit.aes = FALSE) +
    scale_fill_manual(values = ALG_COLOURS) +
    facet_wrap(~ Algorithm, nrow = 2, scales = "free") +
    labs(
      title    = paste("ATE Estimates at NON-CAUSAL SNPs —", sim_label),
      subtitle = wrap_subtitle(paste0(
        scenario_brief(scen_label), "  |  ",
        count_label(df, causal_only = FALSE), "  |  ",
        "True ATE = 0 for all  |  Dashed = 0  |  ",
        "Tight distribution = better  |  Clipped to 2–98th pct"
      )),
      x = "Estimated ATE", y = "Density"
    ) +
    base_theme()

  save_plot(p, file.path(PLOT_DIR, paste0(prefix, "_est_noncausal.png")), w = 16, h = 8)
}

# =============================================================================
# PLOT 3: Bias density — CAUSAL SNPs only, one panel per estimator.
#         Non-causal SNPs are excluded: bias there equals the estimate
#         (true ATE = 0), making it redundant with Plot 2, and the point
#         mass at 0 prevents geom_density from rendering meaningfully.
# =============================================================================
plot_bias <- function(df, sim_label, prefix, scen_label) {
  res <- df %>%
    filter(IS_CAUSAL == TRUE, !is.na(Bias)) %>%
    group_by(Algorithm) %>%
    mutate(
      lo   = quantile(Bias, 0.02, na.rm = TRUE),
      hi   = quantile(Bias, 0.98, na.rm = TRUE),
      Bias = ifelse(Bias < lo | Bias > hi, NA_real_, Bias)
    ) %>%
    ungroup() %>%
    filter(!is.na(Bias)) %>%
    mutate(Algorithm = factor(Algorithm, levels = ALG_ORDER))

  p <- ggplot(res, aes(x = Bias, fill = Algorithm)) +
    geom_density(alpha = 0.6, colour = "white", na.rm = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "black", linewidth = 0.6) +
    scale_fill_manual(values = ALG_COLOURS) +
    facet_wrap(~ Algorithm, nrow = 2, scales = "free") +
    labs(
      title    = paste("Bias at CAUSAL SNPs —", sim_label),
      subtitle = wrap_subtitle(paste0(
        scenario_brief(scen_label), "  |  ",
        count_label(df, causal_only = TRUE), "  |  ",
        "Bias = Estimate − true ATE  |  Dashed = 0  |  ",
        "Centred = unbiased  |  Clipped to 2–98th pct"
      )),
      x = "Bias", y = "Density"
    ) +
    base_theme() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

  save_plot(p, file.path(PLOT_DIR, paste0(prefix, "_bias.png")), w = 16, h = 6)
}

# =============================================================================
# PLOT 4: Ridge plot — all estimators in one panel, causal SNPs only.
#         Each estimator gets its own horizontal lane. Much cleaner than
#         overlapping densities.
# =============================================================================
plot_ridge <- function(df, sim_label, prefix, scen_label) {
  res <- df %>%
    filter(IS_CAUSAL == TRUE, !is.na(Estimate)) %>%
    mutate(Algorithm = factor(Algorithm, levels = rev(ALG_ORDER)))  # rev so top = first

  xlim <- quantile(res$Estimate, c(0.01, 0.99), na.rm = TRUE)

  p <- ggplot(res, aes(x = Estimate, y = Algorithm, fill = Algorithm)) +
    geom_density_ridges(alpha = 0.7, colour = "white",
                        quantile_lines = TRUE, quantiles = 2,  # median line
                        na.rm = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "black", linewidth = 0.6) +
    scale_fill_manual(values = ALG_COLOURS) +
    scale_x_continuous(limits = xlim) +
    labs(
      title    = paste("ATE Estimate Distributions — CAUSAL SNPs —", sim_label),
      subtitle = wrap_subtitle(paste0(
        scenario_brief(scen_label), "  |  ",
        count_label(df, causal_only = TRUE), "  |  ",
        "Each estimator in its own lane  |  Median estimate shown  |  Dashed = 0"
      )),
      x = "Estimated ATE", y = NULL
    ) +
    base_theme()

  save_plot(p, file.path(PLOT_DIR, paste0(prefix, "_ridge_causal.png")), w = 10, h = 7)
}

# =============================================================================
# PLOT 5: Forest plot — mean bias ± 1.96 × SE of mean bias.
#
# Shows mean(Estimate − true_ATE) per estimator, split by causal/non-causal.
# Using bias rather than the raw estimate is appropriate here because
# averaging estimates over SNPs with different randomly-drawn true ATEs is
# not informative; bias is comparable across SNPs regardless of effect size.
# Reference line at 0 = unbiased. Ordered by |mean bias| on causal SNPs.
# =============================================================================
plot_forest_mean_bias <- function(df, sim_label, prefix, scen_label) {

  has_noncausal <- any(df$IS_CAUSAL == FALSE, na.rm = TRUE)
  n_reps        <- length(unique(df$bootstrap_id))
  n_causal      <- df %>% filter(IS_CAUSAL)  %>% pull(SNP) %>% n_distinct()
  n_noncausal   <- df %>% filter(!IS_CAUSAL) %>% pull(SNP) %>% n_distinct()

  summ <- df %>%
    filter(!is.na(Bias)) %>%
    mutate(SNP_type = ifelse(IS_CAUSAL, "Causal SNPs", "Non-causal SNPs")) %>%
    group_by(Algorithm, SNP_type) %>%
    summarise(
      Mean_Bias = mean(Bias, na.rm = TRUE),
      SE_Bias   = sd(Bias,   na.rm = TRUE) / sqrt(sum(!is.na(Bias))),
      N         = sum(!is.na(Bias)),
      .groups   = "drop"
    ) %>%
    mutate(
      bar_lo = Mean_Bias - 1.96 * SE_Bias,
      bar_hi = Mean_Bias + 1.96 * SE_Bias
    )

  alg_order <- summ %>%
    filter(SNP_type == "Causal SNPs") %>%
    arrange(desc(abs(Mean_Bias))) %>%
    pull(Algorithm)

  summ <- summ %>%
    mutate(Algorithm = factor(Algorithm, levels = alg_order))

  caption_txt <- if (has_noncausal) {
    paste0(
      "Mean over ", n_reps, " replicates \u00d7 ", n_causal,
      " causal SNPs (left) / ", n_noncausal, " non-causal SNPs (right)  |  ",
      "DGM: n_causal=5",
      if (grepl("Interaction", sim_label)) ", interaction order 1\u20134" else ", linear (order 1)"
    )
  } else {
    paste0(
      "Mean over ", n_reps, " replicates \u00d7 1 fixed treatment SNP  |  ",
      "True ATE = 1 by construction  |  DGM: n_causal=5",
      if (grepl("Interaction", sim_label)) ", interaction order 1\u20134" else ", linear (order 1)"
    )
  }

  p <- ggplot(summ, aes(x = Mean_Bias, y = Algorithm, colour = Algorithm)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = bar_lo, xmax = bar_hi),
                   height = 0.35, linewidth = 0.8) +
    geom_point(size = 3) +
    scale_colour_manual(values = ALG_COLOURS) +
    labs(
      title    = paste("Mean Bias \u2014", sim_label),
      subtitle = wrap_subtitle(paste0(
        scenario_brief(scen_label), "  |  ",
        "Point = mean bias (Estimate \u2212 true ATE)  |  ",
        "Bar = \u00b11.96 \u00d7 SE of mean bias  |  ",
        "Ordered by |mean bias|, smallest at top"
      )),
      x       = "Mean Bias (Estimate \u2212 true ATE)", y = NULL,
      caption = caption_txt
    ) +
    base_theme()

  # Only facet when there are genuinely two SNP types
  if (has_noncausal) {
    p <- p + facet_wrap(~ SNP_type, scales = "free_x")
  }

  save_plot(p, file.path(PLOT_DIR, paste0(prefix, "_forest_mean_bias.png")), w = 14, h = 6)
}

# =============================================================================
# PLOT 7: Forest plot — mean ATE estimate with CI of the mean.
#
# Point = mean(Estimate), Bar = ±1.96 × SE(mean Estimate), where
# SE(mean Estimate) = sd(Estimate) / sqrt(N).
# Reference line = target effect (causal) / 0 (non-causal).
# =============================================================================
plot_forest_mean_estimate <- function(df, sim_label, prefix, scen_label) {

  has_noncausal <- any(df$IS_CAUSAL == FALSE, na.rm = TRUE)
  n_reps        <- length(unique(df$bootstrap_id))
  n_causal      <- df %>% filter(IS_CAUSAL)  %>% pull(SNP) %>% n_distinct()
  n_noncausal   <- df %>% filter(!IS_CAUSAL) %>% pull(SNP) %>% n_distinct()

  pick_first_numeric <- function(x) {
    vals <- suppressWarnings(as.numeric(unique(x)))
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) return(NA_real_)
    vals[1]
  }

  target_effect <- NA_real_
  if ("beta" %in% colnames(df)) {
    target_effect <- pick_first_numeric(df$beta)
  }
  if (!is.finite(target_effect)) {
    target_effect <- df %>%
      filter(IS_CAUSAL == TRUE, !is.na(true_ATE)) %>%
      summarise(v = mean(true_ATE, na.rm = TRUE)) %>%
      pull(v)
  }

  # Reference line(s): sim2 gets one vertical line at the true ATE;
  # sim1 gets two panels so ref_lines needs both SNP_type levels.
  ref_lines <- if (has_noncausal) {
    data.frame(
      SNP_type = c("Causal SNPs", "Non-causal SNPs"),
      ref_x    = c(target_effect, 0)
    )
  } else {
    data.frame(SNP_type = "Causal SNPs", ref_x = target_effect)
  }

  summ <- df %>%
    filter(!is.na(Estimate)) %>%
    mutate(SNP_type = ifelse(IS_CAUSAL, "Causal SNPs", "Non-causal SNPs")) %>%
    group_by(Algorithm, SNP_type) %>%
    summarise(
      Mean_Estimate = mean(Estimate, na.rm = TRUE),
      SE_Mean       = sd(Estimate,   na.rm = TRUE) / sqrt(sum(!is.na(Estimate))),
      N             = sum(!is.na(Estimate)),
      .groups       = "drop"
    ) %>%
    mutate(
      CI_Low  = Mean_Estimate - 1.96 * SE_Mean,
      CI_High = Mean_Estimate + 1.96 * SE_Mean
    )

  alg_order <- summ %>%
    filter(SNP_type == "Causal SNPs") %>%
    arrange(desc(abs(Mean_Estimate))) %>%
    pull(Algorithm)

  summ <- summ %>%
    mutate(Algorithm = factor(Algorithm, levels = alg_order))

  caption_txt <- if (has_noncausal) {
    paste0(
      "Mean over ", n_reps, " replicates \u00d7 ", n_causal,
      " causal SNPs (left) / ", n_noncausal, " non-causal SNPs (right)  |  ",
      "Variance estimated via EIF"
    )
  } else {
    paste0(
      "Mean over ", n_reps, " replicates \u00d7 1 fixed treatment SNP  |  ",
      "Dashed = true ATE (= 1 by construction)  |  Variance estimated via EIF"
    )
  }

  p <- ggplot(summ, aes(x = Mean_Estimate, y = Algorithm, colour = Algorithm)) +
    geom_vline(data = ref_lines, aes(xintercept = ref_x),
               linetype = "dashed", colour = "grey50", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = CI_Low, xmax = CI_High),
                   height = 0.35, linewidth = 0.8) +
    geom_point(size = 3) +
    scale_colour_manual(values = ALG_COLOURS) +
    labs(
      title    = paste("Mean ATE Estimate \u2014", sim_label),
      subtitle = wrap_subtitle(paste0(
        scenario_brief(scen_label), "  |  ",
        "Point = mean estimate  |  Bar = \u00b11.96 \u00d7 SE(mean estimate)  |  ",
        "Dashed = true ATE  |  Ordered by |mean estimate|"
      )),
      x       = "Mean Estimated ATE", y = NULL,
      caption = caption_txt
    ) +
    base_theme()

  if (has_noncausal) {
    p <- p + facet_wrap(~ SNP_type, scales = "free_x")
  }

  save_plot(p, file.path(PLOT_DIR, paste0(prefix, "_forest_mean_estimate.png")), w = 14, h = 6)
}

# =============================================================================
# PLOT 8: CI width distribution — how wide are the confidence intervals?
#         Wide = uncertain. Narrow = confident (but may be overconfident).
#         Split by causal / non-causal, one panel per estimator.
# =============================================================================
plot_ci_width <- function(df, sim_label, prefix, scen_label) {
  res <- df %>%
    filter(!is.na(CI_Lower), !is.na(CI_Upper)) %>%
    mutate(
      CI_width  = CI_Upper - CI_Lower,
      Algorithm = factor(Algorithm, levels = ALG_ORDER),
      SNP_type  = ifelse(IS_CAUSAL, "Causal", "Non-causal")
    ) %>%
    filter(CI_width < quantile(CI_width, 0.99, na.rm = TRUE))  # trim extreme outliers

  p <- ggplot(res, aes(x = CI_width, fill = Algorithm, colour = SNP_type)) +
    geom_density(alpha = 0.5, linewidth = 0.6, na.rm = TRUE) +
    scale_fill_manual(values  = ALG_COLOURS) +
    scale_colour_manual(values = c("Causal" = "#C0392B", "Non-causal" = "#2980B9"),
                        name   = "SNP type") +
    facet_wrap(~ Algorithm, nrow = 2, scales = "free_y") +
    labs(
      title    = paste("CI Width Distribution —", sim_label),
      subtitle = wrap_subtitle(paste0(
        scenario_brief(scen_label), "  |  ",
        count_label(df), "  |  ",
        "Red = causal SNPs  |  Blue = non-causal SNPs  |  Wider = more uncertain"
      )),
      x = "CI Width (Upper − Lower)", y = "Density"
    ) +
    base_theme() +
    theme(legend.position = "bottom")

  save_plot(p, file.path(PLOT_DIR, paste0(prefix, "_ci_width.png")), w = 16, h = 8)
}

# =============================================================================
# PLOT 9: Summary heatmaps — Coverage and MSE in compact table-style figures.
#         Rows = estimators, Columns = SNP type (Causal / Non-causal).
#         Saved as two separate files: one for Coverage, one for MSE.
# =============================================================================
plot_summary_heatmap <- function(metrics, sim_label, prefix, n_label = "", scen_label = "") {
  has_noncausal <- any(metrics$IS_CAUSAL == FALSE, na.rm = TRUE)

  long <- metrics %>%
    mutate(
      Algorithm = factor(Algorithm, levels = ALG_ORDER),
      SNP_type  = ifelse(IS_CAUSAL, "Causal", "Non-causal")
    ) %>%
    select(Algorithm, SNP_type, Coverage, MSE, MeanBias) %>%
    pivot_longer(cols = c(Coverage, MSE, MeanBias),
                 names_to = "Metric", values_to = "Value") %>%
    mutate(
      Label = case_when(
        Metric == "Coverage" ~ sprintf("%.0f%%", Value * 100),
        Metric == "MSE"      ~ sprintf("%.3f",   Value),
        Metric == "MeanBias" ~ sprintf("%+.3f",  Value)
      )
    )

  cov_data <- long %>% filter(Metric == "Coverage")
  mse_data <- long %>% filter(Metric == "MSE")

  # Narrower plot when there is only one SNP type column (sim2)
  plot_w <- if (has_noncausal) 7 else 5

  p_cov <- ggplot(cov_data, aes(x = SNP_type, y = Algorithm, fill = Value)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = Label), size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                         midpoint = 0.7, limits = c(0, 1),
                         name = "Coverage") +
    scale_y_discrete(limits = rev(ALG_ORDER)) +
    labs(
      title    = paste("Coverage (nominal 95%) \u2014", sim_label),
      subtitle = wrap_subtitle(paste0(scenario_brief(scen_label), "  |  ", n_label)),
      caption  = "Variance estimated via EIF",
      x = NULL, y = NULL
    ) +
    base_theme() +
    theme(legend.position = "right",
          axis.text.x = element_text(size = 11))

  p_mse <- ggplot(mse_data, aes(x = SNP_type, y = Algorithm, fill = log1p(Value))) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = Label), size = 4, fontface = "bold") +
    scale_fill_gradient(low = "#1a9850", high = "#d73027",
                        name = "log(1+MSE)") +
    scale_y_discrete(limits = rev(ALG_ORDER)) +
    labs(
      title    = paste("MSE \u2014", sim_label),
      subtitle = wrap_subtitle(paste0(scenario_brief(scen_label), "  |  ", n_label)),
      caption  = "Variance estimated via EIF",
      x = NULL, y = NULL
    ) +
    base_theme() +
    theme(legend.position = "right",
          axis.text.x = element_text(size = 11))

  save_plot(p_cov, file.path(PLOT_DIR, paste0(prefix, "_heatmap_coverage.png")), w = plot_w, h = 6)
  save_plot(p_mse, file.path(PLOT_DIR, paste0(prefix, "_heatmap_mse.png")),      w = plot_w, h = 6)
}

# =============================================================================
# Main: run all plots for both simulation types
# =============================================================================
for (sim_type in c("linear", "interaction")) {

  all_csv <- file.path(OUT_DIR, sprintf("%s_all.csv", sim_type))
  if (!file.exists(all_csv)) {
    message(sprintf("Missing: %s  — run aggregate_results.R first", all_csv))
    next
  }

  message(sprintf("\n========== %s ==========", toupper(sim_type)))
  df <- read.csv(all_csv, stringsAsFactors = FALSE)

  if ("Algorithm" %in% colnames(df)) {
    n_before <- nrow(df)
    df <- df %>% filter(Algorithm %in% ALLOWED_ALGORITHMS)
    n_dropped <- n_before - nrow(df)
    if (n_dropped > 0) {
      message(sprintf("  Dropped %d rows from unsupported algorithms", n_dropped))
    }
  }

  label      <- if (sim_type == "linear") "Linear (order 1)" else "Interaction (order 1–4)"
  scen_label <- scenario_label(df)
  metrics    <- compute_metrics(df)

  message("  [1/6] Estimate density — causal SNPs")
  plot_estimates_causal(df, metrics, label, sim_type, scen_label)

  message("  [2/6] Bias distributions — causal SNPs")
  plot_bias(df, label, sim_type, scen_label)

  message("  [3/6] Ridge plot — causal SNPs")
  plot_ridge(df, label, sim_type, scen_label)

  message("  [4/6] Forest plot — mean bias")
  plot_forest_mean_bias(df, label, sim_type, scen_label)

  message("  [5/6] Forest plot — mean estimate")
  plot_forest_mean_estimate(df, label, sim_type, scen_label)

  message("  [6/6] Summary heatmaps")
  plot_summary_heatmap(metrics, label, sim_type,
                       n_label = count_label(df), scen_label = scen_label)
}

message("\n========== Done ==========")
message("Outputs in: ", PLOT_DIR)