# ============================================================
#  Grid Violin Plot  —  ATE Estimate
#  Sim2: Fixed treatment SNP, fixed beta = 1, fixed FVE = 0.4
#
#  Facet structure:
#    Rows : LD category  (Strong LD | Moderate LD | Minimal LD)
#    Cols : treatment SNP within each category (3 SNPs per row)
#
#  Each panel shows per-estimator violin + jitter of ATE estimates
#  across 200 replicates.  True ATE = 1 throughout.
#
#  Run from the AnchorGene_sim2 project root:
#    Rscript violin_grid_sim2.R
#
#  Requires: results/<snp>/linear_all.csv  (and interaction_all.csv)
#  for all 9 SNP scenarios — run aggregate_results.R for each first.
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)        

# ── 1. CONFIGURATION ─────────────────────────────────────────

BASE_DIR <- file.path(getwd(), "results")

# Nine treatment SNPs grouped by LD category
# Order within each group controls left-to-right column order
LD_GROUPS <- list(
  "Strong LD"   = c("rs0254", "rs0329", "rs0361"),
  "Moderate LD" = c("rs0293", "rs0314", "rs0387"),
  "Minimal LD"  = c("rs0201", "rs0211", "rs0225")
)

# Flat ordered vector of all 9 SNPs (used for loading and factor levels)
SCENARIOS <- unlist(LD_GROUPS, use.names = FALSE)

# Map each SNP to its LD category (for facet_nested row label)
SNP_TO_CATEGORY <- setNames(
  rep(names(LD_GROUPS), times = lengths(LD_GROUPS)),
  SCENARIOS
)

# Human-readable column labels for each SNP
# Position info added for context
SCENARIO_LABELS <- c(
  rs0254 = "rs0254\n(1,396 kb)",
  rs0329 = "rs0329\n(1,682 kb)",
  rs0361 = "rs0361\n(1,858 kb)",
  rs0293 = "rs0293\n(1,505 kb)",
  rs0314 = "rs0314\n(1,616 kb)",
  rs0387 = "rs0387\n(1,949 kb)",
  rs0201 = "rs0201\n(1,060 kb)",
  rs0211 = "rs0211\n(1,145 kb)",
  rs0225 = "rs0225\n(1,245 kb)"
)

# Estimators to include (C-TMLE0* excluded)
ESTIMATORS <- c("TMLE", "C-TMLE0","C-TMLE0*" ,"C-TMLE1", "Scalable-CTMLE", "Wang-TMLE", "Wang-CTMLE","C-TMLE1-OAL")

TRUE_ATE <- 1.0

PALETTE <- c(
  "TMLE"            = "#4E79A7",
  "C-TMLE0"         = "#59A14F",
  "C-TMLE0*"        = "#70373a",
  "C-TMLE1"         = "#F28E2B",
  "Scalable-CTMLE"  = "#B07AA1",
  "Wang-TMLE"       = "#E15759",
  "Wang-CTMLE"      = "#76B7B2",
  "C-TMLE1-OAL"     = "#FFC0CB"
)

message("Looking for results in: ", BASE_DIR)

# ── 2. LOAD FUNCTION ─────────────────────────────────────────
# Reads <scenario>/linear_all.csv or <scenario>/interaction_all.csv
# for all 9 LD scenarios, tags each with the scenario name and
# its LD category, then binds into one data frame.

load_all <- function(sim_type = c("interaction", "linear")) {
  sim_type <- match.arg(sim_type)
  filename <- paste0(sim_type, "_all.csv")
  rows     <- list()

  for (snp in SCENARIOS) {
    filepath <- file.path(BASE_DIR, snp, filename)

    if (!file.exists(filepath)) {
      message("Not found (skipping): ", filepath)
      next
    }

    df              <- read_csv(filepath, show_col_types = FALSE)
    df$scenario     <- snp
    df$ld_category  <- SNP_TO_CATEGORY[snp]
    rows[[length(rows) + 1]] <- df
  }

  if (length(rows) == 0)
    stop(
      "No files found for sim_type = '", sim_type, "'.\n",
      "Expected paths like: results/rs0254/", filename, "\n",
      "Run aggregate_results.R for each scenario first."
    )

  bind_rows(rows)
}

# ── 3. PLOT FUNCTION ─────────────────────────────────────────
# Called once for "linear" and once for "interaction".
# Returns a ggplot object; saving is done by the caller.

make_violin_plot <- function(sim_type) {

  dat <- load_all(sim_type)

  sim_label <- if (sim_type == "linear") "Linear (order 1)" else "Interaction (order 1\u20134)"

  # ── 3a. Prepare plot data ─────────────────────────────────
  plot_df <- dat |>
    filter(
      IS_CAUSAL == TRUE,
      Algorithm != "C-TMLE0*",
      Algorithm %in% ESTIMATORS
    ) |>
    mutate(
      Algorithm     = factor(Algorithm, levels = ESTIMATORS),
      ld_category   = factor(ld_category,   levels = names(LD_GROUPS)),
      scenario_label = factor(
        SCENARIO_LABELS[scenario],
        levels = SCENARIO_LABELS[SCENARIOS]
      )
    ) |>
    drop_na(Estimate)

  if (nrow(plot_df) == 0)
    stop("No rows remain after filtering for sim_type = '", sim_type, "'.")

  # ── 3b. Per-violin summary stats ──────────────────────────
  stats_df <- plot_df |>
    mutate(covered = (CI_Lower <= TRUE_ATE) & (TRUE_ATE <= CI_Upper)) |>
    drop_na(Estimate, CI_Lower, CI_Upper) |>
    group_by(Algorithm, scenario_label, ld_category) |>
    summarise(
      mean_bias = mean(Bias,    na.rm = TRUE),
      sd_est    = sd(Estimate,  na.rm = TRUE),
      coverage  = mean(covered, na.rm = TRUE),
      .groups   = "drop"
    ) |>
    mutate(
      # Values only — "Bias / SD / Cov" labels appear once as a header
      label = sprintf(
        "%+.2f\n%.2f\n%.0f%%",
        mean_bias, sd_est, coverage * 100
      )
    )

  # One header per panel: "Bias / SD / Cov" sitting left of the first violin
  header_df <- plot_df |>
    distinct(scenario_label, ld_category) |>
    mutate(
      Algorithm = factor(ESTIMATORS[1], levels = ESTIMATORS),
      label     = "Bias:\nSD:   \nCov: "
    )

  # ── 3c. Build plot ────────────────────────────────────────
  ggplot(plot_df, aes(x = Algorithm, y = Estimate, fill = Algorithm)) +

    # Violin
    geom_violin(
      trim      = TRUE,
      scale     = "width",
      alpha     = 0.55,
      linewidth = 0.35,
      colour    = "grey25"
    ) +

    # Jittered individual replicates
    geom_jitter(
      width  = 0.15,
      size   = 0.5,
      alpha  = 0.20,
      colour = "grey20"
    ) +

    # Per-violin annotation: values only (Bias, SD, Coverage)
    geom_text(
      data        = stats_df,
      aes(x       = Algorithm,
          y       = Inf,
          label   = label),
      inherit.aes = FALSE,
      vjust       = 1.2,
      fontface    = "plain",
      size        = 4.2,
      lineheight  = 0.85,
      family      = "Verdana",
      colour      = "grey20"
    ) +

    # Header label "Bias / SD / Cov" — once per panel, left of first violin
    geom_text(
      data        = header_df,
      aes(x       = Algorithm,
          y       = Inf,
          label   = label),
      inherit.aes = FALSE,
      hjust       = 1.5,        # sit just to the left of the first violin
      vjust       = 1.2,
      fontface    = "bold",
      size        = 4.2,
      lineheight  = 0.85,
      family      = "Verdana",
      colour      = "grey30"
    ) +

    # True ATE = 1 (dashed horizontal line)
    geom_hline(
      yintercept = TRUE_ATE,
      colour     = "black",
      linetype   = "dashed",
      linewidth  = 0.6
    ) +

    # 3x3 grid: rows = LD category, cols = SNP within group
    facet_nested_wrap(
      ~ ld_category + scenario_label,
      ncol        = 3,
      scales      = "free_y",
      strip.position = "top",
      nest_line   = element_line(linewidth = 0.5)
    ) +

    # Clip y to data range + small padding; no forced zero
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +

    scale_fill_manual(values = PALETTE, name = "Estimator") +

    labs(
      title    = paste0(sim_label, " - ATE Estimate by Estimator and LD Scenario"),
      subtitle = paste0(
        "Sim2  |  Fixed treatment SNP  |  True ATE = 1  |  FVE = 0.4  |  ",
        "200 replicates per scenario  |  Dashed = true ATE  |  ",
        "Rows: LD category  |  Cols: treatment SNP"
      ),
      x = NULL,
      y = "ATE Estimate"
    ) +

    theme_bw(base_size = 16) +
    theme(
      text               = element_text(family = "Verdana"),
      # Nested strip: outer (LD category) and inner (SNP) rows
      strip.background   = element_rect(fill = "#F0EDE8", colour = "grey75"),
      strip.text         = element_text(face = "bold", size = 16, family = "Verdana"),
      strip.text.y.left  = element_text(angle = 90),
      axis.text.x        = element_text(angle = 35, hjust = 1, size = 14, family = "Verdana"),
      axis.text.y        = element_text(size = 14, family = "Verdana"),
      axis.title.y       = element_text(size = 16, family = "Verdana"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
      legend.position    = "bottom",
      legend.title       = element_text(face = "bold", size = 18, family = "Verdana"),
      legend.text        = element_text(size = 17, family = "Verdana"),
      legend.key.size    = unit(1.3, "lines"),
      plot.title         = element_text(face = "bold", size = 24, hjust = 0.5, family = "Verdana"),
      plot.subtitle      = element_text(size = 18, colour = "grey45", hjust = 0.5, family = "Verdana"),
      plot.margin        = margin(12, 10, 6, 10)
    ) +

    coord_cartesian(clip = "off") +

    guides(fill = guide_legend(nrow = 1))
}

# ── 4. LOOP OVER BOTH SIMULATION TYPES ───────────────────────

for (sim_type in c("linear", "interaction")) {

  message("\n\u2500\u2500 Processing: ", sim_type, " \u2500\u2500")

  p <- make_violin_plot(sim_type)

  # Width: 3 SNP columns × 5 inches each = 15 inches
  # Height: 3 LD category rows × 5 inches each = 15 inches
  outfile <- paste0("violin_grid_sim2_", sim_type, ".png")
  ggsave(
    filename = outfile,
    plot     = p,
    width    = 3 * 7,
    height   = 3 * 6,
    dpi      = 300,
    bg       = "white"
  )
  message("Saved \u2192 ", outfile)

}