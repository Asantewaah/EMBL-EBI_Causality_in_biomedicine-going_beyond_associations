# =============================================================================
# ctmle_mc_analysis.R — SIM2: Fixed Treatment SNP Simulation
# =============================================================================
# Monte Carlo analysis for SIM2: fixed treatment SNP with general causal model.
#
# Key design:
#   - FIXED TREATMENT SNP (one of: rs0020, rs0044, rs0024)
#   - Binary treatment: A = 1{X[treatment_SNP] > 0}
#   - DGM: Y = beta*A + mu_r(X) where beta = 1 (fixed), FVU = 0.4 (fixed)
#   - mu_r(X) = 5 causal terms drawn from the 49 OTHER SNPs (excluding treatment)
#   - True ATE = beta = 1 by construction for every replicate
#   - 200 replicates per treatment SNP (20 batches x 10 replicates)
#   - Two simulation types per batch:
#       linear      : min_order = 1, max_order = 1
#       interaction : min_order = 1, max_order = 4 (nests the linear case)
#
# USAGE (local test):
#   TREATMENT_SNP=rs0020 Rscript ctmle_mc_analysis.R --batch_id=1 --seed=42 --test
#
# USAGE (production batch):
#   TREATMENT_SNP=rs0020 BATCH_ID=1 SEED=42 NSLOTS=10 Rscript ctmle_mc_analysis.R
#
# USAGE (EDDIE array job):
#   TREATMENT_SNP=rs0020 qsub scripts/eddie_mc_array.qsub.sh
#
# OUTPUT:
#   results/<treatment_snp>/linear/linear_batch<ID>_seed<SEED>.csv
#   results/<treatment_snp>/interaction/interaction_batch<ID>_seed<SEED>.csv

suppressPackageStartupMessages({
  library(future.apply)
  library(parallelly)
  library(dplyr)
})

source("src/estimators.R")   # ctmle1, ctmle0, ctmle0star, ctmle1_oal_pkg, scalable_ctmle, tmle, wang_tmle, wang_ctmle, ctmle1_oal

# =============================================================================
# Parse arguments
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(args, flag, default) {
  m <- regmatches(args, regexpr(paste0("--", flag, "=([^ ]+)"), args, perl = TRUE))
  if (length(m) == 0) return(default)
  as.character(sub(paste0("--", flag, "="), "", m))
}

TEST_MODE     <- "--test" %in% args
BATCH_ID      <- as.integer(Sys.getenv("BATCH_ID",       parse_arg(args, "batch_id",      "1")))
SEED          <- as.integer(Sys.getenv("SEED",           parse_arg(args, "seed",          "42")))
N_WORKERS     <- as.integer(Sys.getenv("NSLOTS",         parse_arg(args, "workers",       "1")))
TREATMENT_SNP <- Sys.getenv("TREATMENT_SNP",             parse_arg(args, "treatment_snp", "rs0254"))
OUT_ROOT_BASE <- Sys.getenv("OUT_ROOT",                  parse_arg(args, "out_root",      "results"))

# =============================================================================
# Fixed settings — not varied in sim2
# =============================================================================
GENO_PATH   <- "data/geno.csv"
N_CAUSAL    <- 5L
ALPHA       <- 0.05
BETA        <- 1.0    # fixed treatment effect; true ATE = 1 for every replicate
FVU         <- 0.4    # fixed fraction of variance unexplained

# Batch size: 3 replicates in test mode, 10 in production
B_PER_BATCH <- if (TEST_MODE) 3L else 10L

# Validate treatment SNP
valid_snps <- c("rs0254", "rs0329", "rs0361",   # strong LD
                "rs0293", "rs0314", "rs0387",   # moderate LD
                "rs0201", "rs0211", "rs0225")   # minimal LD
if (!(TREATMENT_SNP %in% valid_snps)) {
  stop("TREATMENT_SNP must be one of: ", paste(valid_snps, collapse = ", "),
       ". Got: '", TREATMENT_SNP, "'")
}

# Output directories — one per sim type, nested under the treatment SNP
OUT_DIR_LIN <- file.path(OUT_ROOT_BASE, TREATMENT_SNP, "linear")
OUT_DIR_INT <- file.path(OUT_ROOT_BASE, TREATMENT_SNP, "interaction")
dir.create(OUT_DIR_LIN, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_INT, recursive = TRUE, showWarnings = FALSE)

message(sprintf(
  "[treatment_snp=%s | batch %d | seed %d | B=%d | workers=%d | test=%s]",
  TREATMENT_SNP, BATCH_ID, SEED, B_PER_BATCH, N_WORKERS, TEST_MODE
))

# =============================================================================
# Parallelisation
# =============================================================================
if (N_WORKERS > 1 && parallelly::supportsMulticore()) {
  plan(multicore, workers = N_WORKERS)
  message(sprintf("Parallel: %d multicore workers (fork-based, shared memory)", N_WORKERS))
} else if (N_WORKERS > 1) {
  plan(multisession, workers = N_WORKERS)
  message(sprintf("Parallel: %d multisession workers", N_WORKERS))
} else {
  plan(sequential)
  message("Sequential execution (set NSLOTS or --workers=N for parallel)")
}

# =============================================================================
# Load genotype data ONCE — shared across all workers via closure
# =============================================================================
message("Loading genotype data...")
X_full  <- as.matrix(read.csv(GENO_PATH, row.names = 1))
snp_ids <- colnames(X_full)

# Restrict to SNPs 201-400 — a 200-SNP window (1,060-2,030 kb) with
# well-characterised LD block structure spanning all three LD categories.
# This reduces the confounder dimension from 999 to ~199 and ensures
# all treatment SNPs and confounders share the same genomic region.
SNP_WINDOW <- 201:400
X_full  <- X_full[, SNP_WINDOW, drop = FALSE]
snp_ids <- colnames(X_full)
message(sprintf("Restricted to SNPs 201-400: %d SNPs spanning %.1f-%.1f kb",
                ncol(X_full),
                as.numeric(sub("kb", "", colnames(X_full)[1])),
                as.numeric(sub("kb", "", colnames(X_full)[ncol(X_full)]))))

if (!(TREATMENT_SNP %in% snp_ids)) {
  stop("TREATMENT_SNP '", TREATMENT_SNP, "' not found in genotype matrix. ",
       "First 10 SNPs: ", paste(head(snp_ids, 10), collapse = ", "))
}

# Mean-impute any missing values (should be none in geno.csv but safe to check)
for (j in seq_len(ncol(X_full))) {
  nas <- is.na(X_full[, j])
  if (any(nas)) X_full[nas, j] <- mean(X_full[, j], na.rm = TRUE)
}

message(sprintf("Loaded: %d samples x %d SNPs", nrow(X_full), length(snp_ids)))
message(sprintf("Treatment SNP: %s (column %d of %d)",
                TREATMENT_SNP, match(TREATMENT_SNP, snp_ids), length(snp_ids)))

# LD block structure is no longer used (cluster-robust inference removed)

# =============================================================================
# DGM helpers
# =============================================================================

# Sample n_causal interaction terms from snp_pool with N(0,1) coefficients.
# snp_pool must already exclude the treatment SNP.
sample_causal_model <- function(snp_pool, n_causal, min_order, max_order) {
  remaining <- snp_pool
  model     <- vector("list", n_causal)
  for (k in seq_len(n_causal)) {
    order    <- sample(min_order:max_order, 1)
    order    <- min(order, length(remaining))
    variants <- sample(remaining, order, replace = FALSE)
    remaining <- setdiff(remaining, variants)
    model[[k]] <- list(variants = variants, coef = rnorm(1))
  }
  model
}

# Evaluate the background signal mu_r(X) — does NOT include the treatment term
compute_mu_y <- function(model, X) {
  mu <- numeric(nrow(X))
  for (term in model) {
    vals <- X[, term$variants, drop = FALSE]
    mu   <- mu + term$coef * apply(vals, 1, prod)
  }
  mu
}

# =============================================================================
# run_one_replicate()
#
# For one replicate:
#   1. Binarise the fixed treatment SNP -> A
#   2. Sample mu_r(X) from the remaining 49 SNPs
#   3. Build full signal: beta*A + mu_r(X)
#   4. Add noise scaled so Var(noise)/Var(Y_total) = fvu
#   5. Run all estimators; record true_ATE = beta = 1
# =============================================================================
run_one_replicate <- function(b_local, global_seed, batch_id, treatment_snp,
                              sim_type, n_causal, min_order, max_order,
                              X, snp_ids, fvu, alpha, beta) {

  # Deterministic per-replicate seed — no collisions across batches or seeds
  rep_seed <- global_seed * 10000L + (batch_id - 1L) * 1000L + b_local
  set.seed(rep_seed)

  # Fixed binary treatment
  A <- ifelse(X[, treatment_snp] > 0, 1.0, 0.0)

  # Background causal model drawn from the 49 SNPs excluding the treatment SNP
  snp_pool <- setdiff(snp_ids, treatment_snp)
  model    <- sample_causal_model(snp_pool, n_causal, min_order, max_order)
  mu_r     <- compute_mu_y(model, X)

  # Full signal: treatment effect + background
  signal     <- beta * A + mu_r

  # Noise scaled so that Var(noise) / Var(Y) = fvu
  var_signal <- var(signal)
  sigma_y    <- if (var_signal < 1e-10) 1 else sqrt(var_signal * fvu / (1 - fvu))
  Y          <- signal + rnorm(nrow(X), 0, sigma_y)

  # Identifiers
  run_id       <- sprintf("%s_%s_b%02d_r%03d", treatment_snp, sim_type, batch_id, b_local)
  bootstrap_id <- (batch_id - 1L) * 10L + b_local

  # Estimator input: W = all SNPs except the treatment SNP
  W      <- X[, snp_pool, drop = FALSE]
  colnames(W) <- paste0("W", seq_len(ncol(W)))
  data_i <- data.frame(A = A, W, Y = Y, stringsAsFactors = FALSE)

  # True ATE = beta = 1 by construction — no analytical computation needed
  res <- tryCatch(
    collect_results_snp(data_i, true_ate = beta, alpha = alpha),
    error = function(e) {
      message("  Error rep=", b_local, " [", treatment_snp, " | ", sim_type, "]: ",
              conditionMessage(e))
      NULL
    }
  )

  if (!is.null(res)) {
    res$run_id        <- run_id
    res$seed          <- global_seed
    res$batch_id      <- batch_id
    res$bootstrap_id  <- bootstrap_id
    res$sim_type      <- sim_type
    res$treatment_snp <- treatment_snp
    res$SNP           <- treatment_snp
    res$IS_CAUSAL     <- TRUE
    res$beta          <- beta
    res$fvu           <- fvu
  }

  res
}

# =============================================================================
# run_batch()
#
# Runs B_PER_BATCH replicates for one (treatment_snp, sim_type) combination.
# Skips safely if the output file already exists (SGE re-submission safe).
# =============================================================================
run_batch <- function(treatment_snp, sim_type, out_dir,
                      n_causal, min_order, max_order) {

  out_file <- file.path(out_dir,
    sprintf("%s_batch%02d_seed%d.csv", sim_type, BATCH_ID, SEED))

  if (file.exists(out_file)) {
    message(sprintf("[%s | %s] Batch %d already done — skipping (%s)",
                    treatment_snp, sim_type, BATCH_ID, out_file))
    return(invisible(out_file))
  }

  message(sprintf("[%s | %s] Starting batch %d (seed=%d, B=%d, workers=%d)",
                  treatment_snp, sim_type, BATCH_ID, SEED, B_PER_BATCH, N_WORKERS))
  t0 <- proc.time()

  results_list <- future_lapply(
    seq_len(B_PER_BATCH),
    function(b_local) {
      run_one_replicate(
        b_local       = b_local,
        global_seed   = SEED,
        batch_id      = BATCH_ID,
        treatment_snp = treatment_snp,
        sim_type      = sim_type,
        n_causal      = n_causal,
        min_order     = min_order,
        max_order     = max_order,
        X             = X_full,
        snp_ids       = colnames(X_full),
        fvu           = FVU,
        alpha         = ALPHA,
        beta          = BETA
      )
    },
    future.seed = TRUE
  )

  combined <- do.call(rbind, Filter(Negate(is.null), results_list))

  if (is.null(combined) || nrow(combined) == 0) {
    message(sprintf("[%s | %s] Batch %d produced no results — check errors above.",
                    treatment_snp, sim_type, BATCH_ID))
    return(invisible(NULL))
  }

  # Canonical column order
  col_order <- c(
    "run_id", "seed", "batch_id", "bootstrap_id", "sim_type",
    "treatment_snp", "SNP", "IS_CAUSAL",
    "beta", "fvu",
    "Algorithm", "Estimate", "Standard.Error",
    "CI_Lower", "CI_Upper", "true_ATE",
    "Bias", "Bias.to.se", "pvalue"
  )
  present  <- intersect(col_order, colnames(combined))
  combined <- combined[, present, drop = FALSE]

  write.csv(combined, out_file, row.names = FALSE)
  elapsed <- (proc.time() - t0)[["elapsed"]]
  message(sprintf("[%s | %s] Batch %d done: %d rows -> %s  (%.1f min)",
                  treatment_snp, sim_type, BATCH_ID, nrow(combined),
                  out_file, elapsed / 60))

  invisible(out_file)
}

# =============================================================================
# Main — run linear and interaction for this treatment SNP
# =============================================================================
message(sprintf("\n========== LINEAR      [treatment_snp=%s] ==========", TREATMENT_SNP))
run_batch(TREATMENT_SNP, "linear",      OUT_DIR_LIN, N_CAUSAL, min_order = 1, max_order = 1)

message(sprintf("\n========== INTERACTION [treatment_snp=%s] ==========", TREATMENT_SNP))
run_batch(TREATMENT_SNP, "interaction", OUT_DIR_INT, N_CAUSAL, min_order = 1, max_order = 4)

message("\n========== Batch complete ==========")
message(sprintf("Outputs in: %s", file.path(OUT_ROOT_BASE, TREATMENT_SNP)))