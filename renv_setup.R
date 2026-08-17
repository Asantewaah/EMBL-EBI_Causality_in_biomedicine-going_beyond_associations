# ============================================================
# renv_setup.R
# ============================================================
# Run ONCE.
# It sets up renv, installs every package the simulation pipeline needs, and
# writes a real renv.lock file (with correct, verified versions/hashes)
# by actually querying CRAN -- not a hand-written approximation.
#
# USAGE:
#   Rscript renv_setup.R
# ============================================================

# ── 1. Install renv itself if not already available ──────────
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# ── 2. Initialize renv for this project ───────────────────────
# bare = TRUE: don't let renv auto-scan/guess dependencies from source
# files (it can miss packages only used inside string-built formulas,
# etc.) -- we supply the full list explicitly below instead.
if (!file.exists("renv.lock")) {
  renv::init(bare = TRUE)
}

# ── 3. Full package list used across the simulation pipeline + course session:
#
# CRAN packages (installed via renv::install(name)):
cran_packages <- c(
  "tidyverse",     # dplyr, tidyr, tibble, ggplot2, purrr, stringr, readr, forcats
  "dplyr",         # explicit, though covered by tidyverse
  "ggplot2",       # explicit, though covered by tidyverse
  "ggrepel",       # non-overlapping plot label placement
  "gridExtra",     # arranging multiple grid-based plots
  "MASS",          # mvrnorm() in simulate_genotypes_CEU.R
  "scales",        # scales::comma / scales::percent in plotting code
  "boot",          # bootstrap resampling utilities
  "viridis",       # colorblind-friendly palettes
  # -- Machine learning / ensemble methods --
  "SuperLearner",  # Super Learner ensemble framework
  "glmnet",        # LASSO / cv.glmnet, used throughout the estimators
  "xgboost",       # extreme gradient boosting
  "ranger",        # fast random forests
  "earth",         # multivariate adaptive regression splines (MARS)
  # -- Causal inference --
  "tmle",          # targeted maximum likelihood estimation
  # -- Statistical modeling --
  "glm2",          # robust glm fitting, used in collaborative CV loss
  "gam",           # generalized additive models
  "caret",         # createFolds() for cross-fitting
  "ctmle",         # C-TMLE package
  "future.apply",  # future_lapply() for parallel replicate execution
  "parallelly",    # supportsMulticore(), worker detection
  # -- I/O --
  "writexl"        # write results to Excel
)

required_packages <- cran_packages

# NOTE: `grid` (low-level grid graphics) is a base-R package that ships
# with every R installation -- it cannot be installed via install.packages()
# and isn't included in the install call above. renv still records its
# version correctly during snapshot() automatically.

# ── 4. Install everything ──────────────────────────────────────
message(sprintf("Installing %d packages into the renv library...",
                length(required_packages)))
renv::install(required_packages)

# ── 5. Snapshot -- writes renv.lock with real, verified versions ──
renv::snapshot(prompt = FALSE)

message("\nDone. renv.lock has been written.")
message("Next steps:")
message("  1. Check that renv.lock now exists and lists all packages above.")
message("  2. Commit renv.lock, .Rprofile, renv/activate.R, and renv/settings.json to git.")
message("  3. Do NOT commit renv/library/ (the actual installed packages) --")
message("     add it to .gitignore if renv::init() didn't already.")
