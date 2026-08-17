# ============================================================
# restore_renv.R
# ============================================================
# For EMBL-EBI training participants.
#
# This installs the exact package versions used to build/test this
# repository, into a project-local library (won't touch or conflict
# with anything else on your machine).
#
# USAGE:
#   Rscript restore_renv.R
#
# (Or, if you're working in RStudio with this repo open as a Project,
# renv activates automatically when you open it -- just run
# renv::restore() at the console instead.)
# ============================================================

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

renv::restore(prompt = FALSE)

message("\nDone. All packages restored to the versions in renv.lock.")
