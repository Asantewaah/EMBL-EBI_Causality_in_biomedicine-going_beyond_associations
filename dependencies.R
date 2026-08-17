# ============================================================
# dependencies.R
# ============================================================
# This file is NOT meant to be run as a script.
#
# renv's default ("implicit") snapshot mode only records packages it can
# detect being used, by scanning the repo's R files for literal
# library()/require()/pkg::fun() patterns. Since this repo's actual
# analysis scripts (which use glmnet, glm2, etc.) live in a separate
# project, renv had no evidence these packages were needed and silently
# dropped some of them from renv.lock during the last snapshot.
#
# This file exists purely to give renv's scanner that evidence, so every
# package in the required list gets captured correctly on the next
# renv::snapshot() -- regardless of what analysis code is or isn't
# checked into this repo yet.
# ============================================================

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(MASS)
library(scales)
library(boot)
library(viridis)

# -- Machine learning / ensemble methods --
library(SuperLearner)
library(glmnet)
library(xgboost)
library(ranger)
library(earth)

# -- Causal inference --
library(tmle)

# -- Statistical modeling --
library(glm2)
library(gam)
library(caret)
library(ctmle)
library(future.apply)
library(parallelly)

# -- I/O --
library(writexl)
