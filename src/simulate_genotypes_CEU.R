# ============================================================
# Simulate Realistic European (CEU) Genotypes
# ============================================================
# LD structure is modelled using empirically calibrated parameters
# derived from HapMap CEU / 1000 Genomes Phase 3 Europeans:
#
#   - LD half-decay distance : ~30–50 kb (we use 40 kb default)
#   - Background LD (>1 Mb)  : r² ≈ 0.01–0.02
#   - Recombination hotspots : ~1 in every 50–100 kb on average
#   - MAF distribution       : Beta(0.5, 0.5) L-shaped (allele-freq spectrum)
#   - FST (CEU vs. other EUR): ~0.004–0.010
#
# Decay model: Hill-Weir (1988)
#   E[r²] ≈ (1 - c)^d  where c = recombination fraction and d = distance
# We approximate this with a two-component model:
#   r(d) = (1 - bg) * exp(-d / lambda) + bg
# where lambda = half_decay_kb / log(2), bg = background LD floor.

library(tidyverse)
library(MASS)   # mvrnorm

set.seed(41)

# ── 0. Parameters ────────────────────────────────────────────

n_samples       <- 10000          # number of individuals to simulate
n_snps          <- 50          # number of SNPs
snp_spacing_kb  <- 10           # physical distance between adjacent SNPs (kb)
                                # 10 kb spacing is typical for dense GWAS arrays

# CEU-calibrated LD decay parameters (empirical, 1000 Genomes / HapMap)
half_decay_kb   <- 80 #80           # distance at which r drops to 0.5 (CEU: 30–50 kb)
background_ld   <- 0.04 #0.04         # residual r at large distance (CEU: ~0.01–0.02)

# Recombination hotspots: randomly placed, each locally boosts decay rate
add_hotspots    <- TRUE
hotspot_prob    <- 0.02 #0.06         # probability any SNP interval is a hotspot
hotspot_factor  <- 10           # hotspots increase local decay rate 10×

# MAF distribution: Beta(0.5, 0.5) gives a realistic U-shaped SFS for Europeans
maf_shape1      <- 0.5          # Beta parameter α
maf_shape2      <- 0.5          # Beta parameter β
maf_min         <- 0.05         # exclude very rare variants
maf_max         <- 0.50

# ── 1. SNP positions ─────────────────────────────────────────
# Small jitter around regular spacing to simulate non-uniform array coverage
base_pos_kb <- cumsum(c(0, rep(snp_spacing_kb, n_snps - 1)))
jitter_kb   <- runif(n_snps, -snp_spacing_kb * 0.3, snp_spacing_kb * 0.3)
pos_kb      <- sort(base_pos_kb + jitter_kb)
pos_kb      <- pos_kb - min(pos_kb)   # start at 0

snp_names <- paste0("rs", formatC(seq_len(n_snps), width = 4, flag = "0"))

# ── 2. CEU allele frequencies ────────────────────────────────
# Beta(0.5, 0.5) approximates the U-shaped site frequency spectrum of Europeans
raw_maf <- rbeta(n_snps, maf_shape1, maf_shape2)
# Rescale to [maf_min, maf_max]
p_anc <- maf_min + raw_maf * (maf_max - maf_min)

message(sprintf("MAF distribution — mean: %.3f  median: %.3f",
                mean(p_anc), median(p_anc)))

# ── 3. Recombination hotspot map ─────────────────────────────
# Each inter-SNP interval independently has a hotspot_prob chance of
# being a hotspot, which compresses the effective genetic distance.
is_hotspot <- c(FALSE, runif(n_snps - 1) < hotspot_prob)
n_hot      <- sum(is_hotspot)
if (add_hotspots) message(sprintf("Recombination hotspots: %d placed", n_hot))

# ── 4. Build CEU LD correlation matrix ───────────────────────
# Two-component exponential decay model calibrated to CEU:
#   r(d) = (1 - bg) * exp(-d / lambda) + bg
# where lambda = half_decay_kb / log(2)
lambda <- half_decay_kb / log(2)

build_ld_matrix <- function(pos_kb, lambda, background_ld,
                            is_hotspot, hotspot_factor, add_hotspots) {
  n <- length(pos_kb)
  R <- matrix(1, n, n)

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Physical distance
      d <- pos_kb[j] - pos_kb[i]

      # Effective decay rate: any hotspot between i and j speeds up decay
      if (add_hotspots && any(is_hotspot[(i + 1):j])) {
        eff_lambda <- lambda / hotspot_factor
      } else {
        eff_lambda <- lambda
      }

      # Hill-Weir-inspired two-component decay
      r_ij <- (1 - background_ld) * exp(-d / eff_lambda) + background_ld
      R[i, j] <- R[j, i] <- r_ij
    }
  }

  # Ensure positive definiteness (small numerical corrections)
  eig    <- eigen(R, symmetric = TRUE)
  eig$values <- pmax(eig$values, 1e-6)
  R_pd   <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  # Rescale diagonal back to 1
  D      <- diag(1 / sqrt(diag(R_pd)))
  R_pd   <- D %*% R_pd %*% D
  R_pd
}

message("Building CEU LD correlation matrix...")
R <- build_ld_matrix(pos_kb, lambda, background_ld,
                     is_hotspot, hotspot_factor, add_hotspots)

# ── 5. Cholesky factorisation for correlated sampling ────────
L <- chol(R)

# ── 6. Generate genotypes under HWE with CEU LD structure ───
message("Simulating genotypes...")

# Draw correlated standard normals via Cholesky
Z <- matrix(rnorm(n_samples * n_snps), nrow = n_samples) %*% L

# Map to correlated uniform via probit (preserves rank correlations)
U <- pnorm(Z)

# Convert to 0/1/2 dosage under Hardy-Weinberg equilibrium per SNP
geno <- matrix(NA_integer_, nrow = n_samples, ncol = n_snps)
for (s in seq_len(n_snps)) {
  p  <- p_anc[s]
  q  <- 1 - p
  p0 <- q^2          # P(AA)
  p1 <- p0 + 2*p*q   # P(Aa)
  geno[, s] <- as.integer(
    cut(U[, s], breaks = c(0, p0, p1, 1),
        include.lowest = TRUE)
  ) - 1L
}

sample_names <- paste0("EUR_", formatC(seq_len(n_samples), width = 4, flag = "0"))
dimnames(geno) <- list(sample_names, snp_names)

# ── 7. Metadata ───────────────────────────────────────────────
metadata <- tibble(
  sample_id  = sample_names,
  population = "CEU",
  ancestry   = "European"
)

# ── 8. Quality checks ────────────────────────────────────────
message("\n── Quality checks ─────────────────────────────────────")

# Observed MAF
obs_p   <- colMeans(geno, na.rm = TRUE) / 2
obs_maf <- pmin(obs_p, 1 - obs_p)
message(sprintf("Observed MAF — mean: %.3f  median: %.3f  range: %.3f–%.3f",
                mean(obs_maf), median(obs_maf), min(obs_maf), max(obs_maf)))

# HWE test per SNP
hwe_p <- function(x) {
  x   <- x[!is.na(x)]; n <- length(x)
  p   <- sum(x) / (2 * n); q <- 1 - p
  obs <- c(sum(x == 0), sum(x == 1), sum(x == 2))
  exp <- c(q^2, 2*p*q, p^2) * n
  if (any(exp < 5)) return(NA_real_)
  pchisq(sum((obs - exp)^2 / exp), df = 1, lower.tail = FALSE)
}
hwe_pvals  <- apply(geno, 2, hwe_p)
n_hwe_fail <- sum(hwe_pvals < 0.05, na.rm = TRUE)
message(sprintf("HWE failures (p<0.05): %d / %d  (expected ~%d by chance)",
                n_hwe_fail, n_snps, round(n_snps * 0.05)))

# Observed LD decay: mean r² in distance bins
r2_mat   <- cor(geno, use = "pairwise.complete.obs")^2
dist_mat <- outer(pos_kb, pos_kb, `-`) %>% abs()

ld_df <- tibble(
  dist_kb = as.vector(dist_mat[upper.tri(dist_mat)]),
  r2      = as.vector(r2_mat[upper.tri(r2_mat)])
) %>%
  mutate(bin = cut(dist_kb,
                   breaks = c(0, 10, 25, 50, 100, 200, 500, Inf),
                   labels = c("<10", "10-25", "25-50", "50-100",
                              "100-200", "200-500", ">500"))) %>%
  group_by(bin) %>%
  summarise(mean_r2 = mean(r2), n_pairs = n(), .groups = "drop")

message("\nLD decay by distance (should match CEU empirical values):")
message(sprintf("  %-10s  %s  %s", "Dist (kb)", "Mean r²", "N pairs"))
walk2(ld_df$bin, seq_len(nrow(ld_df)), function(b, i) {
  message(sprintf("  %-10s  %.4f   %d", b, ld_df$mean_r2[i], ld_df$n_pairs[i]))
})

# Reference CEU empirical r² values (from HapMap/1000G literature):
#   <10 kb  : r² ≈ 0.35–0.50
#   25-50 kb: r² ≈ 0.20–0.30
#   >100 kb : r² ≈ 0.05–0.10
#   >1 Mb   : r² ≈ 0.01–0.02

# ── 9. SNP position annotation ───────────────────────────────
snp_info <- tibble(
  snp_id    = snp_names,
  pos_kb    = round(pos_kb, 2),
  maf       = round(obs_maf, 4),
  is_hotspot_boundary = is_hotspot
)

# ── 10. LD matrix & plots ────────────────────────────────────
# r2_mat was computed above (section 8). Label rows/cols with position in kb
# so axes are interpretable as genomic coordinates.
pos_labels <- paste0(round(pos_kb), "kb")
dimnames(r2_mat) <- list(pos_labels, pos_labels)

# ── 10a. Triangular LD heatmap ───────────────────────────────
r2_upper <- r2_mat
r2_upper[lower.tri(r2_upper, diag = TRUE)] <- NA

ld_long <- r2_upper %>%
  as.data.frame() %>%
  rownames_to_column("SNP1") %>%
  pivot_longer(-SNP1, names_to = "SNP2", values_to = "r2") %>%
  filter(!is.na(r2)) %>%
  mutate(
    SNP1 = factor(SNP1, levels = pos_labels),
    SNP2 = factor(SNP2, levels = pos_labels),
    i    = as.integer(SNP1),
    j    = as.integer(SNP2),
    x    = (i + j) / 2,      # rotated x: midpoint
    y    = (j - i)            # rotated y: SNP separation
  )

# Annotate hotspot boundaries as vertical lines on x-axis
hotspot_kb <- pos_kb[is_hotspot]

p_triangle <- ggplot(ld_long, aes(x = x, y = y, fill = r2)) +
  geom_tile(width = 1, height = 1, colour = NA) +
  # Hotspot positions shown as tick marks at the base
  { if (length(hotspot_kb) > 0)
      geom_vline(
        xintercept = which(is_hotspot),
        colour = "steelblue", linewidth = 0.4, linetype = "dashed", alpha = 0.6
      )
  } +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#FFF3CD", "#FFAB40", "#E64A19", "#B71C1C"),
    values  = c(0, 0.25, 0.5, 0.75, 1),
    limits  = c(0, 1),
    name    = expression(r^2)
  ) +
  scale_x_continuous(
    breaks = seq_along(pos_labels),
    labels = pos_labels,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  coord_fixed() +
  labs(
    title    = "CEU LD matrix (r²) — triangular view",
    subtitle = sprintf("%d SNPs  |  dashed lines = recombination hotspots", n_snps),
    x        = "SNP position",
    y        = "SNP separation"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(colour = "grey50", size = 9),
    axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size  = ifelse(n_snps > 60, 5, 7)),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    panel.grid      = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin     = margin(10, 20, 10, 10)
  )

# ── 10b. LD decay curve (observed vs. theoretical CEU) ───────
# Fine-grained binning for a smooth curve
decay_raw <- tibble(
  dist_kb = as.vector(dist_mat[upper.tri(dist_mat)]),
  r2      = as.vector(r2_mat[upper.tri(r2_mat)])
) %>%
  mutate(bin_kb = cut(dist_kb,
                      breaks = c(seq(0, 200, by = 5), 500, Inf),
                      include.lowest = TRUE)) %>%
  group_by(bin_kb) %>%
  summarise(mid_kb  = mean(dist_kb),
            mean_r2 = mean(r2),
            .groups = "drop")

# Theoretical CEU decay curve
d_seq  <- seq(0, max(pos_kb), length.out = 500)
theory <- tibble(
  dist_kb = d_seq,
  r2      = (1 - background_ld) * exp(-d_seq / lambda) + background_ld
)

# CEU empirical reference points (HapMap / 1000G literature)
ceu_ref <- tibble(
  dist_kb  = c(5,  17,  40,  75,  150, 350),
  r2       = c(0.43, 0.30, 0.20, 0.12, 0.07, 0.04),
  label    = c("<10kb", "10-25kb", "25-50kb", "50-100kb", "100-200kb", "200-500kb")
)

p_decay <- ggplot() +
  geom_line(data = theory, aes(x = dist_kb, y = r2, colour = "Theoretical (CEU)"),
            linewidth = 0.9, linetype = "dashed") +
  geom_line(data = decay_raw, aes(x = mid_kb, y = mean_r2, colour = "Simulated (observed)"),
            linewidth = 1.0) +
  geom_point(data = ceu_ref, aes(x = dist_kb, y = r2),
             shape = 21, fill = "#E64A19", colour = "white", size = 3) +
  geom_text(data = ceu_ref, aes(x = dist_kb, y = r2, label = label),
            vjust = -0.8, hjust = 0, size = 2.8, colour = "grey40") +
  scale_colour_manual(
    values = c("Theoretical (CEU)" = "#1565C0",
               "Simulated (observed)" = "#B71C1C"),
    name   = NULL
  ) +
  scale_x_continuous(limits = c(0, min(max(pos_kb) * 1.05, 500)),
                     labels = scales::comma) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  labs(
    title    = "LD decay — simulated vs. CEU reference",
    subtitle = sprintf("lambda = %.0f kb  |  background r² = %.2f  |  %d hotspots",
                       lambda, background_ld, n_hot),
    x        = "Physical distance (kb)",
    y        = expression(Mean ~ r^2)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(colour = "grey50", size = 9),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

# ── 10c. Save plots ──────────────────────────────────────────
fig_w <- max(8, n_snps * 0.18)

#ggsave("data/ld_triangle_CEU.png", plot = p_triangle,
#       width = fig_w, height = fig_w * 0.55, dpi = 300, bg = "white")

#ggsave("data/ld_decay_CEU.png", plot = p_decay,
#       width = 8, height = 5, dpi = 300, bg = "white")

#write.csv(as.data.frame(geno), "data/geno.csv", row.names = TRUE)
#write.csv(snp_info, "data/snp_info.csv", row.names = FALSE)
#write.csv(as.data.frame(r2_mat), "data/r2_mat.csv", row.names = TRUE)


print(p_decay)
print(p_triangle)

# ── 11. Ready-to-use objects ──────────────────────────────────
message("\nDone! Key objects:")
message("  geno      — integer matrix (samples × SNPs), 0/1/2 coded, CEU LD structure")
message("  metadata  — tibble: sample_id, population, ancestry")
message("  snp_info  — tibble: SNP positions, MAF, hotspot flags")
message("  R         — theoretical LD correlation matrix (SNPs × SNPs)")
message("  r2_mat    — observed r² matrix from simulated genotypes")
message("  ld_df     — observed r² by distance bin (for validation)")
