require(kernlab)
require(ComBatFamily)
require(tidyverse)
require(ggplot2)

setwd("C:\\Users\\janic\\OneDrive - Nanyang Technological University\\projects\\reverse_brain_aging\\")
output_dir <- "results\\v3\\lm_55-100_ebfalse_again\\"

comfam_refs  <- readRDS(paste0(output_dir, "comfam_joint_trainLFC_refs.rds"))
final_models <- readRDS(paste0(output_dir, "final_models_jointcombat_harmonized.rds"))

# ---------------------------------------------------------------------------
# RBIF subject info -- EDIT these for your subject
# ---------------------------------------------------------------------------
rbif_sex <- 1     # <-- EDIT: 0 = M, 1 = F (must match train_sex coding)
rbif_age <- 67.75   # <-- EDIT: this subject's known chronological age

if (is.na(rbif_age)) stop("Set rbif_age before running -- needed for the ComBat covariate and for computing brain-age gap.")

sex_label <- ifelse(rbif_sex == 0, "M", "F")
ref <- comfam_refs[[sex_label]]
mod <- final_models[[sex_label]]

# ---------------------------------------------------------------------------
# Load RBIF movie timeseries + compute sliding-window dFC
# (your extraction code, unchanged)
# ---------------------------------------------------------------------------
movie <- as.data.frame(readRDS('MRI\\RBIF_ses-1_task-movie_fslr32k_219_dtseries.rds'))

# --- Setup ---
roi_data <- movie[, -1]          # drop the subject ID column
roi_data <- as.matrix(roi_data)  # 615 TRs x 219 ROIs
n_roi <- ncol(roi_data)          # 219

# Indices for the upper triangle (excluding diagonal)
upper_idx <- which(upper.tri(matrix(0, n_roi, n_roi)), arr.ind = FALSE)

# =========================================================
# OPTION A: Single static FC matrix (whole scan, one matrix)
# =========================================================
FC_static <- cor(roi_data)              # 219 x 219
FC_static_upper <- FC_static[upper.tri(FC_static)]  # vector, length = 219*218/2 = 23871

# =========================================================
# OPTION B: Dynamic FC time series (sliding window)
# =========================================================
window_len <- 30    # TRs per window 
step <- 1            # slide by 1 TR 

n_tr <- nrow(roi_data)
starts <- seq(1, n_tr - window_len + 1, by = step)

n_edges <- n_roi * (n_roi - 1) / 2
dFC <- matrix(NA, nrow = length(starts), ncol = n_edges)

for (i in seq_along(starts)) {
  s <- starts[i]
  e <- s + window_len - 1
  win_data <- roi_data[s:e, ]
  fc_mat <- cor(win_data)
  dFC[i, ] <- fc_mat[upper.tri(fc_mat)]
}

# dFC is now: (number of windows) x (23871 upper-triangle edges)
dim(dFC)


# ---------------------------------------------------------------------------
# Feature alignment -- POSITIONAL assumption, same caveat as earlier scripts:
# dFC's column order (from fc_mat[upper.tri(fc_mat)]) must match the order
# used when training features were built. Verify before trusting downstream
# results -- see the "FEATURE ALIGNMENT" note in 02_predict_RBIF_slidingwindow.R
# ---------------------------------------------------------------------------
if (ncol(dFC) != length(ref$feature_names)) {
  stop(sprintf("dFC has %d columns, reference expects %d -- ROI count or column order mismatch.",
               ncol(dFC), length(ref$feature_names)))
}
colnames(dFC) <- ref$feature_names

# ---------------------------------------------------------------------------
# Harmonize RBIF windows out-of-sample via predict.comfam
# ---------------------------------------------------------------------------
# verify exact argument names for your installed ComBatFamily version:
args(ComBatFamily:::predict.comfam)

covar_rbif <- data.frame(outcome = rep(rbif_age, nrow(dFC)))

# this is for all windows 
harm_pred <- predict(
  ref$fit,
  newdata  = dFC,
  newbat   = factor(rep("RBIF", nrow(dFC))),
  newcovar = covar_rbif,
  eb       = FALSE,
  model   = "lm",
  formula = y ~ outcome,
)
# 
# # this is for one window
# harm_pred <- predict(
#   ref$fit,
#   newdata  = dFC[1, , drop = FALSE],
#   newbat   = factor(rep("RBIF", nrow(dFC[1, , drop = FALSE]))),
#   newcovar = data.frame(outcome = rep(rbif_age, nrow(dFC[1, , drop = FALSE])))
# )

dFC_harm <- harm_pred$dat.combat

# ---------------------------------------------------------------------------
# Apply trained models
# ---------------------------------------------------------------------------
gpr_preds  <- kernlab::predict(mod$GPR_Linear, dFC_harm)
kqr_preds  <- kernlab::predict(mod$KQR_Linear, dFC_harm)

rbif_results <- data.frame(
  window_start_TR = starts,
  window_mid_TR   = starts + (window_len - 1) / 2,
  GPR_Linear_age  = as.numeric(gpr_preds),
  KQR_Linear_age  = as.numeric(kqr_preds),
  chron_age       = rbif_age
)

rbif_results$GPR_gap  <- rbif_results$GPR_Linear_age - rbif_age
rbif_results$KQR_gap  <- rbif_results$KQR_Linear_age - rbif_age

# ---------------------------------------------------------------------------
# Performance summary -- MAE / gap stats only (see caveat #2 at top of file)
# ---------------------------------------------------------------------------
summarize_gap <- function(gap, label) {
  cat(sprintf("%-12s MAE=%.3f  mean gap=%+.3f  SD gap=%.3f\n",
              label, mean(abs(gap)), mean(gap), sd(gap)))
}
cat("\n--- RBIF single-subject brain-age summary (across", nrow(rbif_results), "windows) ---\n")
summarize_gap(rbif_results$GPR_gap,  "GPR Linear")
summarize_gap(rbif_results$KQR_gap,  "KQR Linear")
cat("(r / bias not reported -- undefined for a single subject's constant true age)\n")

saveRDS(rbif_results, paste0(output_dir, "RBIF_bage_trajectory_jointcombat.rds"))
write.csv(rbif_results, paste0(output_dir, "RBIF_bage_trajectory_jointcombat.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Quick plot
# ---------------------------------------------------------------------------
rbif_results %>%
  pivot_longer(cols = c(GPR_Linear_age, KQR_Linear_age),
               names_to = "model", values_to = "predicted_age") %>%
  ggplot(aes(x = window_start_TR, y = predicted_age, color = model)) +
  geom_line() +
  geom_hline(yintercept = rbif_age, linetype = "dashed", color = "grey40") +
  labs(title = "RBIF brain-age trajectory (movie task, joint-ComBat harmonized)",
       x = "Window number (TR)", y = "Predicted brain age") +
  theme_minimal()


# ---------------------------------------------------------------------------
# Determine windows
# ---------------------------------------------------------------------------

RBIF_results <- readRDS(paste0(output_dir, "RBIF_bage_trajectory_jointcombat.rds"))

find_top_n <- function(df, value_col, n = 3, min_sep = 30, decreasing = TRUE) {
  d <- df[order(df[[value_col]], decreasing = decreasing), ]
  selected <- data.frame()
  for (i in seq_len(nrow(d))) {
    if (nrow(selected) == 0 || all(abs(d$window_mid_TR[i] - selected$window_mid_TR) >= min_sep)) {
      selected <- rbind(selected, d[i, ])
    }
    if (nrow(selected) == n) break
  }
  selected
}

top3_highest_KQR <- find_top_n(RBIF_results, "KQR_gap", n = 3, min_sep = 30, decreasing = TRUE)
top3_lowest_KQR  <- find_top_n(RBIF_results, "KQR_gap", n = 3, min_sep = 30, decreasing = FALSE)

# bracket the peak (Shifted 2 seconds in advance to account for hemodynamic lags)
top3_highest_KQR$range_low  <- top3_highest_KQR$window_mid_TR - 16.5
top3_highest_KQR$range_high <- top3_highest_KQR$window_mid_TR + 12.5

top3_lowest_KQR$range_low  <- top3_lowest_KQR$window_mid_TR - 16.5
top3_lowest_KQR$range_high <- top3_lowest_KQR$window_mid_TR + 12.5

# in terms of mins and seconds
TR <- 0.719

top3_highest_KQR$range_low_sec  <- top3_highest_KQR$range_low  * TR
top3_highest_KQR$range_high_sec <- top3_highest_KQR$range_high * TR

top3_lowest_KQR$range_low_sec  <- top3_lowest_KQR$range_low  * TR
top3_lowest_KQR$range_high_sec <- top3_lowest_KQR$range_high * TR

# convert seconds to MM:SS for readability
format_mmss <- function(sec) {
  m <- floor(sec / 60)
  s <- round(sec %% 60)
  sprintf("%d:%02d", m, s)
}

top3_highest_KQR$range_low_mmss  <- format_mmss(top3_highest_KQR$range_low_sec)
top3_highest_KQR$range_high_mmss <- format_mmss(top3_highest_KQR$range_high_sec)

top3_lowest_KQR$range_low_mmss  <- format_mmss(top3_lowest_KQR$range_low_sec)
top3_lowest_KQR$range_high_mmss <- format_mmss(top3_lowest_KQR$range_high_sec)

top3_highest_KQR
top3_lowest_KQR




# ---------------------------------------------------------------------------
# Apply trained models -- HARMONIZED features
# ---------------------------------------------------------------------------
gpr_preds  <- kernlab::predict(mod$GPR_Linear, dFC_harm)
kqr_preds  <- kernlab::predict(mod$KQR_Linear, dFC_harm)

# ---------------------------------------------------------------------------
# Apply the SAME trained models to the RAW (unharmonized) windows
# to check whether harmonization preserves the movie-related dynamic signal
# ---------------------------------------------------------------------------
gpr_preds_raw <- kernlab::predict(mod$GPR_Linear, dFC)
kqr_preds_raw <- kernlab::predict(mod$KQR_Linear, dFC)

rbif_results <- data.frame(
  window_start_TR    = starts,
  window_mid_TR       = starts + (window_len - 1) / 2,
  GPR_Linear_age       = as.numeric(gpr_preds),
  KQR_Linear_age       = as.numeric(kqr_preds),
  GPR_Linear_age_raw   = as.numeric(gpr_preds_raw),
  KQR_Linear_age_raw   = as.numeric(kqr_preds_raw),
  chron_age            = rbif_age
)

rbif_results$GPR_gap  <- rbif_results$GPR_Linear_age - rbif_age
rbif_results$KQR_gap  <- rbif_results$KQR_Linear_age - rbif_age

rbif_results$GPR_AE   <- abs(rbif_results$GPR_gap)   # per-window absolute error ("MAE" for n=1)
rbif_results$KQR_AE   <- abs(rbif_results$KQR_gap)

rbif_results$KQR_raw_gap  <- rbif_results$KQR_Linear_age_raw - rbif_age
mean(abs(rbif_results$KQR_raw_gap))

# ---------------------------------------------------------------------------
# Performance summary -- MAE / gap stats only 
# ---------------------------------------------------------------------------
summarize_gap <- function(gap, label) {
  cat(sprintf("%-12s MAE=%.3f  mean gap=%+.3f  SD gap=%.3f\n",
              label, mean(abs(gap)), mean(gap), sd(gap)))
}
cat("\n--- RBIF single-subject brain-age summary (across", nrow(rbif_results), "windows) ---\n")
summarize_gap(rbif_results$GPR_gap,  "GPR Linear")
summarize_gap(rbif_results$KQR_gap,  "KQR Linear")
cat("(r / bias not reported -- undefined for a single subject's constant true age)\n")

# ---------------------------------------------------------------------------
# Harmonized vs Unharmonized trajectory agreement
#
# High r  -> the dynamic (movie-related) signal survives harmonization;
#            harmonization is mainly shifting/rescaling the baseline, not
#            distorting the shape.
# Low r   -> harmonization is altering the window-to-window pattern itself,
#            which would be a bigger concern than a simple level shift.
# ---------------------------------------------------------------------------
r_GPR <- cor(rbif_results$GPR_Linear_age, rbif_results$GPR_Linear_age_raw)
r_KQR <- cor(rbif_results$KQR_Linear_age, rbif_results$KQR_Linear_age_raw)

cat("\n--- Harmonized vs Unharmonized trajectory shape agreement ---\n")
cat(sprintf("GPR Linear: r = %.3f\n", r_GPR))
cat(sprintf("KQR Linear: r = %.3f\n", r_KQR))

saveRDS(rbif_results, paste0(output_dir, "RBIF_bage_trajectory_jointcombat.rds"))
write.csv(rbif_results, paste0(output_dir, "RBIF_bage_trajectory_jointcombat.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Comparison plot: harmonized vs unharmonized trajectory, per model
# ---------------------------------------------------------------------------
rbif_results %>%
  select(window_mid_TR, GPR_Linear_age, GPR_Linear_age_raw, KQR_Linear_age, KQR_Linear_age_raw) %>%
  pivot_longer(cols = -window_mid_TR, names_to = "series", values_to = "predicted_age") %>%
  mutate(
    model         = ifelse(grepl("GPR", series), "GPR Linear", "KQR Linear"),
    harmonization = ifelse(grepl("raw", series), "Unharmonized", "Harmonized")
  ) %>%
  ggplot(aes(x = window_mid_TR, y = predicted_age, color = harmonization)) +
  geom_line() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  labs(title = "RBIF brain-age trajectory: harmonized vs unharmonized",
       subtitle = sprintf("Shape agreement -- GPR r=%.2f | KQR r=%.2f", r_GPR, r_KQR),
       x = "Window midpoint (TR)", y = "Predicted brain age") +
  theme_minimal()