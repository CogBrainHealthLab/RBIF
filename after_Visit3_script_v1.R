require(tidyverse)
require(ggplot2)
require(kernlab)
require(ComBatFamily)

setwd("C:\\Users\\janic\\OneDrive - Nanyang Technological University\\projects\\reverse_brain_aging\\")

harm_output_dir <- "results\\after_Visit1\\v3\\lm_55-100_ebfalse_again\\"     # where Script 05's saved objects live

# read in task and rest datasets
task <- as.data.frame(readRDS('MRI\\V3\\RBIF_ses-2_task-agemanipulation_fslr32k_219_dtseries.rds'))
rest <- as.data.frame(readRDS('MRI\\V3\\RBIF_ses-2_task-rest_fslr32k_219_dtseries.rds'))

# setup output directory
output_dir <- "results\\after_Visit3\\v1"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define some variables
TR <- 0.719
counterbalance_order <- 1

##############################
### PREPROCESS RAW TASK FC ###
##############################

# Remove the first column (subject ID)
task <- task[, -1]

# Convert a time in seconds to the corresponding TR (row) index.
# Row 1 = t=0s, so index = floor(time / TR) + 1
sec_to_tr <- function(sec, TR = 0.719) {
  round(sec / TR) + 1
}

# --- Define segment boundaries in TR indices ---
set1_start <- sec_to_tr(0)     # 0
set1_end   <- sec_to_tr(270)   # 376
set2_start <- sec_to_tr(300)   # 417
set2_end   <- sec_to_tr(570)   # 793

cat("image_set_1: TR", set1_start, "to", set1_end,
    "(", (set1_end-set1_start+1), "timepoints )\n")
cat("image_set_2: TR", set2_start, "to", set2_end,
    "(", (set2_end-set2_start+1), "timepoints )\n")

# --- Segment the timeseries ---
image_set_1 <- task[set1_start:set1_end, ]
image_set_2 <- task[set2_start:set2_end, ]

# --- Rename based on counterbalancing order ---
# counterbalance_order should be 1 or 2 (per subject)
segment_by_counterbalance <- function(image_set_1, image_set_2, counterbalance_order) {
  if (counterbalance_order == 1) {
    young_set <- image_set_1
    old_set   <- image_set_2
  } else if (counterbalance_order == 2) {
    young_set <- image_set_2
    old_set   <- image_set_1
  } else {
    stop("counterbalance_order must be 1 or 2")
  }
  list(young_set = young_set, old_set = old_set)
}

result <- segment_by_counterbalance(image_set_1, image_set_2, counterbalance_order)
young_set <- result$young_set
old_set   <- result$old_set


# compute FC

fc_young <- cor(young_set, use = "pairwise.complete.obs")
fc_old   <- cor(old_set,   use = "pairwise.complete.obs")

upper_young <- fc_young[upper.tri(fc_young)]
upper_old   <- fc_old[upper.tri(fc_old)]


##############################
### PREPROCESS RAW REST FC ###
##############################

# Remove the first column (subject ID)
rest <- rest[, -1]

fc_rest <- cor(rest, use = "pairwise.complete.obs")
upper_rest <- fc_rest[upper.tri(fc_rest)]


################
### ANALYSIS ###
################


comfam_refs  <- readRDS(paste0(harm_output_dir, "comfam_joint_trainLFC_refs.rds"))
final_models <- readRDS(paste0(harm_output_dir, "final_models_jointcombat_harmonized.rds"))

rbif_sex <- 1     # <-- EDIT: 0 = M, 1 = F (must match train_sex coding)
rbif_age <- 67.75   # <-- EDIT: this subject's known chronological age

if (is.na(rbif_age)) stop("Set rbif_age before running.")

sex_label <- ifelse(rbif_sex == 0, "M", "F")
ref <- comfam_refs[[sex_label]]
mod <- final_models[[sex_label]]


static_FC <- rbind(upper_young, upper_old, upper_rest)
rownames(static_FC) <- c("young", "old", "rest")

if (ncol(static_FC) != length(ref$feature_names)) {
  stop(sprintf("static_FC has %d columns, reference expects %d -- ROI count or column order mismatch.",
               ncol(static_FC), length(ref$feature_names)))
}
colnames(static_FC) <- ref$feature_names

# ---------------------------------------------------------------------------
# Harmonize all 3 conditions in a single predict.comfam call
# ---------------------------------------------------------------------------
covar_static <- data.frame(outcome = rep(rbif_age, nrow(static_FC)))

harm_pred <- predict(
  ref$fit,
  newdata  = static_FC,
  newbat   = factor(rep("RBIF_visit3", nrow(static_FC))),
  newcovar = covar_static,
  eb       = FALSE      # <-- see caveat above; TRUE recommended given n=3
)

static_FC_harm <- harm_pred$dat.combat

# ---------------------------------------------------------------------------
# Apply trained models
# ---------------------------------------------------------------------------
gpr_preds <- kernlab::predict(mod$GPR_Linear, static_FC_harm)
kqr_preds <- kernlab::predict(mod$KQR_Linear, static_FC_harm)

visit3_results <- data.frame(
  condition      = rownames(static_FC),
  GPR_Linear_age = as.numeric(gpr_preds),
  KQR_Linear_age = as.numeric(kqr_preds),
  chron_age      = rbif_age
)

visit3_results$GPR_gap <- visit3_results$GPR_Linear_age - rbif_age
visit3_results$KQR_gap <- visit3_results$KQR_Linear_age - rbif_age

visit3_results

saveRDS(visit3_results, paste0(output_dir, "\\RBIF_visit3_static_bage_EBFALSE.rds"))
write.csv(visit3_results, paste0(output_dir, "\\RBIF_visit3_static_bage_EBFALSE.csv"), row.names = FALSE)












