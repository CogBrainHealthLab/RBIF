#############################################################################
# 05_combat_joint_harmonization.R
#
# Fits sex-specific ComBat jointly on the 4 training sites + LFC combined.
# Unlike earlier versions, the fitted comfam object is SAVED -- it becomes
# the fixed reference used later to harmonize RBIF out-of-sample via
# predict.comfam() in Script 06.
#############################################################################

require(tidyverse)
require(ggplot2)
require(ComBatFamily)
require(kernlab)
require(pls)

setwd("C:\\Users\\janic\\OneDrive - Nanyang Technological University\\projects\\reverse_brain_aging\\")
output_dir <- "results\\v3\\lm_55-100_ebfalse_again\\"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
train_mri <- as.data.frame(readRDS('MRI\\for_training\\brainageFCtrain.rds'))
train_mri_agelimited <- train_mri[train_mri$age >= 55 & train_mri$age <= 100, ]

test_mri <- as.data.frame(readRDS('MRI\\for_training\\LFC\\LFC_FC219.rds'))
test_beh <- read.csv('MRI\\for_training\\LFC\\LFC_behdata_cleaned.csv')
test_beh_cleaned <- test_beh[test_beh$rsFC == 'yes',]
MCI_idx <- test_beh_cleaned$Diag == "Normal Ageing"
test_beh_cleaned_normal <- test_beh_cleaned[MCI_idx, ]
test_mri_normal <- test_mri[MCI_idx, ]

train_outcome <- train_mri_agelimited$age
train_sex     <- train_mri_agelimited$sex     # 0 = M, 1 = F
train_site    <- train_mri_agelimited$site
train_feat    <- as.matrix(train_mri_agelimited[, !(names(train_mri_agelimited) %in% c("age","sex","site"))])

test_outcome <- test_beh_cleaned_normal$age
test_sex     <- ifelse(test_beh_cleaned_normal$Sex == 1, 0, 1)
test_feat    <- test_mri_normal
colnames(test_feat) <- colnames(train_feat)
test_feat    <- as.matrix(test_feat)

# ---------------------------------------------------------------------------
# Zero-column check (hard stop instead of silent drop)
# ---------------------------------------------------------------------------
check_zero_cols <- function(feat, label) {
  zero_cols <- which(colSums(feat) == 0)
  if (length(zero_cols) > 0) {
    stop(sprintf("%s contains %d all-zero column(s): %s",
                 label, length(zero_cols), paste(head(colnames(feat)[zero_cols], 10), collapse = ", ")))
  }
  cat(sprintf("%s: OK, %d features checked.\n", label, ncol(feat)))
}
check_zero_cols(train_feat, "train_feat")
check_zero_cols(test_feat,  "test_feat")
stopifnot(identical(colnames(train_feat), colnames(test_feat)))

# =============================================================================
# STEP A: Fit sex-specific ComBat JOINTLY on train (4 sites) + LFC
#         -> saved as the fixed reference for RBIF's predict.comfam later
# =============================================================================
sex_levels  <- list(M = 0, F = 1)
comfam_refs <- list()
harmonized_list <- list()

for (s in names(sex_levels)) {
  
  train_idx <- which(train_sex == sex_levels[[s]])
  test_idx  <- which(test_sex  == sex_levels[[s]])
  
  feat_combined    <- rbind(train_feat[train_idx, , drop = FALSE], test_feat[test_idx, , drop = FALSE])
  outcome_combined <- c(train_outcome[train_idx], test_outcome[test_idx])
  site_combined    <- c(as.character(train_site[train_idx]), rep("LFC", length(test_idx)))
  role_combined    <- c(rep("train", length(train_idx)), rep("test", length(test_idx)))
  
  covar_s <- data.frame(outcome = outcome_combined)
  
  cat(sprintf("\nFitting joint ComBat: sex=%s, n=%d (train=%d, LFC=%d), sites=%s\n",
              s, nrow(feat_combined), length(train_idx), length(test_idx),
              paste(unique(site_combined), collapse = ", ")))
  
  comfam_fit <- ComBatFamily::comfam(
    data    = feat_combined,
    bat     = as.factor(site_combined),
    covar   = covar_s,
    model   = "lm",
    formula = y ~ outcome,
    eb      = FALSE       # <-- default true
  )
  
  comfam_refs[[s]] <- list(fit = comfam_fit, feature_names = colnames(feat_combined))
  
  harmonized_list[[s]] <- data.frame(
    comfam_fit$dat.combat,
    age  = outcome_combined,
    sex  = sex_levels[[s]],
    site = site_combined,
    role = role_combined
  )
}

harmonized_all <- dplyr::bind_rows(harmonized_list)

# This is the object Script 06 needs -- DO NOT discard it (unlike pred.allmodels.bysex,
# which fits-and-throws-away comfam objects internally)
saveRDS(comfam_refs, paste0(output_dir, "comfam_joint_trainLFC_refs.rds"))
saveRDS(harmonized_all, paste0(output_dir, "harmonized_trainLFC.rds"))

train_harmonized <- harmonized_all[harmonized_all$role == "train", ]
test_harmonized  <- harmonized_all[harmonized_all$role == "test", ]

# =============================================================================
# STEP B: Benchmark all models on this harmonized train/LFC split
#         (harm = FALSE since harmonization is already done above)
# =============================================================================
source("https://github.com/CogBrainHealthLab/MLtools/blob/main/allregmodelsbysex2.R?raw=TRUE")

feat_cols <- comfam_refs$M$feature_names
train_feat_harm <- as.matrix(train_harmonized[, feat_cols])
test_feat_harm  <- as.matrix(test_harmonized[,  feat_cols])

results_jointcombat <- pred.allmodels.bysex(
  train_outcome = train_harmonized$age,
  train_feat    = train_feat_harm,
  train_sex     = train_harmonized$sex,
  test_outcome  = test_harmonized$age,
  test_feat     = test_feat_harm,
  test_sex      = test_harmonized$sex,
  harm          = FALSE,
  xgb           = FALSE     # set TRUE once xgboost/ParBayesianOptimization are working
)

saveRDS(results_jointcombat, paste0(output_dir, "results_jointcombat.rds"))
results_jointcombat$predmetrics.all
plot.metrics(results_jointcombat$predmetrics.all)


# =============================================================================
# STEP C: Refit + save final GPR-Linear / KQR-Linear models on the
#         harmonized TRAIN data (per sex) -- these are what Script 06 applies
#         to RBIF's sliding windows
# =============================================================================
final_models <- list()

for (s in names(sex_levels)) {
  
  idx       <- which(train_harmonized$sex == sex_levels[[s]])
  feat_s    <- train_feat_harm[idx, , drop = FALSE]
  outcome_s <- train_harmonized$age[idx]
  
  gpr_mod <- kernlab::gausspr(feat_s, outcome_s, kernel = "vanilladot")
  kqr_mod <- kernlab::kqr(feat_s, outcome_s, kernel = "vanilladot")
  
  cat(sprintf("Final models (sex=%s): n=%d\n", s, nrow(feat_s)))
  
  final_models[[s]] <- list(
    GPR_Linear    = gpr_mod,
    KQR_Linear    = kqr_mod,
    feature_names = feat_cols
  )
}

saveRDS(final_models, paste0(output_dir, "final_models_jointcombat_harmonized.rds"))
cat("\nSaved: comfam_joint_trainLFC_refs.rds, final_models_jointcombat_harmonized.rds\n")

