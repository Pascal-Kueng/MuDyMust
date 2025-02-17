


df_kamp <- readRDS('Datasets/kamp.rds')
df_tt <- readRDS('Datasets/t&t.rds')



library(brms)

# Outcome model: mvpa is modeled using a non-linear formula that decomposes social_support
bf_mvpa <- bf(
  mvpa ~ beta0 + beta_w * (social_support - mu_person) + beta_b * mu_person + beta_rel * rel_type +
    (1 | studyID) + (1 | studyID:personID),
  beta0 + beta_w + beta_b + beta_rel ~ 1,
  nl = TRUE
)

# Support model: we model social_support as the sum of the latent person mean and a day-specific deviation,
# with the constraint that the deviation has mean zero.
bf_support <- bf(
  social_support ~ mu_person + delta,
  mu_person ~ 1 + (1 | studyID) + (1 | studyID:personID),
  delta ~ 0,  # This forces the average deviation within each person to be zero.
  nl = TRUE
)

joint_model <- bf_mvpa + bf_support


##################################### Easy alternative ####################

bf_mvpa_control <- bf(
  mvpa ~ (support_cw + support_cb) * rel_type + 
    (1 + support_cw | studyId/personId)
)

bf_mvpa_support <- bf(
  mvpa ~ (control_cw + control_cb) * rel_type +
    (1 + control_cw | studyId/personId)
)

