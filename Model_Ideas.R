


# Modelling examples

## Support
formula_support_emo <- bf(
  mvpa ~ (support_emo_cw + support_emo_cb) * rel_type + 
    (1 + support_emo_cw | studyID/personID)
)

model_support_emo <- brm(
  formula_support_emo, 
  #prior = prior_support_emo,
  data = df_combined, 
  family = gaussian(), 
  chains = 4, 
  iter = 2000, 
  warmup = 1000, 
  cores = 4, 
  backend = 'cmdstan',
  file = file.path('models_cache', 'support_emo')
)

summary(model_support_emo)
check_brms(model_support_emo)


## Control 
formula_control <- bf(
  mvpa ~ (psc_cw + psc_cb + nsc_cw + nsc_cb) * rel_type +
    (1 + psc_cw + nsc_cw | studyID/personID)
)

model_control <- brm(
  formula_control, 
  #prior = prior_control,
  data = df_combined, 
  family = gaussian(), 
  chains = 4, 
  iter = 2000, 
  warmup = 1000, 
  cores = 4, 
  backend = 'cmdstan',
  file = file.path('models_cache', 'control')
)



