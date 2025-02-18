library(brms)
library(tidyverse)
library(bmlm)
#install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
#cmdstanr::check_cmdstan_toolchain(fix = TRUE)
#cmdstanr::install_cmdstan()
library(cmdstanr)

source(file.path('Functions', 'PrettyTables.R'))
source(file.path('Functions', 'ReportModels.R'))
source(file.path('Functions', 'ReportMeasures.R'))
source(file.path('Functions', 'HurdleGaussianFamily.R'))



df_kamp <- readRDS('Datasets/kamp.rds')
df_tt <- readRDS('Datasets/t&t.rds')

# creating necesssary variables

df_kamp <- df_kamp %>%
  rename(
    personID = pid, 
    #day = day, 
    mvpa_10h = mvpa, # !todo: use accelerometer data or self-report? What is more comparable? I think self-report unless all use the same acceleromters and algorightms. 
    mvpa = TotalMVPA,
    support_emo = soz_emo,
    support_inst = soz_prak
  ) %>%
  mutate(
    studyID = 1,
    rel_type = 1, 
    psc = NA, # !todo: set to NA as a test of how models handle it
    nsc = NA, # !todo: set to NA as a test of how models handle it
    muPerson = NA
  ) %>%
  select(
    personID, studyID, rel_type, day, mvpa, support_emo, support_inst, psc, nsc
  ) 

df_tt <- df_tt %>%
  rename(
    personID = userID, # !todo: how should we deal with dyadic datasets?? use coupleId or select one participant per couple?
    #day = day, # !todo: should we center day for each study, so that it is between 0 and 1? Probably not becuase they have different time-spans. 
    psc = ss_psc_more,
    nsc = ss_nsc_more,
    support_emo = ss_emo_pleasure,
    support_inst = ss_inst_prov # !todo: which variable?
  ) %>% 
  mutate(
    studyID = 2,
    rel_type = 2,
    mvpa = ss_pa_min_solo + ss_pa_min_collab,
    muPerson = NA
  ) %>%
  select(
    personID, studyID, rel_type, day, mvpa, support_emo, support_inst, psc, nsc
  )

# combine datasets
df_combined <- df_kamp %>%
  rbind(df_tt) %>%
  bmlm::isolate(
    by = 'personID', 
    value = c('support_emo', 'support_inst', 'psc', 'nsc'),
    which = 'both'
  )


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



