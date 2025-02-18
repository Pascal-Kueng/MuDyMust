# Load required packages
library(brms)

# Define the custom family
hurdle_gaussian <- custom_family(
  "hurdle_gaussian",
  dpars = c("mu", "sigma", "hu"),
  links = c("identity", "log", "logit"),
  lb = c(NA, 0, NA),
  type = "real"
)


# Stan functions for the custom family
stan_funs <- "
  real hurdle_gaussian_lpdf(real y, real mu, real sigma, real hu) {
    if (y == 0) {
      return bernoulli_lpmf(1 | hu);
    } else {
      return bernoulli_lpmf(0 | hu) + normal_lpdf(y | mu, sigma);
    }
  }
  
  real hurdle_gaussian_rng(real mu, real sigma, real hu) {
    if (bernoulli_rng(hu) == 1) {
      return 0;
    } else {
      return normal_rng(mu, sigma);
    }
  }
"

# Include the Stan functions in the model
stanvars <- stanvar(scode = stan_funs, block = "functions")

# Posterior predictions
posterior_predict_hurdle_gaussian <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  hu <- brms::get_dpar(prep, "hu", i = i)
  draws <- length(mu)
  
  # Draw from uniform to decide if zero or from normal
  hu_draws <- runif(draws)
  ifelse(hu_draws < hu, 0, rnorm(draws, mu, sigma))
}

# Posterior expected values
posterior_epred_hurdle_gaussian <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  hu <- brms::get_dpar(prep, "hu")
  # Expected value: E(Y) = (1 - hu)*mu
  mu * (1 - hu)
}

# Log-likelihood calculations (for loo)
log_lik_hurdle_gaussian <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  hu <- brms::get_dpar(prep, "hu", i = i)
  y <- prep$data$Y[i]
  
  if (y == 0) {
    log(hu)
  } else {
    log(1 - hu) + dnorm(y, mean = mu, sd = sigma, log = TRUE)
  }
}

# Wrapper around bf()
bf <- function(..., family = NULL, custom_stanvars = stanvars) {
  # If family is our custom hurdle_gaussian, attach the stanvars to the formula
  if (!is.null(family) && identical(family, hurdle_gaussian)) {
    f <- brms::bf(..., family = family)
    attr(f, "stanvars") <- custom_stanvars
    f
  } else {
    # Otherwise just call bf as normal
    brms::bf(..., family = family)
  }
}

# Wrapper around brm()
brm <- function(formula, data, family = NULL, ...) {
  # If the formula is a brmsformula with stanvars attached, use them
  attached_stanvars <- attr(formula, "stanvars")
  if (!is.null(attached_stanvars)) {
    brms::brm(
      formula = formula,
      data = data,
      family = family, 
      stanvars = attached_stanvars,
      ...
    )
  } else {
    # Otherwise call brms::brm normally
    brms::brm(
      formula = formula,
      data = data,
      family = family,
      ...
    )
  }
}
