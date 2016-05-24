# Fit a three-level model with the Laplace approximation to the likelihood
(mod_Laplace <- glmm(response ~ covariate + (1 | cluster) + (1 | group),
                     data = three_level, family = binomial,
                     method = "Laplace"))

# if we try to fit with adaptive Gaussian quadrature, we get an error
\dontrun{
  (mod_AGQ <- glmm(response ~ covariate + (1 | cluster) + (1 | group),
                   data = three_level, family = binomial, method = "AGQ",
                   control = list(nAGQ = 15)))
}

# We can fit with the Sequential Reduction approximation
\dontrun{
  (mod_SR <- glmm(response ~ covariate + (1 | cluster) + (1 | group),
                  data = three_level, family = binomial, method = "SR",
                  control = list(nSL = 3)))
}
# the estimates of the random effects standard deviations
# are larger than those using the Laplace approximation
