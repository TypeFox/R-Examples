bcMer <-
function(mer, method = NULL, B = NULL, sigma.estimated = TRUE, analytic = TRUE) {
  # A function that calls the bias correction functions.
  #
  # Args: 
  #   mer    = Object of class lmerMod or glmerMod. Obtained by lmer() or glmer()
  #   method = How the bias correction should be evaluated. If NULL than method 
  #            is chosen by family, i.e. analytical if family is Poisson or 
  #            Gaussian and with parametric bootstrap for other. Method may also
  #            be specified before, either "steinian" or "conditionalBootstrap".
  #            "steinian" only available for Gaussian, Poisson and Bernoulli.
  #   B      = Number of Bootstrap replications. Default is NULL then it is 
  #            chosen as maximum of the number of observations and 100.
  #
  # Returns:
  #   bc = Bias correction for a mixed model.
  #
  if (is.null(method)) {
    switch(family(mer)$family,
          binomial = {
            if (is.null(B)) {
              B <- max(length(getME(mer, "y")), 100)
            }
            bc <- conditionalBootstrap(mer, B)
          },
          poisson = {
            bc <- biasCorrectionPoisson(mer)
          },
          gaussian = {
            bc <- biasCorrectionGaussian(mer, sigma.estimated, analytic)
          },
          {
            cat("For this family no bias correction is currently available \n")
            bc <- NA
          }
          )
  } else {
    if(method == "steinian") {
      switch(family(mer)$family,
            binomial = {
              bc <- biasCorrectionBernoulli(mer)
            },
            poisson = {
              bc <- biasCorrectionPoisson(mer)
            },
            gaussian = {
              bc <- biasCorrectionGaussian(mer, sigma.estimated, analytic)
            },
            {
              cat("For this family no bias correction is currently available \n")
              bc <- NA
            }
            )
    }
    if(method == "conditionalBootstrap") {
      if (is.null(B)) {
        B <- max(length(getME(mer, "y")), 100)
      }
      bc <- conditionalBootstrap(mer, B)
    }
  }
  return(bc)
}
