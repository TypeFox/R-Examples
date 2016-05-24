getcondLL <-
function(object) {
  # A function that calls the bias correction functions.
  #
  # Args: 
  #   object = Object of class lmerMod or glmerMod. Obtained by lmer() or glmer()
  #
  # Returns:
  #   cll    = The conditional log-likelihood.
  #
  switch(family(object)$family,
    binomial = {
      cll <- sum(dbinom(x    = getME(object, "y"), 
                        size = length(unique(getME(object, "y"))), 
                        prob = getME(object, "mu"), log = TRUE))
    },
    poisson  = {
      cll <- sum(dpois(x      = getME(object, "y"), 
                       lambda = getME(object, "mu"), log = TRUE))
    },
    gaussian = {
      cll <- sum(dnorm(x    = getME(object, "y"), 
                       mean = getME(object, "mu"), 
                       sd   = sigma(object), log = TRUE))
    },
    {
      cat("For this family no bias correction is currently available \n")
      cll <- NA
    }
    )
  return(cll)
}
