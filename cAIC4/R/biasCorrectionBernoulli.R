biasCorrectionBernoulli <-
function(object) {
  # A function that calculates the bias correction for a generalized linear 
  # mixed models with binary(!) data similar to the centralized Steinian method
  # in Efron (2004).
  #
  # Args: 
  #   object = Object of class lmerMod or glmerMod. Obtained by lmer() or 
  #            glmer(). Needs binary data.
  #
  # Returns:
  #   BC     = (Asymptotic) bias correction (i.e. degrees of freedom) for a 
  #            (generalized) linear mixed model with binary response.
  #  
  y                   <- object@resp$y
	signCor             <- - 2 * y + 1
  mu                  <- object@resp$mu
	eta                 <- qlogis(mu)
  workingMatrix       <- matrix(rep(y, length(y)), ncol = length(y))
  diag(workingMatrix) <- 1 - diag(workingMatrix)
  workingEta          <- diag(apply(workingMatrix, 2, function(x) qlogis(refit(object, newresp = x)@resp$mu) - eta))
	return(sum(mu * (1 - mu) * signCor * workingEta))
}
