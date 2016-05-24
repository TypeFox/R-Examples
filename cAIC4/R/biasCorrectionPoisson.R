biasCorrectionPoisson <-
function(object) {
  # A function that calculates the bias correction for a generalized linear 
  # mixed models with Poisson data, see Lian (2012) & Saefken et al. (2014).
  #
  # Args: 
  #   object = Object of class lmerMod or glmerMod. Obtained by glmer(). With 
  #            family = "poisson".
  # Returns:
  #   BC     = Bias correction (i.e. degrees of freedom) for a (generalized) 
  #            linear mixed model with Poisson response.
  #
  y                   <- object@resp$y
  ind                 <- which(y != 0)
	workingMatrix       <- matrix(rep(y, length(y)), ncol = length(y))
	diag(workingMatrix) <- diag(workingMatrix) - 1
	workingMatrix       <- workingMatrix[, ind]
	workingEta          <- diag(apply(workingMatrix, 2, function(x) refit(object,
                                    newresp = x)@resp$eta)[ind,])
	return(sum(y[ind] * (object@resp$eta[ind] - workingEta)))
}
