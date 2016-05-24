biasCorrectionGaussian <-
function(m, sigma.estimated, analytic) {
  # A function that calls the bias correction functions.
  #
  # Args: 
  #   mer    = Object of class lmerMod. Obtained by lmer().
  #   sigma.estimated = If sigma is estimated. This only is used for the 
  #                     analytical version of Gaussian responses.
  #   analytic = FALSE if the numeric hessian of the (restricted) marginal log-
  #              likelihood from the lmer optimization procedure should be used.
  #              Otherwise (default) TRUE, i.e. use a analytical version that 
  #              has to be computed.
  #
  # Returns:
  #   bc = Bias correction for a mixed model.
  #
  zeroLessModel <- deleteZeroComponents(m)
  if (inherits(zeroLessModel, "lm")) {
    return(ncol(model.matrix(zeroLessModel)))
  }
  model <- getModelComponents(zeroLessModel, analytic)
  if (identical(m, zeroLessModel)) {
    bc       <- calculateGaussianBc(model, sigma.estimated, analytic)
    newModel <- NULL
    new      <- FALSE
  } else {
    bc       <- calculateGaussianBc(model, sigma.estimated, analytic)
    newModel <- zeroLessModel
    new      <- TRUE
  }
  return(list(bc = bc, newModel = newModel, new = new))
}
