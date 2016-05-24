conditionalBootstrap <-
function(object, BootStrRep) {
  # A function that calculates the bias correction for a (generalized) linear 
  # mixed models by the methods in Efron (2004).
  #
  # Args: 
  #   object     = Object of class lmerMod or glmerMod. Obtained by lmer() or 
  #                glmer().
  #   BootStrRep = Number of bootstrap replications.
  #
  # Returns:
  #   bootBC = Bias correction (i.e. degrees of freedom) for a (generalized) 
  #            linear mixed model.
  #  
	dataMatrix    <- simulate(object, nsim = BootStrRep, use.u = TRUE)
  workingEta    <- apply(dataMatrix, 2, function(x) predict(refit(object, newresp = x)))
	dataMatrix    <- dataMatrix - rowMeans(dataMatrix)
	bootBC        <- sum(workingEta * dataMatrix) / ((BootStrRep - 1) * sigma(object)^2)
	return(bootBC)
}
