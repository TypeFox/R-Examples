require(coda)


################################################################################
# 
# @brief General model adequacy test using posterior predictive testing. Prior
#        predictive testing can be achieved by providing samples from the prior
#        instead.
#
# @date Last modified: 2015-05-28
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-11-19, version 1.0
#
# @param    simulationFunction     function      the simulation function
# @param    parameters             matrix        set of parameter samples
# @param    burnin                 scalar        the fraction of samples to burn
# @return                          list          the simlated trees
#
################################################################################


tess.PosteriorPrediction <- function(simulationFunction,parameters,burnin=0.25) {


  samples <- list()
  b <- length(parameters[,1]) * burnin
  for ( i in b:length(parameters[,1])) {

    # get the current set of parameter values
    theta <- parameters[i,]

    # simulate a new observation under the current parameter values
    samples[[i-b+1]] <- simulationFunction(theta)

  }

  return (samples)
}
