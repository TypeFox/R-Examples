

################################################################################
# 
# @brief General model adequacy test using posterior predictive testing. Prior
#        predictive testing can be achieved by providing samples from the prior
#        instead.
#
# @date Last modified: 2012-12-17
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-12-17, version 1.0
#
# @param    samples                list          Samples from the posterior predictive distribution
# @param    observation            any           the observed value
# @param    statistic              function      the function computing the statistic
# @return                          scalar        the upper quantile of observing such a value
#
################################################################################


tess.PosteriorPredictiveTest <- function(samples,observation,statistic) {

  obs <- statistic(observation)

  sampled_statistics <- c()
  count <- 0
  for ( i in 1:length(samples)) {

    # compute the statistic for the i-th sample
    tmp <- statistic(samples[[i]])

    if ( is.finite(tmp) ) {
      count <- count + 1
      sampled_statistics[count] <- tmp
    }
  }

  p <- length(sampled_statistics[sampled_statistics < obs]) / length(sampled_statistics)

  return (list(samples=sampled_statistics,pvalue=p))
}
