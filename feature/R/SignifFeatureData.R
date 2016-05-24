########## R function: SignifFeatureData ##########

# For determining data points that lie in the
# region of significant
# feature for a particular bandwidth and
# significance level.
#
# Last changed: 08 AUG 2011

SignifFeatureData <- function(x, d, dest, SignifFeature)
{
  n <- nrow(x)
  x.ind <- matrix(0, ncol=d, nrow=n)
 
  for (j in 1:d) x.ind[,j] <- findInterval(x[,j], dest$x.grid[[j]])
  return(SignifFeature[x.ind])
}

########## End of SignifFeatureData ##########
