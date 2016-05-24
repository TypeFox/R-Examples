simulHeterogeneous <- function(n=100, nbCluster = 3, nbVariable=1, models="gaussian_pk_sjk", pmiss =0.1, z = NULL)
{
  if (length(nbVariable) != length(models))
  { stop("nbVariable and models must have the same size")}
  if (is.null((z)))
  { z <- rmultinom(n, 2, prob=rep(1/nbCluster, length.out=nbCluster))+1;}
  else
  {  n <- length(z);}

  res <- list(z);
  for (i in 1:length(nbVariable))
  {
    model <- matrix(nrow = n, ncol = nbVariable[i]);
    if (clusterValidCategoricalNames(models[i]))
    {
      # simulate Categorical

    }
    else if (clusterValidDiagGaussianNames(models[i]))
    {
      # simulate Gaussian
    }
    else if (clusterValidGammaNames(models[i]))
    {
      # simulate Gamma
    }
    else if (clusterValidPoissonNames(models[i]))
    {
      # simulate Poisson
    }
    else
    {
      stop("Invalid model name")
    }
  }
  res
}
