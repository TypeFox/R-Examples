parcorplot <- function(mcmcout, parms=NULL, regex=NULL, random=NULL, col=gray(11:0/11), ...){
  if (!(is.mcmc(mcmcout) | is.mcmc.list(mcmcout)))
      mcmcout <- as.mcmc(mcmcout)
  if (!is.mcmc.list(mcmcout))
      mcmcout <- mcmc.list(mcmcout)
  parnames <- parms2plot(varnames(mcmcout), parms, regex, random)
  mat <- cor(as.matrix(mcmcout[, parnames, drop=FALSE]))
  corplot(mat, col=col, ...)
  }
