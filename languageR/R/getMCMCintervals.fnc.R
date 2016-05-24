`getMCMCintervals.fnc` <-
function(fixf, mcmcMatrix, m) {
  stop("MCMC sampling is no longer supported by lme4")
  #nfixed = length(fixf)
  #nsamples = nrow(mcmcMatrix)
  #mat = matrix(0, nrow(m), nsamples)
  #for (i in 1:nsamples) {
  #  mat[,i] = m %*% as.numeric(mcmcMatrix[i,1:nfixed])
  #}
  #hpd = as.data.frame(HPDinterval(t(mat)))
  #return(hpd)
}

