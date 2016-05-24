# lme4 : fixef, VarCorr

overdisp.glmer <- function(model) {
  ## From http://glmm.wikidot.com/faq
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(lme4::VarCorr(model),vpars))+length(lme4::fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model)
  dev <- sum(rp^2)
  prat <- dev/rdf
  cat(paste("Residual deviance: ",round(dev,3)," on ",rdf," degrees of freedom",
    " (ratio: ",round(prat,3),")\n",sep=""))
}
