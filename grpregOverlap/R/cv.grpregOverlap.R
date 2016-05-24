## function: cross-validation
# ------------------------------------------------------------------------------
cv.grpregOverlap <- function(X, y, group, ..., nfolds=10, seed, trace=FALSE) {
  fit <- grpregOverlap(X=X, y=y, group=group, returnX = TRUE, ...)
  cvfit <- cv.grpreg(X = fit$X.latent, y = y, group = fit$grp.vec, ...,
                     nfolds = nfolds, seed = seed, 
                     trace = trace)
  cvfit$fit <- fit
  val <- structure(cvfit, class = c('cv.grpregOverlap', 'cv.grpreg'))
  val
}
# ------------------------------------------------------------------------------
