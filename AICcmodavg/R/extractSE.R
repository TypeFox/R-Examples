##create SE extractor function that includes estimate labels
##generic
extractSE <- function(mod, ...) {
  UseMethod("extractSE", mod)
}


extractSE.default <- function(mod, ...) {
  stop("\nFunction not yet defined for this object class\n")
}


##methods
##coxme objects
extractSE.coxme <- function(mod, ...){
  ##extract vcov matrix
  vcov.mat <- as.matrix(vcov(mod))
  se <- sqrt(diag(vcov.mat))
  fixed.labels <- names(fixef(mod))
  names(se) <- fixed.labels
  return(se)
}


##lmekin objects
extractSE.lmekin <- function(mod, ...){
  ##extract vcov matrix
  vcov.mat <- as.matrix(vcov(mod))
  se <- sqrt(diag(vcov.mat))
  fixed.labels <- names(fixef(mod))
  names(se) <- fixed.labels
  return(se)
}


##mer objects
extractSE.mer <- function(mod, ...){
  ##extract vcov matrix
  vcov.mat <- as.matrix(vcov(mod))
  se <- sqrt(diag(vcov.mat))
  fixed.labels <- names(fixef(mod))
  names(se) <- fixed.labels
  return(se)
}


##merMod objects
extractSE.merMod <- function(mod, ...){
  ##extract vcov matrix
  vcov.mat <- as.matrix(vcov(mod))
  se <- sqrt(diag(vcov.mat))
  fixed.labels <- names(fixef(mod))
  names(se) <- fixed.labels
  return(se)
}
