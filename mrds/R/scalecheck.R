## Code taken from optimx-package.R in package optimx
## Code by John C Nash [aut, cre], Ravi Varadhan [aut], Gabor Grothendieck [ctb]
## See http://cran.r-project.org/web/packages/optimx/index.html
##################################################################
scalecheck<-function(par, lower=lower, upper=upper,dowarn){
  # a function to check the initial parameters and bounds for inputs to optimization codes
  # Arguments:
  # par -- starting parameters supplied 
  #   lower, upper -- lower and upper bounds supplied
  #
  # Returns:
  # list(lpratio, lbratio) -- the log of the ratio of largest to smallest parameters
  #   and bounds intervals (upper-lower) in absolute value (ignoring Inf, NULL, NA)
  ######################################
  if (is.null(par)) { stop("Null parameter vector") }
  npar<-length(par)
  if (is.null(lower)) {
    if (dowarn) warning("Null lower bounds vector")
    lower<-rep(-Inf,npar)
  }
  if (is.null(upper)) {
    if (dowarn) warning("Null upper bounds vector")
    upper<-rep(Inf,npar)
  }
  newpar<-abs(par[which(is.finite(par))])
  logpar<-log10(newpar[which(newpar>0)]) # Change 20100711
  newlower<-abs(lower[which(is.finite(lower))])
  loglower<-log10(newlower[which(newlower>0)]) # Change 20100711
  newupper<-abs(upper[which(is.finite(upper))])
  logupper<-log10(newupper[which(newupper>0)]) # Change 20100711
  bddiff<-upper-lower
  bddiff<-bddiff[which(is.finite(bddiff))]
  lbd<-log10(bddiff[which(bddiff>0)]) # Change 20100711
  lpratio<-max(logpar) - min(logpar)
  if (length(lbd) > 0) {
    lbratio<-max(lbd)-min(lbd)
  } else {
    lbratio<-NA
  }
  ratios<-list(lpratio=lpratio,lbratio=lbratio)
  # return(ratios)
}
# -------------- end scalecheck ----------------- #
#################################################################
