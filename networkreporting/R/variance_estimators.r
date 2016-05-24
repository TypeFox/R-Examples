#####################################################
## variance_estimators.R
##
## this file has the functions that produce
## variance estimates

#####################################################
##' killworth.se
##'
##' compute standard errors for scale-up estimates
##' based on the Killworth estimator
##'
##' note that this is provided for comparison, but
##' that we do not generally recommend using this
##' strategy for estimating variance
##'
##' @param estimates TODO
##' @param d.hat TODO
##' @param total.popn.size TODO
##' @param total TODO
##' @param missing TODO
##' @return the estimated standard error
##' @keywords internal
killworth.se <- function(estimates,
                         d.hat,
                         total.popn.size=NULL,
                         total=TRUE,
                         missing="ignore") {

  stop("killworth.se is not yet implemented.")

  ## TODO -- code below this point not yet altered...

  if (total & is.null(total.popn.size)) {
    stop("you must pass in total.popn.size to get the Killworth variance estimates of the total (rather than proportions).")
  }

  na.rm <- ifelse(missing == "ignore", TRUE, FALSE)

  sum.d.hat <- sum(d.hat, na.rm=na.rm)

  est.props <- ifelse(rep(total, length(estimates)),
                      estimates / total.popn.size,
                      estimates)

  res <- plyr::aaply(est.props,
               1,
               function(est) {
                 return(sqrt((est*(1-est))/sum.d.hat))
               })

  names(res) <- names(estimates)

  if (total) {
    res <- res*total.popn.size
  }

  return(res)
}


