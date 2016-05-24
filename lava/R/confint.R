##' Calculate Wald og Likelihood based (profile likelihood) confidence intervals
##'
##' Calculates either Wald confidence limits: \deqn{\hat{\theta} \pm
##' z_{\alpha/2}*\hat\sigma_{\hat\theta}} or profile likelihood confidence
##' limits, defined as the set of value \eqn{\tau}:
##' \deqn{logLik(\hat\theta_{\tau},\tau)-logLik(\hat\theta)< q_{\alpha}/2}
##'
##' where \eqn{q_{\alpha}} is the \eqn{\alpha} fractile of the \eqn{\chi^2_1}
##' distribution, and \eqn{\hat\theta_{\tau}} are obtained by maximizing the
##' log-likelihood with tau being fixed.
##'
##' @title Calculate confidence limits for parameters
##' @param object \code{lvm}-object.
##' @param parm Index of which parameters to calculate confidence limits for.
##' @param level Confidence level
##' @param profile Logical expression defining whether to calculate confidence
##' limits via the profile log likelihood
##' @param curve if FALSE and profile is TRUE, confidence limits are
##' returned. Otherwise, the profile curve is returned.
##' @param n Number of points to evaluate profile log-likelihood in
##' over the interval defined by \code{interval}
##' @param interval Interval over which the profiling is done
##' @param lower If FALSE the lower limit will not be estimated (profile intervals only)
##' @param upper If FALSE the upper limit will not be estimated (profile intervals only)
##' @param \dots Additional arguments to be passed to the low level functions
##' @return A 2xp matrix with columns of lower and upper confidence limits
##' @author Klaus K. Holst
##' @seealso \code{\link{bootstrap}{lvm}}
##' @keywords models regression
##' @examples
##'
##' m <- lvm(y~x)
##' d <- sim(m,100)
##' e <- estimate(y~x, d)
##' confint(e,3,profile=TRUE)
##' confint(e,3)
##' \donttest{ ## Reduce Ex.timings
##' B <- bootstrap(e,R=50)
##' B
##' }
##' @aliases confint.multigroupfit
##' @export
##' @method confint lvmfit
confint.lvmfit <- function(object,parm=seq_len(length(coef(object))),level=0.95,profile=FALSE,curve=FALSE,n=20,interval=NULL,lower=TRUE,upper=TRUE,...) {
  if (is.character(parm)) {
    parm <- parpos(Model(object),p=parm)
    parm <- parm[attributes(parm)$ord]
  }
  if (!profile) {
    return(confint.default(object,parm=parm,level=level,...))
  }
  res <- c()
  for (i in parm) {
    res <- rbind(res, profci.lvmfit(object,parm=i,level=level,profile=profile,n=n,curve=curve,interval=interval,lower=lower,upper=upper,...))
    if (curve) return(res)
  }
  rownames(res) <- names(coef(object))[parm]
  colnames(res) <- paste((c(0,1) + c(1,-1)*(1-level)/2)*100,"%")
  return(res)
}


##' @export
confint.multigroupfit <- function(object,parm=seq_along(pars(object)),level=0.95,
                                  estimates=TRUE,...) {
  p <- 1-(1-level)/2
  res <- cbind(pars(object),pars(object)) + qnorm(p)*cbind(-1,1)%x%diag(vcov(object))^0.5
  colnames(res) <- paste0(c(1-p,p)*100,"%")
  rownames(res) <- parpos(object); rownames(res)[is.na(rownames(res))] <- ""
  if (estimates) res <- cbind(coef(object,level=0)[,c(1,2,4)],res)
  res[parm,,drop=FALSE]
}
