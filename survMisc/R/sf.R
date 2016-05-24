#' @name sf
#' @title \bold{s}urvival (or hazard) \bold{f}unction
#' based on \eqn{e} and \eqn{n}.
#'
#' @include ten.R
#'
#' @param x One of the following:
#'  \describe{
#'   \item{default}{A numeric vector of events status (assumed sorted by time).}
#'   \item{numeric}{Vectors of events and numbers at risk (assumed sorted by time).}
#'   \item{}{A \code{ten} object.}
#'   \item{}{A \code{stratTen} object.}
#' }
#' @param ... Additional arguments (not implemented).
#' @param n Number at risk.
#' @param what See return, below.
#' @param SCV Include the \bold{S}quared \bold{C}oefficient of
#' \bold{V}ariation, which is calcluated using
#' the mean \eqn{\bar{x}}{mean(x)} and
#' the variance \eqn{\sigma_x^2}{var(x)}:
#'  \deqn{SCV_x = \frac{\sigma_x^2}{\bar{x}^2}}{
#'        SCV[x] = var(x) / mean(x)^2}
#' This measure of \emph{dispersion} is also referred to as
#' the 'standardized variance' or the 'noise'.
#' @param times Times for which to calculate the function.
#'  \cr
#' If \code{times=NULL} (the default), times are used for 
#' which at least one event occurred in at least one covariate group.
#' @param reCalc Recalcuate the values?
#'  \cr
#' If \code{reCalc=FALSE} (the default) and the \code{ten} object already has
#' the calculated values stored as an \code{attribute},
#' the value of the \code{attribute} is returned directly.
#'
#' @return
#' A {data.table} which is stored as an attribute of
#' the \code{ten} object. 
#' \cr
#' If \code{what="s"}, the \bold{s}urvival is returned, based on the
#' Kaplan-Meier or product-limit estimator.
#' This is \eqn{1} at \eqn{t=0} and thereafter is given by:
#' \deqn{\hat{S}(t) = \prod_{t \leq t_i} (1-\frac{e_i}{n_i} )}{
#'       S[t] = prod (1 - e[t]) / n[t] }
#'
#' If \code{what="sv"}, the \bold{s}urvival \bold{v}ariance is returned.
#' \cr
#' Greenwoods estimtor of the variance of the
#' Kaplan-Meier (product-limit) estimator is:
#' \deqn{Var[\hat{S}(t)] = [\hat{S}(t)]^2 \sum_{t_i \leq t}
#'                         \frac{e_i}{n_i (n_i - e_i)}}{
#'       Var(S[t]) = S[t]^2 sum e[t] / (n[t] * (n[t] - e[t]))}
#'
#' If \code{what="h"}, the \bold{h}azard is returned,
#' based on the the Nelson-Aalen estimator.
#' This has a value of \eqn{\hat{H}=0}{H=0} at \eqn{t=0}
#' and thereafter is given by:
#' \deqn{\hat{H}(t) = \sum_{t \leq t_i} \frac{e_i}{n_i}}{
#'       H[t] = sum(e[t] / n[t])}
#'
#' If \code{what="hv"}, the \bold{h}azard \bold{v}ariance is returned.
#' \cr
#' The variance of the Nelson-Aalen estimator is given by:
#' \deqn{Var[\hat{H}(t)] = \sum_{t_i \leq t} \frac{e_i}{n_i^2}}{
#'       Var(H[t]) = sum(e / n^2)}
#'
#' If \code{what="all"} (the default), \emph{all} of the above
#' are returned in a \code{data.table}, along with:
#' \cr
#' Survival, based on the Nelson-Aalen hazard estimator \eqn{H}, 
#' which is:
#'  \deqn{\hat{S_{na}}=e^{H}}{
#'        S[t] = exp(H[t])}
#' Hazard, based on the Kaplan-Meier survival estimator \eqn{S},
#' which is:
#'  \deqn{\hat{H_{km}} = -\log{S}}{
#'        H[t] = -log(S[t])}
#' 
#' @keywords survival 
#'
#' @rdname sf
#' @export
#' 
sf <- function(x, ...) UseMethod("sf")
#'
#' @rdname sf
#' @export
#' 
sf.default <- function(x, ...,
                       what=c("S", "H"),
                       SCV=FALSE,
                       times=NULL){
  stopifnot(all(x >= 0 && x <=1))
  what <- match.arg(what)
  t1 <- ten(x)
  return(sf.ten(t1, what=what, SCV=SCV, times=times))
}
#'
#' @rdname sf
#' @method sf ten
#' @aliases sf.ten
#' @export
#' 
#' @examples
#' data("kidney", package="KMsurv")
#' k1 <- ten(Surv(time=time, event=delta) ~ type, data=kidney)
#' sf(k1)
#' sf(k1, times=1:10, reCalc=TRUE)
#' k2 <- ten(with(kidney, Surv(time=time, event=delta)))
#' sf(k2)
#' ## K&M. Table 4.1A, pg 93.
#' ## 6MP patients
#' data("drug6mp", package="KMsurv")
#' d1 <- with(drug6mp, Surv(time=t2, event=relapse))
#' (d1 <- ten(d1))
#' sf(x=d1$e, n=d1$n, what="S")
#' 
sf.ten <- function(x, ...,
                   what=c("S", "H"),
                   SCV=FALSE,
                   times=NULL,
                   reCalc=FALSE){
  stopifnot(inherits(x, "ten"))
  if (!reCalc & !is.null(attr(x, "sf"))) return (attr(x, "sf"))
  what <- match.arg(what)
  ## name of variance
  nv1 <- paste0(what, "v")
  ## functions to use
  if (what=="S"){
    fun1 <- km
    fun1v <- kmv
  } else {
    fun1 <- na
    fun1v <- nav
  }
  if (attr(x, "ncg")==0){
    cg <- quote(NULL)
  } else {
    cg <- quote(cg)
  }
  if (is.null(times)){
    res1 <- x[, t, by=eval(cg)]
  } else {
    res1 <- data.table::data.table(
      data.table::rbindlist(
        lapply(times, function(i) x[t <= i, list("t1"=i, "t"=max(t)),
                                    by=eval(cg)])))
  }
  res1[, (what) := x[which(t %in% res1$t),
                  fun1(e, n), by=eval(cg)]$V1]
  res1[, (nv1) := x[which(t %in% res1$t),
                 fun1v(e, n), by=eval(cg)]$V1]
  if (SCV) res1[, "SCV" := Sv / S^2]
  data.table::setattr(x, "sf", res1)
  return(attr(x, "sf"))
}
#' 
#' @rdname sf
#' @method sf stratTen
#' @aliases strat.Ten
#' @export
#' @examples
#' data("pbc", package="survival")
#' t1 <- ten(Surv(time, status==2) ~ log(bili) + age + strata(edema), data=pbc)
#' sf(t1)
#' 
sf.stratTen <- function(x, ...,
                        what=c("S", "H"),
                        SCV=FALSE,
                        times=NULL,
                        reCalc=FALSE){
    lapply(x, function (i) sf(i, what=what, times=times, reCalc=reCalc))
    return(lapply(x, attr, "sf"))
}
#' 
#' @rdname sf
#' @method sf numeric
#' @aliases sf.numeric
#' @export
#' @examples
#' ## K&M. Table 4.2, pg 94.
#' data("bmt", package="KMsurv")
#' b1 <- bmt[bmt$group==1, ] # ALL patients
#' t2 <- ten(Surv(time=b1$t2, event=b1$d3))
#' with(t2, sf(x=e, n=n, what="Hv"))
#' ## K&M. Table 4.3, pg 97.
#' sf(x=t2$e, n=t2$n, what="all")
#' 
sf.numeric <- function(x, ...,
                       n=NULL,
                       what=c("all", "S", "Sv", "H", "Hv"),
                       SCV=FALSE,
                       times=NULL){
  what <- match.arg(what)
  if (is.null(n)) {
    what <- substring(what, 1L, 1L)
    if (what=="a") what <- "S"
    return(sf.default(x, what=what, SCV=SCV, time=times))
  }
  stopifnot(is.numeric(x) & is.numeric(n))
  stopifnot(length(x)==length(n))
  if(what == "all") {
    res1 <- data.table::data.table(
      "S"=km(x, n),
      "Sv"=kmv(x, n),
      "H"=na(x, n),
      "Hv"=nav(x, n))
    res1[, "sna" := exp(-H)]
    res1[, "hkm" := -log(S)]
    return(res1)
  }
  res1 <- switch(what,
                 "S"=km(x, n),
                 "Sv"=kmv(x, n),
                 "H"=na(x, n),
                 "Hv"=nav(x, n))
  return(res1)
}
### helper functions
km <- function(e, n) cumprod(1 - (e / n))
kmv <- function(e, n) cumprod(1 - (e / n))^2 * cumsum(e / (n * (n - e)))
na <- function(e, n) cumsum(e / n)
nav <- function(e, n) cumsum(e / (n^2))
## for R CMD check
H <- S <- Hx <- Sv <- NULL
