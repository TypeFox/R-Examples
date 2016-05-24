#' @name ci
#' @title \bold{c}onfidence \bold{i}ntervals for survival curves.
#' 
#' @include ten.R
#' @include sf.R
#'
#' @param x An object of class \code{ten}.
#' @param CI Confidence intervals. As the function currently relies on lookup
#' tables, currently only 90\%, 95\% (the default) and 99\% are supported.
#' @param how Method to use for confidence interval.
#'  \cr
#' \code{point} (the default) uses pointwise confirence intervals.
#'  \cr
#' The alternatives use confidence \emph{bands} (see details).
#' @param trans Transformation to use.
#' \cr
#' The default is \code{trans="log"}.
#' \cr
#' Also supported are linear and arcsine-square root transformations.
#' @param tL \bold{L}ower time point. Used in construction of confidence bands.
#' @param tU \bold{U}pper time point. Used in construction of confidence bands.
#' @inheritParams sf.ten
#' 
#' @return The \code{ten} object is modified in place by the additional of a 
#' \code{data.table} as an \code{attribute}.
#'  \cr
#' \code{attr(x, "ci")} is printed.
#'  \cr
#' This A \code{survfit} object. The \code{upper} and \code{lower}
#' elements in the list (representing confidence intervals)
#' are modified from the original.
#' \cr
#' Other elements will also be shortened if the time range under consideration has been
#' reduced from the original.
#'
#' @details
#' In the equations below
#'  \deqn{\sigma^2_s(t) = \frac{\hat{V}[\hat{S}(t)]}{\hat{S}^2(t)} }{
#'         sigma^2(t) = V[S(t)]/[S(t)]^2}
#' Where \eqn{\hat{S}(t) }{S(t)} is the Kaplan-Meier survival estimate and
#' \eqn{\hat{V}[\hat{S}(t)]}{V[S(t)]} is Greenwood's estimate of its
#' variance.
#'  \cr
#' The \bold{pointwise} confidence intervals are valid for \emph{individual}
#' times, e.g. \code{median} and \code{\link{quantile}} values.
#' When plotted and joined for multiple points they tend to
#' be narrower than the \emph{bands} described below.
#' Thus they tend to exaggerate the impression of certainty
#' when used to plot confidence intervals for a time range.
#' They should not be interpreted as giving the intervals
#' within which the \emph{entire} survival function lies.
#'  \cr
#' For a given significance level \eqn{\alpha}{alpha},
#' they are calculated using the standard normal distribution \eqn{Z}
#' as follows:
#'
#' \itemize{
#'
#' \item linear
#'       \deqn{\hat{S}(t) \pm Z_{1- \alpha} \sigma (t) \hat{S}(t)}{
#'             S(t)+- Z(1-alpha) sigma(t) S(t)}
#'
#' \item log transform
#'       \deqn{ [ \hat{S}(t)^{\frac{1}{\theta}}, \hat{S}(t)^{\theta} ] }{
#'              [S(t)^(1/theta), S(t)^theta]}
#' where
#' \deqn{ \theta = \exp{ \frac{Z_{1- \alpha} \sigma (t)}{ \log{\hat{S}(t)}}} }{
#'        theta = exp ( Z(1-alpha)sigma(t) / log(S(t)) )}
#'
#' \item arcsine-square root transform
#' \cr
#' upper:
#' \cr
#' \deqn{ \sin^2(\max[0, \arcsin{\sqrt{\hat{S}(t)}} -
#'   \frac{Z_{1- \alpha}\sigma(t)}{2}
#'   \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
#'       sin^2(max[0, arcsin S(t)^0.5 - Z(1-alpha)sigma(t)/2 (S(t)/1-S(t))^0.5])}
#' lower:
#' \deqn{ \sin^2(\min[\frac{\pi}{2}, \arcsin{\sqrt{\hat{S}(t)}} +
#'   \frac{Z_{1- \alpha}\sigma(t)}{2}
#'   \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
#'       sin^2(min[pi/2, arcsin S(t)^0.5  + Z(1-alpha)sigma(t)/2 (S(t)/1-S(t))^0.5])}
#'
#' }
#'
#' Confidence \bold{bands} give the values within which the survival function
#' falls within a \emph{range} of timepoints.
#'  \cr \cr
#' The time range under consideration is given so that
#' \eqn{t_l \geq t_{min}}{tL >= min(t)}, the minimum or lowest event time and
#' \eqn{t_u \leq t_{max}}{tU <= max(t)}, the maximum or largest event time.
#'  \cr
#' For a sample size \eqn{n} and \eqn{0 < a_l < a_u <1}:
#' \deqn{a_l = \frac{n\sigma^2_s(t_l)}{1+n\sigma^2_s(t_l)}}{
#' a_l = n*sigma^2(t_l) / [1+n*sigma^2(t_l)]}
#' \deqn{a_u = \frac{n\sigma^2_s(t_u)}{1+n\sigma^2_s(t_u)}}{
#' a_u = n*sigma^2(t_u) / [1+n*sigma^2(t_u)]}
#'
#' For the \bold{Nair} or \bold{equal precision} (\bold{EP}) confidence bands,
#' we begin by obtaining the relevant
#' confidence coefficient \eqn{c_{\alpha}}{c[alpha]}. This is obtained from
#' the upper \eqn{\alpha}{a}-th fractile of the random variable
#' \deqn{U = \sup{|W^o(x)\sqrt{[x(1-x)]}|, \quad a_l \leq x \leq a_u}}{
#' U = sup{ |W(x)[x(1-x)]^0.5|, a_l <= x <= a_u} }
#' Where \eqn{W^o}{W} is a standard Brownian bridge.
#'  \cr
#' The intervals are:
#'
#' \itemize{
#'  \item linear
#'   \deqn{\hat{S}(t) \pm c_{\alpha} \sigma_s(t) \hat{S}(t)}{
#'         S(t)+- c[alpha] sigma(t) S(t)}
#'
#'  \item log transform (the default)
#'    \cr
#'   This uses \eqn{\theta}{theta} as below:
#'    \deqn{\theta = \exp{ \frac{c_{\alpha} \sigma_s(t)}{ \log{\hat{S}(t)}}}}{
#'          theta = exp (c[alpha] * sigma(t) / log(S(t)))}
#'   And is given by:
#'    \deqn{[\hat{S}(t)^{\frac{1}{\theta}}, \hat{S}(t)^{\theta}]}{
#'          [S(t)^(1/theta), S(t)^theta]}
#'
#'  \item arcsine-square root transform
#'    \cr
#'   upper:
#'   \deqn{\sin^2(\max[0, \arcsin{\sqrt{\hat{S}(t)}}
#'         - \frac{c_{\alpha}\sigma_s(t)}{2}
#'         \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}])}{
#'         sin^2(max[0, arcsin S(t)^0.5 - c[alpha]*sigma(t)/2 (S(t)/1-S(t))^0.5])}
#'   lower:
#'   \deqn{\sin^2(\min[\frac{\pi}{2}, \arcsin{\sqrt{\hat{S}(t)}}
#'         + \frac{c_{\alpha}\sigma_s(t)}{2}
#'         \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
#'         sin^2(min[pi/2, arcsin S(t)^0.5 - c[alpha]*sigma(t)/2 (S(t)/1-S(t))^0.5])}
#'
#' }
#'
#' For the \bold{Hall-Wellner} bands the confidence coefficient
#' \eqn{k_{\alpha}}{k[alpha]}
#' is obtained from the upper \eqn{\alpha}{a}-th fractile of a
#' Brownian bridge.
#'  \cr
#' In this case \eqn{t_l} can be \eqn{=0}.
#'  \cr
#' The intervals are:
#'
#' \itemize{
#'
#' \item linear
#'  \deqn{\hat{S}(t) \pm
#'        k_{\alpha} \frac{1+n\sigma^2_s(t)}{\sqrt{n}} \hat{S}(t)}{
#'        S(t)+- k[alpha] [1+n*sigma^2(t)]*S(t) / n^0.5 }
#'
#' \item log transform
#'  \deqn{[\hat{S}(t)^{\frac{1}{\theta}}, \hat{S}(t)^{\theta}]}{
#'        [S(t)^(1/theta), S(t)^theta]}
#'  where
#'  \deqn{\theta = \exp{ \frac{k_{\alpha}[1+n\sigma^2_s(t)]}{
#'        \sqrt{n}\log{\hat{S}(t)}}} }{
#'        theta = exp(k[alpha] * [1 + n * sigma^2(t)] / n^0.5 * log(S(t)))}
#'
#' \item arcsine-square root transform
#' \cr
#' upper:
#' \deqn{ \sin^2(\max[0, \arcsin{\sqrt{\hat{S}(t)}}
#'   - \frac{k_{\alpha}[1+n\sigma_s(t)]}{2\sqrt{n}}
#'   \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
#' sin^2( max[0, arcsin S(t)^0.5 - k[alpha]*[1+n*sigma^2(t)]/(2*n^0.5) (S(t)/1-S(t))^0.5])}
#' lower:
#' \deqn{ \sin^2(\min[\frac{\pi}{2}, \arcsin{\sqrt{\hat{S}(t)}}
#'   + \frac{k_{\alpha}[1+n\sigma^2_s(t)]}{2\sqrt{n}}
#'   \sqrt{ \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{
#' sin^2( min[pi/2, arcsin S(t)^0.5 - k[alpha]*[1+n*sigma^2(t)]/(2*n^0.5) (S(t)/1-S(t))^0.5])}
#'
#' }
#'
#' @source The function is loosely based on \code{km.ci::km.ci}.
#' 
#' @note
#' \itemize{
#'  \item For the Nair and Hall-Wellner bands, the function currently
#'        relies on the lookup tables in \code{package:km.ci}.
#'  \item Generally, the arcsin-square root transform has the best coverage properties.
#'  \item All bands have good coverage properties for samples as small as \eqn{n=20},
#'        except for the \bold{Nair} / \bold{EP} bands with a linear transformation,
#'        which perform poorly when \eqn{n < 200}.
#' }
#' 
#' @keywords survival
#' 
#' @references
#' Nair V, 1984.
#' Confidence bands for survival functions with censored data: a comparative study.
#' \emph{Technometrics}. \bold{26}(3):265-75.
#' \href{http://www.jstor.org/stable/1267553}{JSTOR}.
#' @references
#' Hall WJ, Wellner JA, 1980.
#' Confidence bands for a survival curve from censored data.
#' \emph{Biometrika}. \bold{67}(1):133-43.
#' \href{http://www.jstor.org/stable/2335326}{JSTOR}.
#' 
#' @seealso \code{\link{sf}}
#' @seealso \code{\link{quantile}}
#'
#' @rdname ci
#' @export
#' 
ci <- function(x, ...) UseMethod("ci")
#' 
#' @rdname ci
#' @method ci ten
#' @aliases ci.ten
#' @export
#' @examples
#' ## K&M 2nd ed. Section 4.3. Example 4.2, pg 105.
#' data("bmt", package="KMsurv")
#' b1 <- bmt[bmt$group==1, ] # ALL patients
#' ## K&M 2nd ed. Section 4.4. Example 4.2 (cont.), pg 111.
#' ## patients with ALL
#' t1 <- ten(Surv(t2, d3) ~ 1, data=bmt[bmt$group==1, ])
#' ci(t1, how="nair", trans="lin", tL=100, tU=600)
#' ## Table 4.5, pg. 111.
#' lapply(list("lin", "log", "asi"),
#'        function(x) ci(t1, how="nair", trans=x, tL=100, tU=600))
#' ## Table 4.6, pg. 111.
#' lapply(list("lin", "log", "asi"),
#'        function(x) ci(t1, how="hall", trans=x, tL=100, tU=600))
#' t1 <- ten(Surv(t2, d3) ~ group, data=bmt)
#' ci(t1, CI="0.95", how="nair", trans="lin", tL=100, tU=600)
#'
ci.ten <- function(x,
                   ...,
                   CI=c("0.95", "0.9", "0.99"),
                   how=c("point", "nair", "hall"),
                   trans=c("log", "lin", "asi"),
                   tL=NULL,
                   tU=NULL,
                   reCalc=FALSE){
    stopifnot(inherits(x, "ten"))
    if (!reCalc & !is.null(attr(x, "ci"))) return (attr(x, "ci"))
### 
    trans <- match.arg(trans)
    CI <- 100 * as.numeric(match.arg(CI))
    how <- match.arg(how)
###
    sf(x, SCV=TRUE)
    if (is.null(tL)) tL <- min(x[, min(t), by=cg][, V1])
    if (is.null(tU)) tU <- max(x[, max(t), by=cg][, V1])
    s1 <- data.table::copy(attr(x, "sf")[t >= tL & t <= tU, ])
    if (!"cg" %in% names(s1)) s1[, "cg" := 1]
    if (attr(x, "ncg") >= 1) {
        n1 <- x[, max(ncg), by=cg]
    } else {
        n1 <- data.table::data.table(
            cg=1, V1=x[, max(n)])
    }
### get reference value
    if (how=="point") {
        ## get Z (quantile from normal distribution)
        alpha <- (100 - CI) / 100
        ref1 <- stats::qnorm(1 - alpha / 2)
    }
    if (how=="nair" | how=="hall") {
        A1 <- mapply(FUN=function(cg2, n2)
            data.table::rbindlist(
                list(
                    utils::head(s1[cg==cg2, ], 1),
                    utils::tail(s1[cg==cg2, ], 1)))[, genA(SCV=SCV, n=n2)],
            n1[, cg],
            n1[, V1])
        ## get lookup table for confidence coefficient
        d1s <- paste0("critical.value.", how, ".", CI)
        do.call(utils::data, list(eval(substitute(d1s)), package="km.ci"))
        d1 <- NULL
        do.call(assign, list("d1", eval(parse(text=d1s))))
        ## label lookup table
        if (how=="nair") {
            rownames(d1) <- seq(0.1, 0.98, by=0.02)
            colnames(d1) <- seq(0.02, 0.6, by=0.02)
        } else {
            rownames(d1) <- seq(0.1, 1.0, by=0.02)
            colnames(d1) <- seq(0, 0.6, by=0.02)
        }
        ## if a_L and/or a_U are outside the range of lookup table.
        err1 <- "
Confidence coefficients are outside the range available in the lookup table.
Suggest try narrower range for time i.e.
 increase lower time (tL) and/or decrease upper time (tU).
"
        ref1 <- mapply(function(i, j) tryCatch(d1[i, j], error=function(e) NaN),
                       i=as.character(A1[2, ]),
                       j=as.character(A1[1, ]))
        if (any(is.nan(ref1))) message(err1)
    }
    ## transform function
    tf1 <- switch(trans,
                  lin=LIN,
                  log=LOG,
                  asi=ASI)
    s1[, "lower" := with(
             list(bound=-1,
                  how=how),
             unlist(
                 mapply(   
                     FUN=function(cg2, ref2, n2)
                         s1[cg==cg2, 
                            tf1(how=how, S=S, SCV=SCV, ref=ref2, n=n2, bound=bound)],
                     ref2=ref1,
                     cg2=n1[, cg],
                     n2=n1[, V1])))]
    data.table::set(s1, i=s1[, which(lower < 0)], j="lower", value=0)
    s1[, "upper" := with(
             list(bound=1,
                  how=how),
             unlist(
                 mapply(   
                     FUN=function(cg2, ref2, n2)
                         s1[cg==cg2, 
                            tf1(how=how, S=S, SCV=SCV, ref=ref2, n=n2, bound=bound)],
                     cg2=n1[, cg],
                     ref2=ref1,
                     n2=n1[, V1])))]
    data.table::set(s1, i=s1[, which(upper > 1)], j="upper", value=1)
    data.table::setattr(s1, "CI", CI)
    data.table::setattr(s1, "how", how)
    data.table::setattr(s1, "trans", trans)
    data.table::setattr(s1, "class", c("tenAttr", class(s1)))
    data.table::setattr(x, "ci", s1)
    return(attr(x, "ci"))
}
#' 
#' @rdname ci
#' @method ci stratTen
#' @aliases ci.stratTen
#' @export
#' @examples
#' ## stratified model
#' data("pbc", package="survival")
#' t1 <- ten(coxph(Surv(time, status==2) ~ log(bili) + age + strata(edema), data=pbc))
#' ci(t1)
#' 
ci.stratTen <- function(x,
                        ...,
                        CI=c("0.95", "0.9", "0.99"),
                        how=c("point", "nair", "hall"),
                        trans=c("log", "lin", "asi"),
                        tL=NULL,
                        tU=NULL) {
    lapply(x, ci, CI=CI, how=how, trans=trans, tL=tL, tU=tU)
    lapply(x, attr, "ci")
}
### Helper functions
## generate 'a', lower and upper
genA <- function(SCV, n) round((SCV * n) / (1 + SCV * n), 1)
## linear transform
LIN <- function(how, S, SCV, ref=NULL, n=NULL, bound=NULL){
    if (how=="hall") {
        S + sign(bound) * S * (ref * (1 + SCV * n)) / sqrt(n)
    } else {
        S + sign(bound) * ref * S * sqrt(SCV)
    }
}
## log transform
LOG <- function(how, S, SCV, ref=NULL, n=NULL, bound=NULL){
    if (how=="hall"){
        theta <-  exp(ref * (1 + SCV * n)  / (log(S) * sqrt(n)))
    } else {
        theta <- exp(ref * sqrt(SCV) / log(S))
    }
    S^(theta^sign(bound))
}
## acrsin-sqrt transform
ASI <- function(how, S, SCV, ref=NULL, n=NULL, bound=NULL){
    if (how=="hall") {
        sin(
            asin(sqrt(S)) + sign(bound) * (ref / 2) * (1 + SCV * n) / sqrt(n) *
            sqrt(S / (1 - S)) 
        )^2
    } else {
        sin(
            asin(sqrt(S)) + sign(bound) * (ref / 2) * sqrt(SCV) * 
            sqrt(S / (1 - S))
        )^2
    }
}
## R CMD check
V1 <- SCV <- NULL
