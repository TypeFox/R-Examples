#' @title Simple Bond Functions
#'
#' @description Basic yield, pricing, duration and convexity
#'     calculations. These functions perform simple present value
#'     calculations assuming that all periods between payments are the
#'     same length. Unlike bond functions in Excel, for example,
#'     settlement and maturity dates are not used. By default,
#'     duration is Macaulay duration.
#' 
#' @name bondsimple
#' @aliases bondpv bondyield duration convexity
#'
#' @return Return price, yield, or duration/convexity.
#'
#' @usage
#' bondpv(coupon, mat, yield, principal, freq)
#' bondyield(price, coupon, mat, principal, freq)
#' duration(price, coupon, mat, principal, freq, modified)
#' convexity(price, coupon, mat, principal, freq)
#' 
#' @param coupon annual coupon
#' @param mat maturity in years
#' @param yield annual yield to maturity. If freq > 1, the yield is
#'     freq times the per period yield.
#' @param price price of the bond
#' @param principal maturity payment of the bond, in addition to the
#'     final coupon. Default value is $1,000. If the instrument is an
#'     annuity, set principal to zero.
#' @param freq number of payments per year.
#' @param modified If true, compute modified duration, otherwise
#'     compute Macaulay duration. FALSE by default.
#'
#' @examples
#' coupon <- 6; mat <- 20; freq <- 2; principal <- 100; yield <- 0.045;
#'
#' price <- bondpv(coupon, mat, yield, principal, freq) # 119.7263
#' bondyield(coupon, mat, price=price, principal, freq) # 0.045
#' duration(price, coupon, mat, principal, freq, modified=FALSE) # 12.5043
#' duration(price, coupon, mat, principal, freq, modified=TRUE) # 12.3928
#' convexity(price, coupon, mat, principal, freq) # 205.3245
#'
#' @importFrom stats uniroot
#' 

#' @export
bondpv <- function(coupon, mat, yield, principal=1000, freq=1) {
    couponpd <- coupon/freq
    numpd <- mat*freq
    yieldpd <- yield/freq
    cf <- c(rep(couponpd, numpd-1), principal + couponpd)
    pv <- sum(cf*(1+yieldpd)^(-(1:numpd)))
    return(pv)
}

#' @export
bondyield <- function(price, coupon, mat, principal=1000, freq=1) {
    ## returns annual IRR
    pricediff <- function(x) bondpv(coupon, mat, x, principal, freq) - price
    irr <- uniroot(pricediff, c(-.25, 20), tol=1e-08)
    return(irr$root)
}

#' @export
duration <- function(price, coupon, mat, principal=1000, freq=1,
                     modified=FALSE) {
    ## macaulay duration by default
    yield <- bondyield(price, coupon, mat, principal, freq)
    yieldpd <- yield/freq
    couponpd <- coupon/freq
    numpd <- mat*freq
    cf <- c(rep(couponpd, numpd-1), principal+couponpd)
    dur <- sum((1:(numpd))*cf/(1+yieldpd)^(1:(numpd)))/
        price/freq/(1+yieldpd*modified)
    return(dur)
}

#' @export
convexity <- function(price, coupon, mat, principal=1000, freq=1) {
    ## mat in years, coupon in years, freq=times/year
    numpd <- mat*freq
    couponpd <- coupon/freq
    yield <- bondyield(price, coupon, mat, principal, freq)
    yieldpd <- yield/freq
    cf <- c(rep(couponpd, (numpd)-1), principal+couponpd)
    pv <- sum((1:numpd + (1:numpd)^2)*cf/(1+yieldpd)^(1:numpd))
    conv <- pv/price/(1+yieldpd)^2/freq^2
    return(conv)
}
