##' Bond pricing using yield to maturity.
##' 
##' bond.price computes the price given the yield to maturity
##' bond.duration computes the duration given the yield to maturity
##' bond.yield computes the yield to maturity given the price
##' bond.prices, bond.durations and bond.yields are wrapper functions
##' that use mapply to vectorize bond.price, bond.duration and bond.yield
##' All arguments to bond.prices, bond.durations and bond.yields
##' can be vectors.
##' On the other hand, bond.price, bond.duration and bond.yield do not allow vectors 
##' Standard compounding and day count conventions are supported for all functions.  
##' 
##'
##' @rdname bonds
##' @name bonds
##' @aliases bond.price bond.TCF
##' bond.yield bond.duration bonds
##' bond.prices bond.durations bond.yields
##' @param settle The settlement date for which the bond is traded. Can be a
##' character string or any object that can be converted into date using
##' \code{\link{as.Date}}.
##' @param mature The maturity date of the bond. Can be a character string or
##' any object that can be converted into date using \code{\link{as.Date}}
##' @param coupon The coupon rate in decimal (0.10 or 10e-2 for 10\%)
##' @param freq The frequency of coupon payments: 1 for annual, 2 for
##' semi-annual, 12 for monthly.
##' @param comp.freq The frequency of compounding of the bond yield: 1 for
##' annual, 2 for semi-annual, 12 for monthly. Usually same as freq.
##' @param price The clean price of the bond.
##' @param yield The yield to maturity of the bond
##' @param convention The daycount convention
##' @param modified A logical value used in duration. \code{TRUE} to return Modified
##' Duration, \code{FALSE} otherwise
##' @return \code{bond.TCF} returns a list of three components
##' \item{t}{A vector of cash flow dates in number of years}
##' \item{cf}{A vector of cash flows}
##' \item{accrued}{The accrued interest}
##' @author Prof. Jayanth R. Varma \email{jrvarma@@iimahd.ernet.in}

NULL


##' @rdname bonds
##' @export
bond.price <- function(settle, mature, coupon, freq=2, yield, 
                  convention=c("30/360", "ACT/ACT", "ACT/360", "30/360E"),
                  comp.freq=freq){
  settle <- as.Date(settle)
  mature <- as.Date(mature)  
  TCF <- bond.TCF(settle, mature, coupon, freq, convention)
  if (length(TCF$t) == 1){
      ## short government bond convention = simple interest
    return (TCF$cf[1] / (1 + yield * TCF$t[1]) - TCF$accrued)
  }
  if (comp.freq != freq){
    yield <- equiv.rate(yield, from.freq=comp.freq, to.freq=freq)
  }
  df <- (1+yield/freq)^(-freq*TCF$t)
  sum(TCF$cf*df) - TCF$accrued
}

##' @rdname bonds
##' @export
bond.yield <- function(settle, mature, coupon, freq=2, price,
                  convention=c("30/360", "ACT/ACT", "ACT/360",
                               "30/360E"),
                  comp.freq=freq){
  settle <- as.Date(settle)
  mature <- as.Date(mature)
  TCF <- bond.TCF(settle, mature, coupon, freq, convention)
  if (length(TCF$t) == 1){
      ## short government bond convention = simple interest
      return ( (TCF$cf[1] / (price + TCF$accrued) - 1) / TCF$t[1] )
  }
  fn <- function(r){
    df <- (1+r/freq)^(-freq*TCF$t)
    dirty <- crossprod(TCF$cf, df)
    list(value = dirty - TCF$accrued - price,
         gradient = -as.numeric(crossprod(TCF$cf*df, TCF$t))/
             (1+r/freq))
  }
  equiv.rate(irr.solve(f=fn), from.freq=freq, to.freq=comp.freq)
}

##' @rdname bonds
##' @export
bond.duration <- function(settle, mature, coupon, freq=2, yield,
                     convention=c("30/360", "ACT/ACT", "ACT/360",
                                  "30/360E"),
                     modified=FALSE, comp.freq=freq){
  settle <- as.Date(settle)
  mature <- as.Date(mature)
  TCF <- bond.TCF(settle, mature, coupon, freq, convention)
  if (length(TCF$t) == 1){
      ## short government bond convention = simple interest
      return(yearFraction(settle, mature))
  }
  duration(cf=TCF$cf, cf.t=TCF$t, rate=yield, cf.freq=freq,
           comp.freq=comp.freq, immediate.start=FALSE,
           modified=modified)
}

##' @rdname bonds
##' @export
bond.TCF <- function(settle, mature, coupon, freq=2,
                 convention=c("30/360", "ACT/ACT", "ACT/360", "30/360E")){
  settle <- as.Date(settle)
  mature <- as.Date(mature)
  nextC <- coupons.next(settle, mature, freq)
  prevC <- coupons.prev(settle, mature, freq)
  accrued <- 100 * coupon * yearFraction(prevC, settle, prevC,
                                         nextC, freq, convention)
  start <- yearFraction(settle, nextC, prevC, nextC, freq, convention)
  t <- seq(from=start, by=1/freq, length.out=coupons.n(settle, mature,
                                                       freq))
  cf <- rep(coupon*100/freq, length(t))
  cf[length(cf)] <- 100 + coupon*100/freq
  list(t=t,cf=cf, accrued=accrued)
}

##' Bond pricing using yield to maturity.
##'
##' 
##' Convenience functions for finding coupon dates and number of coupons of a bond.
##' 
##'
##' @name coupons
##' @aliases coupons.n coupons.dates coupons.next coupons.prev coupons
##' @author Prof. Jayanth R. Varma \email{jrvarma@@iimahd.ernet.in}
NULL


##' @rdname coupons
##' @inheritParams bond.price
##' @export
coupons.dates <- function(settle, mature, freq=2){
  settle <- as.Date(settle)
  mature <- as.Date(mature)
  m <- -12/freq
  n <- coupons.n(settle, mature, freq)
  edate(mature, m * (n:1 - 1))
}

##' @rdname coupons
##' @inheritParams bond.price
##' @export
coupons.n <- function(settle, mature, freq=2){
  mature <- as.Date(mature)
  settle <- as.Date(settle)
  n <- as.integer(freq * (mature - settle) / 365.25)
  m <- -12/freq
  while(edate(mature, n * m) <= settle) n <- n - 1
  while(edate(mature, (n + 1) * m) > settle) n <- n + 1
  n+1
}

##' @rdname coupons
##' @inheritParams bond.price
##' @export
coupons.next <- function(settle, mature, freq=2){
  settle <- as.Date(settle)
  mature <- as.Date(mature)
  m <- -12/freq
  n <- coupons.n(settle, mature, freq)
  edate(mature, m * (n-1))
}

##' @rdname coupons
##' @inheritParams bond.price
##' @export
coupons.prev <- function(settle, mature, freq=2){
  settle <- as.Date(settle)
  mature <- as.Date(mature)
  m <- -12/freq
  n <- coupons.n(settle, mature, freq)
  edate(mature, m * n)
}

##' @rdname bonds
##' @export
bond.prices <- function(settle, mature, coupon, freq=2, yield, 
                  convention=c("30/360", "ACT/ACT", "ACT/360", "30/360E"),
                  comp.freq=freq){
  convention <- match.arg(convention)
  mapply(bond.price, settle, mature, coupon, freq, yield, convention,
         comp.freq, USE.NAMES=FALSE)
}

##' @rdname bonds
##' @export
bond.yields <- function(settle, mature, coupon, freq=2, price,
                  convention=c("30/360", "ACT/ACT", "ACT/360",
                               "30/360E"),
                  comp.freq=freq){
  convention <- match.arg(convention)
  mapply(bond.yield, settle, mature, coupon, freq, price, convention,
         comp.freq, USE.NAMES=FALSE)
}

##' @rdname bonds
##' @export
bond.durations <- function(settle, mature, coupon, freq=2, yield,
                           convention=c("30/360", "ACT/ACT",
                                        "ACT/360", "30/360E"),
                     modified=FALSE, comp.freq=freq){
  convention <- match.arg(convention)
  mapply(bond.duration, settle, mature, coupon, freq, yield, convention,
         modified, comp.freq, USE.NAMES=FALSE)
}
