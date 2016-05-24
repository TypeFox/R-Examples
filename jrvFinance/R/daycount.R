##' Day count and year fraction for bond pricing
##' 
##' Implements 30/360, ACT/360, ACT/360 and 30/360E day count conventions. 
##' 
##'
##' @rdname daycount
##' @name daycount
##' @aliases daycount.30.360 daycount.actual yearFraction
##' 
##' @param d1 The starting date of period for day counts
##' @param d2 The ending date of period for day counts
##' @param r1 The starting date of reference period for ACT/ACT day counts
##' @param r2 The ending date of reference period for ACT/ACT day counts
##' @param freq The frequency of coupon payments: 1 for annual, 2 for
##' semi-annual, 12 for monthly.
##' @param convention The daycount convention
##' @param variant Three variants of the 30/360 convention are implemented, but
##' only one variant of ACT/ACT is currently implemented
##' @author Prof. Jayanth R. Varma \email{jrvarma@@iimahd.ernet.in}
##' @references The 30/360 day count was converted from C++ code in the
##' QuantLib library
NULL


##' @rdname daycount
##' @export
yearFraction <- function(d1, d2, r1, r2, freq=2,
                         convention=c("30/360", "ACT/ACT",
                                      "ACT/360", "30/360E")){
  convention <- match.arg(convention)
  if (convention == "ACT/ACT") (1/freq) * daycount.actual(d1,d2) /
                                   daycount.actual(r1,r2)
  else if (convention == "ACT/360") daycount.actual(d1,d2) / 360
  else if (convention == "30/360") daycount.30.360(d1,d2) / 360
  else if (convention == "30/360E") daycount.30.360(d1,d2, "E") / 360
}

##' @rdname daycount
##' @export
daycount.actual <- function(d1, d2, variant = c("bond")){
  as.integer(as.Date(d2) - as.Date(d1))
}

##' @rdname daycount
##' @export
daycount.30.360 <- function(d1, d2, variant = c("US", "EU", "IT")){
  ## The algorithm is taken from the QuantLib source code
  D1 <- as.POSIXlt(d1)
  D2 <- as.POSIXlt(d2)
  dd1 <- D1$mday
  dd2 <- D2$mday
  mm1 <- D1$mon
  mm2 <- D2$mon
  yy1 <- D1$year
  yy2 <- D2$year
  variant <- match.arg(variant)
  if (variant == "US" && dd2 == 31 && dd1 < 30) {
    dd2 <- 1
    mm2 <- mm2 + 1
  }
  if (variant == "IT" && mm1 == 2 && dd1 > 27) dd1 = 30
  if (variant == "IT" && mm2 == 2 && dd2 > 27) dd2 = 30
  360*(yy2-yy1) + 30*(mm2-mm1-1) + max(0,30-dd1) + min(30,dd2)
}

##' Shift date by a number of months
##' 
##' Convenience function for finding the same date in different months. Used for example
##' to find coupon dates of bonds given the maturity date. See \code{\link{coupons}}
##'
##' @param from starting date - a character string or any object that can be
##' converted into date using \code{\link{as.Date}}.
##' @param months Number of months (can be negative)
##' @export
edate <- function(from, months=1){
  from.lt <- as.POSIXlt(from)
  d <- from.lt$mday
  m <- from.lt$mon
  y <- from.lt$year + 1900
  M <- (m + months) %% 12
  Y <- y + (m + months - M) %/% 12
  ld <- last.day.of.month(m+1, y)
  LD <- last.day.of.month(M+1, Y)
  D <- if (d == ld || d > LD) LD else d
  as.Date(ISOdate(Y, M+1, D))
}

## Used internally not exported
last.day.of.month <- function(m, y){
  c(31,28,31, 30,31,30, 31,31,30, 31,30,31)[m] +
      if (is.leap.year(y) && m == 2) 1 else 0
}

## Used internally not exported
is.leap.year <- function(y){
  y %% 4 == 0 && (y %% 100 !=0 || y %% 400 == 0)
}

