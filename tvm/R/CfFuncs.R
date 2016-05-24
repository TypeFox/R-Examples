#' @title Adjusts the discount factors by a spread
#' 
#' @param fd vector of discount factors used to discount cashflows in \code{1:length(fd)} periods
#' @param spread effective spread
#' @examples
#' adjust_disc(fd = c(0.99, 0.98), spread = 0.01)
#' @export
adjust_disc <- function(fd,spread) {
  zeros <- (1 / fd) ^ (1 / seq(along.with = fd))
  zeros_adj <- zeros + spread
  1/(zeros_adj ^ (seq(along.with = zeros_adj)))
}

#' @title Calculates the Total Financial Cost (CFT)
#' 
#' @description This is the IRR of the loan's cashflow, after adding all the extra costs
#' 
#' @details It is assumed that the loan has monthly payments
#' The CFT is returned as an effective rate of periodicty equal to that of the maturity and the rate
#' The interest is calculated over amt + fee
#' 
#' @param amt The amount of the loan
#' @param maturity The maturity of the loan
#' @param rate The loan rate, in  effective rate
#' @param up_fee The fee that the loan taker pays upfront
#' @param per_fee The fee that the loan payer pays every period
#' @examples
#' cft(amt = 100, maturity = 10, rate = 0.05, up_fee = 1, per_fee = 0.1)
#' @export
cft <- function(amt, maturity, rate, up_fee = 0, per_fee = 0) {
  full_amt <- amt + up_fee
  p <- pmt(amt = full_amt, maturity = maturity, rate = rate)
  rate(amt = amt, maturity = maturity, pmt = p + per_fee)
}

#' @title Net Present Value of a periodic cashflow (NPV)
#'  
#' @param i The rate used to discount the cashflow. It must be effective and with a periodicity that matches that of the cashflow
#' @param cf The cashflow
#' @param ts The times on which the cashflow ocurrs. It is assumed that \code{cf[idx]} happens at moment \code{ts[idx]}. If empty, assumes that \code{cf[idx]} happens at period \code{idx - 1}
#' 
#' @return The net present value at
#' 
#' @examples
#' npv(i = 0.01, cf = c(-1, 0.5, 0.9), ts = c(0, 1, 3))
#' @export
npv <- function(i, cf, ts = seq(from = 0, by = 1, along.with = cf)) sum(cf / (1 + i) ^ ts)

#' @title Net Present Value of an irregular cashflow (NPV)
#'  
#' @param i The rate used to discount the cashflow. It must be an effective anual rate (EAR)
#' @param cf The cashflow
#' @param d The dates when each cashflow occurs. Same length as the cashflow
#' @examples
#' xnpv(i = 0.01, cf = c(-1, 0.5, 0.9), d = as.Date(c("2015-01-01", "2015-02-15", "2015-04-10")))
#' @export
xnpv <- function(i, cf, d) sum(cf / ((1 + i) ^ (as.integer(d - d[1]) / 365)))

#' Internal Rate of Return of a periodic cashflow (IRR)
#' 
#' @title The IRR is returned as an effective rate with periodicity equal to that of the cashflow
#' 
#' @param cf The cashflow
#' @param ts The times on which the cashflow ocurrs. It is assumed that \code{cf[idx]} happens at moment \code{ts[idx]}
#' @param interval A length 2 vector that indicates the root finding algorithm where to search for the irr
#' @param ... Other arguments to be passed on to uniroot
#' @examples
#' irr(cf = c(-1, 0.5, 0.9), ts = c(0, 1, 3))
#' @export
irr <- function(cf, ts = seq(from = 0, by = 1, along.with = cf), interval = c(-1, 10), ...) { uniroot(npv, interval = interval, cf = cf, ts = ts, extendInt = "yes", ...)$root }

#' Internal Rate of Return of an irregular cashflow (IRR)
#' 
#' @title The IRR is returned as an effective anual rate
#' 
#' @param cf The cashflow
#' @param d The dates when each cashflow occurs. Same length as the cashflow
#' @param interval A length 2 vector that indicates the root finding algorithm where to search for the irr
#' @param ... Other arguments to be passed on to uniroot
#' @examples
#' xirr(cf = c(-1, 1.5), d = Sys.Date() + c(0, 365))
#' @export
xirr <- function(cf, d, interval = c(-1, 10), ...) { uniroot(xnpv, interval = interval, cf = cf, d = d, extendInt = "yes", ...)$root }

#' @title The value of the payment of a loan with constant payments (french type amortization)
#' 
#' @details The periodicity of the maturity and the rate must match, and this will be the periodicity of the payments
#'
#' @param amt The amount of the loan
#' @param maturity The maturity of the loan
#' @param rate The rate of the loan
#' @examples
#' pmt(amt = 100, maturity = 10, rate = 0.05)
#' @export
pmt <- function(amt, maturity, rate) {  
  return(amt*rate/(1 - (1 + rate) ^ (-maturity)))
}

#' @title The rate of a loan with constant payments (french type amortization)
#' 
#' @details The periodicity of the maturity and the payment must match, and this will be the periodicity of the rate (which is returned as an effective rate)
#' 
#' @param amt The amount of the loan
#' @param maturity The maturity of the loan
#' @param pmt The payments of the loan
#' @param extrema Vector of length 2 that has the minimum and maximum value to search for the rate
#' @param tol The tolerance to use in the root finding algorithm
#' @examples
#' rate(amt = 100, maturity = 10, pmt = 15)
#' @export
rate <- function(amt, maturity, pmt, extrema=c(1e-4,1e9), tol=1e-4) {   
  zerome <- function(r) amt/pmt - (1 - 1 / (1 + r) ^ maturity) / r
  if (zerome(extrema[1]) > 0) return(0)
  if (zerome(extrema[2]) < 0) return(extrema[2])
  return(uniroot(zerome, interval = extrema, tol = tol)$root)
}

rate <- Vectorize(FUN = rate,vectorize.args = c("amt","maturity","pmt"))

#' @title Creates an instance of a loan class
#'
#' @param rate The periodic effective rate of the loan
#' @param maturity The maturity of the loan, measured in the same units as the periodicity of the rate
#' @param amt The amount loaned
#' @param type The type of loan. Available types are \code{c("bullet","french","german")}
#' @param grace_int The number of periods that the loan doesn't pay interest and capitalizes it. Leave in 0 for zero loans
#' @param grace_amort The number of periods that the loan doesn't amortize
#' @examples
#' loan(rate = 0.05, maturity = 10, amt = 100, type = "bullet")
#' @export
loan <- function(rate, maturity, amt, type, grace_int = 0, grace_amort = grace_int) {
  stopifnot(grace_int <= grace_amort)
  stopifnot(grace_amort < maturity)
  l <- structure(list(rate = rate, maturity = maturity, amt = amt, type = type, grace_amort = 0, grace_int = 0), class = c(type,"loan"))
  if (grace_amort > 0 || grace_int > 0) {
    sl <- loan(rate = rate,maturity = maturity - grace_amort, amt = 1, type = type)
    l$cf <- amt*(1 + rate) ^ grace_int * c(rep_len(0, grace_int), rep_len(rate, grace_amort - grace_int), sl$cf)  
  } else {
    l$cf <- cashflow(l)  
  }
  l
}

#' @title Get the cashflow for a loan
#'
#' @description Returns the cashflow for the loan, excluding the initial inflow for the loan taker
#'
#' @param l The loan
#' @examples
#' l <- loan(rate = 0.05, maturity = 10, amt = 100, type = "bullet")
#' cashflow(l)
#' @export
cashflow <- function(l) {
  UseMethod("cashflow")
}

#' @method cashflow loan
#' @export
cashflow.loan <- function(l) {
  stop("Can't get cashflow for a loan without the proper type")
}

#' @method cashflow bullet
#' @export
cashflow.bullet <- function(l) {
  f <- rep_len(l$rate,l$maturity)
  f[l$maturity] <- 1 + f[l$maturity]
  f*l$amt
}

#' @method cashflow german
#' @export
cashflow.german <- function(l) {  
  k <- rep_len(1/l$maturity,l$maturity)
  krem <- c(1,1 - cumsum(k))
  i <- head(krem,l$maturity)*l$rate
  (k + i) * l$amt
}

#' @method cashflow french
#' @export
cashflow.french <- function(l) {  
  rep_len(l$rate / (1 - (1 + l$rate) ^ (-l$maturity)), l$maturity) * l$amt 
}

#' @title Value of a discounted cashflow
#' 
#' @param fd The discount factor vector
#' @param cf The cashflow
#' @examples
#' disc_cf(fd = c(1, 0.99, 0.98, 0.97), cf = c(1, -0.3, -0.4, -0.6))
#' @export
disc_cf <- function(fd, cf) {
  sum(fd * cf)
}

#' @title Remaining capital in a loan
#' 
#' @description The amount that has to be repayed at each moment in a loan, at the end of the period
#' 
#' @param cf The cashflow of the loan, not including the initial inflow for the loan taker
#' @param amt The original amount of the loan
#' @param r The periodic rate of the loan
#' @examples
#' rem(cf = rep_len(0.4, 4), amt = 1, r = 0.2)
#' @export
rem <- function(cf,amt,r) {
  s <- function(t) amt*(1 + r) ^ t - sum(cf[1:t] * (1 + r) ^ (t - (1:t)))
  vapply(X = seq_along(cf), FUN = s, FUN.VALUE = 1)
}

#' @title Find the rate for a loan given the discount factors
#' 
#' @description Thru a root finding process, this function finds the rate that corresponds to a
#' given set of discount factors, as for the loan to have the same present value discounted with the
#' discount factors or with that constant rate
#' 
#' @param m The maturity of the loan
#' @param d The discount factor vector
#' @param loan_type One of the loan types
#' @param interval The interval for the root finding process
#' @param tol The tolerance for the root finding process
#' @examples
#' find_rate(m = 3, d = c(0.99, 0.98, 0.97), loan_type = "bullet")
#' @export
find_rate <- function(m, d, loan_type, interval = c(1e-6, 2), tol = 1e-8) {
  if (length(d) < m) {
    stop("There aren't enough discount factors to discount the entire loan")
  }
  zerome <- function(r) {
    l <- loan(rate = r, maturity = m, amt = 1, type = loan_type)
    sum(d[seq_len(m)] * l$cf) - sum(l$cf / ((1 + r) ^ (seq_along(l$cf))))
  }
  uniroot(f = zerome, interval = interval, tol = tol)$root
}