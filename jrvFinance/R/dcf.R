##' Net Present Value
##' 
##' Computes NPV (Net Present Value) for cash flows with different cash flow and
##' compounding conventions. Cash flows need not be evenly spaced.
##' 
##'
##' @param cf Vector of cash flows
##' @param cf.t Optional vector of timing (in years) of cash flows. If omitted
##' regular sequence of years is assumed.
##' @param rate The interest rate in decimal (0.10 or 10e-2 for 10\%)
##' @param cf.freq Frequency of annuity payments: 1 for annual, 2 for
##' semi-annual, 12 for monthly.
##' @param comp.freq Frequency of compounding of interest rates: 1 for annual,
##' 2 for semi-annual, 12 for monthly, Inf for continuous compounding.
##' @param immediate.start Logical variable which is \code{TRUE} when the first cash
##' flows is at the beginning of the first period (for example, immediate
##' annuities) and \code{FALSE} when the first cash flows is at the end of the first
##' period (for example, deferred annuities)
##' @export
npv <- function(cf, rate, cf.freq=1, comp.freq=1,
                cf.t = seq(from = if (immediate.start) 0 else 1/cf.freq,
                  by = 1/cf.freq, along.with = cf),
                immediate.start=FALSE){
  cc.rate <- equiv.rate(rate, from.freq=comp.freq, to.freq=Inf)
  df <- exp(-cc.rate * cf.t)
  as.numeric(crossprod(cf, df))
}


##' Internal Rate of Return
##' 
##' Computes IRR (Internal Rate of Return) for cash flows with different cash flow and
##' compounding conventions. Cash flows need not be evenly spaced.
##' 
##'
##' @inheritParams npv
##' @param interval the interval c(lower, upper) within which to
##'     search for the IRR
##' @param r.guess the starting value (guess) from which the solver
##'     starts searching for the IRR
##' @param toler the argument \code{toler} for
##'     \code{\link{irr.solve}}.  The IRR is regarded as correct if
##'     abs(NPV) is less than toler.  Otherwise the \code{irr}
##'     function returns \code{NA}
##' @param convergence the argument \code{convergence} for
##'     \code{\link{irr.solve}}
##' @param max.iter the argument \code{max.iter} for
##'     \code{\link{irr.solve}}
##' @param method The root finding method to be used. The
##'     \code{default} is to try Newton-Raphson method
##'     (\code{\link{newton.raphson.root}}) and if that fails to try
##'     bisection (\code{\link{bisection.root}}). The other two
##'     choices (\code{newton} and \code{bisection} force only one of
##'     the methods to be tried.
##' @export
irr <- function(cf, interval = NULL, cf.freq = 1, comp.freq = 1,
                cf.t = seq(from = 0, by = 1/cf.freq, along.with = cf),
                r.guess = NULL, toler = 1e-6, convergence = 1e-8,
                max.iter = 100,
                method = c('default', 'newton', 'bisection')){

    if (sign(max(cf)) == sign(min(cf)))
        return(NA)
    fn <- function(r){
        df <- exp(-r*cf.t)
        list(value=as.numeric(crossprod(cf, df)),
             gradient=-as.numeric(crossprod(cf*df,cf.t)))
    }
    method <- match.arg(method)
    res <- NA
    tryCatch(
        res <- irr.solve(f=fn, interval=interval, r.guess=r.guess, 
                         toler = toler, convergence = convergence,
                         max.iter = max.iter, method = method),
        warning = function(war){}
    )
    if(is.na(res)){
        warning(.irr.warning.msg(method))
    }
    equiv.rate(res, from.freq=Inf, to.freq=comp.freq)
}


##' Equivalent Rates under different Compounding Conventions
##' 
##' Converts an interest rate from one compounding convention to
##' another (for example from semi-annual to monthly compounding or
##' from annual to continuous compounding)
##' 
##'
##' @inheritParams npv
##' @param from.freq Frequency of compounding of the given interest
##'     rate: 1 for annual, 2 for semi-annual, 12 for monthly, Inf for
##'     continuous compounding.
##' @param to.freq Frequency of compounding to which the given
##'     interest rate is to be converted: 1 for annual, 2 for
##'     semi-annual, 12 for monthly, \code{Inf} for continuous
##'     compounding.
##' @export
equiv.rate <- function(rate, from.freq = 1, to.freq = 1){
  cc.rate <- ifelse(from.freq == Inf, rate,
                    log(1+rate/from.freq)*from.freq)
  if (to.freq == Inf) cc.rate else (exp(cc.rate/to.freq)-1)*to.freq
}

##' Duration and Modified Duration
##' 
##' Computes Duration and Modified Duration for cash flows with
##' different cash flow and compounding conventions. Cash flows need
##' not be evenly spaced.
##' 
##'
##' @inheritParams npv
##' @param modified in function duration, \code{TRUE} if modified
##'     duration is desired. \code{FALSE} otherwise.
##' @export
duration <- function(cf, rate, cf.freq=1, comp.freq=1,
                     cf.t = seq(from = ifelse(immediate.start, 0,
                                              1/cf.freq),
                       by = 1/cf.freq, along.with = cf),
                     immediate.start=FALSE, modified=FALSE){
  cc.rate <- equiv.rate(rate, from.freq=comp.freq, to.freq=Inf)
  df <- exp(-cc.rate * cf.t)
  D <- as.numeric(crossprod(cf*df, cf.t)/crossprod(cf, df))
  D / if (modified)  1 + rate/comp.freq else 1
}

