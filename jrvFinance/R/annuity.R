##' Present Value of Annuity and Related Functions
##' 
##' Functions to compute present value and future value of annuities, to find
##' instalment given the present value or future value. Can also find the rate
##' or the number of periods given other parameters.
##' 
##' These functions are based on the Present Value relationship:
##' \deqn{pv = fv \cdot df = 
##' terminal.payment \cdot df + \frac{instalment (1 - df)}{r}}{pv = fv * df = 
##' terminal.payment * df + instalment * (1 - df) / r}
##' where
##' \eqn{df = (1 + r)^{-n.periods}} is the \eqn{n.periods} discount factor and
##' \eqn{r} is the per period interest rate computed using
##' rate, cf.freq and comp.freq.
##' 
##' It is intended that only one of \eqn{pv} or \eqn{fv} is used in any function call, but
##' internally the functions use \eqn{pv + fv \cdot df}{pv + fv * df} as the
##' LHS of the present value relationship under the assumption that only of the
##' two is non zero. 
##' 
##' The function annuity.instalment.breakup regards the annuity as a repayment
##' of a loan equal to pv plus the present value of terminal.payment. The
##' instalment paid in period period.no is broken up into the principal
##' repayment (amortization) and interest components.
##'
##' @rdname annuity
##' @name annuity
##' @aliases annuity.pv annuity.fv annuity.instalment annuity.periods
##' annuity.rate annuity.instalment.breakup
##' @param rate The interest rate in decimal (0.10 or 10e-2 for 10\%)
##' @param n.periods The number of periods in the annuity.
##' @param instalment The instalment (cash flow) per period.
##' @param pv The present value of all the cash flows including the terminal
##' payment.
##' @param fv The future value (at the end of the annuity) of all the cash
##' flows including the terminal payment.
##' @param terminal.payment Any cash flow at the end of the annuity. For
##' example, a bullet repayment at maturity of the unamortized principal.
##' @param immediate.start Logical variable which is \code{TRUE} for immediate
##' annuities (the first instalment is due immediately) and \code{FALSE} for deferred
##' annuities (the first instalment is due at the end of the first period).
##' @param cf.freq Frequency of annuity payments: 1 for annual, 2 for
##' semi-annual, 12 for monthly.
##' @param comp.freq Frequency of compounding of interest rates: 1 for annual,
##' 2 for semi-annual, 12 for monthly, Inf for continuous compounding.
##' @param period.no Used only in \code{annuity.instalment.breakup}. This is the
##' period for which the instalment needs to be broken up into principal and
##' interest parts.
##' @param round2int.digits Used only in \code{annuity.periods}. If the computed
##' number of periods is an integer when rounded to round2int.digits, then the
##' rounded integer value is returned. With the default value of 3, 9.9996 is
##' returned as 10, but 9.9994 and 9.39999999 are returned without any
##' rounding.
##' @return For most functions, the return value is one of the arguments
##' described above. For example \code{annuity.pv} returns \code{pv}. The only exception is
##' \code{annuity.instalment.breakup}. This returns a list with the following
##' components:
##' \item{opening.principal}{The principal balance at the beginning
##' of the period}
##' \item{closing.principal}{The principal balance at the end of
##' the period}
##' \item{interest.part}{The portion of the instalment which
##' represents interest}
##' \item{principal.part}{The portion of the instalment
##' which represents principal repayment}
##' @author Prof. Jayanth R. Varma \email{jrvarma@@iimahd.ernet.in}
NULL

##' @rdname annuity
##' @export
annuity.pv <- function(rate, n.periods=Inf, instalment=1, terminal.payment=0,
                       immediate.start=FALSE, cf.freq=1, comp.freq=1){
  if (rate == 0) return(n.periods * instalment + terminal.payment)
  r = equiv.rate(rate, comp.freq, cf.freq)/cf.freq
  df <- (1 + r)^-n.periods
  adjust <- if (immediate.start) 1 + r else 1
  adjust * (instalment * (1 - df) / r + terminal.payment * df)
}

##' @rdname annuity
##' @export
annuity.fv <- function(rate, n.periods=Inf, instalment=1, terminal.payment=0,
                       immediate.start=FALSE, cf.freq=1, comp.freq=1){
  if (rate == 0) return(n.periods * instalment)
  r = equiv.rate(rate, comp.freq, cf.freq)/cf.freq
  df <- (1 + r)^n.periods
  adjust <- if (immediate.start) 1 + r else 1
  adjust * instalment * (df - 1) / r + terminal.payment
}

##' @rdname annuity
##' @export
annuity.instalment <- function(rate, n.periods=Inf, pv=if(missing(fv)) 1 else 0, fv=0,
                               terminal.payment=0, immediate.start=FALSE, cf.freq=1, comp.freq=1){
  r = equiv.rate(rate, comp.freq, cf.freq)/cf.freq
  df <- (1 + r)^-n.periods
  annuity.pv = pv + (fv - terminal.payment) * df * if (immediate.start) 1 + r else 1
  if (rate == 0) return(annuity.pv /n.periods)
  adjust <- if (immediate.start) 1 + r else 1
  r * annuity.pv / (adjust * (1 - df))
}

##' @rdname annuity
##' @export
annuity.periods <- function(rate, instalment=1, pv=if(missing(fv)) 1 else 0,
                            fv=0, terminal.payment=0,
                            immediate.start=FALSE, cf.freq=1, comp.freq=1,
                            round2int.digits=3){
  if (rate == 0) return((pv + fv - terminal.payment) / instalment)
  r = equiv.rate(rate, comp.freq, cf.freq)/cf.freq
  ## if immediate start we remove the first instalment and
  ## consider the remaining deferred annuity with one less period
  pv = pv - if (immediate.start) instalment else 0
  ## pv + fv * df = terminal.payment * df + instalment * (1 - df) / r
  ## pv + fv * df = terminal.payment * df + instalment / r  - instalment * df / r
  ## fv * df - terminal.payment * df + instalment * df / r =  instalment / r  - pv
  ## df * (fv - terminal.payment  + instalment / r) =  instalment / r  - pv
  df <- (instalment / r  - pv) / (fv - terminal.payment  + instalment / r)
  n <- -log(df) / log(1 + r) + if (immediate.start) 1 else 0
  ## if the result is close to an integer, then round to the integer
  ## the tolerance for this is given by round2int.digits
  as.integer(n) + zapsmall(n - as.integer(n), round2int.digits)
}

##' @rdname annuity
##' @export
annuity.rate <- function(n.periods=Inf, instalment=1, pv=if(missing(fv)) 1 else 0,
                         fv=0, terminal.payment=0,
                         immediate.start=FALSE, cf.freq=1, comp.freq=1){
  if (n.periods == Inf)
    return (pv/instalment)
  f <- function(r){
    df <- (1 + r)^-n.periods            # n period discount factor
    dfg <- -n.periods * df / (1 + r)    # gradient of above
    adjust <- if (immediate.start) 1 + r else 1
    if (r == 0){
      af <- n.periods                          # n period annuity factor
      afg <- -n.periods * (n.periods + 1) / 2  # gradient of above
      ## print(data.frame(df, dfg, af, afg, adjust, instalment, terminal.payment, fv, pv,
      ##                  value=adjust * (instalment * af + terminal.payment * df - fv * df) - pv))
    }else{
      af <- (1 - df) / r                  # n period annuity factor
      afg <- (-(1 - df) - r * dfg)/ r^2       # gradient of above
    }
    adjust <- if (immediate.start) 1 + r else 1
    list(
      value= adjust * (instalment * af + terminal.payment * df - fv * df) - pv,
      gradient = adjust * (instalment * afg + terminal.payment * dfg - fv * dfg))
  }
  equiv.rate(irr.solve(f) * cf.freq, cf.freq, comp.freq)
}

##' @rdname annuity
##' @export
annuity.instalment.breakup <- function(rate, n.periods=Inf, pv=1, immediate.start=FALSE,
                                       cf.freq=1, comp.freq=1, period.no=1){
  instalment <- annuity.instalment(rate=rate, n.periods=n.periods, pv=pv,
                                 immediate.start=immediate.start,
                                 cf.freq=cf.freq, comp.freq=comp.freq)
  r = equiv.rate(rate, comp.freq, cf.freq)/cf.freq
  df <- (1 + r)^(period.no-1)
  opening.principal <- pv * df - annuity.fv(rate=rate, n.periods=period.no-1, 
                                            instalment=instalment,
                                            immediate.start=immediate.start,
                                            cf.freq=cf.freq, comp.freq=comp.freq,)
  list(opening.principal = opening.principal,
       interest.part = opening.principal * r,
       principal.part = instalment - opening.principal * r,
       closing.principal = opening.principal + opening.principal * r - instalment)
}
