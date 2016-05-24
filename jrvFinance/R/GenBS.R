##' Generalized Black Scholes model for pricing vanilla European options
##'
##' Compute values of call and put options as well as the Greeks -
##' the sensitivities of the option price to various input arguments using the
##' Generalized Black Scholes model. "Generalized" means that the asset can
##' have a continuous dividend yield.
##' 
##' The Generalized Black Scholes formula for call options is \cr
##' \eqn{e^{-r t} (s \;  e^{g t} \; Nd1 - X \; Nd2)}{exp(-r * t) * (s * exp(g * t) * Nd1 - X * Nd2)} \cr
##' where \cr
##' \eqn{g = r - div\_yield} \cr
##' \eqn{Nd1 = N(d1)} and \eqn{Nd2 = N(d2)} \cr
##' \eqn{d1 = \frac{log(s / X) + (g + Sigma^2/ 2)  t}{Sigma \sqrt{t}}}{d1 = (log(s / X) + (g + Sigma^2/ 2) * t) / (Sigma * sqrt(t))} \cr
##' \eqn{d2 = d1 - Sigma \sqrt{t}}{d2 = d1 - Sigma * sqrt(t)} \cr
##' N denotes the normal CDF (\link{pnorm})\cr
##' For put options, the formula is \cr
##' \eqn{e^{-r t}  (-s \; e^{g  t} \; Nminusd1 + X \;  Nminusd2)}{exp(-r * t) * (-s * exp(g * t) * Nminusd1 + X * Nminusd2)}\cr
##' where \cr
##' \eqn{Nminusd1 = N(-d1)} and \eqn{Nminusd2 = N(-d2)} \cr
##'
##' @param s the spot price of the asset (the stock price for options on stocks)
##' @param X the exercise or strike price of the option
##' @param r the continuously compounded rate of interest in decimal (0.10 or 10e-2 for 10\%)
##' (use \code{\link{equiv.rate}} to convert to a continuously compounded rate)
##' @param Sigma the volatility of the asset price in decimal  (0.20 or 20e-2 for 20\%)
##' @param t the maturity of the option in years
##' @param div_yield the continuously compounded dividend yield (0.05 or 5e-2 for 5\%)
##' (use \code{\link{equiv.rate}} to convert to a continuously compounded rate)
##' @return A list of the following elements
##' \item{call}{the value of a call option}
##' \item{put}{the value of a put option}
##' \item{Greeks}{a list of the following elements}
##' \item{Greeks$callDelta}{the delta of a call option - the sensitivity to the spot price of the asset}
##' \item{Greeks$putDelta}{the delta of a put option - the sensitivity to the spot price of the asset}
##' \item{Greeks$callTheta}{the theta of a call option - the time decay of the option value
##' with passage of time. Note that time is measured in years. To find a daily theta divided by 365.}
##' \item{Greeks$putTheta}{the theta of a put option}
##' \item{Greeks$Gamma}{the gamma of a call or put option - the second derivative with respect to the spot price
##' or the sensitivity of delta to the spot price}
##' \item{Greeks$Vega}{the vega of a call or put option - the sensitivity to the volatility}
##' \item{Greeks$callRho}{the rho of a call option - the sensitivity to the interest rate}
##' \item{Greeks$putRho}{the rho of a put option - the sensitivity to the interest rate}
##' \item{extra}{a list of the following elements}
##' \item{extra$d1}{the d1 of the Generalized Black Schole formula}
##' \item{extra$d2}{the d2 of the Generalized Black Schole formula}
##' \item{extra$Nd1}{is \link{pnorm}(d1)}
##' \item{extra$Nd2}{is \link{pnorm}(d2)}
##' \item{extra$Nminusd1}{is \link{pnorm}(-d1)}
##' \item{extra$Nminusd2}{is \link{pnorm}(-d2)}
##' \item{extra$callProb}{the (risk neutral) probability that the call will be exercised = Nd2}
##' \item{extra$putProb}{the (risk neutral) probability that the put will be exercised = Nminusd2}
##'
##' @export
##' @importFrom stats pnorm dnorm
GenBS <- function(s, X, r, Sigma, t, div_yield = 0){
  g = r - div_yield ## the growth rate of the asset price
  d1 <- (log(s / X) + (g + Sigma*Sigma / 2) * t) / (Sigma * sqrt(t))
  ## the division above can be 0/0=NaN if s=X and either t=0 or g=Sigma=0
  ## we set d1 to 0 in this case since numerator is of higher order in t and Sigma
  d1 <- ifelse(is.nan(d1), 0, d1)
  d2 <- d1 - Sigma * sqrt(t)
  Nd1 <- pnorm(d1)
  Nd2 <- pnorm(d2)
  Nminusd1 <- pnorm(-d1)
  Nminusd2 <- pnorm(-d2)
  ## Black Scholes call price
  call <- exp(-r * t) * (s * exp(g * t) * Nd1 - X * Nd2)
  put <- exp(-r * t) * (-s * exp(g * t) * Nminusd1 + X * Nminusd2)
  callDelta <- exp(-div_yield*t)*Nd1
  putDelta <- -exp(-div_yield*t)*Nminusd1
  a <- -s * dnorm(d1) * Sigma * exp(-div_yield*t) / (2 * sqrt(t))
  ## if t is 0, the above is 0/0 = NaN but dnorm(d1) term goes to 0 faster
  ## so we set the result to 0 in this case
  a <- ifelse(is.nan(a), 0,  a)
  b = r * X * exp(-r * t) * Nd2
  c = div_yield * s * Nd1* exp(-div_yield*t)
  callTheta <- a - b + c
  b = r * X * exp(-r * t) * Nminusd2
  c = div_yield * s * Nminusd1 * exp(-div_yield*t)
  putTheta <- a + b - c
  Gamma <- dnorm(d1) * exp(-div_yield*t)/ (s * Sigma * sqrt(t))
  ## if t is 0, the above is 0/0 = NaN but dnorm(d1) term goes to 0 faster
  ## so we set the result to 0 in this case
  Gamma <- ifelse(is.nan(Gamma), 0,  Gamma)
  Vega <- s * sqrt(t) * dnorm(d1) * exp(-div_yield*t)
  callRho <- X * t * exp(-r * t) * Nd2
  putRho <- -X * t * exp(-r * t) * Nminusd2
  callProb <- Nd2
  putProb <- Nminusd2
  Greeks <- list(callDelta=callDelta, putDelta=putDelta, callTheta=callTheta,
                 putTheta=putTheta, Gamma=Gamma, Vega=Vega, callRho=callRho, putRho=putRho)
  extra <- list(d1=d1, d2=d2, Nd1=Nd1, Nd2=Nd2, Nminusd1=Nminusd1,
                Nminusd2=Nminusd2, callProb=callProb, putProb=putProb)
  list(call=call, put=put, Greeks=Greeks, extra=extra)
}

##' Generalized Black Scholes model implied volatility
##'
##' Find implied volatility given the option price using the generalized Black Scholes model.
##' "Generalized" means that the asset can have a continuous dividend yield.
##' 
##' \code{GenBSImplied} calls \code{\link{newton.raphson.root}} and
##' if that fails \code{\link{uniroot}}
##'
##' @param price the price of the option
##' @param PutOpt \code{TRUE} for put options, \code{FALSE} for call options
##' @param toler passed on to \code{\link{newton.raphson.root}}
##' The implied volatility is regarded as correct if the solver is able to
##' match the option price to within less than toler. Otherwise the function returns \code{NA}
##' @param max.iter passed on to \code{\link{newton.raphson.root}}
##' @param convergence passed on to \code{\link{newton.raphson.root}}
##' @inheritParams GenBS
##'
##' @export
##' @importFrom stats uniroot
GenBSImplied <- function(s, X, r, price, t, div_yield, PutOpt=FALSE,
                         toler=1e-6, max.iter=100, convergence=1e-8){
  ## discount the exercise price to eliminate r
  X = X * exp(-r * t)
  ## adjust the spot price to eliminate div_yield
  s = s * exp(-div_yield * t)
  ## use put call parity to convert put option into call option
  if (PutOpt) price = price + (s - X)
  SminusX = s - X
  SplusX = s + X
  if (price < SminusX || price < 0 || price > s){
    warning("Implied volatility is undefined because price is outside theoretical bounds")
    return(NA)
  }
  if (price == SminusX || price == 0) 
    ## if price equals intrinsic value, volatility is zero
    return(0)
  if (X == 0 && price != s) {
    ## if x is 0, option price must equal stock price
    warning("Implied volatility is undefined because price is outside theoretical bounds")
    return(NA)
  }
  ## we use a Taylor series approximation for an initial guess
  guess <- GenBSImpliedGuess(price, SminusX, SplusX, t)
  ## since price exceeds intrinsic value, 0 is a lower bound for the volatility
  ## we now seek an upper bound by doubling 100% volatility until the price is exceeded
  upper <- 100e-2
  while (GenBS(s, X, 0, upper, t, 0)$call < price && upper < .Machine$integer.max) upper <-  upper * 2
  if (upper > .Machine$integer.max) upper <-  Inf
  f <- function(sigma){
    temp <- GenBS(s, X, 0, sigma, t, 0)
    list(value = temp$call - price, gradient = temp$Greeks$Vega)
  }
  res <- newton.raphson.root(f=f, guess=guess, lower=0, upper=upper, max.iter=max.iter,  
                           toler=toler, convergence=convergence)
  if (! is.na(res))
    return (res)
  ## if newton.raphson.root fails, we try bisection
  res2 <- uniroot(function(x){f(x)$value}, c(0, upper), tol=convergence)
  if (abs(res2$f.root > toler)){
    warning("Error finding implied volatility")
    return (NA)                         # uniroot also failed so return NA
  }
  warning("Found implied volatility using uniroot after newton.raphson.root failed")
  ## Above warning is to prevent confusion from warning message in newton.raphson.root
  return (res2$root)
}


 GenBSImpliedGuess <- function(price, SminusX, SplusX, t){
   ## We use a Taylor series approximation similar to
   ## Brenner, H. and Subrahmanyam, M. G. (1994), "A simple approximation to option valuation
   ## and hedging in the Black Scholes model", Financial Analysts Journal, Mar-Apr 1994, 25-28
   ## Assume a call option and that r has been eliminated by discounting the exercise price 
   ## Approximate log(s/x) as 2(s-x)/(s+x) = SminusX(s-x)/H
   ## where SminusX = s-x and H=(s+x)/2.
   ## We then approximate the normal integral by a Taylor series
   ## if price is the call price, this gives the approximation
   ## price/H = SigmaRootT(1/Root2Pi +SminusX/(2*H*SigmaRootT)
   ##                            +SminusX^2/(2*H^2*SigmaRootT^2*Root2Pi)
   ## Letting z=SigmaRootT*H/Root2Pi, this yields a quadratic equation for z:
   ## z^2 - (price-SminusX/2)*z +SminusX^2/(4*Pi) = 0
   ## The solution is
   ## z = (price - SminusX/2)/2 + Radical
   ## where Radical is the square root of (price-SminusX/2)^2 - SminusX^2/Pi
   ## The linear approximation is obtained by dropping the constant term in the quadratic
   ## to give the linear equation:
   ## z = (price-SminusX/2)

  root2pi = sqrt(2 * pi)
  h = 0.5 * SplusX
  temp = (price - 0.5 * SminusX)
  radical = temp*temp - SminusX*SminusX / pi
  if (radical < 0)                   # Try Linear Approximation
    sigmaRootT = (root2pi / h) * temp
  else{                              # Try Quadratic Approximation
    radical = sqrt(radical)
    sigmaRootT = (root2pi / h) * (temp / 2 + radical)
  }
   sigmaRootT / sqrt(t)
}

