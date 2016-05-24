#' @title Option pricing with jumps
#'
#' @description The functions \code{cashjump}, \code{assetjump}, and
#'     \code{mertonjump} return call and put prices, as vectors named
#'     "Call" and "Put", or "Call1", "Call2", etc. in case inputs are
#'     vectors. The pricing model is the Merton jump model, in which
#'     jumps are lognormally distributed. 
#'
#'
#' @seealso McDonald, Robert L., \emph{Derivatives Markets}, 3rd Edition
#'     (2013) Chapter 24
#'
#' @name jumps
#'
#' @aliases assetjump cashjump  mertonjump
#'
#' @return A vector of call and put prices computed using the Merton
#'     lognormal jump formula.
#'
#' @usage
#' assetjump(s, k, v, r, tt, d, lambda, alphaj, vj)
#' cashjump(s, k, v, r, tt, d, lambda, alphaj, vj)
#' mertonjump(s, k, v, r, tt, d, lambda, alphaj, vj)
#'
#'
#' @param s Stock price
#' @param k Strike price of the option
#' @param v Volatility of the stock, defined as the annualized
#'     standard deviation of the continuously-compounded return
#' @param r Annual continuously-compounded risk-free interest rate
#' @param tt Time to maturity in years
#' @param d Dividend yield, annualized, continuously-compounded
#' @param lambda Poisson intensity: expected number of jumps per year
#' @param alphaj Mean change in log price conditional on a jump
#' @param vj Standard deviation of change in log price conditional on
#'     a jump
#'
#' @details Returns a scalar or vector of option prices, depending on
#' the inputs
#' @importFrom stats dpois ppois
#' 
#' @seealso bscall bsput
#'
#' @examples
#' s <- 40; k <- 40; v <- 0.30; r <- 0.08; tt <- 2; d <- 0;
#' lambda <- 0.75; alphaj <- -0.05; vj <- .35;
#' bscall(s, k, v, r, tt, d)
#' bsput(s, k, v, r, tt, d)
#' mertonjump(s, k, v, r, tt, d, 0, 0, 0)
#' mertonjump(s, k, v, r, tt, d, lambda, alphaj, vj)
#'
#' ## following returns the same price as previous
#' c(1, -1)*(assetjump(s, k, v, r, tt, d, lambda, alphaj, vj) -
#' k*cashjump(s, k, v, r, tt, d, lambda, alphaj, vj))
#'
#' ## return call prices for different strikes
#' kseq <- 20:60
#' cp <- mertonjump(s, kseq, v, r, tt, d, lambda, alphaj,
#'     vj)[paste0('Call', 1:length(kseq))]
#'
#' ## Implied volatilities: Compute Black-Scholes implied volatilities
#' ## for options priced using the Merton jump model
#' vimp <- sapply(1:length(kseq), function(i) bscallimpvol(s, kseq[i],
#'     r, tt, d, cp[i]))
#' plot(kseq, vimp, main='Implied volatilities', xlab='Strike',
#'     ylab='Implied volatility', ylim=c(0.30, 0.50))


#' @export
cashjump <- function(s, k, v, r, tt, d, lambda, alphaj, vj) {
    Call <- .jumpprice(s, k, v, r, tt, d, lambda, alphaj, vj, cashcall)
    Put <- .jumpprice(s, k, v, r, tt, d, lambda, alphaj, vj, cashput)
    return(c(Call=Call, Put=Put))
}

#' @export
assetjump <- function(s, k, v, r, tt, d, lambda, alphaj, vj) {
    Call <- .jumpprice(s, k, v, r, tt, d, lambda, alphaj, vj,
               assetcall)
    Put <- .jumpprice(s, k, v, r, tt, d, lambda, alphaj, vj,
               assetput)
    return(c(Call=Call, Put=Put))
}

#' @export
mertonjump <- function(s, k, v, r, tt, d, lambda, alphaj, vj) {
    Call <- .jumpprice(s, k, v, r, tt, d, lambda, alphaj, vj,
               bscall)
    Put <- .jumpprice(s, k, v, r, tt, d, lambda, alphaj, vj,
               bsput)
    return(c(Call=Call, Put=Put))

}


.jumpprice <- function(s, k, v, r, tt, d, lambda, alphaj, vj,
                       pricingfunc) {
    eps <- 1e-10
    lambdap <- lambda * exp(alphaj)
    if (ppois(0, lambdap*tt) > 1-eps) {
        price <- pricingfunc(s, k, v, r, tt, d)
    } else {
        cum <- 0.0
        i <- 0
        kappa  <-  exp(alphaj) - 1
        while (ppois(i,lambdap*tt) <= 1-eps) {
            p <- dpois(i,lambdap * tt)
            vp <- (v^2 + i*vj^2 / tt)^(0.5)
            rp <- r - lambda * kappa + i * alphaj / tt
            cum <- cum + p * pricingfunc(s, k, vp, rp, tt, d)
            i <- i + 1
        }
        price <- cum
    }
    return(price)
}
