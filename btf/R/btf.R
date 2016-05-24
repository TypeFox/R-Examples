##' Bayesian trend filtering via Eigen
##'
##' Fits Bayesian trend filtering hierarchical model to univariate function.
##' Two conditional priors are available: double exponential or generalized double Pareto. 
##'
##' @param y response vector
##' @param x inputs corresponding to y observations
##' @param k degree of polynomial fit
##' @param iter number of samples to draw from posterior
##' @param cond.prior choose the conditional prior on f|sigma
##' @param alpha shape parameter for prior on lambda
##' @param rho rate parameter for prior on lambda
##' @param debug boolean telling btf to check for NaNs or not
##' @aliases btf
##' @author Edward A. Roualdes
##' @seealso \code{\link[genlasso]{trendfilter}}
##' @references R. J. Tibshirani. Adaptive piecewise polynomial estimation via trend filtering. The Annals of Statistics, 42(1):285-323, 2014.
##' @examples
##' # Cubic trend filtering
##' # from genlasso::trendfilter
##' \dontrun{n <- 100
##' beta0 = numeric(100)
##' beta0[1:40] <- (1:40-20)^3
##' beta0[40:50] <- -60*(40:50-50)^2 + 60*100+20^3
##' beta0[50:70] <- -20*(50:70-50)^2 + 60*100+20^3
##' beta0[70:100] <- -1/6*(70:100-110)^3 + -1/6*40^3 + 6000
##' beta0 <- -beta0
##' beta0 <- (beta0-min(beta0))*10/diff(range(beta0))
##' y <- beta0 + rnorm(n)
##' bfit <- btf(y=y, k=3)
##' plot(bfit, col='grey70')}
##' @export
btf <- function(y='vector', x=NULL, k='int', iter=1e4, cond.prior=c('gdp', 'dexp'), alpha=NULL, rho=NULL, debug=FALSE) {

    ## checks
    if ( missing(y) ) stop('must provide response vector y.')
    if ( !is.numeric(y) ) stop('repsonse vector y must be numeric.')
    if ( any(diff(x) == 0) ) stop("elements of x must be unique.")
    if ( k < 0 || round(k) != k ) stop("order k must be nonnegative integer.")
    if ( is.unsorted(x) ) stop("x must be in increasing order.")
    if ( k>3 ) warning(paste("For numerical stability, do not run Bayesian trend filtering with a polynomial order larger than 3."))
    n <- length(y)
    if ( missing(x) ) x <- seq_len(n)/n
    nx <- length(x)
    if ( nx != n ) stop("length of x and y differ.")

    D <- genDelta(n, k, x)
    
    ## which conditional prior?
    cond.prior <- match.arg(cond.prior)
    if ( cond.prior == 'gdp' ) {
        if ( missing(alpha) ) alpha <- -1.0 else alpha <- alpha  
        if ( missing(rho) ) rho <- 0.01 else rho <- rho
        
        ## run sampler
        chain <- gdPBTF(iter, y, k, D, alpha, rho, debug)

    } else if ( cond.prior == 'dexp' ) {
        if ( missing(alpha) ) alpha <- 1 else alpha <- alpha
        if ( missing(rho) ) rho <- 1e-4 else rho <- rho

        ## run sampler
        chain <- dexpBTF(iter, y, k, D, alpha, rho, debug)

    } else {
        stop("specified value of cond.prior not understood.")
    }

    ## tidying
    chain <- as.mcmc(chain)
    varnames(chain) <- c(paste('beta', seq_len(n), sep=''),
                         's2', 'lambda',
                         paste('omega', seq_len(n-k-1), sep=''), 'alpha')

    ## append some shit for plotting
    attr(chain, 'y') <- y
    attr(chain, 'x') <- x
    class(chain) <- c('btf', class(chain))
    chain
}
