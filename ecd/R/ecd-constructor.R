#' Constructor of ecd class
#' 
#' Construct an ecd class by providing the required parameters.
#' The default is the standard cusp distribution. 
#' Cusp is validated by \code{eps = max(.Machine$double.eps*1000, 1e-28)}.
#'
#' @param alpha numeric, the flatness parameter. Default: 0.
#' @param gamma numeric, the sharpness parameter. Default: 0.
#' @param sigma numeric, the scale parameter. Must be positive. Default: 1.
#' @param beta  numeric, the skewness parameter. Default: 0.
#' @param mu    numeric, the location parameter. Default: 0.
#' @param cusp  logical, indicate type of cusp. The singular points in cusp 
#'              requires special handling.
#'              Default: 0, not a cusp.
#'              1: cusp with alpha specified.
#'              2: cusp with gamma specified.
#' @param lambda numeric, the leading exponent for the special model. Default: 3.
#' @param with.stats logical, also calculate statistics, default is \code{TRUE}.
#' @param with.quantile logical, also calculate quantile data, default is \code{FALSE}.
#' @param bare.bone logical, skip both const and stats calculation, default is \code{FALSE}.
#'                  This for debug purpose for issues on integrating \eqn{e^y(x)}.
#' @param verbose logical, display timing information, for debugging purpose, default is \code{FALSE}.
#'
#' @return An object of ecd class 
#'
#' @keywords constructor
#'
#' @author Stephen H. Lihn
#'
#' @export ecd
#' 
#' @examples
#' d <- ecd()
#' d <- ecd(1,1)
#' d <- ecd(alpha=1, gamma=1)

### <======================================================================>
"ecd" <- function(alpha = 0, gamma = 0, sigma = 1, beta = 0, mu = 0, 
                  cusp = 0, lambda = 3,
                  with.stats = TRUE, with.quantile = FALSE,
                  bare.bone = FALSE, verbose=FALSE)
{
    call <- match.call()

    if(!is.numeric(alpha)){
      stop("Parameter 'alpha' must be numeric!\n")
    }
    if(!is.numeric(gamma)){
        stop("Parameter 'gamma' must be numeric!\n")
    }
    if(!is.numericMpfr(sigma)){
        stop("Parameter 'sigma' must be numericMpfr!\n")
    }
    if(sigma <= 0){
      stop("Parameter 'sigma' must be positive!\n")
    }
    if(!is.numeric(beta)){
        stop("Parameter 'beta' must be numeric!\n")
    }
    if(!is.numericMpfr(mu)){
        stop("Parameter 'mu' must be numericMpfr!\n")
    }
    if(!is.numeric(lambda)){
        stop("Parameter 'lambda' must be numeric!\n")
    }

    # assertion for lambda
    if( lambda != 3){
        if ( alpha != 0 | gamma != 0 ) {
            stop("When lambda is not 3, alpha and gamma must be zero!\n")
        }
    }
    
    
    # handle cusp
    if(alpha == 0 & gamma == 0 & beta == 0 & cusp == 0 & lambda == 3) {
        cusp <- 1
    }

    gamma2 <- ecd.adj_gamma(gamma)
    R <- sqrt(alpha^2+ gamma2^2)        
    theta2 <- ifelse(R==0, 0, acos(alpha/R))
    theta <- ifelse(gamma2>=0, theta2, 2*pi-theta2)
    
    model <- .ecd.model(alpha, gamma, sigma, beta, mu, cusp, lambda)

    use.mpfr <- function() {
    	sum <- alpha+gamma+sigma+beta+mu+lambda
    	ifelse(class(sum)=="mpfr", TRUE, FALSE)
    }

    d <- new("ecd", call = call,
               alpha = unname(alpha),
               gamma = unname(gamma),
               sigma = unname(sigma),
               beta  = unname(beta),
               mu    = unname(mu),
               cusp  = unname(cusp),
               lambda = unname(lambda),
               R     = unname(R),
               theta = unname(theta),
               use.mpfr = unname(use.mpfr()),
               const = NaN,
               stats = list(),
               quantile = new("ecdq"),
               model = model)
    
    if (bare.bone) return(d)
    
    # validate cusp
    if (cusp > 0) {
        if(! (d@alpha >= 0 & d@gamma <= 0)) {
            stop(paste("Failed to validate cusp by a>=0 and r<=0, alpha=",
                       d@alpha, "gamma=", d@gamma))
        }
        D <- discr(d)
        Dabs <- 16*(4*abs(d@gamma)^3 + 27*abs(d@alpha)^2)
        if (Dabs > 0) D <- D/Dabs
        eps <- max(.Machine$double.eps*1000, 1e-28)
        if (abs(D) > eps) {
            stop(paste("Failed to validate cusp by discr:", D))
        }
    }
    # -------------
    if (verbose) print(paste(Sys.time(), "ecd constructor: setup_const"))

    const(d) <- ecd.setup_const(d, verbose=verbose)
    
    # statistics
    if (with.stats) {
        if (verbose) print(paste(Sys.time(), "ecd constructor: stats"))
        stats(d) <- ecd.stats(d, verbose=verbose)
    }
    if (with.quantile) {
        if (verbose) print(paste(Sys.time(), "ecd constructor: quantile"))
        quantile(d) <- ecdq(d, verbose=verbose)
    }
    if (verbose) print(paste(Sys.time(), "ecd constructor: done"))
   
    invisible(d)
}
### <---------------------------------------------------------------------->
