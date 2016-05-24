#' Polar constructor of ecd class
#' 
#' Construct an ecd class by specifying R and theta. They are converted to 
#' \code{alpha} and \code{gamma}, then passed onto the ecd constructor.
#'
#' @param R numeric, the radius parameter. Default is \code{NaN}.
#' @param theta numeric, the angle parameter. Default: \code{NaN}.
#' @param sigma numeric, the scale parameter. Must be positive. Default: 1.
#' @param beta  numeric, the skewness parameter. Default: 0.
#' @param mu    numeric, the location parameter. Default: 0.
#' @param cusp  logical, indicate type of cusp (0,1,2).
#' @param with.stats logical, also calculate statistics, default is \code{TRUE}.
#' @param with.quantile logical, also calculate quantile data, default is \code{FALSE}.
#' @param bare.bone logical, skip both const and stats calculation, default is \code{FALSE}.
#'                  This for debug purpose for issues on integrating \eqn{e^y(x)}.
#' @param verbose logical, display timing information, for debugging purpose, default is \code{FALSE}.
#'
#' @return The ecd class 
#'
#' @keywords ecd constructor
#'
#' @author Stephen H. Lihn
#'
#' @export
#' 
#' @examples
#' d <- ecd.polar(R=1, theta=0.5*pi)

### <======================================================================>
"ecd.polar" <- function(R = NaN, theta = NaN, sigma = 1, beta = 0, mu = 0, cusp = 0,
                        with.stats = TRUE, with.quantile = FALSE, 
                        bare.bone = FALSE, verbose = FALSE)
{
    if (R < 0) {
        stop("R must be positive")
    }
    
    # cusp override
    if( abs(theta/pi - 7/4) <= .Machine$double.eps) {
        cusp <- 1
    }

    if (cusp > 0) {
        if(! is.na(theta)) {
            if( abs(theta/pi - 7/4) > .Machine$double.eps) {
                stop("Cusp constructor requires theta be either NaN or exactly 7/4 pi within machine double eps")
            }
        }
        
        th45 <- 45/180*pi
        
        if (cusp==1) return(ecd.cusp(
        	alpha= R*cos(th45),
        	sigma = sigma, mu = mu,
        	with.stats = with.stats, 
        	with.quantile = with.quantile,
        	bare.bone = bare.bone,
            verbose = verbose
        	))
        	
        if (cusp==2) return(ecd.cusp(
        	gamma= ecd.adj2gamma(-R*sin(th45)),
        	sigma = sigma, mu = mu,
        	with.stats = with.stats, 
        	with.quantile = with.quantile,
        	bare.bone = bare.bone,
            verbose = verbose
        	))
        
        stop(paste("Unknown cusp value:", cusp))
    }

    if (cusp != 0) {
        stop("Failed to determine the type of cusp")
    }
    
    ecd(alpha = R*cos(theta),
        gamma = ecd.adj2gamma(R*sin(theta)),
        sigma = sigma, beta = beta, mu = mu,
        cusp  = cusp, 
        with.stats = with.stats, 
        with.quantile = with.quantile,
        bare.bone = bare.bone,
        verbose = verbose)
}
### <---------------------------------------------------------------------->
