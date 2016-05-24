#' Cusp constructor of ecd class
#' 
#' Construct an ecd class for cusp distribution by specifying either alpha or gamma,
#' but not both. At the moment, it can't handle beta. 
#'
#' @param alpha numeric, the flatness parameter. Default: NaN.
#' @param gamma numeric, the sharpness parameter. Default: NaN.
#' @param sigma numeric, the scale parameter. Must be positive. Default 1.
#' @param mu    numeric, the location parameter. Default: 0.
#' @param with.stats logical, also calculate statistics, default is \code{TRUE}.
#' @param with.quantile logical, also calculate quantile data, default is \code{FALSE}.
#' @param bare.bone logical, skip both const and stats calculation, default is \code{FALSE}.
#'                  This for debug purpose for issues on integrating \eqn{e^y(x)}.
#' @param verbose logical, display timing information, for debugging purpose, default is \code{FALSE}.
#'
#' @return The ecd class 
#'
#' @keywords ecd cusp constructor
#'
#' @author Stephen H. Lihn
#'
#' @export
#' 
#' @examples
#' d <- ecd.cusp(alpha=1)
#' d <- ecd.cusp(gamma=-1)

### <======================================================================>
"ecd.cusp" <- function(alpha = NaN, gamma = NaN, sigma = 1, mu = 0,
                       with.stats = TRUE, with.quantile = FALSE,
                       bare.bone = FALSE, verbose = FALSE)
{
    cusp <- 0
    # type 1
    if (! is.na(alpha)) {
        if(! is.na(gamma)) {
            stop("Type 1 cusp constructor requires gamma=NaN")
        }
        cusp <- 1
        gamma <- ecd.cusp_a2r(alpha)
    }
    else if (! is.na(gamma)) {
        if(! is.na(alpha)) {
            stop("Type 2 cusp constructor requires alpha=NaN")
        }
        cusp <- 2
        alpha <- ecd.cusp_r2a(gamma)
    }
    if (cusp == 0) {
        stop("Failed to determine the type of cusp")
    }
    
    ecd(alpha = alpha,
        gamma = gamma,
        sigma = sigma,
        mu    = mu,
        cusp  = cusp,
        with.stats = with.stats, 
        with.quantile = with.quantile,
        bare.bone = bare.bone,
        verbose = verbose)
}
### <---------------------------------------------------------------------->
