#' Constructor of \code{ecdattr} class for the Elliptic Database (ECDB)
#' 
#' Construct an \code{ecdattr} class by providing the required parameters.
#' This object has one-to-one correspondence to the rows in \code{ECDATTR} table.
#' This is used primarily as object wrapper for safe update to ECDB.
#'
#' @param alpha numeric, must be an integer after multiplied by 1000000.
#' @param gamma numeric, must be an integer after multiplied by 1000000.
#'              \code{NaN} if cusp is 1.
#' @param cusp numeric, representing type of cusp. Only 0 (default) and 1 are allowed.
#' @param use.mpfr logical, whether to use mpfr for ecd object, default is \code{FALSE}.
#'
#' @return an object of \code{ecdattr} class
#'
#' @keywords ecdattr constructor
#'
#' @export ecdattr
#' @include ecdattr-class.R
#' 
#' @examples
#' a <- ecdattr(1,1)
#' b <- ecdattr(alpha=1, cusp=1)
#'
### <======================================================================>
"ecdattr" <- function(alpha, gamma=NaN, cusp=0, use.mpfr=FALSE)
{
    if (length(alpha) > 1) {
        stop("ecdattr constructor can only handle one object at a time")
    }
    if (length(gamma) > 1) {
        stop("ecdattr constructor can only handle one object at a time")
    }

    C <- 1000000
    alpha_m <- round(alpha * C)
    gamma_m <- round(gamma * C)

    if (cusp == 1) {
        if (! is.nan(gamma)) stop("gamma must be NaN when cusp=1")
        gamma <- ecd.cusp_a2r(alpha) # cusp is alpha driven
        gamma_m <- ceiling(gamma * C)   
    }
    # we also want to reverse validate the cusp setting for very few collisions
    else if (cusp == 0) {
        gamma2 <- ecd.cusp_a2r(alpha) # cusp is alpha driven
        gamma2_m <- ceiling(gamma2 * C)
        if (alpha >= 0 & gamma_m == gamma2_m) {
            cusp <- 1
        }
    } else {
        stop("Cusp must be 0 or 1")
    }
    
    call <- match.call()
    attr <- new("ecdattr", call = call,
                alpha = alpha,
                gamma = gamma,
                cusp = cusp,
                use.mpfr = use.mpfr,
                enriched = FALSE,
                ecd = new("ecd"),
                alpha_m = alpha_m,
                gamma_m = gamma_m,
                attr = list()
                )
    
    return(attr)
}
### <---------------------------------------------------------------------->
