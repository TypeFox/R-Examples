#' The moments, characteristic function (CF), and moment generating function (MGF)
#' of standard cusp distribution.
#' 
#' The moments of standard cusp distribution are calculated via Gamma function.
#' The CF and MGF are calculated as sum of moment terms. The CF is a complex number.
#' Since the terms in MGF is ultimately diverging, the sum is truncated before
#' the terms are increasing.
#'
#' @param n integer vector specifying the n-th moments
#' @param t numeric vector for CF and MGF
#' @param mu length-one numeric, specifying mean for CF and MGF
#' @param sigma length-one numeric, specifying volatility for CF and MGF
#' @param rel.tol relative tolerance
#' @param show.warning logical, to show warning or not.
#'
#' @return the values of the moments, CF, MGF
#'
#' @keywords ecd cusp
#'
#' @export ecd.cusp_std_moment
#' @export ecd.cusp_std_cf
#' @export ecd.cusp_std_mgf
#'
#' @importFrom utils tail
#'
#' @examples
#' ecd.cusp_std_moment(c(2,4))
#'
### <======================================================================>
"ecd.cusp_std_moment" <- function(n)
{
    n2 <- ifelse(floor(n/2)*2==n & n>=0, n, NaN)
    ifelse(!is.nan(n2), 2/sqrt(pi)*gamma(3/2*(n2+1)), 0)
}
### <---------------------------------------------------------------------->
#' @rdname ecd.cusp_std_moment
"ecd.cusp_std_cf" <- function(t, mu=0, sigma=1, rel.tol=1e-8, show.warning=FALSE)
{
    #
    if (length(t)>1) {
        return(sapply(t, ecd.cusp_std_cf,
                      mu=mu, sigma=sigma, rel.tol=rel.tol))
    }

    chr <- function(n) {
        mnt <- ecd.cusp_std_moment(n)
        fac <- gamma(n+1)
        pow <- (-1)^(n/2) * (sigma*t)^n
        mnt*pow/fac
    }
    
    NMAX <- 50
    BLOCK <- 10
    ttl <- 1
    err <- NULL
    for(n in seq(2, NMAX, by=BLOCK)) {
        c <- chr(seq(n, n+BLOCK-2, by=2))
        abs_c <- abs(c)
        i <- which(abs_c == min(abs_c))
        c <- c[1:i]
        ttl <- ttl + sum(c)
        err <- abs(tail(c,1)/ttl)
        if (err < rel.tol) break
        if (i < length(abs_c)) {
            if (show.warning & err > rel.tol) {
                warning(paste("std cusp cf truncation error larger than rel.tol:", 
                               sprintf("%.1e",err), "vs", rel.tol))
            }
            break
        }
    }
    exp(1i*t*mu) * ttl
}
### <---------------------------------------------------------------------->
#' @rdname ecd.cusp_std_moment
"ecd.cusp_std_mgf" <- function(t, mu=0, sigma=1, rel.tol=10e-8, show.warning=FALSE)
{
    #
    M <- ecd.cusp_std_cf(-1i*t, mu=mu, sigma=sigma, rel.tol=rel.tol,
                         show.warning=show.warning)
    if (Im(M) != 0) stop("ERROR: cusp_std_mgf must be a real number!")
    return(Re(M))
}
### <---------------------------------------------------------------------->


