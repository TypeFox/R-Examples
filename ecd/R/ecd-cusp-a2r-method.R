#' Conversion between alpha and gamma for cusp distribution
#' 
#' \code{ecd.cusp_a2r} converts from alpha to gamma.
#' \code{ecd.cusp_r2a} converts from gamma to alpha.
#'
#' @param alpha numeric
#' @param gamma numeric
#'
#' @return gamma for \code{a2r}; alpha for \code{r2a}.
#'
#' @keywords cusp
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecd.cusp_a2r
#' @export ecd.cusp_r2a
#'
#' @examples
#' gamma <- ecd.cusp_a2r(alpha=1)
#' alpha <- ecd.cusp_r2a(gamma=1)

### <======================================================================>
"ecd.cusp_a2r" <- function(alpha)
{
    -(27*abs(alpha)^2/4)^(1/3)
}
### <---------------------------------------------------------------------->
#' @rdname ecd.cusp_a2r
"ecd.cusp_r2a" <- function(gamma)
{
    sqrt(4*abs(gamma)^3/27)
}
