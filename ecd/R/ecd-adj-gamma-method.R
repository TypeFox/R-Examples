#' Discriminant-adjusted gamma
#' 
#' Adjust gamma by diescriminant conversion formula so that
#' the critical line is a straight 45-degree line. 
#' The inverse adjustment is also provided.
#'
#' @param gamma numeric, the gamma paramter
#' @param adj_gamma numeric, the disciminant-adjusted gamma
#'
#' @return adjusted gamma (or the reverse of adjustment)
#'
#' @keywords discriminant
#'
#' @export ecd.adj_gamma
#' @export ecd.adj2gamma
#'
#' @examples
#' gamma2 <- ecd.adj_gamma(c(1,2))
#' gamma <- ecd.adj2gamma(c(1,2))

### <======================================================================>
"ecd.adj_gamma" <- function(gamma)
{
    sign(gamma)*abs(ecd.cusp_r2a(gamma))
}
### <---------------------------------------------------------------------->
#' @rdname ecd.adj_gamma
"ecd.adj2gamma" <- function(adj_gamma)
{
    sign(adj_gamma)*abs(ecd.cusp_a2r(adj_gamma))
}
### <---------------------------------------------------------------------->

