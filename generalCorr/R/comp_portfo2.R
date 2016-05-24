#' Compares two vectors (portfolios) using stochastic dominance of orders 1 to 4.
#'
#' Given two vectors of portfolio returns this function calls the internal function wtdpapb
#' to report the simple means of four sophsiticated measures of stochastic dominance.
#'
#' @param xa {data on returns for portfolio A in the form of a T by 1 vector}
#' @param xb {data on returns for portfolio B in the form of a T by 1 vector}
#' @return returns four numbers as averages of four sophisticated measures of stochastic
#' dominance measurements called SD1 to SD4.
#' @note It is possible to modify this function to report the median or standard
#' deviation or any other descriptive statistic by changing the line in the
#' code \code{oumean = apply(outb, 2, mean)} toward the end of the function.
#' A trimmed mean may be of interest when outliers are suspected.
#' @note require(np)
#' @note Make sure that functions wtdpapb, bigfp, stochdom2 are in the memory.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
## @seealso  \code{\link{help}},
#'## @references
#' @keywords stochastic dominance SD1 SD2 SD3 SD4 wtdpapb bigfp
#' @examples
#'
#' set.seed(30)
#' xa=sample(20:30)
#' xb=sample(32:40)
#' gp = comp_portfo2(xa, xb)
#' ##positive SD1 to SD4 means xb dominates xa as it should
#'
#' @export

comp_portfo2 <-
function(xa, xb){  #simplified:NAs already out
gp=wtdpapb(xa,xb)  #function above
stdo2=stochdom2(dj=gp$dj, wpa=gp$wpa, wpb=gp$wpb)
outb=cbind(stdo2$sd1b,stdo2$sd2b, stdo2$sd3b, stdo2$sd4b)
#print(outb)
#column sums for 4 orders of stochastic dominance
oumean=apply(outb,2,mean)
return(oumean)}
