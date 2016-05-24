#' Scan Statistic
#' 
#' \code{scan.stat} calculates the scan statistic for various distributions.
#' 
#' @param yin The sum of the response values inside the window.  Generally, the sum of the cases.
#' @param ein The expected value of the response in the window.  Generally, the estimated overall risk for all regions combined, multiplied by the population size of the window.
#' @param eout The expected value of the response outside the window.
#' @param ty The sum of all responses in the study area.  Generally, the total number of cases.
#' @param type The type of scan statistic to implement.  Currently, only "poisson" is implemented.
#' @return A vector of scan statistics. 
#' @author Joshua French
#' @export
#' @references Kulldorff, M. (1997) A spatial scan statistic. Communications in Statistics -- Theory and Methods 26, 1481-1496.
#' @examples 
#' # statistic for most likely cluster of New York leukemia data
#' scan.stat(106, 62.13, 552 - 62.13, 552)

scan.stat = function(yin, ein, eout, ty, type = "poisson")
{
  if(length(yin) != length(ein)) 
    stop("length(yin) != length(ein)")
  if(length(yin) != length(eout)) 
    stop("length(yin) != length(eout)")
  if(type == "poisson")
  {
    yout = ty - yin
    tall = yin * (log(yin) - log(ein)) + yout * (log(yout) - log(eout))
    # correct test statistics 
    tall[yin/ein <= yout/eout] = 0
  }
  return(tall)
}