#' diluted Earnings Per Share
#'
#' @param ni     net income
#' @param pd     preferred dividends
#' @param cpd    dividends on convertible preferred stock
#' @param cdi    interest on convertible debt
#' @param tax    tax rate
#' @param w      weighted average number of common shares outstanding
#' @param cps    shares from conversion of convertible preferred stock
#' @param cds    shares from conversion of convertible debt
#' @param iss    shares issuable from stock options
#' @seealso \code{\link{EPS}}
#' @seealso \code{\link{iss}}
#' @seealso \code{\link{was}}
#' @export
#' @examples
#' diluted.EPS(ni=115600,pd=10000,cdi=42000,tax=0.4,w=200000,cds=60000)
#' 
#' diluted.EPS(ni=115600,pd=10000,cpd=10000,w=200000,cps=40000)
#' 
#' diluted.EPS(ni=115600,pd=10000,w=200000,iss=2500)
#' 
#' diluted.EPS(ni=115600,pd=10000,cpd=10000,cdi=42000,tax=0.4,w=200000,cps=40000,cds=60000,iss=2500)
diluted.EPS <- function(ni, pd, cpd=0,cdi=0,tax=0,w,cps=0,cds=0,iss=0){
  basic = (ni-pd)/w
  diluted = (ni-pd+cpd+cdi*(1-tax))/(w+cps+cds+iss)
  if(diluted > basic){
    diluted = (ni-pd+cpd)/(w+cps+iss)
  }
  return(diluted)
}
