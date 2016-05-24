#' Extended confidence interval
#' 
#' As in \code{\link{ci}},calculates the confidence interval around mean using
#' \code{\link{t.test}}, but also returns mean, sd, CV (Coefficient of Variation), 2 given Quantiles, min and max
#' 
#' @return A dataframe with the lower and upper confidence interval, as well as
#'        the level used, and mean, sd, CV (Coefficient of Variation), 2 given Quantiles, min and max
#' @note Since the discovery of \code{\link{summary}}, I don't really use this anymore.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2010
#' @seealso \code{\link{ci}}, \code{\link{t.test}}, \code{\link{summary}}
#' @keywords htest
#' @export
#' @examples
#' 
#' yourdata <- c(5:8,3,14)
#' cie(yourdata)          # confidence interval with the default confidence level (95%)
#' cie(yourdata, lev=0.99)# specified with a different confidence level 
#' cie(yourdata, 0.99, 4) # returns 4 decimal places
#' cie(yourdata, digits=2)# rounds to 2 decimal places with default level
#' cie                    # shows the function itself
#' 
#' @param dat vector with the data to use.
#' @param lev numeric. confidence level. DEFAULT: 0.95
#' @param digits integer. Number of digits rounded to in output. DEFAULT: 3
#' @param p1 numeric. lower percentile passed to \code{\link{quantile}}. DEFAULT: 0.05
#' @param p2 numeric. upper percentile passed to \code{\link{quantile}}. DEFAULT: 0.95
#' 
cie <- function(
  dat,
  lev=.95,
  digits=3,
  p1=0.05,
  p2=0.95)
{
t(round(data.frame(  CI.lower = t.test(dat, conf.level=lev)$conf.int[1] ,
                     CI.upper = t.test(dat, conf.level=lev)$conf.int[2] ,
                     level    = lev ,
                     mean     = mean(dat, na.rm=TRUE) ,
                     sd       = sd(dat, na.rm=TRUE) ,
                     CV       = sd(dat, na.rm=TRUE)/mean(dat, na.rm=TRUE) ,
                     median   = median(dat, na.rm=TRUE) ,
                     Quant.p1 = quantile(dat, prob=p1, na.rm=TRUE) ,
                     Quant.p2 = quantile(dat, prob=p2, na.rm=TRUE) ,
                     p1       = p1 ,
                     p2       = p2 ,
                     min      = min(dat, na.rm=TRUE) ,
                     max      = max(dat, na.rm=TRUE)
       ),digits))
}
