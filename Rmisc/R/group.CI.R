#' Group Confidence Interval
#' 
#' Calculates the confidence interval of grouped data
#' 
#' @keywords multivariate
#' 
#' @param x an `aggregate` compatible formula
#' @param data a data frame (or list) from which the variables in formula should be taken
#' @param ci the confidence interval to be calculated
#' 
#' @return A data frame consisting of one column for each grouping factor plus 
#' three columns for the upper bound, mean and lower bound of the confidence interval 
#' for each level of the grouping factor
#' 
#' @export
#' 
#' @examples
#' require(latticeExtra)
#' with(group.CI(weight~feed,chickwts),
#'  segplot(feed~weight.lower+weight.upper,center=weight.mean)
#' )
#' 
#' require(Hmisc)
#' with(group.CI(Temp~Month,airquality),
#'  xYplot(Cbind(Temp.mean,Temp.lower,Temp.upper)~numericScale(Month),type="b",ylim=c(60,90))
#' )
#' 
group.CI <-
function(x,data,ci=.95) {
 return(group.UCL(x,data,FUN=CI,ci=ci))
}
