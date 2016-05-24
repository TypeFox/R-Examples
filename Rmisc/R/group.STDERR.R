#' Group Standard Error Interval
#' 
#' Calculates the standard error interval of grouped data.
#' 
#' @keywords multivariate
#' 
#' @param x an `aggregate` compatible formula
#' @param data a data frame (or list) from which the variables in formula should be taken.
#' 
#' @return A data frame consisting of one column for each grouping factor plus 
#' three columns for the upper bound, mean and lower bound of the standard error 
#' interval for each level of the grouping factor.
#' 
#' @export
#' 
#' @examples
#' require(latticeExtra)
#' with(group.STDERR(weight~feed,chickwts),
#'  segplot(feed~weight.lower+weight.upper,center=weight.mean)
#' )
#' 
#' require(Hmisc)
#' with(group.STDERR(Temp~Month,airquality),
#'  xYplot(Cbind(Temp.mean,Temp.lower,Temp.upper)~numericScale(Month),type="b",ylim=c(60,90))
#' )
#' 
group.STDERR <-
function(x,data) {
 return(group.UCL(x,data,FUN=STDERR))
}
