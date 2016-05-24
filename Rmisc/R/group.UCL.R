#' Group Upper-Center-Lower
#' 
#' Applies a function which calculates a parameter with lower/uper 
#' bounds to groups of data.
#' 
#' @keywords multivariate
#' 
#' @param x an `aggregate` compatible formula
#' @param data a data frame (or list) from which the variables in formula should be taken.
#' @param FUN the function to apply to each group
#' @param ... extra params passed on to aggregate
#' 
#' @return A data frame consisting of one column for each grouping factor plus 
#' three columns for the upper bound, mean and lower bound of the standard error 
#' interval for each level of the grouping factor.
#' 
#' @export
#' 
#' @examples
#' require(latticeExtra)
#' with(group.UCL(weight~feed,chickwts,FUN=CI),
#'  segplot(feed~weight.lower+weight.upper,center=weight.mean)
#' )
#' 
#' require(Hmisc)
#' with(group.UCL(Temp~Month,airquality,FUN=STDERR),
#'  xYplot(Cbind(Temp.mean,Temp.lower,Temp.upper)~numericScale(Month),type="b",ylim=c(60,90))
#' )
#'
group.UCL <-
function(x,data,FUN,...) {
 d <- aggregate(x,data,FUN=FUN,...)
 y <- colnames(d)[ncol(d)]
 n <- as.data.frame(d[,y])
 colnames(n) <- sapply(list("upper","mean","lower"),function(l){paste(y,l,sep=".")})
 d[ncol(d)] <- NULL
 return(cbind(d,n))
}
