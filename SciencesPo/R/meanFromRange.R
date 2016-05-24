#' @title  Estimates Mean and Standard Deviation from the Median and the Range

#' @description When conductig a meta-analysis study, it is not always possible to recover from reports the mean and standard deviation values, but rather the medians and range of values. This function provides an approach to convert the median/range into a mean and a variance.
#'
#' @references
#' Hozo1, Stela P.; et al (2005) Estimating the mean and variance from the median, range, and the size of a sample. \emph{BMC Medical Research Methodology}, 5:13.
#'
#' @param low The min of the data.
#' @param med The median of the data.
#' @param high The max of the data
#' @param n The size of the sample.
#' @export
#' @examples
#' meanFromRange(5,8,12,10)
#'
`meanFromRange` <-function(low,med,high,n)
{
  mn<-(low+2*med+high)/4+(low-2*med+high)/(4*n)
  s=sqrt((low*low+med*med+high*high+(n-3)*((low+med)^2+(med+high)^2)/8-n*mn*mn)/(n-1))
  c(mn,s)
}
