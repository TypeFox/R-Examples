#' @name ddmvnorm
#' @title Normal depth versus depth plot
#' @export
#' 
#' @param x The data sample for DD plot.
#' @param size size of theoretical set
#' @param robust Logical. Dafault \code{FALSE}. If \code{TRUE}, robust measures are used to specify the parameters of theoretical distribution.
#' @param alpha cutoff point for robust measure of covariance.
#' @param title title of a plot.
#' @param ... Parameters passed to \code{depth}
#' 
#' @description
#' Produces a  normal DD plot of a multivirate dataset.
#'
#' @details
#' 
#' In the first step the location and scale of x are estimated and theoretical sample from normal distribution with those parameters is generated. The plot presents the depth of empirical points with respect to dataset x and with respect to the theoretical sample.
#' 
#' @return
#' Returns the normal depth versus depth plot of multivariate dataset \code{x}. 
#'
#'
#' @references
#' 
#' Liu, R.Y., Parelius, J.M. and Singh, K. (1999), Multivariate analysis by data depth: Descriptive statistics, graphics and inference (with discussion),  Ann. Statist.,  27, 783-858.
#' 
#' Liu, R.Y., Singh K. (1993), A Quality Index Based on Data Depth and Multivariate Rank Test, \emph{Journal of the American Statistical Association} vol. 88.
#'
#'  @author Daniel Kosiorowski, Mateusz Bocian, Anna Wegrzynkiewicz and Zygmunt Zawadzki from Cracow University of Economics.
#'  
#' @seealso \code{\link{ddPlot}} to generate ddPlot to compare to datasets or to compare a dataset with other distributions.
#'
#' @examples
#' # EXAMPLE 1
#' norm = mvrnorm(1000, c(0,0,0), diag(3))
#' con = mvrnorm(100, c(1,2,5), 3*diag(3))
#' sample = rbind(norm, con)
#' ddMvnorm(sample, robust=TRUE)
#' 
#' # EXAMPLE 2
#' data(under5.mort,inf.mort,maesles.imm)
#' data1990=na.omit(cbind(under5.mort[,1],inf.mort[,1],maesles.imm[,1]))
#' ddMvnorm(data1990, robust=FALSE)
#' 
ddMvnorm <-function(x, size = nrow(x), robust=FALSE, alpha=0.05, title = "ddMvnorm", ...)     
{
  depth_sample <- depth(x, x,  ...)  
  
  if(robust == TRUE) 
  {
  	varcov <- cov(x[depth_sample>=quantile(depth_sample, alpha),])
  	location <- depthMedian(x, ...)
  } 
  else
  { 
  	location <- apply(x, 2, mean)
  	varcov  <- cov(x) 
  }
  theoretical <- mvrnorm(size, location, varcov)  
  depth_theoretical <- depth(x, theoretical, ...)
  ddplot = new("DDPlot",X = depth_sample, Y = depth_theoretical, title = title)
  
  return(ddplot)
}
