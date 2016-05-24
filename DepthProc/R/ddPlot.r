#' @name ddPlot
#' @title Depth versus depth plot
#' @export
#' @description
#' Produces a DD plot which allows to compare two multivirate datasets or to compare a subject dataset with theoretical distribution.
#'
#' @param x The first or only data sample for ddPlot.
#' @param y The second data sample. \code{x} and \code{y} must be of the same space.
#' @param scale logical. determines whether the dispersion is to be aligned.
#' @param location determines whether the location is to be aligned to 0 vector with depth median.
#' @param name_x name for data set x. It will be passed to drawing function.
#' @param name_y name for data set y.
#' @param title title of the plot.
#' @param ... Parameters passed to depth function
#'
#' @details
#'  
#' For two probability distributions  \eqn{ F }  and  \eqn{ G } , both in  \eqn{ {{{R}}^{d}} } , we can define {depth vs. depth} plot being very useful generalization of the one dimensional quantile-quantile plot:    \deqn{ DD(F,G)=\left\{ \left( D({z},F),D({z},G) \right),{z}\in {{{R}}^{d}} \right\} }    
#' Its sample counterpart calculated for two samples  \eqn{ {{{X}}^{n}}=\{{{X}_{1}},.,{{X}_{n}}\} }  from  \eqn{ F } , and  \eqn{ {{Y}^{m}}=\{{{Y}_{1}},...,{{Y}_{m}}\} }  from  \eqn{ G }  is defined as
#' \deqn{ DD({{F}_{n}},{{G}_{m}})=\left\{ \left( D({z},{{F}_{n}}),D({z},{{G}_{m}}) \right),{z}\in \{{{{X}}^{n}}\cup {{{Y}}^{m}}\} \right\}}  
#'
#' @references
#' Liu, R.Y., Parelius, J.M. and Singh, K. (1999), Multivariate analysis by data depth: Descriptive statistics, graphics and inference (with discussion), \emph{Ann. Statist.}, \bold{27}, 822--831.
#'           
#'            Liu, R.Y., Singh K. (1993), A Quality Index Based on Data Depth and Multivariate Rank Test, \emph{Journal of the American Statistical Association} vol. 88.
#'  
#'  @author Daniel Kosiorowski, Mateusz Bocian, Anna Wegrzynkiewicz and Zygmunt Zawadzki from Cracow University of Economics.
#' 
#' @examples
#' require(sn)
#' require(mvtnorm)
#' 
#' # EXAMPLE 1: Location difference
#' standard = mvrnorm(1000, c(0,0), diag(2))
#' shift    =  mvrnorm(1000, c(0.5, 0), diag(2))
#' ddPlot(x = standard, y = shift, title = "Difference in position")
#' ddPlot(x = standard, y = shift, location = TRUE, title = "Location aligned")
#' 
#' ## EXAMPLE 2: Scale difference
#' standard <- mvrnorm(1000, c(0,0), diag(2))
#' scale <- mvrnorm(1000, c(0,0), 4*diag(2))
#' ddPlot(x=standard, y=scale)
#' ddPlot(x=standard, y=scale, scale=TRUE)
#'   
ddPlot <- function (x, y, scale = FALSE, location = FALSE, name_x = "X", name_y = "Y", title = "Depth vs. depth plot", ...) 
{
    if (ncol(x) != ncol(y)) {
      print("Wrong dimensions of the datasets! ncol(x)!=ncol(y)")
    }
    else {
      if (scale == TRUE) 
      {
        depth_sample_x <- depth(x, x, name = name_x, ...)
        depth_sample_y <- depth(y, y, name = name_y, ...)
        varcovx <- cov(x[which(depth_sample_x >= median(depth_sample_x)), ])
        varcovy <- cov(y[which(depth_sample_y >= median(depth_sample_y)), ])
        x_new <- t(solve(chol(varcovx)) %*% t(x))
        y_new <- t(solve(chol(varcovy)) %*% t(y))
      }
      else 
      {
        x_new <- x
        y_new <- y
      }
    }
    if (location == TRUE) 
    {
        medx <- depthMedian(x_new,  ...)
        medy <- depthMedian(y_new,  ...)
        x_new <- sweep(x_new, 2, medx, "-")
         y_new <- sweep(y_new, 2, medy, "-")
      
    }
    data <- rbind(x_new, y_new)
    depth_x <- depth(data, x_new,  name = name_x, ...)
    depth_y <- depth(data, y_new,  name = name_y, ...)
    
    ddplot = new("DDPlot", X = depth_x, Y = depth_y, title = title)

  
  return(ddplot)
}
