#' Median polish multidimensional.
#'
#' Fits an additive model for multidimensional array, using Tukey's median polish procedure.
#' 
#' @param data class \code{\link{ConstructMPst}}.
#' @param eps real number greater than \code{0}, default 0.01. A tolerance for convergence: see Details
#' @param maxiter the maximum number of iterations. Default 10.
#' @param na.rm logical. If the data contains NA's. Default TRUE.
#' @param \dots ignored.
#' @details The model fitted is additive \eqn{\mu +\alpha _{a} + \beta _{b} + \xi _{c} + \tau _{t}}, where \eqn{\mu} is an overall mean, \eqn{\alpha_{a}} is the \eqn{a}-th row effect, \eqn{\beta_{b}} is the effect \eqn{b}-th column effect, \eqn{\xi_{c}} is the \eqn{c}-th layer effect, \eqn{\tau _{t}} is the \eqn{t}-th time effect. The algorithm works by alternately removing medians of every spatio - temporal dimensions, and continues until the proportional reduction in the sum of absolute residuals is less than eps or until there have been maxiter iterations. If na.rm is FALSE the presence of any NA value in x will cause an error, otherwise NA values are ignored. MedianPolishM returns an object of class MedianPolishM (see below). There is a plotting method for this class, \code{\link{plot.MedianPolishM}}.
#' @return An object of class medpolish with the following named components in a list:
#' @return \item{residuals}{the residuals.}
#' @return \item{overall}{the fitted constant term.}
#' @return \item{effects}{the fitted every space - time effects.}
#' @return \item{iter}{number of iterations used in the range maxiter.}
#' @references Hoaglin, D. C., Mosteller, F., & Tukey, J. W. (Eds.). (2011). \emph{Exploring data tables, trends, and shapes} (Vol. 101). John Wiley & Sons.\href{http://www.wiley.com/WileyCDA/WileyTitle/productCd-047004005X.html}{[link]}
#' @importFrom reshape2 melt  
#' @export 
#' 
MedianPolishM.ConstructMPst <-
  function(data, eps, maxiter, na.rm,...){   
    
    stopifnot(inherits(data, "ConstructMPst"))    

    MedPolish<-
      MedianPolishM.default(data$Value, eps, maxiter, na.rm)
    MedPolish$data<-data$Value
    MedPolish$eps<-eps
    MedPolish$Gr<-2
    return(MedPolish)
    
  }

