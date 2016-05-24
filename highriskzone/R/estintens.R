#' Estimates the intensity of the point pattern.
#'
#' Estimates the intensity of the point pattern by a kernel method
#' (See \code{\link[spatstat]{density.ppp}}).
#' @param ppdata  data of class ppp
#' @param covmatrix (Optional) Covariance matrix of the kernel of a normal distribution
#' @param weights (Optional) vector of weights attached to each observation
#' @importFrom stats density
#' @export
#' @return A list of
#'    \item{ intensest }{ Estimated intensity (object of class "im", see \code{\link[spatstat]{density.ppp}}). }
#'    \item{ covmatrix }{ Covariance matrix. If \code{covmatrix = NULL}, the matrix is estimated by \code{\link[ks]{Hscv}}. }
#' @seealso \code{\link[spatstat]{density.ppp}}, \code{\link[ks]{Hscv}}, \code{\link[spatstat]{eval.im}}
#' @examples
#' \dontrun{
#' data(craterA)
#' #change npixel = 50 with 1000 and you get a nicer picture
#' spatstat::spatstat.options(npixel=50)
#' int <- est_intens(craterA)
#' plot(int$intensest, main = "pixel image of intensity")
#' plot(craterA$window, main = "contour plot of intensity"); contour(int$intensest, add =TRUE)
#' }



est_intens <- function(ppdata, covmatrix=NULL, weights=NULL){


  #check if input arguments have correct values
  if ( !is.ppp(ppdata) ) {
    stop("data is not of class ppp")
  }
  if ( !is.null(covmatrix) && (!isSymmetric(covmatrix) | !all(eigen(covmatrix)$values > 0)) ){
    stop("covmatrix has to be symmetric and positive semidefinit")
  }

  # if covmatrix is not given, estimate it
  if(is.null(covmatrix)){
    covmatrix <- Hscv(cbind(ppdata$x, ppdata$y))
  }

  #estimate intensity
  lambda.hat <- density(ppdata, varcov=covmatrix, weights=weights)
  lambda.hat <- eval.im(ifelse(lambda.hat < 0, 0, lambda.hat))

  result <- list(intensest=lambda.hat, covmatrix=covmatrix)
  return(result)
}