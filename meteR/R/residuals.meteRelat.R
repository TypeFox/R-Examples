#' @title Compute residuals between METE predictions and date of a meteRelat object
#'
#' @description
#' \code{residuals.meteRelat} computes residuals between METE predictions and 
#' data of a meteDist object
#'
#' @details
#' See Examples. Typically not called directly by the user and rather used for 
#' calculating the mean square error with \code{mse.meteRelat}.
#' 
#' @param object a \code{meteRelat} object
#' @param ... arguments to be passed
#' @export
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(spp=arth$spp,
#'                 abund=arth$count,
#'                 power=arth$mass^(.75),
#'                 minE=min(arth$mass^(.75)))
#' ebar1 <- ebar(esf1)
#' residuals(ebar1)
#' 
#' @return a numeic vector giving residuals for each data point
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso \code{mse.meteDist}
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.
#' @importFrom stats approxfun

residuals.meteRelat <- function(object, ...) {
  xy <- object$pred
  names(xy) <- c('x', 'y')
  fun <- approxfun(xy)
  
  out <- object$obs[[2]] - fun(object$obs[[1]])
  
  return(out)
} 

#==========================================================================
# #' @title Computes mean squared error for meteRelat
# #'
# #' @description
# #' \code{mse.meteRelat} computes mean squared error for rank or cdf between METE prediction and data
# #'
# #' @details
# #' See Examples.
# #' 
# #' @param x a \code{meteRelat} object
# #' @export
# #' 
# #' @examples
# #' data(arth)
# #' esf1 <- meteESF(spp=arth$spp,
# #'                 abund=arth$count,
# #'                 power=arth$mass^(.75),
# #'                 minE=min(arth$mass^(.75)))
# #' ebar1 <- ebar(esf1)
# #' mse(ebar1)
# #' 
# #' @return numeric; the value of the mean squared error.
# #'
# #' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# #' @seealso mseZ.meteDist
# #' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# # @aliases - a list of additional topic names that will be mapped to
# # this documentation when the user looks them up from the command
# # line.
# # @family - a family name. All functions that have the same family tag will be linked in the documentation.
# 
# 
# mse.meteRelat <- function(x) {
#   resid <- residuals(x)
#   
#   return(mean(resid^2))
# }
