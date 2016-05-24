#' Format and visualise functions used to generate virtual species with \code{\link{generateSpFromFun}}
#' 
#' @description
#' This function is a helper function to simplify the formatting of functions
#' for generateSpFromFun.
#' @details
#' This function formats the \code{parameters} argument of \code{\link{generateSpFromFun}}.  
#' For each environmental variable, provide a vector containing the function name, and its arguments.  
#' 
#' 
#' For example, assume we want to generate a species responding to two environmental variables bio1 and bio2.
#' \itemize{
#' \item{The response to bio1 is a normal response (\code{\link{dnorm}}), of mean 1 and standard deviation 0.5.}
#' \item{The response to bio2 is a linear response (\code{\link{linearFun}}), of slope (a) 2 and intercept (b) 5.}
#' }
#' The correct writing is:
#' 
#' \code{formatFunctions(
#' bio1 = c(fun = "dnorm", mean = 1, sd = 0.5),
#' bio2 = c(fun = "linearFun", a = 2, b = 5))}
#' 
#' 
#' 
#' @section Warning:
#' Do not use 'x' as a name for your environmental variables.
#' @param x NULL or a \code{RasterStack}. If you want to visualise the functions,
#' provide your \code{RasterStack} here.
#' @param rescale \code{TRUE} or \code{FALSE}. If \code{TRUE}, individual response
#' plots are rescaled between 0 and 1 with the formula (val - min) / (max - min).
#' @param ... the parameters to be formatted. See details.
#' @export
#' @import raster
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' with help from C. N. Meynard, C. Bellard & F. Courchamp
#' @examples
#' my.parameters <- formatFunctions(variable1 = c(fun = 'dnorm',
#'                                             mean = 0.00012, sd = 0.0001),
#'                               variable2 = c(fun = 'linearFun', a = 1, b = 0))
#' 
#' 
#' my.parameters <- formatFunctions(bio1 = c(fun = "logisticFun", 
#'                                          alpha = -12.7, beta = 68),
#'                                  bio2 = c(fun = "linearFun", 
#'                                           a = -0.03, b = 191.2),
#'                                  bio3 = c(fun = "dnorm", 
#'                                           mean = 86.4, sd = 19.1),
#'                                  bio4 = c(fun = "logisticFun", 
#'                                           alpha = 2198.5, beta = 11381.4))
#' \dontrun{
#' # An example using worldclim data
#' bio1.4 <- getData('worldclim', var='bio', res=10)[[1:4]]
#' my.parameters <- formatFunctions(x = bio1.4,
#'                                  bio1 = c(fun = "logisticFun", 
#'                                           alpha = -12.7, beta = 68),
#'                                  bio2 = c(fun = "linearFun", 
#'                                           a = -0.03, b = 191.2),
#'                                  bio3 = c(fun = "dnorm", 
#'                                           mean = 86.4, sd = 19.1),
#'                                  bio4 = c(fun = "logisticFun", 
#'                                           alpha = 2198.5, beta = 11381.4))
#' }
          
formatFunctions <- function(x = NULL, rescale = TRUE, ...)
{
  details <- list()
  args <- list(...)
  for (i in names(args))
  {
    if(!("fun" %in% names(args[[i]])))
    {
      stop(paste("No function was correctly provided for variable", i))
    }
    details[[i]]$fun <- args[[i]]["fun"]
    args[[i]] <- args[[i]][-(names(args[[i]]) %in% "fun")]
    details[[i]]$args <- as.list(args[[i]])
    a <- sapply(args[[i]], as.numeric)
    if(any(is.na(a)))
    {
      details[[i]]$args[[which(!is.na(a))]] <- sapply(details[[i]]$args[[which(!is.na(a))]], as.numeric)
    } else
    {
      details[[i]]$args <- a
    }
  }
  if (is(x, "Raster"))
  {
    
    plotResponse(x = x, parameters = details, rescale = rescale, approach = "response")
  } else if (!is.null(x))
  {
    stop("x must either be NULL or a raster stack of environmental variables")
  }
  return(details)
}

