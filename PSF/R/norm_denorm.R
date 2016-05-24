#' Normalize the data
#'
#' takes a series of data and normalize it
#' @param data_in is a series of data
#' @return the normalized data series
#' @export
#' @examples
#' #iris data as input
#' x <- na.omit(iris[2])
#' y <- data_norm(x)
#' y
data_norm <- function(data_in)
{
  x <- (data_in - min(data_in))/(max(data_in)-min(data_in))
  dmax <- max(data_in)
  dmin <- min(data_in)
  output <- list("x"=x, "dmax"=dmax, "dmin"=dmin)
  return(output)
}

#' Denormalize tha data
#'
#' takes a series of normalized data along with min and max values and denormalize it
#' @param dnorm is the series of normalized data
#' @param dmax is the maximum value in the original (Data before normalization) series
#' @param dmin is the minimum value in the original (Data before normalization) series
#' @return the original data string after denormalization of normalized data
#' @export
data_denorm <- function(dnorm,dmax,dmin)
{
  y <- dnorm * (dmax - dmin) + dmin
  return(y)
}
