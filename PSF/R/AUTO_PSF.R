#' AUTO_PSF() takes input data and autogenerate optimum Window size (W) and cluster size (K)
#'
#' Predicts the "next_val" numbers of future values
#' @param data_in as Input data, in any format (data matrix data frame or vector). All variables should be numeric and NA values will get removed while execution.
#' @param next_val as Integer number. It states the number of predicted values to be obtained.
#' @return Function returns "next_val" numbers of future values and corresponding plot.
#' @export
#' @examples
#' ## Generate 100 random numbers within some limits
#' x <- sample(1:7, 100, replace = TRUE)
#' AUTO_PSF(x, 2)      # next_val = 5
#'
#' ## iris data as input data
#' AUTO_PSF (iris[2], 4) #next_val = 4

AUTO_PSF <- function(data_in, next_val)
{
  k <- optimum_k(data_in)
  w <- optimum_w(data_in,next_val)[1]
  w1 <- unlist(w)
  pred <- pred_for_w_plot(data_in,w1,k,next_val)[1]
  #pred <- unlist(pred)
  #pred1 <- data.frame(pred)
  options(warn=-1)
  return(pred)
}
