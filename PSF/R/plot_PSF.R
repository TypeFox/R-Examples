#' Generate the plot for Forecasted Values
#'
#' Takes Input data and Predicted data as input
#' @param data_in as Input data, in any format (data matrix data frame or vector). All variables should be numeric.
#' @param predicted_data as predicted data, in any format (data matrix data frame or vector). All variables should be numeric
#' @return PSF_plot returns the plot of Predicted data
#' @export
#' @importFrom graphics plot points title
#' @examples
#' ## Generate 100 random numbers within some limits
#' x <- sample(1:7, 100, replace = TRUE)
#' y <- sample(1:4, 20, replace = TRUE)
#' plot_PSF(x,y)

plot_PSF <- function(data_in,predicted_data)
{
  if(!is.vector(data_in))
  {
    data_in <- data_in[, 1]
  }

  if(!is.vector(predicted_data))
  {
    predicted_data <- predicted_data[, 1]
  }


  data3 <- c(data_in,predicted_data)
  #PSF_plot <- plot(data_in,xlim=c(minX, maxX), ylim=c(minY, maxY),type = "o",col=c("red"), xlab = "Time Series Data", ylab=NA)
  PSF_plot <- plot(data_in,type = "o", xlim = c(min(data_in),length(data3)),  col=c("red"), xlab = "Time Series Data", ylab=NA)
  title("Forecast with PSF")
  points(length(data_in):(length(data_in)+length(predicted_data)),
         data3[length(data_in):(length(data_in)+length(predicted_data))],
         col="blue", type ="o",pch=22)
  options(warn=-1)
  return(PSF_plot)
}
