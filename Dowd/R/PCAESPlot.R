#' ES plot
#' 
#' Estimates ES plot using principal components analysis
#' 
#' @param Ra Matrix return data set where each row is interpreted as a set of daily observations, and each column as the returns to each position in a portfolio
#' @param position.data Position-size vector, giving amount invested in each position
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Computes PCA ES
#'    Ra <- matrix(rnorm(15*20),15,20)
#'    position.data <- rnorm(20)
#'    PCAESPlot(Ra, position.data)
#'
#' @export
PCAESPlot <- function(Ra, position.data){
  # Check that inputs have correct dimensions
  return.data<-as.matrix(Ra)
  pcaes.95 <- double(10)
  pcaes.99 <- double(10)
  for (i in 1:10) {
    pcaes.95[i] <- PCAES(return.data, position.data, i, .95)
    pcaes.99[i] <- PCAES(return.data, position.data, i, .99)
  }
  t <- 1:10
  par(mfrow=c(2,1))
  plot(t, pcaes.99, type="l")
  plot(t, pcaes.95, type="l")
  
}