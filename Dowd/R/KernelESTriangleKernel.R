#' Calculates ES using triangle kernel approach
#' 
#' The output consists of a scalar ES for specified confidence level.
#' 
#' @param Ra Profit and Loss data set
#' @param cl VaR confidence level
#' @return Scalar VaR
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # VaR for specified confidence level using triangle kernel approach
#'    Ra <- rnorm(30)
#'    KernelESTriangleKernel(Ra, .95)
#'
#' @export
KernelESTriangleKernel <- function(Ra, cl){
  PandL <- as.vector(Ra)
  n <- 1000
  delta.cl <- (1 - cl) / n
  VaR <- double(999)
  for (i in 1:(n - 1)) {
    if(i<(n-1)){
      VaR[i] <- KernelVaRTriangleKernel(PandL, cl + i * delta.cl, FALSE)
    } else if (i == n-1) {
      VaR[i] <- KernelVaRTriangleKernel(PandL, cl + i * delta.cl, TRUE)
    }
  }
  ES <- mean(VaR)
  return(ES)
  
}
