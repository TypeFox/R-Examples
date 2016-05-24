#' Calculates ES using box kernel approach
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
#'    # VaR for specified confidence level using box kernel approach
#'    Ra <- rnorm(30)
#'    KernelESBoxKernel(Ra, .95)
#'
#' @export
KernelESBoxKernel <- function(Ra, cl){
  PandL <- as.vector(Ra)
  n <- 1000
  delta.cl <- (1 - cl) / n
  VaR <- double(999)
  for (i in 1:(n - 1)) {
    if(i<(n-1)){
      VaR[i] <- KernelVaRBoxKernel(PandL, cl + i * delta.cl, FALSE)
    } else if (i == n-1) {
      VaR[i] <- KernelVaRBoxKernel(PandL, cl + i * delta.cl, TRUE)
    }
  }
  
  ES <- mean(VaR)
  return(ES)
  
}
