#' @title VaR for Generalized Pareto
#'
#' @description Estimates the Value at Risk of a portfolio assuming losses are 
#' distributed as a generalised Pareto.
#'
#' @param Ra Vector of daily Profit/Loss data 
#' @param beta Assumed scale parameter
#' @param zeta Assumed tail index
#' @param threshold.prob Threshold probability corresponding to threshold u and 
#' x
#' @param cl VaR confidence level
#' 
#' @return Expected Shortfall
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' McNeil, A., Extreme value theory for risk managers. Mimeo, ETHZ, 1999.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Computes ES assuming generalised Pareto for following parameters
#'    Ra <- 5 * rnorm(100)
#'    beta <- 1.2
#'    zeta <- 1.6
#'    threshold.prob <- .85
#'    cl <- .99
#'    GParetoVaR(Ra, beta, zeta, threshold.prob, cl)
#'
#' @export
GParetoVaR <- function(Ra, beta, zeta, threshold.prob, cl){
  
  if ( max(cl) >= 1){
    stop("Confidence level(s) must be less than 1")
  }
  if ( min(cl) <= 0){
    stop("Confidence level(s) must be greater than 0")
  }
  
  x <- as.vector(Ra)
  n <- length(x)
  x <- sort(x)
  Nu <- threshold.prob * n
  Nu <- ((Nu >= 0) * floor(Nu) + (Nu < 0) * ceiling(Nu))
  u <- x[n - Nu]
  y <- u+(beta/zeta)*((((1/threshold.prob)*(1-cl))^(-zeta))-1)
  
}