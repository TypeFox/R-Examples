#' @title A blocked correlated data simulation.
#' 
#' @description Creates a data simulation of \emph{n} observations with \emph{signal} 
#' groups of (\emph{p0/signal}) signal variables and (\emph{p-p0}) noise variables. Random noise is 
#' added to all columns. The default values, with \emph{n=100} create
#' the simulation of Zou and Hastie (2005).
#' 
#' @param n number of observations
#' @param p number of coordinate directions in the design matrix (default 40)
#' @param p0 number of signal coordinate directions in the design matrix (default 15)
#' @param signal number of signal groups (default 3)
#' @param sigma within group correlation coefficient (default sqrt(0.01))
#' @param beta.true specify the true simulation parameters. (default NULL = generated from other arguments)
#'    
#' @references Zou, H. and Hastie, T. (2005) Regularization and variable selection via the elastic net 
#'    \emph{J. R. Statist. Soc. B}, 67, Part 2, pp. 301-320
#'
#' @return list of 
#' \itemize{
#' \item x simulated design matrix
#' \item y simulated response vector
#' \item beta.true true beta parameters used to create the simulation
#' }
#' 
#' @examples
#' #--------------------------------------------------------------------------
#' # Example: Elastic net simulation
#' #  
#' # For elastic net simulation data, see Zou, H. and Hastie, T. (2005) 
#' # Regularization and variable selection via the elastic net J. R. Statist. Soc. B
#' # , 67, Part 2, pp. 301-320
#'   # Set the RNG seed to create a reproducible simulation
#'   set.seed(432) # Takes an integer argument
#'   
#'   # Creata simulation with 100 observations.
#'   dta <- elasticNetSim(n=100)
#'   
#'   # The simulation contains a design matrix x, and response vector y
#'   dim(dta$x)
#'   length(dta$y)
#'   print(dta$x[1:5,])
#'
#' @export elasticNetSim

elasticNetSim <- function(n, p=40, p0=15, signal=3, sigma=sqrt(0.01), beta.true=NULL){
  if(is.null(beta.true)){
    beta.true <- rep(0, p)
    beta.true[1:p0] <- signal
  }
  x1 <- rnorm(n) + matrix(rnorm(n * 5, sd = sigma), ncol = 5)
  x2 <- rnorm(n) + matrix(rnorm(n * 5, sd = sigma), ncol = 5)
  x3 <- rnorm(n) + matrix(rnorm(n * 5, sd = sigma), ncol = 5)
  x4 <- matrix(rnorm(n * 25), ncol = 25)
  x <- cbind(x1, x2, x3, x4)
  y <- x %*% beta.true + rnorm(n, sd = 15)
  return(list(x = x, y = y, beta.true=beta.true))
}
