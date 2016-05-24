#' SIR model with two types of imports (P 6.6).
#' @description Solves a model with two types of stochastic imports and demographic stochasticit
#' @param pars \code{\link{vector}} with 5 parameters: transmission rate, recovery rate and per capita death rate. The names of these values must be "beta", "gamma", "mu", "epsilon" and "delta"  respectively. All parameters must be positive. The birth rate is assumed to be constant and equal to mu * N, therefore preventing extinction of the host population.
#' @param init \code{\link{vector}} with 3 values: the initial population size that are susceptible, infectious and the total population size. The names of these values must be "X", "Y" and "N", respectively. All initial conditions must be positive and all refer to integer numbers.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 6.6 from page 210 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' @return \code{\link{list}}. The first, second and third elements are the vectors (\code{*$pars}, \code{*$init}, \code{*$time}, containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fourth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. The following columns contain the number of susceptibles, infectious, recovered and boolean for epsilon and delta imports.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals}
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- c(beta = 1, gamma = 0.1, mu = 5e-4,
#'                 epsilon = 2e-5, delta = 0.01)
#' initials <- c(X = 5, Y = 2, N = 50)
#' 
#' # Solve and plot.
#' sir.2types.imports <- SIR2TypesImports(parameters, initials, 2 * 365)
#' PlotMods(sir.2types.imports)
#' 
SIR2TypesImports <- function(pars = NULL, init = NULL, time = NULL, ...) {
  if (is.null(pars)) {
    stop("undefined 'pars'")
  }
  if (is.null(pars)) {
    stop("undefined 'inits'")
  }
  if (is.null(pars)) {
    stop("undefined 'time'")
  }
  
  output <- c(time = 0, init['X'], init['Y'], Z = (init[['N']] - init[['X']] - init[['Y']]), ep = 0, dt = 0)
  res <-  data.frame(time = 0, X = init['X'], Y = init['Y'], Z = (init['N'] - init['X'] - init['Y']), ep = 0, dt = 0)
  N <- init[['N']]
  
  rate <- vector(length = 8, mode = 'numeric')
  change <- matrix(nrow = 8, ncol = 3)
  change[1,] <- c(-1, +1, 0)
  change[2,] <- c(0, -1, +1)
  change[3,] <- c(+1, 0, 0)
  change[4,] <- c(-1, 0, 0)
  change[5,] <- c(0, -1, 0)
  change[6,] <- c(0, 0, -1)
  change[7,] <- c(-1, +1, 0)
  change[8,] <- c(0, +1, 0)
  
  
  while(output[1] < time){
    
    N <- sum(output[2:4])
    rate[1] <- pars[['beta']] * output[2] * output[3] / N
    rate[2] <- pars[['gamma']] * output[3]
    rate[3] <- pars[['mu']] * N
    rate[4] <- pars[['mu']] * output[2]
    rate[5] <- pars[['mu']] * output[3]
    rate[6] <- pars[['mu']] * output[4]
    rate[7] = pars[['epsilon']] * output[2]
    rate[8] = pars[['delta']]
    
    r <- runif(2);
    
    ts <- -log(r[2])/(sum(rate));
    m <- min(which(cumsum(rate) >= r[1] * sum(rate)));
    
    output[1:4] <- output[1:4] + c(ts, change[m, ])
    output[5] = 0
    output[6] = 0
    if (m == 6)
      output[5] = 1
    if (m == 7)
      output[6] = 1
    res <- rbind(res, output)
    
  }
  
  rownames(res) <- NULL
  return(list(pars = pars,
              init = init,
              time = time,
              results = res[ , 1:4]))
}