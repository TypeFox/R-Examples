#' SIR model with demographic stochasticity (P 6.4).
#' @description Solves a SIR model with demographic stochasticity
#' @param pars \code{\link{vector}} with 3 parameters: transmission rate, recovery rate and per capita death rate. The names of these values must be "beta", "gamma" and "mu", respectively. All parameters must be positive and all rates are specified in days. The birth rate is assumed to be constant and equal to mu * N, therefore preventing extinction of the host population.
#' @param init \code{\link{vector}} with 3 values: the initial population size that are susceptible, infectious and the total population size. The names of these values must be "X", "Y" and "N", respectively. All initial conditions must be positive.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 6.4 from page 203 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' @return \code{\link{list}}. The first, second and third elements are the vectors (\code{*$pars}, \code{*$init}, \code{*$time}, containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fourth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. The following columns contain the proportion of susceptibles, infectious and recovered.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals}
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- c(beta = 1, gamma = 1 / 10, mu = 5e-4)
#' initials <- c(X = 500, Y = 25, N = 5e3)
#' 
#' # Solve and plot.
#' sir.demog.stoch <- SIRDemogStoch(pars = parameters, 
#'                                  init = initials, time = 2 * 365)
#' PlotMods(sir.demog.stoch)
#' 
SIRDemogStoch <- function(pars = NULL, init = NULL, time = NULL, ...) {
  if (is.null(pars)) {
    stop("undefined 'pars'")
  }
  if (is.null(pars)) {
    stop("undefined 'inits'")
  }
  if (is.null(pars)) {
    stop("undefined 'time'")
  }
  
  output <- c(time = 0, init['X'], init['Y'], Z = (init[['N']] - init[['X']] - init[['Y']]))
  res <-  data.frame(time = 0, X = init['X'], Y = init['Y'], Z = (init['N'] - init['X'] - init['Y']))
  N <- init[['N']]
  
  rate <- vector(length = 6, mode = 'numeric')
  change <- matrix(nrow = 6, ncol = 3)
  change[1,] <- c(-1, +1, 0)
  change[2,] <- c(0, -1, +1)
  change[3,] <- c(+1, 0, 0)
  change[4,] <- c(-1, 0, 0)
  change[5,] <- c(0, -1, 0)
  change[6,] <- c(0, 0, -1)
  
  while(output[1] < time){
    
    N <- sum(output[2:4])
    rate[1] <- pars[['beta']] * output[2] * output[3] / N
    rate[2] <- pars[['gamma']] * output[3]
    rate[3] <- pars[['mu']] * N
    rate[4] <- pars[['mu']] * output[2]
    rate[5] <- pars[['mu']] * output[3]
    rate[6] <- pars[['mu']] * output[4]
    
    r <- runif(2);
    
    ts <- -log(r[2])/(sum(rate));
    m <- min(which(cumsum(rate) >= r[1] * sum(rate)));
    
    output <- output + c(ts, change[m, ])
    res <- rbind(res, output)
    
  }
  
  rownames(res) <- NULL
  return(list(pars = pars,
              init = init,
              time = time,
              results = res))
}