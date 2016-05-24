#' SIR model with Scaled additive noise (P 6.2).
#' @description Solves a SIR model with scaled additive noise.
#' @param pars \code{\link{vector}} with 5 parameters: transmission rate, recovery rate, per capita death rate, the total population size and the number of steps that will change noise term. The names of these values must be "beta", "gamma", "mu", "N" and "step", respectively. All parameters must be positive and all rates are specified in days. The birth rate is assumed to be constant and equal to mu * N, therefore preventing extinction of the host population. Noise terms are generated as a function of the time step and its magnitude is a function of the rate of each process.
#' @param init \code{\link{vector}} with 2 values: the initial population size that are susceptible and infectious. The names of these values must be "X" and "Y", respectively. All initial conditions must be positive.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 6.2 from page 197 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second, third and fourth elements are the vectors (\code{*$pars}, \code{*$init}, \code{*$time}, containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. The following columns contain the proportion of susceptibles, infectious and recovered.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals} 
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- c(beta = 1, gamma = 1 / 10, mu = 1 / (50 * 365),
#'                 N = 1e6, step = 1)
#' initials <- c(X = 1e5, Y = 500)
#' 
#' # Solve and plot.
#' sir.scaled.additive.noise <- 
#' SIRScaledAdditiveNoise(pars = parameters, 
#'                        init = initials, time = 5 * 365)
#' PlotMods(sir.scaled.additive.noise)
#' 
SIRScaledAdditiveNoise <- function(pars = NULL, init = NULL, time = NULL, ...) {
  if (is.null(pars)) {
    stop("undefined 'pars'")
  }
  if (is.null(pars)) {
    stop("undefined 'inits'")
  }
  if (is.null(pars)) {
    stop("undefined 'time'")
  }
  init2 <- NULL
  init2 <- init
  function1 <- function(pars = NULL, init = NULL, time = NULL) {
    function2 <- function(time, init, pars) {
      with(as.list(c(init, pars)), {
        
        dX = (mu * N + sqrt(mu * N) * xi1) -
          (beta * X * Y / N + sqrt(beta * X * Y / N) * xi2) - 
          (mu * X + sqrt(mu * X) * xi3);
        
        # python equations same as matlab
        #         dX = (mu * N + sqrt(mu * N) * xi1) -
        #           (beta * X * Y / N + sqrt(beta * X * Y /N) * xi2)  -
        #           (mu * V[1] + np.sqrt(mu*V[1]) * xi3); # problema detectado aqui
        
        dY = (beta * X * Y / N + sqrt(beta * X * Y / N) * xi2) -
          (gamma * Y + sqrt(gamma * Y) * xi4) - (mu * Y + sqrt(mu * Y) * xi5)
        
        list(c(dX, dY))
      })
    }
    
    output <- ode(times = time, 
                  func = function2, 
                  y = init, parms = pars, ...)
    return(output)
  }
  
  T0 <- 0;
  MaxTime <- time;
  output <- data.frame(time = numeric(), X = numeric(), Y = numeric())
  sqrtStep <- sqrt(pars['step'])
  
  while (T0 <= MaxTime & init['X'] > 0 & init['Y'] > 0) {
    
    pars <- c(pars[1:5], xi = rnorm(5)/sqrtStep)
    
    output <- rbind.data.frame(output,
                               function1(pars = pars, init = init,
                                         time = seq(T0, T0 + pars["step"], 1)))
    init['X'] <- tail(output[,'X'],1)
    init['Y'] <- tail(output[,'Y'],1)
    T0 <- T0 + pars['step']
  }
  
  return(list(model = function1,
              pars = pars,
              init = init2,
              time = time,
              results = as.data.frame(unique(output))))
}