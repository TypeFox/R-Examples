#' SIR model with constant additive noise (P 6.1).
#' @description Solves a SIR model with constant additive noise added to the transmission rate.
#' @param pars \code{\link{vector}} with 5 values: the transmission rate, the recovery rate, the birth (deadth) rate, the amount of noise experienced in the transmission rate and the population size assumed to be constant. The names of these values must be "beta", "gamma", "mu", "noise", and "N" respectively.
#' @param init \code{\link{vector}} with 2 values: the initial number of susceptibles and infectious. The names of these values must be "X", and "Y", respectively. "X" and "Y" must be positive and are numbers not proportions.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param step step size to set the integration step and to scale the noise term.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 6.1 from page 194 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second, third and fourth elements are the vectors \code{*$pars}, \code{*$init} and \code{*$time}, containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. Second and third columns contain the number of susceptibles and recovered.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals} 
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- c(beta = 1, gamma = 1 / 10, mu = 1 / (50 * 365),
#'                 noise = 10, N = 1e6)
#' initials <- c(X = 1e5, Y = 500)
#' 
#' # Solve and plot.
#' sir.additive.noise <- SIRAdditiveNoise(pars = parameters, init = initials,
#'                                        time = 0:(2 * 365), step = 1)
#' PlotMods(sir.additive.noise)
#' 
SIRAdditiveNoise <- function(pars = NULL, init = NULL, time = NULL, step = 1, ...) {
  if (is.null(pars)) {
    stop("undefined 'pars'")
  }
  if (is.null(pars)) {
    stop("undefined 'inits'")
  }
  if (is.null(pars)) {
    stop("undefined 'time'")
  }
  function1 <- function(pars = NULL, init = NULL, time = NULL) {
    function2 <- function(time, init, pars) {
      with(as.list(c(init, pars, step = step)), {
        Noise <- noise * rnorm(1) / sqrt(step)
        dX <- mu * N - beta * X * Y / N - Noise - mu * Y
        dY <- beta * X * Y / N + Noise -  mu * Y - gamma * Y
        list(c(dX, dY))
      })
    }
    
    init <- c(init['X'], init['Y'])
    output <- matrix(0, nrow = ceiling(time[length(time)] / step), ncol = 2)
    output <- rbind(as.numeric(init), output)
    t <- 1
    
    while (t <= time[length(time)] & init[1] > 0 & init[2] > 0) {
      sqrt.step <- sqrt(step)
      output.tmp <- ode(times = time[1]:step, 
                        func = function2, 
                        y = init, parms = pars, ...)
      init <- output.tmp[nrow(output.tmp), -1]
      output[t + 1, ] <- output.tmp[nrow(output.tmp), -1]
      t <- t + 1
    }
    
    output <- cbind(time, output)
    colnames(output) <- c('time', 'X', 'Y')
    return(output)
  }
  
  output <- function1(pars = pars, init = init, time = time)
  return(list(model = function1,
              pars = pars,
              init = init,
              time = time,
              results = as.data.frame(output)))
}