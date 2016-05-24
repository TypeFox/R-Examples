#' Simple SIS model (P 2.5).
#' @description Solves a simple SIS model without births or deaths.
#' @param pars \code{\link{vector}} with 2 values: the transmission and recovery rates. The names of these values must be "beta", and "gamma", respectively.
#' @param init \code{\link{vector}} with 2 values: the initial proportion of proportion of susceptibles and infectious. The names of these values must be "S" and "I", respectively.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 2.5 from page 39 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' 
#' All parameters must be positive and S + I <= 1.
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second, third and fourth elements are the vectors (\code{*$pars}, \code{*$init}, \code{*$time}, containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. Second and third columns contain the proportion of susceptibles and infectious.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals} 
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- c(beta = 1.4247, gamma = 0.14286)
#' initials <- c(S = 1 - 1e-06, I = 1e-06)
#' 
#' # Solve and plot.
#' sis <- SIS(pars = parameters, init = initials, time = 0:70)
#' PlotMods(sis)
#' 
SIS <- function(pars = NULL, init = NULL, time = NULL, ...) {
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
      with(as.list(c(init, pars)), {
        dS <- - beta * S * I + gamma * I
        dI <- beta * S * I - gamma * I
        list(c(dS, dI))
      })
    }
    init <- c(init['S'], init['I'])
    output <- ode(times = time, 
                  func = function2, 
                  y = init, parms = pars, ...)
    return(output)
  }
  
  output <- function1(pars = pars, init = init, time = time)
  return(list(model = function1,
              pars = pars,
              init = init,
              time = time,
              results = as.data.frame(output)))
}