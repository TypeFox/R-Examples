#' SIR model with births and deaths (P 2.2).
#' @description Solves a simple SIR model with equal births and deaths.
#' @param pars \code{\link{vector}} with 3 values: the per capita death rate (equal to the population level birth rate), the transmission rate, and the recovery rate. The names of these values must be "mu", "beta", and "gamma", respectively.
#' @param init \code{\link{vector}} with 3 values: the initial proportion of proportion of susceptibles, infectious and recovered. The names of these values must be "S", "I" and "R", respectively.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 2.2 from page 27 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' 
#' All parameters must be positive and S + I + R <= 1.
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second, third and fourth elements are the vectors \code{*$pars}, \code{*$init}, \code{*$time}, containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. Second, third and fourth columns contain the proportion of susceptibles, infectious and recovered.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals} 
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- c(mu = 1 / (70 * 365),
#'                        beta = 520 / 365, gamma = 1 / 7)
#' initials <- c(S = 0.1, I = 1e-4, R = 1 - 0.1 - 1e-4)
#' 
#' # Solve and plot.
#' sir.birth.death <- SIRBirthDeath(pars = parameters, init = initials, 
#'                                  time = 0:(60 * 365))
#' PlotMods(sir.birth.death)
#' 
SIRBirthDeath <- function(pars = NULL, init = NULL, time = NULL, ...) {
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
        dS <- mu - beta * S * I - mu * S
        dI <- beta * S * I - gamma * I - mu * I
        dR <- gamma * I - mu * R
        list(c(dS, dI, dR))
      })
    }
    init <- c(init['S'], init['I'], init['R'])
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