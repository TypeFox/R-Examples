#' SIR model with disease induced mortality: frequency-dependent transmission (P 2.4).
#' @description Solves a SIR model with a probability of mortality, and unequal births and deaths
#' @param pars \code{\link{vector}} with 5 values: the probability than an infected individual dies from the disease before recovering, the per capita death rate from natural causes,  the population level birth rate, the transmission rate, and the recovery rate. The names of these values must be "rho", "mu", "nu", "beta", and "gamma", respectively. All parameters must be positive.
#' @param init \code{\link{vector}} with 3 values: the initial number of susceptibles, infectious and recovered. The names of these values must be "X", "Y" and "Z", respectively.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 2.4 from page 36 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second, third and fourth elements are the vectors (\code{*$pars}, \code{*$init}, \code{*$time}, containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. Second, third and fourth columns contain the number of susceptibles, infectious and recovered.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals} 
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters<- c(rho = 0.5,mu = 1 / (70 * 365.0),nu= 1 / (70 * 365.0),
#'                       beta = 520 / 365.0, gamma = 1 / 7)
#' initials <- c(X = 0.2, Y = 1e-4, Z = 0)
#' 
#' # Solve and plot.
#' sir.induced.mortality2 <- SIRInducedMortality2(pars = parameters, 
#'                                  init = initials, 
#'                                  time = 0:1e4)
#' PlotMods(sir.induced.mortality2)
SIRInducedMortality2 <- function(pars = NULL, init = NULL, time = NULL, ...) {
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
        dX <- nu - beta * X * Y / sum(X,Y,Z) - mu * X
        dY <- beta * X * Y / sum(X,Y,Z) - ((gamma + mu) / (1 - rho)) * Y
        dZ <- gamma * Y - mu * Z
        list(c(dX, dY, dZ))
      })
    }
    init <- c(init['X'], init['Y'], init['Z'])
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