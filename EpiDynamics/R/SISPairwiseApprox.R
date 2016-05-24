#' Pairwise SIS approximation model (P 7.8).
#' @description Solves a pairwise approximation to the SIS model on a random network of N individuals, each with n contacts.
#' @param pars \code{\link{vector}} with 4 values: the transmission rate across a contact, the recovery rate for infectious individuals, the number of connections per individual in the population and the number of individuals in the population. The names of these values must be "tau", "gamma", "n" and "N" respectively.
#' @param init \code{\link{vector}} with 3 values: the initial number of of susceptibles, infectious and susceptible-infectious pairs. The names of these values must be "X", "Y" and "XY", respectively.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 7.8 from page 285 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' 
#' All parameters must be positive.
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second, third and fourth elements are the vectors \code{*$pars}, \code{*$init} and \code{*$time}, containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. Second, third and fourth columns contain the number of susceptibles, infectious and susceptible-infectious pairs.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals} 
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' n <- 4; N <- 1e4; Y <- 1; X <- N - Y
#' parameters <- c(tau = 0.1, gamma = 0.05, n = n, N = N)
#' initials <- c(X = X, Y = Y, XY = n * Y * X / N)
#' 
#' # Solve and plot.
#' sis.pairwise.approx <- SISPairwiseApprox(pars = parameters,
#'                                          init = initials, time = 0:100)
#' PlotMods(sis.pairwise.approx)
#' 
SISPairwiseApprox <- function(pars = NULL, init = NULL, time = NULL, ...) {
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
        dX <- gamma * (N - X) - tau * XY
        dXY <- tau * (n - 1) * (n * X - XY) * XY / (n * X) +
          gamma * (n * N - n * X - XY) - tau * XY - tau * (n - 1) * XY *
          XY / (n * X) - gamma * XY
        list(c(dX, NA, dXY))
      })
    }
    init <- c(init['X'], init['Y'], init['XY'])
    output <- ode(times = time, 
                  func = function2, 
                  y = init, parms = pars, ...)
    output[ , 'Y'] <- pars['N'] - output[ , 'X']
    return(output)
  }
  
  output <- function1(pars = pars, init = init, time = time)
  return(list(model = function1,
              pars = pars,
              init = init,
              time = time,
              results = as.data.frame(output)))
}