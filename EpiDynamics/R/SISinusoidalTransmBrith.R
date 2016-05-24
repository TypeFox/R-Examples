#' Rabbit Hemorrhagic Disease model with sinusoidal transmission rate and per capita birth rate (P 5.4).
#' @description Solves the Rabbit Hemorrhagic Disease, in which both transmission rate and birth rates can be seasonally forced.
#' @param pars \code{\link{list}} with: the mean transmission rate, the amplitude of sinuoidal forcing (transmission), the mean birth rate, the amplitude of sinuoidal forcing for the birth rate, the frequency of the oscillations, the recovery rate, the per capita death rate,  the mortality rate due to infection, and the carrying capacity. The names of these values must be "beta0", "beta1", "alpha0", "alpha1", "w", "gamma", "mu", "m" and "K", respectively. All parameters must be positive, alpha1, beta1 <= 1.
#' @param init \code{\link{vector}} with 3 values: the initial numbers of susceptible hosts (rabbits), infectious hosts (rabbits) and total population size. The names of these values must be "X", "Y" and "N", respectively. All initial values must be positive and X(0) + Y(0) <= N(0).
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 5.4 from page 186 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second element is a list with the \code{*$pars} argument. The third and fourth elements are the vectors \code{*$init} and \code{*$time}, containing the \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. Second, third and fourth columns contain the number of susceptibles, infectious and recovered.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals}.
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- list(beta0 = 0.936, beta1 = 0.1, alpha0 = 0.02, alpha1 = 0.1,
#'                    w = 2 * pi / 365, gamma = 0.025,  mu = 0.01, m = 0.475,
#'                    K = 10000)
#' initials <- c(X = 0.5, Y = 0.01, N = 0.6)
#' 
#' # Solve and plot.
#' sis.rhdm <- SISinusoidalTransmBrith(pars = parameters,
#'                                     init = initials,
#'                                     time = 0:(60 * 365))
#' PlotMods(sis.rhdm)
#'                         
SISinusoidalTransmBrith <- function(pars = NULL, init = NULL, time = NULL, ...) {
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
        dX = alpha0 * (1 + alpha1 * sin(w * time)) * N - beta0 * 
          (1 + beta1 * sin(w * time)) * X * Y - (mu + N / K) * X # dX/dt
        dY = beta0 * (1 + beta1 * sin(w * time)) * X * 
          Y - (gamma + m + mu + N / K) * Y # dY/dt
        dN = (alpha0 * (1 + alpha1 * sin(w * time)) - mu - N / K) * 
          N - m * Y # dN/dt
        list(c(dX, dY, dN))
      })
    }
    init <- c(init['X'], init['Y'], init['N'])
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