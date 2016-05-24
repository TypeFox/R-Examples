#' SIR model with sinusoidal forcing (P 5.1).
#' @description Solves a SIR model with sinusoidal forcing of the transmission rate.
#' @param pars \code{\link{list}} with 4 values: the death rate, the mean transmission rate, a scalar (or a \code{\link{vector}} to create bifurcations) with the amplitude of sinusoidal forcing, the frequency of oscillations and the recovery rate. The names for these values must be: "mu", "beta0", "beta1", "omega" and "gamma", respectively.
#' @param init \code{\link{vector}} with 3 values: initial proportion of proportion of susceptibles, infectious and recovered. The names of these values must be "S", "I" and "R", respectively. All parameters must be positive and S + I <= 1.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 5.1 from page 160 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' 
#' When beta1 is a vector in \code{pars}, it must be a sequence between 0 and 1.
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second, third and fourth elements are the vectors (\code{*$pars}, \code{*$init}, \code{*$time}, containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. Second, third and fourth columns contain the proportion of susceptibles, infectious and recovered.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals} 
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- list(beta0 = 17 / 13, beta1 = 0.1, gamma = 1 / 13,
#'                    omega = 2 * pi / 365, mu = 1 / (50 * 365))
#' 
#' initials <- c(S = 1 / 17, I = 1e-4, 
#'               R = 1 - 1 / 17 - 1e-4)
#' 
#' # Solve and plot.
#' sir.sinusoidal.forcing <- SIRSinusoidalForcing(pars = parameters, 
#'                                                init = initials, 
#'                                                time = 0:(60 * 365))
#' PlotMods(sir.sinusoidal.forcing)
#'
#' # Solve bifurcation dynamics for 20 years.
#' # If max(time) < 3650, bifurcation dynamics are solved for 3650 time-steps.
#' parameters2 <- list(beta0 = 17 / 13, beta1 = seq(0.001, 0.251, by = 0.001),
#'                    gamma = 1 / 13, omega = 2 * pi / 365, mu = 1 / (50 * 365))
#' # Uncomment the following lines (running it takes more than a few seconds):
#' # bifur <- SIRSinusoidalForcing(pars = parameters2, 
#' #                               init = initials,
#' #                               time = 0:(20 * 365))
#' # PlotMods(bifur, bifur = TRUE)
#' 
SIRSinusoidalForcing <- function(pars = NULL, init = NULL, time = NULL, ...) {
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
        beta <- beta0 * (1 + beta1 * sin(omega * time))
        dS <- mu - beta * S * I - mu * S
        dI <- beta * S * I - mu * I - gamma * I
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
  
  if (length(pars$beta1) == 1) {
    return(list(model = function1,
                pars = pars,
                init = init,
                time = time,
                results = as.data.frame(function1(pars = pars,
                                                  init = init, time = time))))
  } else {
    end <- max(time)
    beta2 <- pars$beta1
    if (end < 3650) {end <- 3650}
    bifur <- matrix(rep(NA, length(beta2) * 10), ncol = 10)
    for (i in 1:length(beta2)) {
      pars$beta1 <- beta2[i]
      output <- function1(pars = pars, init = init, time = time)
      init <- output[nrow(output), -1]
      for (j in 0:9) {
        bifur[i, j + 1] <- output[end - j * 365, 'I']
      }
    }
    return(data.frame(beta1 = beta2, bifur))
  }
}