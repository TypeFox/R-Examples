#' SIR model with tau leap method (P 6.5).
#' @description SIR model with demographic stochasticity approximated using the tau-leap method and assuming Poisson distributions.
#' @param pars \code{\link{vector}} with 5 values: the transmission, recovery and death rates, the population size assumed to be constant and the time step. The names of these values must be "beta", "gamma", "mu", "N" and "tau" respectively.
#' @param init \code{\link{vector}} with 3 values: the initial number of susceptibles, infectious and recovered, respectively. The names of these values must be "X", "Y" and "Z" respectively.
#' @param end.time end time to be simulated.
#' @details This is the R version of program 6.5 from page 204 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' @return \code{\link{list}}. The first three elements are the vectors \code{*$pars}, \code{*$init} and \code{*$time}, containing the \code{pars}, \code{init} and \code{end.time} arguments of the function. The fourth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as time steps determined by the parameters \code{tau} and \code{end.time}. The first column contains the time steps. The second, third and fourth columns contain the number of susceptibles, infectious and recovered.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals}
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- c(beta = 1, gamma = 1 / 10, mu = 5e-4, N = 50, tau = 1)
#' initials <- c(X = 5, Y = 1, Z = 44)
#' end.time <- 2 * 365
#' 
#' # Solve and plot.
#' sir.demog.stoch <- SIRTauLeap(pars = parameters, init = initials,
#'                               end.time = end.time)
#' PlotMods(sir.demog.stoch)
#' 
SIRTauLeap <- function(pars, init, end.time) {
  init2 <- init
  Equations <- function(pars, init, end.time) {
    with(as.list(c(pars, init)), {
      rate <- rep(0, 6)
      change <- matrix(0, nrow = 6, ncol = 3)
      N <- X + Y + Z
      tau <- 1
      rate[1] <- beta * X * Y / N
      change[1, ] <- c(-1, 1, 0)
      rate[2] <- gamma * Y
      change[2, ] <- c(0, -1, 1)
      rate[3] <- mu * N
      change[3, ] <- c(1, 0, 0)
      rate[4] <- mu * X
      change[4, ] <- c(-1, 0, 0)
      rate[5] <- mu * Y
      change[5, ] <- c(0, -1, 0)
      rate[6] <- mu * Z
      change[6, ] <- c(0, 0, -1)
      init <- c(X = X, Y = Y, Z = Z)
      for (i in 1:6) {
        num <- rpois(1, rate[i] * tau)
        num.min <- min(num, init[which(change[i, ] < 0)])
        init <- init + change[i, ] * num.min
      }
      return(init)
    })
  }
  
  X <- Y <- Z <- double()
  t <- 0
  time <- seq(0, end.time, by = pars['tau'])
  for (t in time) {
    tmp <- Equations(pars, init, end.time)
    X <- c(X, init['X'])
    Y <- c(Y, init['Y'])
    Z <- c(Z, init['Z'])
    init <- tmp
  }
  return(list(pars = pars,
              init = init2,
              time = time,
              results = data.frame(time, X, Y, Z)))
}                                    