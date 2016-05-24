#' SIS model with demographic stochasticity (P 6.3).
#' @description SIS model with event-driven or demographic stochasticity.
#' @param pars \code{\link{vector}} with 3 values: the transmission and recovery rates and the population size assumed to be constant. The names of these values must be "beta", "gamma", and "N" respectively.
#' @param init initial number of infectious.
#' @param end.time end time to be simulated.
#' @details This is the R version of program 6.3 from page 202 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' @return \code{\link{list}}. The first three elements are the vectors \code{*$pars}, \code{*$init} and \code{*$time}, containing the \code{pars}, \code{init} and \code{end.time} arguments of the function. The fourth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as time steps created during the stochastic simulations. The second column contains the number of infectious.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals}
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- c(beta = 0.03, gamma = 1 / 100, N = 100)
#' initials <- 70
#' 
#' # Solve and plot.
#' sis.demog.stoch <- SISDemogStoch(pars = parameters,
#'                                  init = initials, end.time = 10 * 365)
#' PlotMods(sis.demog.stoch)
#' 
SISDemogStoch <- function(pars, init, end.time) {
  init2 <- NULL
  init2 <- init
  equations <- function(pars, init, time.step) {
    with(as.list(pars), {
      rate1 <- beta * (N - init) * init / N
      rate2 <- gamma * init
      r1 <- runif(1)
      r2 <- runif(1)
      step <- -log(r2) / (rate1 + rate2)
      if (r1 < (rate1 / (rate1 + rate2))) {
        init <- init + 1
      } else {
        init <- init - 1
      }
      return(c(step, init))
    })
  }  
  lop <- 1
  output <- time.step <- 0
  output.tmp <- c(0, 0)
  while (time.step[lop] < end.time & init > 0) {
    output.tmp <- equations(pars, init, output.tmp[1])
    lop <- lop + 1
    time.step <- c(time.step, (time.step[lop - 1] + output.tmp[1]))
    output <- c(output, init)
    lop <- lop + 1
    time.step <- c(time.step, (time.step[lop - 1]))
    output <- c(output, output.tmp[2])
    init <- output.tmp[2]
  }
  return(list(pars = pars,
              init = init2,
              time = end.time,
              results = data.frame(time = time.step, Y = output)))
}