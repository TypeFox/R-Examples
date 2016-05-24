#' SIR model with 2 age classes (P 3.3).
#' @description Solves a SIR model two different age-groups.
#' @param pars \code{\link{vector}} with 8 values: 4 transmission rates, 1 recovery rate, rate at which children mature, death rate in the childhood group and death rate in the adult group. The names of these values must be "betaCC","betaCA","betaAC", "betaAA", "gamma", "lC","muC" and "muA", respectively. The letters after the word "beta" denote transmission to any group from any group, e.g., "betaCA" represent transmission to children group from adult group. All parameters must be positive. Parameters "nC" na "nu" (proportion of the population that are in the childhood group and birth rate into the childhood class, respectively) are not defined explicitly, but calculated as:  nC = muA/(lC+muA) and nu = (lC+nuA)nC. All rates are specified in years and all parameters must be positive
#' @param init \code{\link{vector}} with 4 values: the initial proportion of the population that are both susceptible and in the childhood group, the initial proportion of the population that are both infectious and in the childhood group, the initial proportion of the population that are both susceptible and in the adult group, and the initial proportion of the population that are both infectious and in the adult group. The names of these values must be "SC", "IC", "SA" and "IA", respectively. Requirements: SC + IC <= nC, and SA + IA <= nA = 1 - nC.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 3.3 from page 79 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second, third and fourth elements are the vectors (\code{*$pars}, \code{*$init}, \code{*$time}, containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. The following columns contain the proportion of susceptibles, infectious and recovered.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals} 
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- c(betaCC = 100, betaCA = 10, betaAC = 10, betaAA = 20,
#'                 gamma = 10, lC = 0.0666667, muC = 0.0, muA = 0.016667)
#' initials <- c(SC = 0.1, IC = 0.0001, SA = 0.1, IA = 0.0001)
#' 
#' # Solve and plot.
#' sir2AgeClasses <- sir2AgeClasses(pars = parameters, 
#'                           init = initials, time = seq(0, 100, 0.01))
#' PlotMods(sir2AgeClasses, variables = c('IA', 'IC'), grid = FALSE)
#' 
sir2AgeClasses <- function(pars = NULL, init = NULL, time = NULL, ...) {
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
        
        dSC = nu - SC * (betaCC * IC + betaCA * IA) - muC * SC - lC * SC
        dIC = SC * (betaCC * IC + betaCA * IA) - gamma * IC - muC * IC - lC * IC
        dSA = lC * SC - SA * (betaAC * IC + betaAA * IA) - muA * SA
        dIA = lC * IC + SA * (betaAC * IC + betaAA * IA) - gamma * IA - muA * IA
        
        list(c(dSC, dIC, dSA, dIA))
      })
    }
    
    init <- c(init['SC'], init['IC'], init['SA'], init['IA'])
    pars <- c(pars['betaCC'], pars['betaCA'], pars['betaAC'],
              pars['betaAA'], pars['gamma'], pars['lC'], pars['muC'],
              pars['muA'], nu = (pars[['lC']] + pars[['muA']]) * 
                pars[['muA']] / (pars[['lC']] + pars[['muA']]),
              nC = pars[['muA']]/ (pars[['lC']] + pars[['muA']]))
    
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