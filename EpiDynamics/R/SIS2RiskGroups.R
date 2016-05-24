#' SIS model with 2 risk groups (P 3.1).
#' @description Solves a SIS model with high-risk (H) and low-risk (L).
#' @param pars \code{\link{vector}} with 6 values: 4 transmission rates, 1 recovery rate and  the proportion of the population that are in the high-risk group. The names of these values must be "betaHH","betaHL","betaLL", "betaLH", "gamma" and "nH", respectively. The letters after the word "beta" denote transmission to any group from any group, e.g., "betaHL" represent transmission to high-risk group from low-risk group. All parameters must be positive.
#' @param init \code{\link{vector}} with 2 values: the initial proportion of infectious in the high-risk group and the intial proportion of infectious in the low-risk group. The names of these values must be "IH" and "IL", respectively.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 3.1 from page 58 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' 
#' All parameters must be positive and  nH <= 1, IH(0) <= nH, IL(0) <= 1-nH.
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second, third and fourth elements are the vectors (\code{*$pars}, \code{*$init}, \code{*$time}, containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. The following columns contain the proportion of susceptibles and infectious.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals} 
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- c(betaHH = 10.0, betaHL = 0.1, betaLH = 0.1,
#'                 betaLL = 1.0, gamma = 1, nH = 0.2)
#' initials <- c(IH = 0.00001, IL = 0.001)
#' 
#' # Solve and plot.
#' sis2risk.groups<- SIS2RiskGroups(pars = parameters,
#'                                  init = initials, time = 0:15)
#' PlotMods(sis2risk.groups, variables = c('IL', 'IH'), grid = FALSE)
#' 
SIS2RiskGroups <- function(pars = NULL, init = NULL, time = NULL, ...) {
  function2 <- function(time, init, pars) {
    with(as.list(c(init, pars)), {
      dS = - beta * S * I
      dI = beta * S * I - gamma * I
      dR = gamma * I
      list(c(dS, dI, dR))
    })
  }
  function1 <- function(pars = NULL, init = NULL, time = NULL) {
    function2 <- function(time, init, pars) {
      with(as.list(c(init, pars)), {
        
        dSH = - betaHH * SH * IH - betaHL * SH * IL + gamma * IH
        dIH = betaHH * SH * IH + betaHL * SH * IL - gamma * IH
        dSL = - betaLL * SL * IL - betaLH * SL * IH + gamma * IL
        dIL = betaLL * SL * IL + betaLH * SL * IH - gamma * IL
        
        list(c(dSH, dIH, dSL, dIL))
      })
    }
    
    init <- c(SH = pars[['nH']] - init[['IH']] ,init['IH'],
              SL = 1 - pars[['nH']] - init[['IL']], init['IL'])
    
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