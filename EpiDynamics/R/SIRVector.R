#' SIR model for mosquito vectors (P 4.4).
#' @description Solves a simple SIR model for mosquito vectors.
#' @param pars \code{\link{vector}} with 8 values: the human mortality rate (muH), the mosquito mortality rate (muM), the human birth rate (vH), the mosquito birth rate (vM), the transmission probability from humans to mosquitos following a bite (betaHM), the transmission probability from mosquitos to humans, following a bite (betaMH), the human recovery rate (gamma) and the rate at which humans are bitten (r). Abbreviations in parenthesis indicate the names which must be given to the values.
#' @param init \code{\link{vector}} with 4 values: the initial numbers of susceptible humans, susceptible mosquitos, infectious humans and infectious mosquitos. The names of this values must be "XH", "XM", "YH" and "YM", respectively.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 4.2 from page 123 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second, third and fourth elements are vectors (\code{*$pars}, \code{*$init}, \code{*$time}, respectively) containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. Second, third and fourth columns contain the proportion of susceptibles, infectious and recovered.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals} 
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- c(muH = 5.5e-5, muM = 0.143,
#'                 vH = 5.5e-2, vM = 1.443e3, 
#'                 betaHM = 0.5, betaMH = 0.8, 
#'                 gamma = 0.033, r = 0.5 / 1e3)
#' 
#' initials <- c(XH = 1e3, XM = 1e4, YH = 1, YM = 1)
#' 
#' # Solve and plot.
#' sir.vector <- SIRVector(pars = parameters,
#'                         init = initials,
#'                         time = 0:1000)                                 
#' PlotMods(sir.vector)
#' 
SIRVector <- function(pars = NULL, init = NULL, time = NULL, ...) {
  if (is.null(pars)) {
    stop("undefined 'pars'")
  }
  if (is.null(pars)) {
    stop("undefined 'inits'")
  }
  if (is.null(pars)) {
    stop("undefined 'time'")
  }
  SolveSIRMosVecfu <- function(pars = NULL, init = NULL, time = NULL) {
    SolveSIRMosVec.fu <- function(time, init, pars) {
      with(as.list(c(init, pars)), {
        dXH = vH - XH * r * (betaHM * YM) - muH * XH
        dXM = vM - XM * r * (betaMH * YH) - muM * XM
        dYH = XH * r * (betaHM * YM) - gamma * YH - muH * YH
        dYM = XM * r * (betaMH * YH) - muM * YM
        list(c(dXH, dXM, dYH, dYM))
      })
    }
    init <- c(init['XH'], init['XM'], init['YH'], init['YM'])
    output <- ode(times = time, 
                  func = SolveSIRMosVec.fu, 
                  y = init, parms = pars, ...)
    return(output)
  }
  
  output <- SolveSIRMosVecfu(pars = pars, init = init, time = time)
  return(list(model = SolveSIRMosVecfu,
              pars = pars,
              init = init,
              time = time,
              results = as.data.frame(output)))
}
