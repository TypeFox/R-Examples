#' SIR model with partial immunity (P 4.1).
#' @description Solves a model with partial immunitiy.
#' @param pars \code{\link{vector}} with 9 values: the death, birth, transmission and recovery rates, the modified susceptiblity to strain i for those individuals recovered from the other strain and the modified transmission rate of strain i from those individuals that have recovered from the other strain. The name of these values must be "mu", "v", "beta1", "beta2", "gamma1", "gamma2", "alpha1", "aplha2", "a1", "a2". The numbers 1 and 2 at the end of parameters names stand for strain 1 and strain 2.
#' @param init \code{\link{vector}} with 8 values: In this formulation NAB refers to the proportion of the population who are A with respect to strain 1 and B with respect to strain 2. Thus, initial values must be named: "NSS", "NIS", "NRS", "NRI", "NSI", "NSR", "NIR" and "NRR". The sum of initial values must be <= 1.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 4.1 from page 118 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second, third and fourth elements are the vectors (\code{*$pars}, \code{*$init}, \code{*$time}, containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. Second, third and fourth columns contain the proportion of susceptibles, infectious and recovered.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals} 
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- c(mu = 1 / (70 * 365), beta1 = 260 / 365, 
#'                 beta2 = 520 / 365, gamma1 = 1 / 7, 
#'                 gamma2 = 1 / 7, alpha1 = 0.5,
#'                 alpha2 = 0.4, a1 = 0.4, a2 = 0.5)
#' 
#' initials <- c(NSS = 0.1, NIS = 1e-4, NRS = 0.02, NRI = 0,
#'               NSI = 1e-4, NSR = 0.5, NIR = 0, NRR = 0.3798)
#' 
#' # Solve and plot.
#' sir.partial.immunity <- SIRPartialImmunity(pars = parameters, 
#'                                            init = initials,
#'                                            time = 0:(100 * 365))
#' PlotMods(sir.partial.immunity, variables = c('NIS', 'NIR'), grid = FALSE)
#'                                       
SIRPartialImmunity <- function(pars = NULL, init = NULL, time = NULL, ...) {
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
        lambda1 <- beta1 * (NIS + a1 * NIR)  
        lambda2 <- beta2 * (NSI + a2 * NRI)
        dNSS <- mu - NSS * (lambda1 + lambda2) - mu * NSS
        dNIS <- NSS * lambda1 - gamma1 * NIS - mu * NIS
        dNRS <- gamma1 * NIS - alpha2 * NRS * lambda2 - mu * NRS
        dNRI <- alpha2 * NRS * lambda2 - gamma2 * NRI - mu * NRI
        dNSI <- NSS * lambda2 - gamma2 * NSI - mu * NSI
        dNSR <- gamma2 * NSI - alpha1 * NSR * lambda1 - mu * NSR
        dNIR <- alpha1 * NSR * lambda1 - gamma1 * NIR - mu * NIR
        dNRR <- gamma1 * NIR + gamma2 *NRI - mu * NRR
        list(c(dNSS, dNIS, dNRS, dNRI, dNSI, dNSR, dNIR, dNRR))
      })
    }
    init <- c(init['NSS'], init['NIS'], init['NRS'], init['NRI'],
              init['NSI'], init['NSR'], init['NIR'], init['NRR'])
    
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