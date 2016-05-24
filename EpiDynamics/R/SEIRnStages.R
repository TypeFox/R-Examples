#' SEIR model with n stages (P 3.5).
#' @description Solves a SEIR model with multiple stages to create gamma-distributed exposed and infectious periods.
#' @param pars \code{\link{vector}} with 5 values: the transmission rate, the removal or recovery rate, the death rate (we assume that nu=mu), the number of stages in the infected period and the number of stages in the exposed period. The names of these elements must be "beta", "gamma", "mu", "n" and "m", respectively, see example. All rates are specified in days and all rates and parameters must be positive, moreover, m < n.
#' @param init \code{\link{vector}} with n + 1 values: initial proportions of the population that are susceptible and infected. The vector must be named, see example. Requirements: S +  all Infected <= 1.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 3.5 from page 94 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second, third and fourth elements are vectors (\code{*$pars}, \code{*$init}, \code{*$time}, respectively) containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. The following columns contain the proportion of susceptibles and infected.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals} 
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' n <- 13
#' parameters <- list(beta = 17 / 5, gamma = 1 / 13, mu = 1 / (55 * 365),
#'                    n = n, m = 8)
#'
#' initials <- c(S = 0.05, I = 0.00001 * rep(1, n) / n)
#' 
#' # Solve and plot.
#' # Uncomment the following lines (running it takes more than a few seconds):
#' # seir.n.stages <- SEIRnStages(pars = parameters, 
#' #                              init = initials, 
#' #                              time = seq(1, 30 * 365, 1))
#' # PlotMods(seir.n.stages, variables = 2)
#' # PlotMods(seir.n.stages, variables = 3:13, grid = F)


SEIRnStages <- function(pars = NULL, init = NULL, time = NULL, ...) {
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
        
        dif.eq <- vector(length = n + 1)
        Ii <- c(eval(parse(text = paste('I', m + 1, sep=''))))
        
        for(g in (m + 2):n){
          Ii <- c(Ii, eval(parse(text = paste('I', g, sep=''))))
        }
        
        infection <- beta * S * sum(Ii)
        dif.eq[1] <- mu - infection - mu * S #dS/dt
        dif.eq[2] <- infection - gamma * n * I1 - mu * I1 #dI1/dt
        
        for (i in 3:(n + 1)){
          
          Ib <- eval(parse(text = paste('I', i - 2, sep='')))
          I <- eval(parse(text = paste('I', i - 1, sep='')))
          
          dif.eq[i] <- gamma * n * Ib - gamma * n * I - mu * I #dI/dt
          
        }
        list(dif.eq)
      })
    }
    
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