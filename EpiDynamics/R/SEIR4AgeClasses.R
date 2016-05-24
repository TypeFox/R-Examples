#' SEIR model with 4 age classes and yearly aging (P 3.4).
#' @description Solves a SEIR model with four different age-groups and yearly "movements" between the groups mimicking the school year
#' @param pars \code{\link{list}} with: a matrix with the transmission rates, the rate at which individuals move from the exposed to the infectious classes, the recovery rate, a \code{\link{vector}} with death rates in each age group, and a \code{\link{vector}} with birth rates into the childhood class. The names of these elements must be "beta", "sigma", "gamma", "mu", and "nu", respectively, see example.
#' @param init \code{\link{vector}} with 16 values: initial proportions of the population that are susceptible, exposed, infectious and recovered in a particular age group. The vector must be named, see example. Requirements: S + E + I <= n  for each age group and all values must be positive.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 3.4 from page 87 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' 
#' All rates are specified in days. Moreover, a \code{\link{vector}} n with the proportion of each age group. All parameters must be positive.
#' @return \code{\link{list}} of class \code{SolveSIR4ACYA}. The first element, \code{*$model}, is the model function. The second, third and fourth elements are vectors (\code{*$pars}, \code{*$init}, \code{*$time}, respectively) containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. The following columns contain the proportion of susceptibles, exposed, infectious and recovered.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals} 
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- list(beta = matrix(c(2.089, 2.089, 2.086, 2.037,
#'                                    2.089, 9.336, 2.086, 2.037,
#'                                    2.086, 2.086, 2.086, 2.037,
#'                                    2.037, 2.037, 2.037, 2.037),
#'                                  nrow = 4, ncol = 4),
#'                     sigma = 0.125, gamma = 0.2,
#'                     mu = c(0, 0, 0, 1) / (55 * 365),
#'                     nu = c(1 / (55 * 365), 0, 0, 0),
#'                     n = c(6, 4, 10, 55) / 75)
#'                           
#' initials <- c(S = c(0.05, 0.01, 0.01, 0.008), 
#'               E = c(0.0001, 0.0001, 0.0001, 0.0001),
#'               I = c(0.0001, 0.0001, 0.0001, 0.0001),
#'               R = c(0.0298, 0.04313333, 0.12313333, 0.72513333))
#' 
#' # Solve and plot.
#' # Uncomment the following lines (running it takes more than a few seconds):
#' # seir4.age.classes <- SEIR4AgeClasses(pars = parameters, 
#' #                                      init = initials,
#' #                                      time = 0:36500)
#' # PlotMods(seir4.age.classes,
#' #          variables = c('I1', 'I2', 'I3', 'I4'), grid = F)
#' 
SEIR4AgeClasses <- function(pars = NULL, init = NULL, time = NULL, ...) {
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
        
        dif.eq <- vector(length = 16)
        Ii <- c(eval(parse(text = 'I1')), eval(parse(text = 'I2')), eval(parse(text = 'I3')), eval(parse(text = 'I4')))
        for (i in 1:4){
          
          S <- eval(parse(text = paste('S', i, sep='')))
          E <- eval(parse(text = paste('E', i, sep='')))
          I <- eval(parse(text = paste('I', i, sep='')))
          R <- eval(parse(text = paste('R', i, sep='')))
          
          infection <- beta[i,] %*% t(t(Ii)) * S
          dif.eq[i] <- nu[i] * n[4] - infection - mu[i] * S #dS/dt
          dif.eq[4 + i] <- infection - sigma * E - mu[i] * E #dE/dt
          dif.eq[2 * 4 + i] <- sigma * E - gamma * I - mu[i] * I #dI/dt
          dif.eq[3 * 4 + i] <- gamma * I - mu[i] * R #dR/dt
          
        }
        list(dif.eq)
      })
    }
    
    output <- ode(times = time, 
                  func = function2, 
                  y = init, parms = pars, ...) 
    return(output)
  }
  
  T0 <- min(time);
  MaxTime <- max(time);
  output <- cbind.data.frame(time = 0, 
                             data.frame(t(vector(16, mode='numeric'))))
  names(output)[2:17] <- names(init)
  
  while (T0 < MaxTime) {
    output <- rbind.data.frame(output,
                               function1(pars = pars, init = init,
                                         time = seq(T0, T0 + 365,
                                                    time[2] - time[1])))
    
    init['S1'] <- output[T0 + 365, 'S1'] - output[T0 + 365, 'S1'] /6
    init['S2'] <- output[T0 + 365, 'S2'] + output[T0 + 365, 'S1'] / 6 -
      output[T0 + 365, 'S2'] / 4
    init['S3'] <- output[T0 + 365, 'S3'] + output[T0 + 365, 'S2'] / 4 -
      output[T0 + 365, 'S3'] / 10
    init['S4'] <- output[T0 + 365, 'S4'] + output[T0 + 365, 'S3'] / 10
    
    init['E1'] <- output[T0 + 365, 'E1'] - output[T0 + 365, 'E1'] / 6
    init['E2'] <- output[T0 + 365, 'E2'] + output[T0 + 365, 'E1'] / 6 -
      output[T0 + 365, 'E2'] / 4
    init['E3'] <- output[T0 + 365, 'E3'] + output[T0 + 365, 'E2'] / 4 -
      output[T0 + 365, 'E3'] / 10
    init['E4'] <- output[T0 + 365, 'E4'] + output[T0 + 365, 'E3'] / 10
    
    init['I1'] <- output[T0 + 365, 'I1'] - output[T0 + 365, 'I1'] / 6
    init['I2'] <- output[T0 + 365, 'I2'] + output[T0 + 365, 'I1'] / 6 -
      output[T0 + 365, 'I2'] / 4
    init['I3'] <- output[T0 + 365, 'I3'] + output[T0 + 365, 'I2'] / 4 -
      output[T0 + 365, 'I3'] / 10
    init['I4'] <- output[T0 + 365, 'I4'] + output[T0 + 365, 'I3'] / 10
    
    init['R1'] <- output[T0 + 365, 'R1'] - output[T0 + 365, 'R1'] / 6
    init['R2'] <- output[T0 + 365, 'R2'] + output[T0 + 365, 'R1'] / 6 -
      output[T0 + 365, 'R2'] / 4
    init['R3'] <- output[T0 + 365, 'R3'] + output[T0 + 365, 'R2'] / 4 -
      output[T0 + 365, 'R3'] / 10
    init['R4'] <- output[T0 + 365, 'R4'] + output[T0 + 365, 'R3'] / 10
    
    T0 <- T0 + 365 + time[2] -time[1]
    
  }
  
  return(list(model = function1,
              pars = pars,
              init = init,
              time = time,
              results = output[-1,]))
}