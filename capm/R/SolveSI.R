#' Modelling of sterilization and immigration of comapnion animals.
#' @description System of ordinary differential equations to simulate the effect of sterilization and immigration on population dynamics.
#' @param pars \code{\link{vector}} of length 4. The values are point estimates of birth rate, death rate, carrying capacity and sterilization rate. The names of this values must be "b", "d", "k" and "s", respectively.
#' @param init \code{\link{vector}} of length 2. The values are initial population size and initial proportion of sterilized animals. The names of this values must be "n" and "q", respectively.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param dd string equal to \code{b} or \code{d} to define if density-dependece act on birth or death rartes respectively.
#' @param im a number representing the total of immigrants per time unit.
#' @param s.range optional sequence (between 0 and 1) of the sterilization rates to be simulated.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details The implemented model is described by Amaku, et. al., 2009 and the function is a wrapper around the defaults of \link[deSolve]{ode} function, whose help page must be consulted for details.
#' @return \code{\link{list}}. The first element, \code{name}, is a string with the name of the function, the second element, \code{model}, is the model function. The third, fourth and fifth elements are vectors (\code{pars}, \code{init}, \code{time}, respectively) containing the \code{pars}, \code{init} and \code{time} arguments of the function. The sisxth element \code{results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time, second column the population size and third column the proportion of sterilized animals. If \code{s.range} is specified, fourth column contains its specific instances.
#' @note Logistic growth models are not intended for scenarios in which population size is greater than carrying capacity and growth rate is negative.
#' @references Amaku M, Dias R and Ferreira F (2009). Dinamica populacional canina: potenciais efeitos de campanhas de esterilizacao. Revista Panamericana de Salud Publica, 25(4), pp. 300-304.
#' 
#' Soetaert K, Cash J and Mazzia F (2012). Solving differential equations in R. Springer.
#' 
#' \url{http://oswaldosantos.github.io/capm}
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions from estimates   
#' # obtained in examples section from svysumm function but
#' # estimating a proportion insted of a total for births.
#' pars.solve.si = c(b = 0.245, d = 0.101, 
#'                  k = 98050, s = 0.048)
#' init.solve.si = c(n = 89137, q = 0.198)
#' 
#' # Solve for a specific sterilization rate.
#' solvesi.pt = SolveSI(pars = pars.solve.si, 
#'                      init = init.solve.si, 
#'                      time = 0:15, dd = 'b',
#'                      im = 100, method = 'rk4')
#' 
#' # Solve for a range of sterilization rates.
#' solvesi.rg = SolveSI(pars = pars.solve.si,
#'                      init = init.solve.si,
#'                      time = 0:15, dd = 'b', im = 100, 
#'                      s.range = seq(0, .4, l = 50),
#'                      method = 'rk4')
#' 
SolveSI <- function(pars = NULL, init = NULL, time = NULL, dd = 'b', im = 0, s.range = NULL, ...) {
  SolveSIfu <- function(pars = NULL, init = NULL, time = NULL) {
    SolveSI.fu <- function(time, init, pars) {
      with(as.list(c(init, pars)), {
        if (dd == 'b') {
          nat = b - (b - d) * (n / k)
          mor = d
        }
        if (dd == 'd') {
          nat = b
          mor = d + (b - d) * (n / k)
        }
        dn = n * (nat * (1 - q) - mor) + im
        dq = (1 - q) * (s - q * nat)
        list(c(dn, dq))
      })
    }
    init <- c(init['n'], init['q'])
    SolveSI.out <- ode(times = time, func = SolveSI.fu, 
                           y = init, parms = pars, ...)
    return(as.data.frame(SolveSI.out))
  }
  
  if (!is.null(s.range)) {
    output <- NULL
    paras <- pars
    for(i in 1:length(s.range)) {
      paras['s'] = s.range[i]
      tmp = SolveSIfu(pars = paras, init = init, time = time)
      output = rbind(output, tmp)
    }   
    s.rate <- rep(s.range , each = length(time))   
    SolveSI <- list(
      name = 'SolveSI',
      model = SolveSIfu,
      pars = pars,
      init = init,
      time = time,
      results = as.data.frame(cbind(output, s.rate))
    )
    class(SolveSI) <- 'capmModels'
    return(SolveSI)
    
  } else {
    output <- SolveSIfu(pars = pars, init = init, time = time)
    SolveSI <- list(
      name = 'SolveSI',
      model = SolveSIfu,
      pars = pars,
      init = init,
      time = time,
      results = as.data.frame(output)
    )
    class(SolveSI) <- 'capmModels'
    return(SolveSI)
  }
}