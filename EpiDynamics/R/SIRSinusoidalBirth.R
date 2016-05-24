#' SIR model with sinusoidal births (P 5.3).
#' @description Solves a SIR model with sinusoidal forcing of the birth rate.
#' @param pars \code{\link{list}} with: the transmission rate, the mean birth rate, a scalar (or a \code{\link{vector}} to create bifurcations) with the amplitude of sinuoidal forcing for the birth rate, the frequency of the oscillations, the per capita death rate and the recovery rate. The names of these values must be "beta", "alpha0", "alpha1", "w", "mu"  and "gamma", respectively. All rates must be specified in days and alpha1 <= 1.
#' @param init \code{\link{vector}} with 3 values: the initial proportion of susceptibles and infectious. The names of these values must be "S" and "I", respectively. All parameters must be positive and S + I <= 1.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 5.3 from page 184 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani. To create bifurcations, \code{alpha1} must be a vector. For bifurcations, if max(time) < 3650), time is defined as c(0:3650).
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second element is a \code{\link{list}} with the the \code{*$pars} argument. The third and fourth elements are the vectors (\code{*$init}, \code{*$time}, containing the \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. Second, third and fourth columns contain the proportion of susceptibles, infectious and recovered.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals}. It is important to note that we wrote equations for the R population based on equations in the website and because of the dynamic S + I + R fluctuates around 1. Then, using R = 1 - S - I solves this inconsistency.
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#'                                   
#' # Parameters and initial conditions (bifurcation plot of infectious)
#' parameters <- list(beta = 17 / 13, alpha0 = 1 / (50 * 365),
#'                    alpha1 = 0.25, w = 2 * pi / 365 ,
#'                    mu = 1 / (50 * 365), gamma =  1 / 13)
#' 
#' parameters2 <- list(beta = 17 / 13, alpha0 = 1 / (50 * 365),
#'                    alpha1 = seq(0, 0.99, 0.01), w = 2 * pi / 365 ,
#'                    mu = 1 / (50 * 365), gamma =  1 / 13)
#' 
#' initials <- c(S = 1 / 17, I = 1e-4, R =  1 - (1 / 17 + 1e-4))
#' 
#' # Solve and plot.
#' sir.sinusoidal.birth <- SIRSinusoidalBirth(pars = parameters,
#'                                            init = initials, 
#'                                            time = 0:(20 * 365))
#' PlotMods(sir.sinusoidal.birth)
#' 
#' # Bifurcations
#' # Uncomment the following lines (running it takes more than a few seconds):
#' # bifurcation <- SIRSinusoidalBirth(pars = parameters2,
#' #                                   init = initials, 
#' #                                   time = 0:(20 * 365))
#' # PlotMods(bifur, bifur = TRUE)
SIRSinusoidalBirth <- function(pars = NULL, init = NULL, time = NULL, ...) {
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
        dS = alpha0 * (1 + alpha1 * sin(w * time)) - beta * S * I - mu * S
        dI = beta * S * I - gamma * I - mu * I
        dR = gamma * I - mu * R
        list(c(dS, dI, dR))
      })
    }
    init <- c(init['S'], init['I'], init['R'])
    output <- ode(times = time, 
                  func = function2, 
                  y = init, parms = pars, ...)
    return(output)
  }
  
  if (length(pars$alpha1) == 1){
    output <- function1(pars = pars, init = init, time = time)    
    return(list(model = function1,
                pars = pars,
                init = init,
                time = time,
                results = as.data.frame(output)))
  }
  
  else{
    if (max(time) < 3650){
      time <- c(0:3650)
    }
    
    alpha2 <- pars$alpha1
    alpha1.range <- length(pars$alpha1)
    bifur.infectious <- matrix(0,nrow=alpha1.range,ncol=10)
    
    for (i in 1:alpha1.range){      
      pars['alpha1'] <- alpha2[i]
      solve.result <- SIRSinusoidalBirth(pars = pars, init = init, time = time)
      init <- c(S = as.numeric(tail(solve.result$results, 1)[2]),
                I = as.numeric(tail(solve.result$results, 1)[3]), 
                R = as.numeric(tail(solve.result$results, 1)[4]))
      for (j in 0:9){
        bifur.infectious[i,j+1] <-
          solve.result$results[((max(time) - j * 365) - 1), 'I']
      }      
    }  
    
    return(data.frame(alpha1 = alpha2, bifur.infectious))
  }
  
}