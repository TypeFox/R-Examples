#' Partial immunity model that cycles (P 4.2).
#' @description Solves multi-strain where the strains are arranged in a circle and each strain offers partial immunity (in terms of reduced transmission) to its neighbours.
#' @param pars \code{\link{vector}} with 2 vectors and 2 values in the following order: beta, gamma, mu and a. The first element of each vector corresponds to the first strain, the second element to the second strain and so on. mu is the per capita death rate and alpha is the modified transmission rate due to partial immunity.
#' @param init \code{\link{vector}} with 3 vectors in the following order: S, P and L. The first element of each vector corresponds to the first strain, the second element to the second strain and so on. Si + Pi + Li = 1 but sum(Si) could be greater than 1.
#' @param time time sequence for which output is wanted; the first value of times must be the initial time.
#' @param ... further arguments passed to \link[deSolve]{ode} function.
#' @details This is the R version of program 4.2 from page 123 of "Modeling Infectious Disease in humans and animals" by Keeling & Rohani.
#' @return \code{\link{list}}. The first element, \code{*$model}, is the model function. The second, third and fourth elements are the vectors (\code{*$pars}, \code{*$init}, \code{*$time}, containing the \code{pars}, \code{init} and \code{time} arguments of the function. The fifth element \code{*$results} is a \code{\link{data.frame}} with up to as many rows as elements in time. First column contains the time. Second, third and fourth columns contain the proportion of susceptibles, infectious and recovered.
#' @references Keeling, Matt J., and Pejman Rohani. Modeling infectious diseases in humans and animals. Princeton University Press, 2008.
#' \href{http://www.modelinginfectiousdiseases.org/}{Modeling Infectious Diseases in Humans and Animals} 
#' @seealso \link[deSolve]{ode}.
#' @export
#' @examples 
#' # Parameters and initial conditions.
#' parameters <- c(beta = rep(40, 4), gamma = rep(9.98, 4),
#'                        mu = 0.02, a = 0.4 )
#' 
#' initials <- c(S = c(0.08, 0.1, 0.1, 0.11),
#'                        P = c(0.4, 0.3, 0.3, 0.29),
#'                        L = c(0.15, 0.02, 0.1, 0.01))
#' 
#' # Solve and plot.
#' mlti.strain.pi <- MultiStrainPartialImmunity(pars = parameters, 
#'                                              init = initials, 
#'                                              time = 0:200)
#' PlotMods(mlti.strain.pi, variables = c('L1', 'L2', 'L3', 'L4'), grid = FALSE)
#'                                  
MultiStrainPartialImmunity <- function(pars = NULL, init = NULL, time = NULL, ...) {
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
      V <- init
      P <- pars
      N <- length(V) / 3
      with(as.list(c(init, pars)), {
        Y <- rep(0, 3 * N)
        for (i in 1:N) {
          r <- (i %% N) + 1
          l <- ((i + N - 2) %% N) + 1
          Y[i] <- mu - V[i] * (V[(2 * N) + i] + V[(2 * N) + l] +
                                 V[(2 * N) + r]) - mu * V[i]
          Y[N + i] <- V[i] * (V[(2 * N) + l] + V[(2 * N) + r]) - 
            V[N+i] * V[(2 * N) + i] - mu * V[N + i]
          Y[(2 * N) + i] <- P[i] * (V[i] + a * V[N + i]) * 
            V[(2 * N) + i] - P[N + i] * V[(2 * N) + i] -
            mu * V[(2 * N) + i]
        }
        list(Y)
      })
    }
    tmp <- NULL
    for (i in 1:length(init)) {
      tmp[i] <- init[i]
    }
    init <- tmp    
    output <- ode(times = time, 
                  func = function2, 
                  y = init, parms = pars, ...)
    return(output)
  }
  
  output <- function1(pars = pars, init = init, time = time)
  colnames(output)[-1] <- names(init)
  return(list(model = function1,
              pars = pars,
              init = init,
              time = time,
              results = as.data.frame(output)))
}