#' \pkg{spsann} annealing schedule
#' 
#' Set the control parameters for the annealing schedule of \pkg{spsann} functions.
#' 
#' @inheritParams spJitter
#' 
#' @param x.max,x.min,y.max,y.min Numeric value defining the minimum and maximum quantity of random noise to 
#' be added to the projected x- and y-coordinates. The units are the same as of the projected x- and 
#' y-coordinates. If missing, they are estimated from \code{candi}, \code{x.min} and \code{y.min} being set 
#' to zero, and \code{x.max} and \code{y.max} being set to half the maximum distance in the x- and 
#' y-coordinates, respectively.
#' 
#' @param initial.acceptance Numeric value between 0 and 1 defining the initial acceptance probability, i.e. 
#' the proportion of proposed system configurations that should be accepted in the first chain. The 
#' optimization is stopped and a warning is issued if this value is not attained. Defaults to 
#' \code{initial.acceptance = 0.95}.
#' 
#' @param initial.temperature Numeric value larger than 0 defining the initial temperature of the system. A 
#' low \code{initial.temperature}, combined with a low \code{initial.acceptance} result in the algorithm to 
#' behave as a greedy algorithm, i.e. only better system configurations are accepted. Defaults to 
#' \code{initial.temperature = 0.001}.
#' 
#' @param temperature.decrease Numeric value between 0 and 1 used as a multiplying factor to decrease the 
#' temperature at the end of each Markov chain. Defaults to \code{temperature.decrease = 0.95}.
#' 
#' @param chains Integer value defining the maximum number of chains, i.e. the number of cycles of jitters at 
#' which the temperature and the size of the neighbourhood should be kept constant. Defaults to 
#' \code{chains = 500}.
#' 
#' @param chain.length Integer value defining the length of each Markov chain relative to the number of 
#' sample points. Defaults to \code{chain.length = 1}, i.e. one time the number of sample points.
#' 
#' @param stopping Integer value defining the maximum allowable number of Markov chains without improvement of 
#' the objective function value. Defaults to \code{stopping = 10}.
#' 
#' @return
#' A list with a set of control parameters of the annealing schedule.
#' 
#' @references 
#' Aarts, E. H. L.; Korst, J. H. M. Boltzmann machines for travelling salesman problems. \emph{European 
#' Journal of Operational Research}, v. 39, p. 79-95, 1989.
#' 
#' Černý, V. Thermodynamical approach to the travelling salesman problem: an efficient simulation algorithm.
#' \emph{Journal of Optimization Theory and Applications}, v. 45, p. 41-51, 1985.
#' 
#' Brus, D. J.; Heuvelink, G. B. M. Optimization of sample patterns for universal kriging of environmental 
#' variables. \emph{Geoderma}, v. 138, p. 86-95, 2007.
#' 
#' Kirkpatrick, S.; Gelatt, C. D.; Vecchi, M. P. Optimization by simulated annealing. \emph{Science}, v. 220,
#' p. 671-680, 1983.
#' 
#' Metropolis, N.; Rosenbluth, A. W.; Rosenbluth, M. N.; Teller, A. H.; Teller, E. Equation of state 
#' calculations by fast computing machines. \emph{The Journal of Chemical Physics}, v. 21, p. 1087-1092, 1953.
#' 
#' van Groenigen, J.-W.; Stein, A. Constrained optimization of spatial sampling using continuous 
#' simulated annealing. \emph{Journal of Environmental Quality}. v. 27, p. 1078-1086, 1998.
#' 
#' Webster, R.; Lark, R. M. \emph{Field sampling for environmental science and management}. London: 
#' Routledge, p. 200, 2013.
#' 
#' @author 
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[spsann]{optimACDC}}, \code{\link[spsann]{optimCORR}},
#' \code{\link[spsann]{optimDIST}}, \code{\link[spsann]{optimMKV}},
#' \code{\link[spsann]{optimMSSD}}, \code{\link[spsann]{optimPPL}},
#' \code{\link[spsann]{optimSPAN}}, \code{\link[spsann]{optimUSER}}.
#' @export
#' @examples
#' schedule <- scheduleSPSANN()
# FUNCTION #####################################################################
scheduleSPSANN <-
  function (initial.acceptance = 0.95, initial.temperature = 0.001,
            temperature.decrease = 0.95, chains = 500, chain.length = 1,
            stopping = 10, x.max, x.min = 0, y.max, y.min = 0, cellsize) {
    
    if (initial.acceptance > 1 || initial.acceptance < 0)
      stop ("'initial.acceptance' must be between 0 and 1")
    if (initial.temperature <= 0)
      stop ("'initial.temperature' must be larger than 0")
    if (temperature.decrease >= 1 || temperature.decrease <= 0)
      stop ("'temperature.decrease' must be between 0 and 1")
    if (chains < 1 || !pedometrics::isNumint(chains))
      stop ("'chains' must be an integer larger than 0")
    if (chain.length < 1 || !pedometrics::isNumint(chain.length))
      stop ("'chain.length' must be an integer larger than 0")
    if (stopping < 1 || !pedometrics::isNumint(stopping))
      stop ("'stopping' must be an integer larger than 0")
    if (missing(x.max)) x.max <- NULL
    # if (missing(x.min)) x.min <- NULL
    if (missing(y.max)) y.max <- NULL
    # if (missing(y.min)) y.min <- NULL
    if (missing(cellsize)) {
      cellsize <- NULL
    } else {
     if (length(cellsize) == 1) {
       cellsize <- rep(cellsize, 2)
     }
    }
    
    # Output
    res <- list(
      initial.acceptance = initial.acceptance, initial.temperature = initial.temperature,
      temperature.decrease = temperature.decrease, chains = chains, chain.length = chain.length, 
      stopping = stopping, x.max = x.max, x.min = x.min, y.max = y.max, y.min = y.min, cellsize = cellsize)
    return (res)
  }
