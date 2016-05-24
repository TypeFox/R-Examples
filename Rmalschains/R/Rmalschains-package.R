#' This package implements an algorithm family for continuous
#' optimization called memetic algorithms with local search chains
#' (MA-LS-Chains). 
#' 
#' One of the main issues to optimize a real-coded function is the capability of the algorithm to
#' realize a good exploration of the search space and, at the same time, to exploit the most
#' promising region to obtain accurate solutions.
#' 
#' Memetic algorithms are hybridizations of genetic algorithms with local search methods. 
#' They are especially suited for continuous optimization, as they combine the power of evolutionary
#' algorithms to explore the search space with a local search method to find the local optimum of a 
#' promising region. In these algorithms, it is recommended to increase the effort invested in the 
#' local search (measured in number of evaluations, called intensity) in the improvement of the 
#' most promising solution. However, it is not easy to decide the right intensity for each solution. 
#' 
#' MA-LS-Chains is a steady-state memetic algorithm, which combines a steady-state genetic algorithm with various
#' different local search methods. In contrast to the generational approach, where all 
#' individuals are substituted in an iteration, in the steady-state genetic algorithm 
#' in each iteration only one solution, the worst one, is subtituted in the population.
#' This makes it possible to not lose the improvement obtained by the local search over the individuals.
#' 
#' For MA-LS-Chains, the current state of the local search algorithm is stored along with the individuals. So, it becomes
#' possible to run the local search a fixed number of iterations, stop it, and possibly later continue the previous 
#' local search over the same individual. In this way, MA-LS-Chains
#' controls the application of the local search to the most promising solutions.
#'
#' The package implements various different local search strategies:
#' 
#' \itemize{
#' \item CMA-ES The Covariance Matrix Adaptation Evolution Strategy
#' \item SW A Solis Wets solver
#' \item SSW Subgrouping Solis Wets 
#' \item Simplex
#' }  
#' 
#' CMA-ES is a very effective local search strategy, but quite complicated, and it does not scale well if the amount of 
#' parameters to optimize is huge. The Solis Wets solver is pretty simple and therewith fast. The SSW strategy is an
#' adapted version of the Solis Wets solver for high dimensional data, so that the algorithm with this type of local search
#' scales well with the dimensionality of the data. It applies the Solis Wets solver to randomly chosen subgroups of 
#' variables (Subgrouping Solis Wets).
#' 
#' All the local search methods can also be used directly, without making use of the evolutionary algorithm. 
#' 
#' The package contains some demos illustrating its use. To get a list of them, type:
#' 
#' \code{library(Rmalschains)}
#' 
#' \code{demo()}
#' 
#' The demos currently available are \code{claw}, \code{rastrigin}, \code{sphere}, \code{rastrigin_highDim}, and \code{rastrigin_inline}. So in order to, 
#' e.g., execute the \code{claw} demo, type
#' 
#' \code{demo(claw)}
#'
#' All algorithms are implemented in C++, and they run pretty fast. A usual processing to speed up optimization is to implement the objective function also in C/C++.
#' However, a bottleneck in this approach is that the function needs to be passed as an R function, so that the optimizer needs to go back from C++ to R to C/C++ in each call 
#' of the target function. The package provides an interface which allows to pass the C/C++ target function directly as a pointer. See the \code{rastrigin_inline} 
#' demo for how to do that. The demo also shows how an environment can in this approach be used to pass additional parameters to the target function. 
#'  
#' For theoretical background of the algorithm, the reader may refer to the cited literature, where the algorithms where originally proposed.
#' 
#' Additional information is also available at the project website: 
#' 
#' \url{http://sci2s.ugr.es/dicits/software/Rmalschains}
#' 
#' @title Getting started with the Rmalschains package
#' @name Rmalschains-package
#' @aliases Rmalschains
#' @docType package
#' @title Getting started with the Rmalschains package
# @encoding UTF-8
# @encoding Latin-1
#' @author Christoph Bergmeir \email{c.bergmeir@@decsai.ugr.es} 
#'
#' Daniel Molina \email{dmolina@@decsai.ugr.es}
#' 
#' José M. Benítez \email{j.m.benitez@@decsai.ugr.es}
#'
#' DiCITS Lab, Sci2s group, DECSAI, University of Granada. \url{http://dicits.ugr.es}.
#' 
#' Additional information is also available at our group's website of continuous optimization:
#' 
#' \url{http://sci2s.ugr.es/EAMHCO/}
#' 
#' @references 
#' 
#' Molina, D., Lozano, M., Sánchez, A.M., Herrera, F.
#' Memetic algorithms based on local search chains for large scale continuous optimisation problems: MA-SSW-Chains
#' (2011) Soft Computing, 15 (11), pp. 2201-2220. 
#'
#' Molina, D., Lozano, M., Herrera, F.
#' MA-SW-Chains: Memetic algorithm based on local search chains for large scale continuous global optimization
#' (2010) 2010 IEEE World Congress on Computational Intelligence, WCCI 2010 - 2010 IEEE Congress on Evolutionary Computation, CEC 2010. 
#'
#' Molina, D., Lozano, M., García-Martínez, C., Herrera, F.
#' Memetic algorithms for continuous optimisation based on local search chains
#' (2010) Evolutionary Computation, 18 (1), pp. 27-63.
#' 
#' @keywords optimization, MA-LS-Chains
#' @seealso \code{\link{malschains}}, \code{\link{malschains.control}}
#' @useDynLib Rmalschains .registration=TRUE
#' @import Rcpp
# @exportPattern "^[[:alpha:]]+"
#' @examples
#' 
#' ##############################################
#' #Example for maximization of the claw function
#' ##############################################
#' 
#' claw <- function(xx) {
#'   x <- xx[1]
#'   y <- (0.46 * (dnorm(x, -1, 2/3) + dnorm(x, 1, 2/3)) +
#'         (1/300) * (dnorm(x, -0.5, 0.01) + dnorm(x, -1,
#'               0.01) + dnorm(x, -1.5, 0.01)) + (7/300) *
#'         (dnorm(x, 0.5, 0.07) + dnorm(x, 1, 0.07) + dnorm(x,
#'               1.5, 0.07)))
#'   return(y)
#' }
#' 
#' #use MA-CMA-Chains
#' res.claw <- malschains(function(x) {-claw(x)}, lower=c(-3), upper=c(3), 
#'                        maxEvals=50000, control=malschains.control(popsize=50, 
#'                        istep=300, ls="cmaes", optimum=-5))
#' 
#' \dontrun{
#' #use only the CMA-ES local search               
#' res.claw2 <- malschains(function(x) {-claw(x)}, lower=c(-3), upper=c(3), verbosity=0,
#'                        maxEvals=50000, control=malschains.control(ls="cmaes", 
#'                            lsOnly=TRUE, optimum=-5))
#' 
#' #use only the Simplex local search               
#' res.claw3 <- malschains(function(x) {-claw(x)}, lower=c(-3), upper=c(3), verbosity=0,
#'                        maxEvals=50000, control=malschains.control(ls="simplex", 
#'                            lsOnly=TRUE, optimum=-5))
#'                     
#' x <- seq(-3, 3,length=1000)
#' claw_x <- NULL
#' for (i in 1:length(x)) claw_x[i] <- claw(x[i])
#' 
#' plot(x,claw_x, type="l")
#' points(res.claw$sol, -res.claw$fitness, col="red")
#' points(res.claw2$sol, pch=3, -res.claw2$fitness, col="blue")
#' points(res.claw3$sol, pch=3, -res.claw3$fitness, col="green")
#' 
#' 
#' ##############################################
#' #Example for the rastrigin function
#' ##############################################
#'  
#' rastrigin <- function(x) {
#'   
#'   dimension <- length(x)
#'   
#'   res <- 0.0
#'   for (i in 1:dimension) {
#'     res <- res + (x[i]*x[i] - 10.0*cos(2.0*pi*x[i]) + 10.0)
#'   }
#' 
#'   res 
#' }
#' 
#' res.rastrigin1 <- malschains(rastrigin, lower=seq(-1.0, -1.0, length=30), 
#'                              upper=seq(1.0, 1.0, length=30), maxEvals=50000, 
#'                              control=malschains.control(effort=0.8, alpha=0.3, 
#'                              popsize=20, istep=100, ls="simplex"))
#' 
#' 
#' res.rastrigin2 <- malschains(rastrigin, lower=seq(-1.0, -1.0, length=30), 
#'                              upper=seq(1.0, 1.0, length=30), maxEvals=50000, 
#'                              initialpop = seq(0.1, 0.1, length=30), 
#'                              control=malschains.control(popsize=50, 
#'                              istep=300, ls="cmaes"))
#' 
#' res.rastrigin1
#' res.rastrigin2
#' }
NULL

