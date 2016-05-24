#' Experimental Differential Evolution - ExpDE
#'
#' Modular implementation of the Differential Evolution Algorithm for 
#' the experimental investigation of the effects of different operators 
#' on the performance of the algorithm.
#'
#' This routine is used to launch a differential evolution algorithm for the 
#' \strong{minimization} of a given problem instance using different variants of 
#' the recombination, mutation and selection operators. The input parameters 
#' that describe those operators receive list objects describing the operator 
#' variants to be used in a given optimization procedure.
#'
#' @section Mutation Parameters:
#' \code{mutpars} is used to inform the routine the type of differential 
#' mutation to use, as well as any mutation-related parameter values. The 
#' current version accepts the following options:
#' 
#' \itemize{
#'    \item \code{\link{mutation_rand}}
#'    \item \code{\link{mutation_best}}
#' }
#' 
#' \code{mutpars} receives a list object with name field \code{mutpars$name} 
#' (containing the name of the function to be called, e.g., 
#' \code{name = "mutation_rand"}) as well as whatever parameters that function 
#' may require/accept (e.g., \code{mutpars$f = 0.7}, \code{mutpars$nvecs = 2}, 
#' etc.). See the specific documentation of each function for details. 
#' 
#' Some examples are provided in the \code{Examples} section below.
#'
#' @section Recombination parameters:
#' As with the mutation parameters, \code{recpars} is used to define the 
#' desired recombination strategy. The current version accepts the following 
#' options:
#' 
#' \itemize{
#'    \item \code{\link{recombination_arith}}   
#'    \item \code{\link{recombination_bin}}
#'    \item \code{\link{recombination_blxAlpha}}
#'    \item \code{\link{recombination_blxAlphaBeta}}
#'    \item \code{\link{recombination_eigen}}
#'    \item \code{\link{recombination_exp}}
#'    \item \code{\link{recombination_flat}}
#'    \item \code{\link{recombination_geo}}
#'    \item \code{\link{recombination_lbga}}
#'    \item \code{\link{recombination_linear}}
#'    \item \code{\link{recombination_mmax}}
#'    \item \code{\link{recombination_npoint}}
#'    \item \code{\link{recombination_onepoint}}
#'    \item \code{\link{recombination_pbest}}
#'    \item \code{\link{recombination_sbx}}
#'    \item \code{\link{recombination_wright}}    
#' }
#' 
#' \code{recpars} receives a list object with name field \code{recpars$name} 
#' (containing the name of the function to be called, e.g., 
#' \code{name = "recombination_bin"}) as well as whatever parameters that 
#' function may require/accept (e.g., \code{recpars$cr = 0.8}, 
#' \code{recpars$minchange = TRUE}, etc.). See the specific documentation of 
#' each function for details. 
#' 
#' Some examples are provided in the \code{Examples} section below.
#' 
#' @section Selection parameters:
#' \code{selpars} follows the same idea as \code{mutpars} and \code{recpars}, 
#' and is used to define the selection operators. Currently, only the standard 
#' DE selection, \code{\link{selection_standard}}, is implemented.
#'
#' @section Stop criteria:
#' \code{stopcrit} is similar to \code{recpar} and the other list arguments, 
#' but with the difference that multiple stop criteria can be defined for the
#' algorithm. The names of the stop criteria to be used are passed in the 
#' \code{stopcrit$names} field, which must contain a character vector. Other 
#' parameters to be used for stopping the algorithm (e.g., the maximum number 
#' of iterations \code{stopcrit$maxiter}) can also be included as 
#' \code{stopcrit} fields. Currently implemented criteria are:
#' 
#' \itemize{
#'    \item \code{"stop_maxiter"} (requires additional field 
#'      \code{stopcrit$maxiter = ?} with the maximum number of iterations).
#'    \item \code{"stop_maxeval"} (requires additional field 
#'      \code{stopcrit$maxevals = ?} with the maximum number of function calls).
#'  }
#'  
#'  See \code{\link{check_stop_criteria}} for details.
#'  
#' @section Problem description:
#' The \code{probpars} parameter receives a list with all definitions related 
#' to the problem instance to be optimized. There are three required fields in 
#' this parameter:
#' 
#' \itemize{
#'    \item \code{probpars$name}, the name of the function that represents the 
#'      problem to be solved.
#'    \item \code{probpars$xmin}, a vector containing the lower bounds of all 
#'      optimization variables (i.e., a vector of length M, where M is the 
#'      dimension of the problem).
#'    \item \code{probpars$xmax}, a vector containing the upper bounds of all 
#'      optimization variables.
#' }
#' 
#' \strong{Important}: the objective function routine must receive a matrix of 
#' row vectors to be evaluated in the form of an input parameter named either 
#' "x" or "X" or "Pop" (any one of the three is allowed).
#' 
#' @section Random Seed:
#' The \code{seed} argument receives the desired seed for the PRNG. This value 
#' can be set for reproducibility purposes. The value of this parameter defaults 
#' to NULL, in which case the seed is arbitrarily set using 
#' \code{as.numeric(Sys.time())}.
#'
#' @section Showpars:
#' \code{showpars} is a list containing parameters that control the printed
#' output of \code{ExpDE}. Parameter \code{showpars} can have the following 
#' fields:
#' \itemize{
#'    \item \code{showpars$show.iters = c("dots", "numbers", "none")}: type of 
#'      output. Defaults to \code{"numbers"}.
#'    \item \code{showpars$showevery}: positive integer that determines how 
#'      frequently the routine echoes something to the terminal. Defaults 
#'      to \code{1}.
#'  }
#' 
#' @param popsize population size
#' @param mutpars list of named mutation parameters.
#'    See \code{Mutation parameters} for details.
#' @param recpars list of named recombination parameters.
#'    See \code{Recombination parameters} for details.
#' @param selpars list of named selection parameters.
#'    See \code{Selection parameters} for details.
#' @param stopcrit list of named stop criteria parameters. 
#'    See \code{Stop criteria} for details.
#' @param probpars list of named problem parameters.
#'    See \code{Problem Description} for details.
#' @param seed seed for the random number generator. 
#'    See \code{Random Seed} for details.
#' @param showpars parameters that regulate the echoing of progress indicators
#'    See \code{Showpars} for details.
#'
#' @return A list object containing the final population (sorted by performance)
#', the performance vector, and some run statistics.
#' @author Felipe Campelo (\email{fcampelo@@ufmg.br}) and Moises Botelho 
#'         (\email{moisesufop@@gmail.com})
#'
#' @examples
#' # DE/rand/1/bin with population 40, F = 0.8 and CR = 0.5
#' popsize  <- 100
#' mutpars  <- list(name = "mutation_rand", f = 0.8)
#' recpars  <- list(name = "recombination_bin", cr = 0.5, minchange = TRUE)
#' selpars  <- list(name = "selection_standard")
#' stopcrit <- list(names = "stop_maxiter", maxiter = 100)
#' probpars <- list(name  = "sphere",
#'                 xmin = rep(-5.12,10), xmax = rep(5.12,10))
#' seed <- NULL
#' showpars <- list(show.iters = "numbers", showevery = 1)
#' ExpDE(popsize, mutpars, recpars, selpars, stopcrit, probpars, seed, showpars)
#'
#' # DE/rand/2/blxAlpha
#' recpars  <- list(name = "recombination_blxAlpha", alpha = 0.1)
#' mutpars  <- list(name = "mutation_rand", f = 0.8, nvecs = 2)
#' ExpDE(popsize, mutpars, recpars, selpars, stopcrit, probpars)
#' 
#' # DE/best/1/sbx
#' recpars  <- list(name = "recombination_sbx", eta = 10)
#' mutpars  <- list(name = "mutation_best", f = 0.6, nvecs = 1)
#' ExpDE(popsize, mutpars, recpars, selpars, stopcrit, probpars)
#'
#' # DE/best/1/eigen/bin
#' recpars  <- list(name = "recombination_eigen", 
#'                  othername = "recombination_bin", 
#'                  cr = 0.5, minchange = TRUE)
#' showpars <- list(show.iters = "dots", showevery = 10)
#' stopcrit <- list(names = "stop_maxeval", maxevals = 10000)
#' ExpDE(popsize, mutpars, recpars, selpars, stopcrit, probpars, seed = 1234)
#' 
#' @references 
#' F. Campelo, M. Botelho, "Experimental Investigation of Recombination 
#' Operators for Differential Evolution", Genetic and Evolutionary 
#' Computation Conference, Denver/CO, July 2016 (Accepted).
#' @export

ExpDE <- function(popsize,
                  mutpars  = list(name = "mutation_rand",
                                  f    = 0.2),
                  recpars  = list(name  = "recombination_bin",
                                  cr    = 0.8, 
                                  nvecs = 1),
                  selpars  = list(name = "standard"),
                  stopcrit,
                  probpars,
                  seed = NULL,
                  showpars = list(show.iters = "none"))
{
  # ========== Error catching and default value definitions
  # Check seed
  stopifnot(is.null(seed) || seed > 0,
            is.null(seed) || is.numeric(seed),
            is.null(seed) || seed == floor(seed))
  
  if (is.null(seed)) {seed <- as.numeric(Sys.time())}
  set.seed(seed)
  
  # ==========
  
  # Generate initial population
  X <- create_population(popsize  = popsize,
                         probpars = probpars)

  # Evaluate the initial population
  J <- evaluate_population(probpars = probpars,
                           Pop      = X)

  # Prepare for iterative cycle:
  keep.running  <- TRUE     # stop criteria flag
  t             <- 0        # counter: iterations
  nfe           <- popsize  # counter: number of function evaluations


  # Iterative cycle
  while(keep.running){
    # Update iteration counter
    t <- t + 1          
    
    # Reset candidate vector performance values
    G <- NA * J

    # Mutation
    M <- do.call(mutpars$name,
                 args = list(X       = X,
                             mutpars = mutpars))


    # Recombination
    U <- do.call(recpars$name,
                 args = list(X       = X,
                             M       = M,
                             recpars = recpars))
    
    # Repair U
    U <- matrix(pmax(0, pmin(U, 1)), 
                byrow = FALSE, 
                nrow  = nrow(U))

    # Evaluate U 
    # Some recombination operators evaluate the 'offspring' solutions, so only
    # the 'unevaluated' ones need to be dealt with here.
    toeval <- is.na(G)
    G[toeval] <- evaluate_population(probpars = probpars,
                                     Pop      = U[toeval, ])
    nfe <- nfe + sum(toeval)

    # Selection
    next.pop <- do.call(selpars$name,
                        args = list(X = X,
                                    U = U,
                                    J = J,
                                    G = G))

    # Stop criteria
    keep.running <- check_stop_criteria()

    # Compose next population
    X <- next.pop$Xsel
    J <- next.pop$Jsel
    
    # Echo progress
    print_progress()
  }

  X <- denormalize_population(probpars, X[order(J), ])
  J <- sort(J)
  return(list(X     = X,
              Fx    = J,
              Xbest = X[1,],
              Fbest = J[1],
              nfe   = nfe,
              iter  = t))
}
