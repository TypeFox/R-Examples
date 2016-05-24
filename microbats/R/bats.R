#' Bat Algorithm
#'
#' @description The function \code{bat_optim} implements a nature-inspired metaheuristic algorithm that deals with both continuous and discrete
#' optimization problems. The algorithm is based on the echolocation behavior of microbats and uses frequency tuning.
#'
#' @param D integer: the dimension of the search variables
#' @param NP integer: the population size, typically between 10 and 40
#' @param N_Gen integer: the number of generations, or iterations
#' @param A numeric; loudness, between 0 and 1, can be constant or decreasing
#' @param r numeric; pulse rate, must be positive, can be constant or decreasing
#' @param Qmin minimum frequency
#' @param Qmax maximum frequency
#' @param Lower lower bound of the search variables
#' @param Upper upper bound of the search variables
#' @param FUN the objective function to optimize, must be supplied by a user
#' @param ... optional arguments to \code{FUN}
#' @author Seong Hyun Hwang, Rachel Myoung Moon
#' @details
#' \code{bat_optim} implements the standard bat algorithm in three robust steps.
#' The first step is to initialize the parameters of algorithm to
#' generate and evaluate the initial population from which to determine the best solution.
#'
#' Secondly, a population of virtual microbats are moved in a \emph{d}-dimensional search or solution space according to
#' the updating rules of the algorithm: each bat is encoded with a velocity and a location at each iteration in the
#' search space. The location is a solution vector, and the current best solution is achieved.
#'
#' Then the current best solution is improved using random walks. The new solution is evaluated and updated.
#' See \emph{References} below for more details.
#'
#' In essence, frequency tuning acts as mutation to vary the solutions locally; hence, increasing the range of
#' frequencies leads to a global search. The mutation, compared with genetic algorithms, has no crossover but
#' depends on loudness and pulse emission. So technically, varying loudness and pulse emission rates can also
#' make the search intensive approaching the global optimality.
#'
#' One of the advantages of the bat algorithm is that it can converge very quickly at the initial stage and
#' can switch from exploration to exploitation when the optimality is approaching.
#' @return
#' Returns a list of four values: minimum fitness, population of solutions, fitness, best solution(s)
#' @examples
#' # find the x-value that gives the minimum of the quadratic function y = x^2 - 3x
#' # x should then be 1.5
#' quad_func <- function(D, sol) {
#'  val = 0
#'  for (i in 1:D) {
#'    val <- val + sol[i] * sol[i] - sol[i] * 3
#'  }
#'  return(val)
#' }
#'
#' # run a simulation using the standard bat algorithm
#' set.seed(5)  # for reproducive results
#' fit <- bat_optim(D = 1, NP = 20, N_Gen = 100, A = 0.5, r = 0.5,
#'                  Qmin = 0, Qmax = 2, Lower = -10, Upper = 10, FUN = quad_func)
#' x <- fit$best_solution
#'
#' @export
#' @importFrom stats runif rnorm
#' @references
#' [1] Yang, X.-S. "A new metaheuristic bat-inspired algorithm." Nature inspired cooperative strategies for optimization (NICSO 2010). Springer Berlin Heidelberg, 2010. 65-74.
#'
#' [2] Fister, I. Jr., Fister, I., Yang, X.-S., Fong, S., Zhuang, Y. "Bat algorithm: Recent advances." IEEE 15th International Symposium on Computational Intelligence and Informatics (CINTI), IEEE, 2014. 163-167.

bat_optim <- function(D, NP, N_Gen, A, r, Qmin, Qmax, Lower, Upper, FUN, ...)
{
  # initialize the parameters
  f_min <- 0
  Lb <- matrix(rep(Lower, D), nrow = 1)
  Ub <- matrix(rep(Upper, D), nrow = 1)
  Q <- matrix(rep(0, NP), nrow = 1)
  v <- matrix(0, nrow = NP, ncol = D)
  Sol <- matrix(0, nrow = NP, ncol = D)
  Fitness <- matrix(rep(0, NP), nrow = 1)
  best <- matrix(rep(0, D), nrow = 1)

  cat("Initializing the virtual microbats...\n")
  for (i in 1:NP) {
    Q[i] <- 0
    for (j in 1:D) {
      rnd <- runif(1, min = 0, max = 1)
      v[i, j] <- 0
      Sol[i, j] <- Lb[j] + (Ub[j] - Lb[j]) * rnd
    }
    Fitness[i] <- FUN(D, Sol[i,])
  }

  cat("Finding the best bat\n")
  i = 1
  j = 1
  for (i in 1:NP) {
    if (Fitness[i] < Fitness[j]) j <- i
  }
  for (i in 1:D) {
    best[i] <- Sol[j, i]
  }
  f_min <- Fitness[j]

  cat("Moving the bats via random walk\n")
  S <- matrix(0, nrow = NP, ncol = D)
  for (t in 1:N_Gen) {
    for (i in 1:NP) {
      rnd <- runif(1, min = 0, max = 1)
      Q[i] <- Qmin + (Qmin - Qmax) * rnd
      for (j in 1:D) {
        v[i, j] <- v[i, j] + (Sol[i, j] - best[j]) * Q[i]
        S[i, j] <- Sol[i, j] + v[i, j]
        if (S[i, j] < Lb[j]) S[i, j] <- Lb[j]
        if (S[i, j] > Ub[j]) S[i, j] <- Ub[j]
      }
      rnd <- runif(1, min = 0, max = 1)
      if (rnd > r) {
        for (j in 1:D) {
          S[i, j] <- best[j] + 0.001 * rnorm(1, mean = 0, sd = 1)
        }
      }
      Fnew <- FUN(D, S[i, ])
      rnd <- runif(1, min = 0, max = 1)
      if (Fnew <= Fitness[i] & rnd < A) {
        for (j in 1:D) {
          Sol[i, j] <- S[i, j]
        }
        Fitness[i] <- Fnew
      }
      if (Fnew <= f_min) {
        for (j in 1:D) {
          best[j] <- S[i, j]
        }
        f_min <- Fnew
      }
    }
  }
  return(list(min_fitness = f_min, population_solutions = Sol, fitness = Fitness, best_solution = best))
}
