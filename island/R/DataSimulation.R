# This is the script for the functions and generics for the data generation
# tool.

# The only function should be envelope or newenvelope. Something like that.

#' Data simulation of colonization-extinction dynamics
#'
#' \code{data_generation} simulates species richness data according to the
#' stochastic model of island biogeography \cr \code{PA_simulation} simulates
#' presence-absence data according to the stochastic model of island
#' biogeography
#'
#' To simulate community assembly, we need an initial vector of
#' presence-absence, from which the subsequent assembly process will be
#' simulated. This initial vector is considered as \code{x[, column]}.
#'
#' @param x A dataframe with the vector of initial absences and presences.
#' @param column A number indicating the column with the initial
#'   presence-absence data.
#' @param transitions A vector with the transition probabilities of the
#'   simulation, in the form (T01, T10).
#' @param iter Number of times that the specified dynamics should be repeated.
#' @param times Number of temporal steps to simulate.
#' @examples data_generation(as.data.frame(rep(0, 100)), 1, c(0.5, 0.5), 5, 25)
#' data_generation(alonso[[1]], 3, c(0.5, 0.5), 5, 25)
#' PA_simulation(as.data.frame(c(rep(0, 163), rep(1, 57))), 1, c(0.13, 0.19),
#' 20)
#'
#' @return A matrix with species richness representing each row consecutive
#'   samples and each column a replica of the specified dynamics  or a matrix
#'   with presence-absence data for the specified dynamics, each row
#'   representing a species and each column consecutive samplings.
#'
#' @seealso \code{\link{cetotrans}} to obtain the transition probabilities
#'   asociated with a colonization-extinction pair.
#'
#' @export
data_generation <- function(x, column, transitions, iter, times) {
  y <- data.frame()
  total <- matrix(NA, ncol = iter, nrow = times)
  for (i in 1:iter) {
    y <- cbind(x[, column])
    for (j in 1:times) {
      T01 <- transitions[1]
      T10 <- transitions[2]
      al <- stats::runif(nrow(y))
      for (k in 1:nrow(y)) {
        if (y[k, 1] == 0) {
          if (al[k] < T10) y[k, 1] <- 1
        }  else {
          if (al[k] < T01) y[k, 1] <- 0
        }
      }
      total[j, i] <- colSums(y)
    }
  }
  total
}

#' @rdname data_generation
#' @export
PA_simulation <- function(x, column, transitions, times) {
  result <- matrix(NA, ncol = times, nrow = nrow(x))
  y <- as.matrix(x[, column])
  result[, 1] <- y
  for (j in 2:times) {
    T01 <- transitions[1]
    T10 <- transitions[2]
    al <- stats::runif(nrow(x))
    for (k in 1:nrow(x)) {
      if (y[k, 1] == 0) {
        if (al[k] < T10) {y[k, 1] <- 1}
      } else {
        if (al[k] < T01) {y[k, 1] <- 0}
      }
    }
    result[, j] <- y
  }
  result
}
