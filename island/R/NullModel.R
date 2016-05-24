### This file is intented to have the functions for the Null model and the
### R-squared of our predictions.

#' Model prediction error
#'
#' \code{r_squared} evaluates \eqn{R^2} for our simulated dynamics. \cr
#' \code{simulated_model} Error of the stochastic model. \cr \code{null_model}
#' Error of the null model. \cr
#'
#' The importance of assesing how well a model predicts new data is paramount.
#' The most used metric to assess this model error is \eqn{R^2}.  \eqn{R^2} is
#' always refered to a null model and is defined as follows: \deqn{R^{2} = 1 -
#' \epsilon^{2} / \epsilon^{2}_0}{R^2 = 1 - \epsilon^2 / \epsilon^2_0} where
#' \eqn{\epsilon^2} is the prediction error defined as the mean squared
#' deviation of model predictions from actual observations, and
#' \eqn{\epsilon^2_0} is a null model error, in example, an average of squared
#' deviations evaluated with a null model.
#'
#' Our null model corresponds with a random species model with no time
#' correlations, in which we draw randomly from a uniform distribution a number
#' of species between 0 and number of species observed in the species pool. The
#' expectation of the sum of squared errors under the null model is evaluated
#' analytically in Alonso et al. (2015).
#'
#' @note The value of \eqn{R^2} depends critically on the definition of the null
#'   model. Note that different definitions of the null model will lead to
#'   different values of \eqn{R^2}.
#'
#' @return \code{r_squared}  gives the value of \eqn{R^2} for the predictions of
#'   the model. \cr \cr \code{null_model} gives the average of squared
#'   deviations of the null model predictions from actual observations,
#'   \eqn{\epsilon^2_0}. \cr \cr \code{simulated_model} gives the average of
#'   squared deviations of the model predictions from the actual observations,
#'   \eqn{\epsilon^2}.
#'
#' @param observed A vector with the actual observed species richness.
#' @param simulated A vector with the simulated species richness.
#' @param sp Number of species in the species pool.
#'
#' @references Alonso, D., Pinyol-Gallemi, A., Alcoverro T. and Arthur, R..
#'   (2015) Fish community reassembly after a coral mass mortality: higher
#'   trophic groups are subject to increased rates of extinction. \emph{Ecology
#'   Letters}, \bold{18}, 451--461.
#'
#' @examples idaho.sim <- data_generation(as.data.frame(c(rep(0, 163),
#' rep(1, 57))), 1, c(0.162599, 0.111252), 250, 20)
#' idaho.me <- c(57, apply(idaho.sim, 1, quantile, 0.5))
#' r_squared(colSums(idaho[[1]][,3:23]), idaho.me, 220)
#'
#' null_model(colSums(idaho[[1]][,3:23]), 220)
#'
#' simulated_model(colSums(idaho[[1]][,3:23]), idaho.me)
#'
#' @export
r_squared <- function(observed, simulated, sp) {
  sol <- 1 - (simulated_model(observed, simulated) / null_model(observed, sp))
  names(sol) <- "R-Squared"
  sol
}

#' @rdname r_squared
#' @export
null_model <- function(observed, sp) {
  total <- 0

  for (i in 1:length(observed)) {
    partial <- 0
    for (j in 0:sp) {
      value <- (j - observed[i]) ^ 2
      partial <- partial + value
    }
    names(total) <- "Null Model Error"
    total <- total + partial / sp
  }

  total <- total / length(observed)
  total
}

#' @rdname r_squared
#' @export
simulated_model <- function(observed, simulated) {
  total <- 0

  for (i in 1:length(observed)) {
    value <- (simulated[i] - observed[i]) ^ 2
    total <- total + value
  }
  names(total) <- "Simulated Model Error"
  total <- total / length(observed)
  total
}
