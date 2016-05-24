# This is the script for the functions and generics for the evaluation of the
# different models.


#' Akaike Information Criterium
#'
#' \code{akaikeic} calculates the Akaike Information Criterium (AIC) of a model.
#' \cr \code{akaikeicc} calculates the corrected Akaike Information Criterium
#' (AICc) for small samples.
#'
#' \deqn{AIC = 2 * k - 2 * lnL} \deqn{AICc = 2 * k - 2 * lnL + 2 * k *
#'   (k + 1) / (n - k - 1)}
#'
#' @param lnL Log-likelihood of the model.
#' @param k Number of parameters of the model.
#' @param n Sample size.
#' @export
#' @return A number with the AIC value for a model with k parameters and
#'   Log-likelihood lnL, or the AICc value for a model with k parameters,
#'   Log-likelihood lnL and sample size n.
#' @seealso \code{\link{weight_of_evidence}}
#' @examples akaikeic(-1485.926, 3)
#' akaikeicc(736.47, 6, 15)
#' akaikeicc(736.47, 6, 100)
akaikeic <- function(lnL, k) {

  2 * k - 2 * lnL
}

#' @rdname akaikeic
#' @export
akaikeicc <- function(lnL, k, n) {

  2 * k - 2 * lnL + ( (2 * k * (k + 1)) / (n - k - 1))
}

#' Weight of evidence
#'
#' \code{weight_of_evidence} calculates the weight of evidence of a set of
#' nested models.
#'
#' Calculates the weight of evidence in favor of model i being the actual
#' Kullback-Leibler best model given a set of models R for your data.
#'
#' @param data A dataframe with the names of the models in the first column and
#'   their AIC values in the second column.
#'
#' @details \deqn{w_i = exp( - 1/2 * \Delta_i) / \Sigma exp( - 1/2 * \Delta_r)}
#'   \deqn{r = 1, R}
#'
#' @return A dataframe with the names of the analysed models, their AIC
#'   differences with respect to the best model and the w_i of each one.
#' @seealso \code{\link{akaikeic}}
#' @references K. P.  Burnham, D. R. Anderson. \emph{Model selection and
#'   multimodel inference: a practical information-theoretic approach} (New
#'   York:Springer, ed. 2, 2002).
#'
#' @examples models <- c("Best_3k", "Best_4k", "Best_5k", "Best_6k", "Best_7k",
#'   "Best_8k", "Best_9k")
#'
#'   aks <- c(2977.852, 2968.568, 2957.384, 2952.618,
#'   2949.128, 2947.038, 2943.480)
#'
#'   weight_of_evidence(cbind(models, aks))
#'
#' @export
weight_of_evidence <- function(data) {
  ## This function returns the associated weight of evidence that a model has
  ## in function of the data.

  ## Input must be a dataframe with the names of the models in column 1 and
  ## their AIC values in the second column.

  ## Initialize the output dataframe
  out <- data.frame()
  ## Initializing variables; x represents the relative likelihood of the
  ## model given the MLEs of model parameters based on the same data, w is
  ## considered as the weight of evidence in favor of the model being the
  ## actual Kullback-Leibler best model for the data given the current set of
  ## models, IncAIC is the difference of the AIC of the model compared with
  ## that with the minimum AIC, and SumIncAIC its the sum of the AIC
  ## differences of all the models in data. See Burnham&Anderson (1998) for
  ## details.
  x <- c(nrow=1, ncol=nrow(data))
  w <- c(nrow=1, ncol=nrow(data))
  IncAIC <- c(nrow=1, ncol=nrow(data))
  SumIncAIC <- 0

  for (i in 1:nrow(data)) {
    IncAIC[i] <- as.numeric(data[i, 2]) - as.numeric(min(data[, 2]))
    x[i] <- exp( - (1 / 2) * IncAIC[i])
    if (i == length(data[, 2])) SumIncAIC <- sum(x)
    }

  for (i in 1:nrow(data)) {
    w[i] <- x[i] / SumIncAIC
    }

  out <- cbind(data[, 1], as.numeric(IncAIC), as.numeric(w))
  colnames(out) <- c("Model", "IncAIC", "w")
  out
}


