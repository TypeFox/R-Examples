#' Guess
#' 
#' Infers the category of a new observation.
#' 
#' Given an already trained and smoothed Probabilistic neural network, the function \code{guess} gives the category with the highest probability, and the probabilities of each category.
#' 
#' @param nn A trained and smoothed Probabilistic neural network.
#' @param X A vector describing a new observation.
#' 
#' @seealso \code{\link{pnn-package}}, \code{\link{learn}}, \code{\link{smooth}}, \code{\link{perf}}, \code{\link{norms}}
#' 
#' @return A \code{list} of the guessed category and the probabilities of each category.
#' 
#' @examples
#' library(pnn)
#' data(norms)
#' pnn <- learn(norms)
#' pnn <- smooth(pnn, sigma=0.8)
#' guess(pnn, c(1,1))
#' guess(pnn, c(1,1))$category
#' guess(pnn, c(1,1))$probabilities
#' guess(pnn, c(2,1))
#' guess(pnn, c(1.5,1))
#' @export
guess <- function(nn, X) {
    X <- matrix(X, ncol=nn$k)
    probs <- guess.probabilities.of.each.category(nn, X)
    if(is.na(probs[1])) return(NA)
    category <- names(probs[probs == max(probs)])
    results <- list(category=category, probabilities=probs)
    return(results)
}
