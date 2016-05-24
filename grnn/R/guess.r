#' Guess
#' 
#' Infers the value of a new observation.
#' 
#' @examples
#' n <- 100
#' set.seed(1)
#' x <- runif(n, -2, 2)
#' y0 <- x^3
#' epsilon <- rnorm(n, 0, .1)
#' y <- y0 + epsilon
#' grnn <- learn(data.frame(y,x))
#' grnn <- smooth(grnn, sigma=0.1)
#' guess(grnn, -2)
#' guess(grnn, -1)
#' guess(grnn, -0.2)
#' guess(grnn, -0.1)
#' guess(grnn, 0)
#' guess(grnn, 0.1)
#' guess(grnn, 0.2)
#' guess(grnn, 1)
#' guess(grnn, 2)
#' @param nn A trained and smoothed General regression neural network.
#' @param X A vector describing a new observation.
#' @seealso \code{\link{grnn-package}}
#' @export
guess <- function(nn, X) {
    results <- Y(nn$Xa, nn$Ya, X, nn$sigma)
    return(results)
}
