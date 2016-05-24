#' Perf
#' 
#' Performance of a Probabilist neural network.
#' 
#' The function \code{perf} uses a hold-out method. This method takes the training set used by the function \code{\link{learn}} and iterate over each observation trying to guess the current observation with a reduced training set (without the current observation).It generates:
#' \itemize{
#'  \item Two lists of \code{observed} and \code{guessed} values.
#'  \item the following statistics: number of \code{success} and \code{fails}, a sucess rate (\code{success_rate}) and a \code{bic} indicator.
#' }
#' 
#' @param nn A trained and smoothed Probabilist neural network.
#' 
#' @return A probabilist neural network updated with its performance.
#' 
#' @seealso \code{\link{pnn-package}}, \code{\link{learn}}, \code{\link{smooth}}, \code{\link{guess}}, \code{\link{norms}}
#' 
#' @examples
#' library(pnn)
#' data(norms)
#' pnn <- learn(norms)
#' pnn <- smooth(pnn, sigma=0.8)
#' pnn <- perf(pnn)
#' pnn$observed
#' pnn$guessed
#' pnn$success
#' pnn$fails
#' pnn$success_rate
#' pnn$bic
#' @export
perf <- function(nn) {
    ts_hdo <- hold_out(nn$set)
    success <- 0
    guessed <- c()
    while(length(ts <- ts_hdo()) > 1) {
        model0 <- learn(set=ts$rest, category.column=nn$category.column)
        model0$sigma <- nn$sigma
        X <- ts$one[-nn$category.column]
        category <- guess(model0, as.matrix(X))$category
        if(!is.na(category) & category == ts$one[nn$category.column]) success <- success + 1
        guessed <- c(guessed, category)
    }
    nn$observed <- nn$set[,nn$category.column]
    nn$guessed <- factor(guessed)
    nn$success <- success
    nn$fails <- nn$n - nn$success
    nn$success_rate <- nn$success / nn$n
    nn$bic <- nn$n * log( (nn$n - nn$success) / (nn$n - 1) ) + nn$k * log(nn$n)
    return(nn)
}
