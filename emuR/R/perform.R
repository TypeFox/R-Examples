##' Performance (hit rate) of a confusion matrix
##' 
##' Performs (hit rate) of a confusion matrix
##' 
##' 
##' @param data A confusion matrix.
##' @return Caluculates the accuracy (total score) of the confusion matrix,
##' returning percentage of correct, and incorrect matches.
##' @seealso confusion
##' @keywords misc
##' @export perform
"perform" <- function(data)
{
  ## calculates total score in a confusion matrix, data
  k <- 0
  for(j in 1:nrow(data)) {
    k <- k + data[j, j]
  }
  total <- sum(data)
  wrong <- total - k
  correct <- (k/total) * 100
  wrong <- wrong/total * 100
  labcol <- c("correct %", "incorrect %")
  m <- cbind(correct, wrong)
  dimnames(m) <- list(NULL, labcol)
  m
}
