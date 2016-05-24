#' Compute Bonferroni Intervals
#'
#' This function will compute apply the Bonferroni correction to the selected
#' populations. 
#' 
#' @param X is a matrix or data frame that contains the responses. Each column
#' represents a different population.
#' @param alpha denotes the significance level of the intervals to be formed.
#' @param k corresponds to the number of populations to be selected.
#'
#' @export
#'
#' @details If there are p populations, then the Bonferroni correction will be
#' applied by using \eqn{\alpha / p} instead of just p.
#'
#' @examples
#' set.seed(18)
#' p <- 10; n <- 10
#' Xmat <- matrix(rnorm(p*n), nrow=n, ncol=p)
#' colnames(Xmat) <- paste("p.", 1:p, sep="")
#' bonferroniIntervals(Xmat, alpha=0.1, k=4)
#'
#' @seealso
#' \code{\link{asymmetricIntervals}}, \code{\link{bootstrapIntervals}}
#'
#' @return The function returns a matrix with k rows and 3 columns. This is
#' similar to the output of the predict.lm function of R.
#'
bonferroniIntervals <- function(X, alpha=0.05, k=2) {
  # Unstack the X data frame
  X.stack <- stack(data.frame(X))
  p <- ncol(X)
  if(k > p) 
    stop("More populations selected than available. Please ensure that k < p.")
  
  # fit the linear model
  mod <- lm(values ~ ind, data=X.stack)
  
  # Set up the prediction matrix
  pop.names <- unique(X.stack$ind)
  pred.df <- data.frame(ind = pop.names)

  
  # call predict.lm and return the appropriate number of rows
  out <- predict(mod, newdata=pred.df, interval="confidence", 
    level=1-alpha/p, type="response")
  pop.order <-  order(out[,"fit"], decreasing=TRUE)
  out <- out[pop.order[1:k], ]
  if(is.matrix(out)) {
    rownames(out) <- pop.names[pop.order[1:k]]
  } else {
    out <- matrix(out, nrow=1)
    rownames(out) <- pop.names[pop.order[1]]
    colnames(out) <- c("fit", "lwr", "upr")
  }
  
  out
}
