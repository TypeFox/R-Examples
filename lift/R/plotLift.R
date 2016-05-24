#' Compute the top decile lift and plot the lift curve
#'
#' \code{plotLift} plots the commonly used lift curve by ordering the data by the predictions, and computing the proportion of positives for each bucket.
#'
#' @param predicted A numeric vector with the classifier's predicted scores / probabilities
#' @param labels An integer vector containing binary labels with values {0,1}
#' @param cumulative boolean. Should the cumulative lift curve be plotted or not?
#' @param n.buckets scalar. How many buckets should be used. One can use more buckets with large datasets
#' @param ... additional parameters to the \code{plot} function
#'
#' @return lift curve
#'
#' @examples
#' data(churn)
#' plotLift(churn$predictions,churn$labels)
#' @author Steven Hoornaert, Michel Ballings, Dirk Van den Poel, Maintainer: \email{Steven.Hoornaert@@UGent.be}
plotLift <- function(predicted,labels,cumulative=TRUE, n.buckets=10,...){
  if(is.factor(labels)) labels <- as.integer(as.character(labels))
  lift <- cbind(predicted,labels)
  lift <- lift[order(-lift[,1]),]

  buckets <- ceiling(seq_along(lift[,2])/floor(length(lift[,2])/n.buckets))

  cap <- floor(length(lift[,2])/n.buckets) * n.buckets
  lift <-  aggregate(lift[1:cap,2], by=list(buckets[1:cap]), mean)
  ylab <- "Lift"
  if (cumulative) {
    lift[,2] <- cumsum(lift[,2])/seq_along(lift[,2])
    ylab <- "Cumulative lift"
  }
  plot(lift[,1],lift[,2]/mean(labels), type="l", ylab=ylab, xlab="bucket",... )
}
