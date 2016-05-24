#' Compute the top decile lift and plot the lift curve
#'
#' \code{TopDecileLift} computes the commonly used top decile lift by ordering the data by the predictions, and computing the proportion of positives in the top 10\%.
#'
#' @param predicted A numeric vector with the classifier's predicted scores / probabilities
#' @param labels An integer vector containing binary labels with values {0,1}
#'
#' @return a scalar denoting the top decile lift
#' @examples
#' data(churn)
#' TopDecileLift(churn$predictions,churn$labels)
#' @author Steven Hoornaert, Michel Ballings, Dirk Van den Poel, Maintainer: \email{Steven.Hoornaert@@UGent.be}
TopDecileLift<- function(predicted,labels){
  if(is.factor(labels)) labels <- as.integer(as.character(labels))
  lift <- cbind(predicted,labels)
  lift <- lift[order(-lift[,1]),]
  lift <- as.numeric(mean(lift[1:round(nrow(lift)/10,0),2]) /
                           mean(lift[,2]))
  round(lift,digits=3)
}
