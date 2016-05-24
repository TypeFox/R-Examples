#' Confusion matrix and overall accuracy of predicted binary response
#' 
#'  Takes in actual binary response, predicted probabilities and cutoff value, and 
#'  returns confusion matrix and overall accuracy
#'  @param y actual binary response variable
#'  @param yhat predicted probabilities corresponding to the actual binary response
#'  @param cutoff threshold value in the range 0 to 1
#'  @details 
#'  When we predict a binary response, first thing that we want to check is accuracy of
#'  the model for a particular cutoff value. This function does just that and provides 
#'  confusion matrix (numbers and percentage) and overall accuracy. Overall accuracy is
#'  calculated as (TP + TN)/(P + N). 
#'  
#'  The output is a list from which the individual elements can be picked as shown in
#'  the example.
#'  @return a three element list: confusion matrix as a table, confusion matrix (percentages)
#'          as a table and overall accuracy value
#'  @author Akash Jain
#'  @seealso \code{\link{ks}}, \code{\link{auc}}, \code{\link{iv}}, \code{\link{splitdata}}
#'  @examples
#'  # A 'data.frame' with y and yhat
#' df <- data.frame(y = c(1, 0, 1, 1, 0),
#'                  yhat = c(0.86, 0.23, 0.65, 0.92, 0.37))
#'
#' # Accuracy tables and overall accuracy figures 
#' ltAccuracy <- accuracy(y = df[, 'y'], yhat = df[, 'yhat'], cutoff = 0.7)
#' accuracyNumber <- ltAccuracy$accuracyNum
#' accuracyPercentage <- ltAccuracy$accuracyPer
#' overallAccuracy <- ltAccuracy$overallAcc
#'  @export
accuracy <- function(y, yhat, cutoff) {
  if(length(unique(y)) != 2 | (class(y) != 'integer' && class(y) != 'numeric' && class(y) != 'factor')) {
    stop('Invalid input: y should be integer or factor vector representing a binary response')
  } else if(class(yhat) != 'numeric' | max(yhat) > 1 | min(yhat) < 0) {
    stop('Invalid input: yhat should be numeric vector of predicted probabilities in the range 0 to 1')
  } else if(length(y) != length(yhat)) {
    stop('Invalid input: vectors y and yhat should have the same length')
  } else if(class(cutoff) != 'numeric' | length(cutoff) != 1 | cutoff > 1 | cutoff < 0) {
    stop('Invalid input: cutoff should be a numeric value between 0 and 1')
  } else {
    ypred <- ifelse(yhat > cutoff, 1, 0)
    accuracyNum <- table(y, ypred)
    accuracyPer <- (accuracyNum/rowSums(accuracyNum)) * 100
    overallAccuracy <- (sum(diag(accuracyNum))/sum(colSums(accuracyNum))) * 100
    lt <- list(accuracyNum = accuracyNum,
               accuracyPer = round(accuracyPer, digits = 2),
               overallAcc = round(overallAccuracy, digits = 2))
    return(lt)
  }
}