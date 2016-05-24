#' Area under curve of predicted binary response
#' 
#'  Takes in actual binary response and predicted probabilities, and returns auc value
#'  @param y actual binary response
#'  @param yhat predicted probabilities corresponding to the actual binary response
#'  @details
#'  Area under the receiver operating characteristic (ROC) curve is the most sought after
#'  criteria for judging how good model predictions are.
#'  
#'  \code{auc} function calculates the true positive rates (TPR) and false positive
#'  rates (FPR) for each cutoff from 0.01 to 1 and calculates the area using trapezoidal
#'  approximation. A ROC curve is also generated.
#'  @return area under the ROC curve
#'  @author Akash Jain
#'  @seealso \code{\link{accuracy}}, \code{\link{ks}}, \code{\link{iv}}, \code{\link{gini}}, \code{\link{splitdata}}
#'  @examples
#'  # A 'data.frame' with y and yhat
#' df <- data.frame(y = c(1, 0, 1, 1, 0, 0, 1, 0, 1, 0),
#'                  yhat = c(0.86, 0.23, 0.65, 0.92, 0.37, 0.45, 0.72, 0.19, 0.92, 0.50))
#'
#' # AUC figure
#' AUC <- auc(y = df[, 'y'], yhat = df[, 'yhat'])
#'  @export
auc <- function(y, yhat) {
  if(length(unique(y)) != 2 | (class(y) != 'integer' && class(y) != 'numeric' && class(y) != 'factor')) {
    stop('Invalid input: y should be integer or factor vector representing a binary response')
  } else if(class(yhat) != 'numeric' | max(yhat) > 1 | min(yhat) < 0) {
    stop('Invalid input: yhat should be numeric vector of predicted probabilities in the range 0 to 1')
  } else if(length(y) != length(yhat)) {
    stop('Invalid input: vectors y and yhat should have the same length')
  } else {
    aucData <- function(y, yhat, cutoff) {
      ypred <- ifelse(yhat > cutoff, 1, 0)
      accuracyNum <- table(y, ypred)
      tpr <- accuracyNum[4]/sum(accuracyNum[2,])
      fpr <- 1 - (accuracyNum[1]/sum(accuracyNum[1,]))
      aucData <- data.frame(tpr = tpr, fpr = fpr)
      return(aucData)
    }
    aucs <- function(data, i) {
      area <- 0.5 *(data$fpr[i] - data$fpr[i-1]) * (data$tpr[i] + data$tpr[i-1])
      return(area)
    }
    ltData <- lapply(seq(0.01, 1, by = 0.01), function(cutoff) aucData(y, yhat, cutoff))
    data <- do.call('rbind', ltData)
    data <- data[order(data[, 'fpr']), ]
    vtAUC <- sapply(2:nrow(data), function(i) aucs(data, i))
    auc <- sum(vtAUC, na.rm = TRUE)
    par(xaxs = 'i', yaxs = 'i')
    plot(x = data$fpr, 
         y = data$tpr, 
         type = 'l', 
         xlab = 'False Positive Rate', 
         ylab = 'True Positive Rate',
         xlim = c(0, 1),
         ylim = c(0, 1),
         col = 'green',
         lwd = 2)
    return(auc)
  }
}