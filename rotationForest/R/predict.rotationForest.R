#' Predict method for rotationForest objects
#'
#' Prediction of new data using rotationForest.
#' 
#' @param object An object of class \code{rotationForest}, as created by the function \code{rotationForest}
#' @param newdata A data frame with the same predictors as in the training data.
#' @param all Return the predictions per tree instead of the average.
#' @param ... Not used currently.
#' @examples 
#' data(iris)
#' y <- as.factor(ifelse(iris$Species[1:100]=="setosa",0,1))
#' x <- iris[1:100,-5]
#' rF <- rotationForest(x,y)
#' predict(object=rF,newdata=x)
#' @references Rodriguez, J.J., Kuncheva, L.I., 2006. Rotation forest: A new classifier ensemble method. IEEE Trans. Pattern Anal. Mach. Intell. 28, 1619-1630. doi:10.1109/TPAMI.2006.211
#' @seealso \code{\link{rotationForest}}
#' @return A vector containing the response scores.
#' @author Authors: Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
#' @method predict rotationForest
#' @keywords classification
predict.rotationForest <- function(object,newdata,all=FALSE,...){
  newdata <- data.frame(sapply(newdata,as.numeric))
  #check if column names are in same order and if are the same number of columns
  if (!identical(colnames(newdata),object$columnnames)) stop("Variable names and/or order of variables in newdata is not identical to training set. Please check if variables are exactly the same in both sets.")
  predicted <- matrix(NA,nrow=nrow(newdata),ncol=length(object$models))
  for (i in 1:length(object$models)){

    #create PCs
    final <- data.frame(as.matrix(newdata) %*% as.matrix(object$loadings[[i]]))
  
    predicted[,i] <- predict(object$models[[i]], final ,type='prob')[,2]
    
    
  }
  
  if (all) {
    predicted 
  } else {
    rowMeans(predicted)
  }
  
}