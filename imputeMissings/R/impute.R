#' Impute missing values with the median/mode or \code{randomForest}
#'
#' When the median/mode method is used: character vectors and factors are imputed with the mode. Numeric and integer vectors are imputed with the median.
#' When the random forest method is used predictors are first imputed with the mean/median and each variable is then predicted and imputed with that value.
#' For predictive contexts there is a \code{compute} and an \code{impute} function. The former is used on a training set to learn the values (or random forest models) to impute (used to predict).
#' The latter is used on both the training and new data to impute the values (or deploy the models) learned by the \code{compute} function.
#'
#' @param data A data frame with dummies or numeric variables. Categorical variables (i.e., non-dummy / indicator) variables work up to a certain number but are not recommended.
#' @param object If NULL \code{impute} will call \code{compute} on the current dataset. Otherwise it will accept the output of a call to \code{compute}
#' @param method Either "median/mode" or "randomForest". Only works if object = NULL
#' @return An imputed data frame.
#' @examples
#' #Example:
#' #create some data
#' data <- data.frame(V1=as.factor(c('yes','no','no',NA,'yes','yes','yes')),
#'                   V2=as.character(c(1,2,3,4,4,4,NA)),
#'                   V3=c(1:6,NA),V4=as.numeric(c(1:6,NA)))
#' #demonstrate function
#' object <- compute(data,method="randomForest")
#' object <- compute(data,method="median/mode")
#'
#' impute(data)
#' impute(data,object=compute(data, method="randomForest"))
#' impute(data,method="randomForest")
#'
#' @seealso \code{\link{compute}}
#' @author Matthijs Meire, Michel Ballings, Dirk Van den Poel, Maintainer: \email{Matthijs.Meire@@UGent.be}

impute <- function(data,  object=NULL, method="median/mode"){
  if (is.null(object)) object <- compute(data, method=method)
  if (!identical(colnames(data),names(object))) stop('Variable names and variable positions need to be identical in compute and impute')
  if (!class(object[[1]])=="randomForest"){
  data <- data.frame(sapply(1:ncol(data), function(i) {
        data[is.na(data[,i]),i] <- object[[i]]
        return(data[,i,drop=FALSE])
    }, simplify = FALSE))
  } else {
    # first impute predictors only with median/mode
    predictorsImputed <- impute(data, method="median/mode")
    # then use that data to predict response
    for (i in 1:ncol(data)){
      predicted <- predict(object[[i]],newdata=predictorsImputed[,-i], type="response")
      NAs <- is.na(data[,i])
      data[NAs,i] <- predicted[NAs]
    }

  }
  data
}




