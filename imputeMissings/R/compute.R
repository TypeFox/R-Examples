#' Compute the missing values to later impute them in another dataset
#'
#' When the median/mode method is used: character vectors and factors are imputed with the mode. Numeric and integer vectors are imputed with the median.
#' When the random forest method is used predictors are first imputed with the mean/median and each variable is then predicted and imputed with that value.
#' For predictive contexts there is a \code{compute} and an \code{impute} function. The former is used on a training set to learn the values (or random forest models) to impute (used to predict).
#' The latter is used on both the training and new data to impute the values (or deploy the models) learned by the \code{compute} function.
#'
#' @param data A data frame with dummies or numeric variables. Categorical variables (i.e., non-dummy / indicator) variables work up to a certain number but are not recommended.
#' @param method Either "median/mode" or "randomForest"
#' @param ... additional arguments for \code{randomForest}
#' @return Values or models used for imputation
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
#' @seealso \code{\link{impute}}
#' @author Matthijs Meire, Michel Ballings, Dirk Van den Poel, Maintainer: \email{Matthijs.Meire@@UGent.be}
compute <- function (data, method="median/mode", ...)
{
    if (method=="median/mode") {
        Mode <- function(x) {
            xtab <- table(x)
            xmode <- names(which(xtab == max(xtab)))
            return(xmode[1])
        }

        values <- sapply(data, function(x) {
            if (class(x) %in% c("character", "factor"))
                Mode(x)
            else if (class(x) %in% c("numeric", "integer"))
                median(x, na.rm = TRUE)

        }, simplify = FALSE)
    } else if (method=="randomForest") {

        data <- impute(data, method="median/mode")

        values <- list()
        for (i in 1:ncol(data)){

          values[[i]] <- randomForest(x=data[,-i,drop=FALSE],y= data[,i], ...)
        }
        names(values) <- colnames(data)
    }
    values
}



