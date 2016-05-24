#' @encoding UTF-8
#' @title Generate dummy variables
#'
#' @description Provides an alternative to generate dummy variables
#'
#' @param x a column position to generate dummies
#' @param data the data object as a data.frame
#' @param drop A logical value. If \code{TRUE}, unused levels will be omitted
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#'
#' @details A matrix object
#'
#' @keywords Models
#'
#' @examples
#' df <- data.frame(y = rnorm(25), x = runif(25,0,1), sex = sample(1:2, 25, rep=TRUE))
#'
#' dummy(df$sex)
#'
#' @export
`dummy` <-
  function (x, data = NULL, drop = TRUE)
  {
    if (is.null(data)) {
      varname <- as.character(sys.call(1))[2]
      varname <- sub("^(.*\\$)", "", varname)
      varname <- sub("\\[.*\\]$", "", varname)
    }
    else {
      if (length(x) > 1)
        stop("More than one variable to create dummies at same  time.")
      varname <- x
      x <- data[, varname]
    }
    if (drop == FALSE && class(x) == "factor") {
      x <- factor(x, levels = levels(x), exclude = NULL)
    }
    else {
      x <- factor(x, exclude = NULL)
    }
    if (length(levels(x)) < 2) {
      warning(varname, " has only 1 dimension. Generating dummy variable anyway.")
      return(matrix(rep(1, length(x)), ncol = 1, dimnames = list(rownames(x),
                                                                 c(paste(varname, "_", x[[1]], sep = "")))))
    }
    mat <- model.matrix(~x - 1, stats::model.frame(~x - 1), contrasts = FALSE)
    colnames.mm <- colnames(mat)
    cat(" ", varname, ":", ncol(mat), "dummy variables generated\n")
    mat <- matrix(as.integer(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = list(NULL,
                                                                                       colnames.mm))
    colnames(mat) <- sub("^x", paste(varname, "_", sep = ""), colnames(mat))
    if (!is.null(row.names(data)))
      rownames(mat) <- rownames(data)
    return(mat)
  }
NULL


#' @encoding UTF-8
#' @title Extraction of categorical values as a preprocessing step for making dummy variables
#'
#' @description  \code{categories} stores all the categorical values that are present in the factors and character vectors of a data frame. Numeric and integer vectors are ignored. It is a preprocessing step for the \code{dummy} function. This function is appropriate for settings in which the user only wants to compute dummies for the categorical values that were present in another data set. This is especially useful in predictive modeling, when the new (test) data has more or other categories than the training data.
#'
#' @param x data frame containing factors or character vectors that need to be transformed to dummies. Numerics, dates and integers will be ignored.
#' @param p select the top p values in terms of frequency. Either "all" (all categories in all variables), an integer scalar (top p categories in all variables), or a vector of integers (number of top categories per variable in order of appearance.
#' @examples
#' #create toy data
#' (traindata <- data.frame(xvar=as.factor(c("a","b","b","c")),
#'                          yvar=as.factor(c(1,1,2,3)),
#'                          var3=c("val1","val2","val3","val3"),
#'                          stringsAsFactors=FALSE))
#' (newdata <- data.frame(xvar=as.factor(c("a","b","b","c","d","d")),
#'                        yvar=as.factor(c(1,1,2,3,4,5)),
#'                        var3=c("val1","val2","val3","val3","val4","val4"),
#'                        stringsAsFactors=FALSE))
#'
#' categories(x=traindata,p="all")
#' categories(x=traindata,p=2)
#' categories(x=traindata,p=c(2,1,3))
#' @seealso \code{\link{dummy}}
#' @return  A list containing the variable names and the categories
#' @author Authors: Michel Ballings, and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
#' @export
`categories` <- function(x, p="all"){
  categoricals <- which(sapply(x,function(x) is.factor(x) || is.character(x)))
  x <- data.frame(x[,categoricals])
  cats <- sapply(1:ncol(x),function(z) {
    cats <- table(x[,z])
    if(is.numeric(p) && length(p) == 1) {
      names(sort(cats,decreasing=TRUE)[1:if(length(cats) <= p) length(cats) else p])
    } else if (is.numeric(p) && length(p) >= 1) {
      names(sort(cats,decreasing=TRUE)[1:if(length(cats) <= p[z]) length(cats) else p[z]])
    } else if (p=="all") {
      names(cats)
    }
  },simplify=FALSE)
  names(cats) <- names(x)
  cats
}
NULL

