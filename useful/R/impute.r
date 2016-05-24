#' @title  simple.impute
#' @description Generic function for simple imputation.
#' @details Provides the ability to simply impute data based on a simple measure such as mean or median.  For more robust imputation see the packages Amelia, mice or mi.
#' @aliases simple.impute
#' @export simple.impute
#' @importFrom stats median
#' @author Jared P. Lander
#' @param x An object to be imputed
#' @param fun The function with which to fill in missing values
#' @param \dots Further arguments
#' @return An object with the missing values imputed.
#' @examples 
#' theDF <- data.frame(A=1:10, B=1:10, C=1:10)
#' theDF[c(1, 4, 6), c(1)] <- NA
#' theDF[c(3, 4, 8), c(3)] <- NA
#' 
#' simple.impute(theDF$A)
#' simple.impute(theDF$A, mean)
#' simple.impute(theDF$A, constant(4))
#' simple.impute(theDF)
#' simple.impute(theDF, mean)
#' simple.impute(theDF, constant(4))
#' 
simple.impute <- function(x, fun=median, ...)
{
    UseMethod('simple.impute')
}

#' @title simple.impute.default
#' @description Function for imputing a vector with missing data.
#' @details Provides the ability to simply impute data based on a simple measure such as mean or median.  For more robust imputation see the packages Amelia, mice or mi.
#' @aliases simple.impute.default
#' @export
#' @export simple.impute.default
#' @importFrom stats median
#' @author Jared P. Lander
#' @param x A numeric or integer vector
#' @param fun The function with which to fill in missing values
#' @param \dots Further arguments
#' @return An object with the missing values imputed.
#' @examples 
#' theDF <- data.frame(A=1:10, B=1:10, C=1:10)
#' theDF[c(1, 4, 6), c(1)] <- NA
#' theDF[c(3, 4, 8), c(3)] <- NA
#' 
#' simple.impute.default(theDF$A)
#' simple.impute.default(theDF$A, mean)
#' simple.impute.default(theDF$A, constant(4))
#' 
simple.impute.default <- function(x, fun=median, ...)
{
    # find missing values
    theNA <- is.na(x)
    # replace with the constructed value
    x[theNA] <- fun(x[!theNA])
    
    return(x)
}


#' @title simple.impute.data.frame
#' @description Function for imputing a data.frame with missing data.
#' @details Provides the ability to simply impute data based on a simple measure such as mean or median.  For more robust imputation see the packages Amelia, mice or mi.
#' 
#' Each column is imputed independently.
#' @aliases simple.impute.data.frame
#' @export
#' @export simple.impute.data.frame
#' @importFrom dplyr mutate_each_ funs
#' @importFrom magrittr "%>%"
#' @importFrom stats median
#' @author Jared P. Lander
#' @param x A data.frame
#' @param fun The function with which to fill in missing values
#' @param \dots Further arguments
#' @return A data.frame with the missing values imputed.
#' @examples 
#' theDF <- data.frame(A=1:10, B=1:10, C=1:10)
#' theDF[c(1, 4, 6), c(1)] <- NA
#' theDF[c(3, 4, 8), c(3)] <- NA
#' 
#' simple.impute.data.frame(theDF)
#' simple.impute.data.frame(theDF, mean)
#' simple.impute.data.frame(theDF, constant(4))
#' 
simple.impute.data.frame <- function(x, fun=median, ...)
{
    . <- NULL
    x %>% mutate_each_(funs(simple.impute(., fun=fun)), vars=names(x))
}

#' @title  simple.impute.tbl_df
#' @description Function for imputing a tbl_df with missing data.
#' @details Provides the ability to simply impute data based on a simple measure such as mean or median.  For more robust imputation see the packages Amelia, mice or mi.
#' 
#' Each column is imputed independently.
#' @aliases simple.impute.tbl_df
#' @export
#' @importFrom stats median
#' @author Jared P. Lander
#' @param x A data.frame
#' @param fun The function with which to fill in missing values
#' @param \dots Further arguments
#' @return A data.frame with the missing values imputed.
#' @examples 
#' theDF <- data.frame(A=1:10, B=1:10, C=1:10)
#' theDF[c(1, 4, 6), c(1)] <- NA
#' theDF[c(3, 4, 8), c(3)] <- NA
#' 
#' simple.impute.data.frame(theDF)
#' simple.impute.data.frame(theDF, mean)
#' simple.impute.data.frame(theDF, constant(4))
#' 
simple.impute.tbl_df <- function(x, fun=median, ...)
{
    simple.impute.data.frame(x=x, fun=fun, ...)
}

#' @title constant
#' @description Helper function for imputing contstants
#' @details Returns a function that always returns the value of n.
#' @export constant
#' @aliases constant
#' @author Jared P. Lander
#' @param n The value to return
#' @return A function that when used simply returns n.
#' @examples 
#' constant(4)(1:10)
#' 
#' theDF <- data.frame(A=1:10, B=1:10, C=1:10)
#' theDF[c(1, 4, 6), c(1)] <- NA
#' theDF[c(3, 4, 8), c(3)] <- NA
#' simple.impute(theDF, constant(4))
#' 
constant <- function(n=1)
{
    function(x, ...) n
}
