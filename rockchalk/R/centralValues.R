##' Central Tendency estimates for variables
##'
##' This is needed for the creation of summaries and predicted values
##' of regression models. It takes a data frame and returns a new data
##' frame with one row in which the mean or mode of the columns is
##' reported.
##'
##' @param x a data frame
##' @return a data frame with the same variables and one row, the summary indicators.
##' @author Paul E. Johnson \email{pauljohn@@ku.edu}
##' @export
##' @examples
##' myDat <- data.frame(x=rnorm(100), y=rpois(100,l=4), z = cut(rnorm(100), c(-10,-1,0,10)))
##' centralValues(myDat)

centralValues <- function(x){
    if( !is.data.frame(x)) stop("represent: x must be a data frame!")
    nc <- NCOL(x)
    nams <- colnames(x)
    represents <- x[1, , drop = FALSE] ## row 1, so we inherit the df's properties, labels
    for (i in 1: nc) {
        xvar <- x[ , i]
        if (is.numeric(xvar)){
            represents[1, i] <- mean(xvar, na.rm = TRUE)
        } else if (is.factor(xvar)) {
            xvartable <- table(xvar)
            xvarmode <- names(xvartable)[which.max(xvartable)]
            represents[1, i] <- xvarmode          
        } else {
            xvar <- factor(xvar)
            xvartable <- table(xvar)
            xvarmode <- names(xvartable)[which.max(xvartable)]
            represents[1, i] <- xvarmode      
        }
    }
    as.data.frame(represents)
    rownames(represents) <- NULL
    represents
}
