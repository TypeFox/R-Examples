##' Create a Data Frame from All Combinations of Factors
##'
##' Simple wrapper of the 'expand.grid' function.  If x is a table
##' then a data frame is returned with one row pr individual
##' observation.
##' @title Create a Data Frame from All Combinations of Factors
##' @param _data Data.frame
##' @param ... vectors, factors or a list containing these
##' @author Klaus K. Holst
##' @export
##' @examples
##' dd <- Expand(iris, Sepal.Length=2:8, Species=c("virginica","setosa"))
##' summary(dd)
##'
##' T <- with(warpbreaks, table(wool, tension))
##' Expand(T)
Expand <- function(`_data`,...) {
    if (missing(`_data`)) {
        return(expand.grid(...))
    }
    if (inherits(`_data`,"table")) {
        M <- as.data.frame(`_data`)
        idx <- rep(seq(nrow(M)),M[,ncol(M)])
        return(M[idx,-ncol(M),drop=FALSE])
    }
    if (!inherits(`_data`,"data.frame")) {
        return(expand.grid(`_data`,...))
    }
    dots <- list(...)
    nn <- names(dots)
    for (n in nn) {
        y <- dots[[n]]
        if (is.factor(`_data`[1,n])) {
            dots[[n]] <- factor(y,levels=levels(`_data`[1,n]))
        }
    }
    do.call("expand.grid",dots)
}
