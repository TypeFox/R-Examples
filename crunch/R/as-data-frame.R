CrunchDataFrame <- function (dataset) {
    ## S3 constructor method for CrunchDataFrame. terms.formula doesn't seem
    ## to like S4 subclass of environment
    stopifnot(is.dataset(dataset))
    out <- new.env()
    out$.crunchDataset <- dataset
    with(out, for (v in aliases(allVariables(dataset))) eval(parse(text=paste0("delayedAssign('", v, "', as.vector(.crunchDataset[['", v, "']]))"))))
    class(out) <- "CrunchDataFrame"
    return(out)
}

setOldClass("CrunchDataFrame")

##' @export
dim.CrunchDataFrame <- function (x) dim(x$.crunchDataset)

##' @export
names.CrunchDataFrame <- function (x) names(x$.crunchDataset)

##' as.data.frame method for CrunchDataset
##'
##' This method is defined principally so that you can use a CrunchDataset as
##' a \code{data} argument to other R functions (such as
##' \code{\link[stats]{lm}}). Unless you give it the \code{force==TRUE}
##' argument, this function does not in fact return a \code{data.frame}: it
##' returns an object with an interface like a data.frame, such that you get
##' R vectors when you access its columns (unlike a CrunchDataset, which
##' returns CrunchVariable objects). This allows modeling functions that
##' require select columns of a dataset to retrieve only those variables from
##' the remote server, rather than pulling the entire dataset into local
##' memory.
##'
##' @param x a CrunchDataset
##' @param row.names part of as.data.frame signature. Ignored.
##' @param optional part of as.data.frame signature. Ignored.
##' @param force logical: actually coerce the dataset to \code{data.frame}, or
##' leave the columns as unevaluated promises. Default is \code{FALSE}.
##' @param ... additional arguments passed to as.data.frame.default
##' @return an object of class \code{CrunchDataFrame} unless \code{force}, in
##' which case the return is a \code{data.frame}.
##' @name dataset-to-R
NULL

##' @rdname dataset-to-R
##' @export
as.data.frame.CrunchDataset <- function (x, row.names = NULL, optional = FALSE,
                                        force=FALSE, ...) {
    # message("CrunchDataset to CrunchDataFrame")
    out <- CrunchDataFrame(x)
    if (force) {
        out <- as.data.frame(out)
    }
    return(out)
}

##' @rdname dataset-to-R
##' @export
as.data.frame.CrunchDataFrame <- function (x, row.names = NULL, optional = FALSE, ...) {
    # message("CrunchDataFrame to data.frame")
    x <- x$.crunchDataset
    default.stringsAsFactors <- function () FALSE
    limit <- min(c(10000, getOption("crunch.data.frame.limit")))
    if (nrow(x) * ncol(x) > limit) {
        halt("Dataset too large to coerce to data.frame. ",
            "Consider subsetting it first")
    }
    out <- lapply(x, as.vector)
    names(out) <- names(x)
    # return(as.data.frame(out, ...))
    return(structure(out, class="data.frame", row.names=c(NA, -nrow(x))))
}
