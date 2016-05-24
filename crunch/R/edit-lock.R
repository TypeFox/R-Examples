##' Lock and unlock a dataset for editing
##'
##' Crunch allows a single active editor. If you have edit privileges but are
##' not currently editing the dataset, you must unlock the dataset before
##' making changes. You may then lock the dataset when you're done editing.
##' @param dataset a \code{CrunchDataset}
##' @return \code{dataset}, invisibly, after having set the current editor.
##' @export
lock <- function (dataset) {
    crPATCH(self(dataset), body=toJSON(list(current_editor=NULL)))
    invisible(dataset)
}

##' @rdname lock
##' @export
unlock <- function (dataset) {
    crPATCH(self(dataset), body=toJSON(list(current_editor=userURL())))
    invisible(dataset)
}
