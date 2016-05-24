forks <- function (dataset) {
    stopifnot(is.dataset(dataset))
    return(ForkCatalog(crGET(shojiURL(dataset, "catalogs", "forks"))))
}

##' Create a fork of a dataset
##'
##' As with many other version control systems, in Crunch you can fork a
##' dataset's revision history, effectively making a copy on which you can work
##' independently of the original dataset. You can then merge those change back
##' to the original dataset or keep working independently.
##' @param dataset The \code{CrunchDataset} to fork
##' @param forkname character name to give the fork. If omitted, one will be
##' provided for you
##' @return The new fork, a \code{CrunchDataset}.
##' @export
forkDataset <- function (dataset, forkname) {
    if (missing(forkname)) {
        nforks <- length(forks(dataset))
        prefix <- ifelse(nforks, paste0("Fork #", nforks + 1, " of"), "Fork of")
        forkname <- paste(prefix, name(dataset))
    }
    fork_url <- crPOST(shojiURL(dataset, "catalogs", "forks"),
        body=toJSON(list(element="shoji:entity", body=list(name=forkname))))
    updateDatasetList()
    invisible(entity(datasetCatalog()[[fork_url]]))
}

##' Merge changes to a dataset from a fork
##'
##' @param dataset The \code{CrunchDataset} to merge to
##' @param fork The \code{CrunchDataset}, perhaps forked from \code{dataset},
##' that is to be merged in.
##' @param autorollback logical If the merge fails, should \code{dataset} be
##' restored to its state prior to the merge, or should it be left in its
##' partially merged state for debugging and manual fixing? Default is
##' \code{TRUE}, i.e. the former.
##' @return \code{dataset} with changes from \code{fork} merged to it.
##' @export
mergeFork <- function (dataset, fork, autorollback=TRUE) {
    m <- crPOST(shojiURL(dataset, "catalogs", "actions"), body=toJSON(list(
        element="shoji:entity",
        body=list(dataset=self(fork), autorollback=autorollback)
    )))
    if (is.character(m)) {
        ## Should only be character if the merge got Jobbed and returned
        ## 202 Location with the progress URL
        pollProgress(m)
    }
    return(refresh(dataset))
}
