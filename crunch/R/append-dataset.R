##' Append one Crunch dataset to another
##'
##' @param dataset1 a CrunchDataset
##' @param dataset2 another CrunchDataset, or possibly a data.frame. If
##' \code{dataset2} is not a Crunch dataset, it will be uploaded as a new
##' dataset before appending.
##' @param cleanup logical: if the append operation fails or is aborted, should
##' the intermediate batch created on \code{dataset1} be deleted? Default is
##' \code{TRUE}; \code{FALSE} may be useful if you want to review the append
##' conflicts in the web application.
##' @return A CrunchDataset with \code{dataset2} appended to \code{dataset1}
##' @export
##' @importFrom httpcache dropCache
appendDataset <- function (dataset1, dataset2, cleanup=TRUE) {

    stopifnot(is.dataset(dataset1))
    batch_url <- addBatchToDataset(dataset1, dataset2)
    dataset <- try(acceptAppendResolutions(batch_url, dataset1), silent=TRUE)
    if (is.error(dataset)) {
        if (cleanup) {
            d <- try(crDELETE(batch_url,
                    drop=dropCache(shojiURL(dataset1, "catalogs", "batches"))),
                silent=TRUE)
            if (is.error(d)) {
                warning("Batch ", batch_url,
                    " could not be deleted. It may still be processing.")
            }
        } else {
            message("Batch URL: ", batch_url) ## So you can fix and retry
        }
        rethrow(dataset)
    }
    invisible(dataset)
}

addBatchToDataset <- function (dataset1, dataset2) {
    if (!is.dataset(dataset2)) {
        ## TODO: compose batch directly, not as dataset?
        temp.ds.name <- paste("Appending to", name(dataset1), now())
        message("Creating ", dQuote(temp.ds.name), " as temporary dataset")
        dataset2 <- newDataset(dataset2, name=temp.ds.name)
    }

    ## Validate
    if (identical(self(dataset1), self(dataset2))) {
        halt("Cannot append dataset to itself")
    }

    batches_url <- shojiURL(dataset1, "catalogs", "batches")
    body <- list(
        element="shoji:entity",
        body=list(
            dataset=self(dataset2)
        )
    )
    invisible(crPOST(batches_url, body=toJSON(body)))
}

acceptAppendResolutions <- function (batch_url, dataset, ...) {

    status <- pollBatchStatus(batch_url, batches(dataset), until=c("ready", "imported"))

    batch <- ShojiObject(crGET(batch_url))
    cflicts <- flattenConflicts(batch@body$conflicts)

    if (status == "conflict") {
        failures <- formatFailures(cflicts)
        for (i in failures) message(i)
        err <- c("There are conflicts that cannot be resolved automatically.",
            "Please manually address them and retry.")
        halt(paste(err, collapse=" "))
    }

    resolutions <- formatConflicts(cflicts)
    ## Report on what was done/will be done
    for (i in resolutions) message(i)

    invisible(refresh(dataset))
}
