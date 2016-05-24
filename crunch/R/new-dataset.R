##' Upload a data.frame to Crunch to make a new dataset
##'
##' @param x a data.frame or other rectangular R object
##' @param name character, the name to give the new Crunch dataset. Default is
##' the name of the R object passed in \code{x}
##' @param ... additional arguments passed to \code{ \link{createDataset}}
##' @return If successful, an object of class CrunchDataset.
##' @export
newDataset <- function (x, name=as.character(substitute(x)), ...) {

    Call <- match.call()
    is.2D <- !is.null(dim(x)) && length(dim(x)) %in% 2
    if (!is.2D) {
        halt("Can only make a Crunch dataset from a two-dimensional data ",
            "structure")
    }

    if (is.data.frame(x) || nrow(x) > 1000000) {
        Call[[1]] <- as.name("newDatasetByCSV")
    } else {
        Call[[1]] <- as.name("newDatasetByColumn")
    }
    ds <- eval.parent(Call)
    invisible(ds)
}


##' Upload a data.frame column-by-column to make a new dataset
##'
##' Use this version if you have lots of variables, under 1M rows, perhaps
##' backed by ff or other memory-mapped files, and time to kill.
##'
##' @param x a data.frame or other rectangular R object
##' @param name character, the name to give the new Crunch dataset. Default is
##' the name of the R object passed in \code{x}
##' @param ... additional arguments passed to \code{ \link{createDataset}}
##' @return If successful, an object of class CrunchDataset.
##' @seealso \code{\link{newDataset}} \code{\link{newDatasetByCSV}}
##' @export
newDatasetByColumn <- function (x, name=as.character(substitute(x)), ...) {

    ds <- createDataset(name=name, ...)
    vardefs <- lapply(names(x),
        function (v) toVariable(x[[v]], name=v, alias=v))
    ds <- addVariables(ds, vardefs)
    invisible(ds)
}

##' Upload a file to Crunch to make a new dataset
##'
##' Use this import method if you have an SPSS data file. Reading such a file
##' into R as a data.frame will result in lost metadata. You can just send it
##' directly to Crunch and let the server process it.
##'
##' @param file character, the path to a file to upload. This should either be
##' a .csv or .sav (SPSS) file.
##' @param name character, the name to give the new Crunch dataset. Default is
##' the file name
##' @param ... additional arguments passed to \code{ \link{createDataset}}
##' @return On success, an object of class \code{CrunchDataset}.
##' @export
newDatasetFromFile <- function (file, name=basename(file), ...) {
    if (!file.exists(file)) {
        halt("File not found")
    }
    ds <- createDataset(name, ...)
    ds <- addSourceToDataset(ds, createSource(file))
    invisible(ds)
}

##' @importFrom httr upload_file
createSource <- function (file, ...) {
    crPOST(sessionURL("sources"),
        body=list(uploaded_file=upload_file(file)), ...)
}

##' Create an empty dataset
##'
##' Use only if you're writing a function to create a Crunch dataset from a
##' custom data structure. If you have a data.frame, just call
##' \code{\link{newDataset}} on it.
##'
##' @param name character, the name to give the new Crunch dataset. This is
##' required.
##' @param ... additional arguments for the POST to create the dataset, such as
##' "description".
##' @return An object of class CrunchDataset.
##' @seealso \code{\link{newDataset}}
##' @keywords internal
##' @export
createDataset <- function (name, ...) {
    dataset_url <- crPOST(sessionURL("datasets"),
        body=toJSON(list(name=name, ...)))
    updateDatasetList()
    ds <- entity(datasetCatalog()[[dataset_url]])
    invisible(ds)
}

addSourceToDataset <- function (dataset, source_url, ...) {
    batches_url <- shojiURL(dataset, "catalogs", "batches")
    body <- list(
        element="shoji:entity",
        body=list(
            source=source_url
        )
    )
    batch_url <- crPOST(batches_url, body=toJSON(body), ...)

    status <- try(pollBatchStatus(batch_url, batches(dataset),
        until=c("ready", "imported")))
    if (is.error(status)) {
        halt("Error importing file")
    } else if (status %in% "ready") {
        crPATCH(batch_url, body=toJSON(list(status="importing")))
        pollBatchStatus(batch_url, refresh(batches(dataset)))
    }

    invisible(refresh(dataset))
}

.delete_all_my_datasets <- function () {
    lapply(urls(datasetCatalog()), crDELETE)
}

##' Upload a data.frame to Crunch to make a new dataset
##'
##' This function uses the CSV+JSON import format, which is faster and more
##' effective for certain dataset sizes and shapes than
##' \code{\link{newDatasetByColumn}}.
##'
##' @param x a data.frame or other rectangular R object
##' @param name character, the name to give the new Crunch dataset. Default is
##' the name of the R object passed in \code{x}
##' @param ... additional arguments passed to \code{ \link{createDataset}}
##' @return If successful, an object of class CrunchDataset.
##' @seealso \code{\link{newDataset}} \code{\link{newDatasetByColumn}}
##' @importFrom utils write.csv
##' @export
newDatasetByCSV <- function (x, name=as.character(substitute(x)), ...) {

    ## Get all the things
    message("Processing the data")
    vars <- lapply(names(x),
        function (i) toVariable(x[[i]], name=i, alias=i))
    names(vars) <- names(x)

    ## Extract the data
    cols <- lapply(vars, function (v) v[["values"]])
    filename <- tempfile()
    gf <- gzfile(filename, "w")
    write.csv(cols, file=gf, na="", row.names=FALSE)
    close(gf)

    ## Drop the columns from the metadata and compose the payload
    vars <- lapply(vars, function (v) {
        v[["values"]] <- NULL
        return(v)
    })
    meta <- shojifyMetadata(vars, name=name, ...)

    ## Send to Crunch
    ds <- createWithMetadataAndFile(meta, filename)
    invisible(ds)
}

##' Make a dataset with metadata and a CSV
##'
##' This function just takes what you give it and POSTs to Crunch. No
##' validation, no automatic wrapping in the Shoji envelope, etc.
##'
##' @param metadata a list representation of the dataset's metadata, which
##' will be JSON-serialized and POSTed.
##' @param file a path to a CSV file, optionally zipped, that corresponds to
##' the above metadata.
##' @param strict logical: must the metadata exactly match the data? Default is
##' TRUE.
##' @param cleanup logical: if the file upload fails, delete the dataset?
##' Default is TRUE.
##' @return On success, a new dataset.
##' @export
##' @keywords internal
createWithMetadataAndFile <- function (metadata, file, strict=TRUE, cleanup=TRUE) {
    message("Uploading metadata")
    dataset_url <- crPOST(sessionURL("datasets"), body=toJSON(metadata))
    updateDatasetList()
    ds <- entity(datasetCatalog()[[dataset_url]])

    message("Uploading data")
    batches_url <- ds@catalogs$batches
    if (!strict) {
        batches_url <- paste0(batches_url, "?strict=0")
    }
    if (substr(file, 1, 5) == "s3://") {
        ## S3 upload
        batch <- try(crPOST(batches_url, body=toJSON(list(
            element="shoji:entity",
            body=list(
                url=file
            )))))
    } else {
        ## Local file. Send it as file upload
        batch <- try(crPOST(batches_url,
            body=list(file=httr::upload_file(file))))
    }
    if (is.error(batch) && cleanup) {
        delete(ds, confirm=FALSE)
        rethrow(batch)
    }
    message("Done!")
    return(refresh(ds))
}

# createWithMetadataAndFile <- function (metadata, file) {
#     ## Overwrite that and dump out the CSV and JSON, for debugging purposes
#     message("Writing metadata")
#     cat(toJSON(metadata), file="example.json")
#
#     message("Writing data")
#     file.copy(file, "example.csv.gz", overwrite=TRUE)
#     message("Done!")
#     print(c("example.json", "example.csv.gz"))
#     halt()
# }

##' Wrap variable metadata inside a dataset entity
##'
##' @param metadata list of variable metadata
##' @param order a valid "order" payload: list containing either aliases or
##' list(group, entities)
##' @param ... dataset entity metadata. "name" is required.
##' @param return list suitiable for JSONing and POSTing to create a dataset
##' @export
##' @keywords internal
shojifyMetadata <- function (metadata, order=list(list(group="ungrouped",
                            entities=I(names(metadata)))), ...) {
    return(list(element="shoji:entity",
                 body=list(...,
                           table=list(element="crunch:table",
                                      metadata=metadata,
                                      order=order))))
}
