pollBatchStatus <- function (batch.url, catalog, until="imported", wait=1) {

    ## Configure polling interval. Will increase by rate (>1) until reaches max
    max.wait <- 30
    increase.by <- 1.2

    starttime <- Sys.time()
    timeout <- crunchTimeout()
    timer <- function (since, units="secs") {
        difftime(Sys.time(), since, units=units)
    }
    status <- catalog[[batch.url]]$status
    while (status %in% c("idle", "importing", "analyzing", "appended") && timer(starttime) < timeout) {
        Sys.sleep(wait)
        catalog <- refresh(catalog)
        status <- catalog[[batch.url]]$status
        wait <- min(max.wait, wait * increase.by)
    }

    if (status %in% "idle") {
        halt("Append process failed to start on the server")
    } else if (status %in% c("analyzing", "importing", "appended")) {
        halt("Timed out. Check back later. Consider also increasing options(crunch.timeout)")
    } else if (status == "error") {
        halt("There was an error appending the datasets. Please contact support@crunch.io")
    } else if (status %in% c(until, "conflict")) {
        return(status)
    } else {
        halt(status)
    }
}

crunchTimeout <- function () {
    opt <- getOption("crunch.timeout")
    if (is.null(opt) || !is.numeric(opt)) opt <- 60
    return(opt)
}

flattenConflicts <- function (x) {
    ## x is list, keys are variable URLs, objects are arrays of objects with keys "message" and "resolution" (in R-speak, s/array/list/, s/object/list/, s/keys/names/)

    ## Filter out variables that have empty conflicts, if any
    x <- Filter(function (a) length(a$conflicts), x)

    if (length(x) == 0) {
        ## Bail, but return the right shape
        return(data.frame(message=c(), resolution=c(), name=c(), url=c()))
    }

    ## flatten object to data.frame with url, message, resolution
    dfconflicts <- function (clist) {
        data.frame(message=clist$message,
            resolution=clist$resolution %||% NA_character_,
            stringsAsFactors=FALSE)
    }

    out <- mapply(function (i, d) {
        df <- do.call(rbind, lapply(d$conflicts, dfconflicts))
        df$url <- i
        meta <- d$source_metadata %||% d$metadata
        df$name <- meta$name %||% meta$references$name
        return(df)
    }, i=names(x), d=x, SIMPLIFY=FALSE)
    return(do.call(rbind, out))
}

formatConflicts <- function (flat) {
    ## flat is the output of flattenConflicts
    if (nrow(flat)) {
        x <- groupConflicts(flat)
        return(vapply(x, formatConflictMessage, character(1), USE.NAMES=FALSE))
    } else {
        return("No conflicts.")
    }
}

groupConflicts <- function (flat) {
    ## reshape conflicts to be by conflict-resolution, not by variable
    ## flat is the output of flattenConflicts

    ## split by message, then by resolution
    return(unlist(lapply(split(flat, flat$message),
        function (x) split(x, x$resolution)), recursive=FALSE))
}

formatConflictMessage <- function (x) {
    ## receives a data.frame with variable URLs and common conflict and resolution
    conflict <- paste("Conflict:", unique(x$message))
    resolution <- paste("Resolution:", unique(x$resolution))
    ## those should be length 1 by construction

    varnames <- x$name %||% x$url
    vars <- paste0(nrow(x), " variable", ifelse(nrow(x) > 1, "s", ""), ": ",
        serialPaste(dQuote(unique(varnames))))
    return(paste(conflict, resolution, vars, sep="; "))
}

formatFailures <- function (flat) {
    ## Keep only the conflicts without resolution:
    flat <- flat[is.na(flat$resolution),,drop=FALSE]
    failures <- split(flat, flat$message)
    return(vapply(failures, formatFailureMessage, character(1),
        USE.NAMES=FALSE))
}

formatFailureMessage <- function (x) {
    ## receives a data.frame with variable URLs and common conflict and resolution (which is NA)
    conflict <- paste("Critical conflict:", unique(x$message))
    varnames <- x$name %||% x$url
    vars <- paste0(nrow(x), " variable", ifelse(nrow(x) > 1, "s", ""), ": ",
        serialPaste(dQuote(unique(varnames))))
    return(paste(conflict, vars, sep="; "))
}
