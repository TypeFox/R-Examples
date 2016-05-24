date_ranges <- function(start, end, by) {
    start <- as.Date(start)
    end <- as.Date(end)
    by <- match.arg(by, c("day", "week", "month", "quarter", "year"))
    dates <- seq.Date(start, end, by = by)
    res <- cbind(start = as.character(dates),
                 end = as.character(c(dates[-1] - 1, end)))
    as.data.frame(res, stringsAsFactors = FALSE)
}

parse_date <- function(x) {
    stopifnot(is.character(x))
    if (grepl("daysAgo", x, fixed = TRUE))
        x <- Sys.Date() - as.numeric(sub("^([0-9]+).*", "\\1", x))
    else if (x == "today")
        x <- Sys.Date()
    else if (x == "yesterday")
        x <- Sys.Date() - 1
    return(as.character(x))
}

#' @importFrom utils txtProgressBar setTxtProgressBar capture.output
#' @importFrom plyr rbind.fill
#' @importFrom stats aggregate.data.frame
#' @include get-data.R
#' @include utils.R
fetch_by <- function(path, query, by, token) {
    query$start.date <- parse_date(query$start.date)
    query$end.date <- parse_date(query$end.date)
    dates <- date_ranges(query$start.date, query$end.date, by)
    n <- nrow(dates)
    message("Batch processing mode enabled.\n",
            sprintf("Fetch data by %s: from %s to %s.", by, query$start.date, query$end.date))
    pages <- vector(mode = "list", length = n)
    pb <- txtProgressBar(min = 0, max = n, initial = 0, style = 3)
    for (i in 1:n) {
        query$start.date <- dates$start[i]
        query$end.date <- dates$end[i]
        suppressMessages(capture.output(pages[[i]] <- get_data(path, query, token)))
        setTxtProgressBar(pb, i)
    }
    cl <- vapply(pages, is.null, logical(1))
    if (all(cl))
        return(NULL)
    else if (any(cl))
        pages <- pages[!cl]
    res <- pages[[1]]
    res$rows <- rbind.fill(lapply(pages, .subset2, "rows"))
    names(res$query) <- rename_params(names(res$query))
    res$query$start.date <- pages[[1]]$query$`start-date`
    res$query$end.date <- pages[[length(pages)]]$query$`end-date`
    res$totalResults <- sum_by(pages, "totalResults")
    res$containsSampledData <- any(unlist(lapply(pages, .subset2, "containsSampledData")))
    if (isTRUE(res$containsSampledData)) {
        res$sampleSize <- sum_by(pages, "sampleSize")
        res$sampleSpace <- sum_by(pages, "sampleSpace")
    }
    if (is.null(query$dimensions))
        res$rows <- as.data.frame(as.list(colSums(res$rows)))
    else if (!is.null(query$dimensions) && !any(grepl("date", query$dimensions))) {
        mets <- parse_params(query$metrics)
        dims <- parse_params(query$dimensions)
        res$rows <- aggregate.data.frame(res$rows[mets], res$rows[dims], sum)
    }
    close(pb)
    if (grepl("ga:users|ga:[0-9]+dayUsers", query$metrics))
        warning("The 'ga:users' or 'ga:NdayUsers' total value for several days is not the sum of values for each single day.", call. = FALSE)
    return(res)
}
