# Get the Google Analytics API data
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom plyr rbind.fill.matrix rbind.fill
#' @include url.R
#' @include request.R
#' @include convert.R
get_data <- function(path = NULL, query = NULL, token) {
    # Set limits
    if (grepl("management", paste(path, collapse = "/"))) {
        limit <- 1000L
        items <- "items"
    } else {
        limit <- 10000L
        items <- "rows"
    }
    if (is.null(query$start.index))
        query$start.index <- 1L
    total <- ifelse(!is.null(query$max.results) && query$max.results > limit, query$max.results, NA)
    if (is.null(query$max.results) || query$max.results > limit) {
        query$max.results <- limit
        pagination <- getOption("rga.pagination", TRUE)
    } else
        pagination <- FALSE
    if (grepl("data/realtime", paste(path, collapse = "/")))
        pagination <- FALSE
    # Make request
    res <- api_request(get_url(path, query), token)
    if (res$totalResults == 0L || is.null(res[[items]]) || length(res[[items]]) == 0L)
        return(NULL)
    # Pagination
    if (isTRUE(pagination) && query$max.results < res$totalResults) {
        message(sprintf("API response contains more then %d items. Batch processing mode enabled.", query$max.results))
        if (is.na(total) || total >= res$totalResults)
            total <- res$totalResults
        sidx <- seq.int(query$start.index, total, query$max.results)
        midx <- diff(c(sidx, total + 1L))
        pages <- vector(mode = "list", length = length(sidx))
        pb <- txtProgressBar(min = 0, max = length(sidx), initial = 1, style = 3)
        for (i in 2L:length(sidx)) {
            query$start.index <- sidx[i]
            query$max.results <- midx[i]
            pages[[i]] <- api_request(get_url(path, query), token)[[items]]
            setTxtProgressBar(pb, i)
        }
        pages[[1L]] <- res[[items]]
        if (is.matrix(pages[[1L]]))
            pages <- rbind.fill.matrix(pages)
        if (is.data.frame(pages[[1L]]))
            pages <- rbind.fill(pages)
        else if (is.list(pages[[1L]]))
            pages <- unlist(pages, recursive = FALSE, use.names = FALSE)
        res[[items]] <- pages
        close(pb)
    }
    res[[items]] <- build_df(res)
    if (nrow(res[[items]]) < res$totalResults)
        warning(sprintf("Only %d observations out of %d were obtained. Set max.results = NULL (default value) to get all results.", nrow(res[[items]]), res$totalResults), call. = FALSE)
    return(res)
}
