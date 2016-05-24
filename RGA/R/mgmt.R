# Remove some field and convert dates
#' @importFrom lubridate ymd_hms
fix_mgmt <- function(x) {
    x <- x[!grepl("(self|parent|child)Link", names(x))]
    if (!is.null(x$created))
        x$created <- ymd_hms(x$created)
    if (!is.null(x$updated))
        x$updated <- ymd_hms(x$updated)
    return(x)
}

# Get the Management API data
#' @include get-data.R
list_mgmt <- function(path, query, token) {
    json_content <- get_data(path, query, token)
    if (is.null(json_content$items) || length(json_content$items) == 0) {
        message("No results were obtained.")
        return(invisible(NULL))
    }
    res <- fix_mgmt(json_content$items)
    attr(res, "username") <- json_content$username
    return(res)
}

# Get the Management API data
#' @include url.R
#' @include request.R
get_mgmt <- function(path, token) {
    res <- api_request(get_url(path), token)
    res <- fix_mgmt(res)
    return(res)
}
