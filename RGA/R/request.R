# Error printing function
#' @importFrom httr http_status
#' @include utils.R
error_reasons <- function(x) {
    code_message <- http_status(as.numeric(x$error$code))$message
    errors <- x$error$errors
    errors$reason <- gsub("^([[:alpha:]])", "\\U\\1", to_separated(errors$reason, sep = " "), perl = TRUE)
    errors$message <- gsub("\n", ". ", errors$message, fixed = TRUE)
    res <- sprintf("%s: %s", errors$reason, errors$message)
    if (!is.null(errors$location)) {
        errors$location <- rename_params(errors$location)
        idx_inv <- grep("Invalid parameter", errors$reason)
        idx_req <- unique(grep("Required", errors$reason), grep("parameter", errors$locationType))
        if (length(idx_inv))
            res[idx_inv] <- sprintf("%s '%s': %s", errors$reason[idx_inv], errors$location[idx_inv], errors$message[idx_inv])
        if (length(idx_req))
            res[idx_req] <- sprintf("%s %s: '%s'", errors$reason[idx_req], errors$locationType[idx_req], errors$location[idx_req])
    }
    paste(c(code_message, res), collapse = "\n")
}

# Process response
#' @importFrom httr status_code http_error content parse_media
#' @importFrom jsonlite fromJSON
#' @include utils.R
process_response <- function(response) {
    stopifnot(inherits(response, "response"))
    status_code <- status_code(response)
    if (status_code == 204L)
        return(NULL)
    if (!http_error(response)) {
        text <- content(response, as = "text")
        if (text == "")
            stop("No output to parse.", call. = FALSE)
        res <- fromJSON(text, flatten = TRUE)
    } else {
        if (status_code == 404L)
            stop(sprintf("The requested URL not found. URL: %s.", strsplit(response$url, "?", fixed = TRUE)[[1L]][1L]), call. = FALSE)
        type <- parse_media(response$headers$`Content-type`)
        if (type$complete == "application/json") {
            res <- fromJSON(content(response, as = "text"))
            stop(error_reasons(res), call. = FALSE)
        } else {
            res <- content(response, as = "text")
            stop(sprintf("HTTP error %s:\n%s.", status_code, res), call. = FALSE)
        }
    }
    return(res)
}

# Get a Google Analytics API response
#' @importFrom httr config GET accept_json
#' @importFrom stats runif
#' @include auth.R
api_request <- function(url, token) {
    if (missing(token) && is.null(get_token()))
        api_request(url, token = authorize(cache = FALSE))
    if (missing(token) && !is.null(get_token()))
        token <- get_token()
    if (validate_token(token))
        config <- config(token = token)
    attempts <- getOption("rga.retry.attempts", 5L) + 1L
    for (i in 0L:attempts) {
        response <- GET(url, config = config, accept_json())
        res <- tryCatch(process_response(response), error = identity)
        if (!inherits(res, "error"))
            break
        else if (grepl("User rate limit exceeded|Quota exceeded", res$message) & i < attempts)
            Sys.sleep(2L^i + runif(1L))
        else
            stop(res)
    }
    return(res)
}
