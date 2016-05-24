
#' @importFrom utils URLencode
url_encode_request <- function(url, ...) {
    validate_arg(url, required = TRUE, class = 'character', length_min = 1, length_max = 1)

    args <- list(...)
    args <- args[!sapply(args, is.null)]
    args <- sapply(args, as.character, simplify = FALSE)
    args <- sapply(args, utils::URLencode, simplify = FALSE)
    args <- paste(mapply(paste, sep = '=', names(args), args), collapse = '&')

    paste(url, args, sep = '?')
}
