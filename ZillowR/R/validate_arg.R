
#' @importFrom methods is
validate_arg <- function(x,
    required = FALSE,
    class = NULL,
    length_min = NULL,
    length_max = NULL,
    inclusion = NULL,
    exclusion = NULL,
    format = NULL,
    value_min = NULL,
    value_max = NULL
) {
    y <- c()

    if (!is.null(required) && length(required) == 1L && required %in% TRUE) {
        if (is.null(x)) y <- c(y, sprintf("'%s' is required", deparse(substitute(x))))
    }

    if (!is.null(class) && length(class) == 1L && is.character(class)) {
        if (!is.null(x) && !methods::is(x, class)) y <- c(y, sprintf("'%s' must be %s", deparse(substitute(x)), class))
    }

    if (!is.null(length_min) && length(length_min) == 1L && is.numeric(length_min)) {
        if (!is.null(x) && length(x) < length_min) y <- c(y, sprintf("'%s' must be at least length %i", deparse(substitute(x)), length_min))
    }

    if (!is.null(length_max) && length(length_max) == 1L && is.numeric(length_max)) {
        if (!is.null(x) && length(x) > length_max) y <- c(y, sprintf("'%s' must be no more than length %i", deparse(substitute(x)), length_max))
    }

    if (!is.null(inclusion)) {
        if (!is.null(x) && !all(x %in% inclusion)) y <- c(y, sprintf("'%s' must be one of: %s", deparse(substitute(x)), paste(inclusion, collapse = ', ')))
    }

    if (!is.null(exclusion)) {
        if (!is.null(x) && any(x %in% exclusion)) y <- c(y, sprintf("'%s' must not be any of: %s", deparse(substitute(x)), paste(exclusion, collapse = ', ')))
    }

    if (!is.null(format) && length(format) == 1L && is.character(format)) {
        if (!is.null(x) && !all(grepl(format, x))) y <- c(y, sprintf("'%s' must match: %s", deparse(substitute(x)), format))
    }

    if (!is.null(value_min) && length(value_min) == 1L && is.numeric(value_min)) {
        if (!is.null(x) && !all(x >= value_min, na.rm = TRUE)) y <- c(y, sprintf("'%s' values must be greater than or equal to %i", deparse(substitute(x)), value_min))
    }

    if (!is.null(value_max) && length(value_max) == 1L && is.numeric(value_max)) {
        if (!is.null(x) && !all(x <= value_max, na.rm = TRUE)) y <- c(y, sprintf("'%s' values must be less than or equal to %i", deparse(substitute(x)), value_max))
    }

    return(y)
}
