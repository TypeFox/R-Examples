is.hypoexp <- function(x) {
    UseMethod('is.hypoexp')
}

is.hypoexp.default <- function(x) {
    if (!is.list(x)) {
        x <- as.list(x)
    }
    'hypoexp' %in% names(x)
}

is.hypoexp.claiminfo <- function(x) {
    'hypoexp' %in% names(x)
}

is.hypoexp.riskproc <- function(x) {
    is.hypoexp(x[['claims']])
}
