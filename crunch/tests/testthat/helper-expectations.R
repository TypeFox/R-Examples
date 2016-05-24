object_sort <- function (x) {
    if (is.list(x)) {
        x <- as.list(x) ## For S4 subclasses
        if (!is.null(names(x))) {
            x <- x[sort(names(x))]
        }
        return(lapply(x, object_sort))
    }
    return(x)
}

expect_json_equivalent <- function (object, expected, ...) {
    expect_equivalent(object_sort(object), object_sort(expected), ...)
}

expect_output <- function (object, ...) {
    testthat::expect_output(print(object), ...)
}
