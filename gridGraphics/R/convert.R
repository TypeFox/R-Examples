
# Handle log axes

tx <- function(x, par) {
    if (par$xlog && !is.null(x)) {
        log10(x)
    } else {
        x
    }
}

ty <- function(x, par) {
    if (par$ylog && !is.null(x)) {
        log10(x)
    } else {
        x
    }
}
