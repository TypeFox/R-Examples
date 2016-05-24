##
##  h u m p s . R  Matlab Test Function
##


humps <- function(x) {
    if (missing(x))
        x <- seq(0.0, 1.0, by=0.05)

    stopifnot(is.numeric(x))

    1/((x-0.3)^2 + 0.01) + 1/((x-0.9)^2 + 0.04) - 6
}

sinc <- function(x) {
    stopifnot(is.numeric(x))

    sin(pi * x) / (pi * x)
}
