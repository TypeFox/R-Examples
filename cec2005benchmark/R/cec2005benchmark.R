cec2005benchmark <- function (i, x) {
    if (is.numeric(i) && i >= 1 && i <= 25) {
        if (is.vector(x)) {
            row <- 1; col <- length(x)
        } else if (is.matrix(x)) {
            row <- nrow(x); col <- ncol(x)
        } else {
            stop("x should be a vector or a matrix")
        }
        if (!(col %in% c(2, 10, 30, 50))) {
            stop("only 2, 10, 30, 50 variables are supported")
        }
        extdatadir <- system.file("extdata", package = "cec2005benchmark")
        f <- .C("cec2005benchmark", extdatadir = as.character(extdatadir), 
                i = as.integer(i), x = as.double(x), row = as.integer(row),
                col = as.integer(col), f = double(row), 
                PACKAGE = "cec2005benchmark")$f
    } else {
        stop("i should be an integer between 1 and 25")
    }
    return(f)
}

cec2005benchmark1 <- function (x) cec2005benchmark(1, x)
cec2005benchmark2 <- function (x) cec2005benchmark(2, x)
cec2005benchmark3 <- function (x) cec2005benchmark(3, x)
cec2005benchmark4 <- function (x) cec2005benchmark(4, x)
cec2005benchmark5 <- function (x) cec2005benchmark(5, x)
cec2005benchmark6 <- function (x) cec2005benchmark(6, x)
cec2005benchmark7 <- function (x) cec2005benchmark(7, x)
cec2005benchmark8 <- function (x) cec2005benchmark(8, x)
cec2005benchmark9 <- function (x) cec2005benchmark(9, x)
cec2005benchmark10 <- function (x) cec2005benchmark(10, x)
cec2005benchmark11 <- function (x) cec2005benchmark(11, x)
cec2005benchmark12 <- function (x) cec2005benchmark(12, x)
cec2005benchmark13 <- function (x) cec2005benchmark(13, x)
cec2005benchmark14 <- function (x) cec2005benchmark(14, x)
cec2005benchmark15 <- function (x) cec2005benchmark(15, x)
cec2005benchmark16 <- function (x) cec2005benchmark(16, x)
cec2005benchmark17 <- function (x) cec2005benchmark(17, x)
cec2005benchmark18 <- function (x) cec2005benchmark(18, x)
cec2005benchmark19 <- function (x) cec2005benchmark(19, x)
cec2005benchmark20 <- function (x) cec2005benchmark(20, x)
cec2005benchmark21 <- function (x) cec2005benchmark(21, x)
cec2005benchmark22 <- function (x) cec2005benchmark(22, x)
cec2005benchmark23 <- function (x) cec2005benchmark(23, x)
cec2005benchmark24 <- function (x) cec2005benchmark(24, x)
cec2005benchmark25 <- function (x) cec2005benchmark(25, x)
