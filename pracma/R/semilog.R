##
##  s e m i l o g . R
##


semilogx <- function(x, y, ...) {
    plot(x, y, log = "x", ...)
    grid()
}


semilogy <- function(x, y, ...) {
    plot(x, y, log = "y", ...)
    grid()
}


loglog <- function(x, y, ...) {
    plot(x, y, log = "xy", ...)
    grid()
}
