##
##  f f t . R  Fourier Transform
##


ifft <- function(x) {
    if (length(x) == 0)
        return(c())
    if ( (!is.vector(x, mode="numeric") && !is.vector(x, mode="complex")))
        stop("Argument 'x' must be real or complex vector.")
    
    fft(x, inverse = TRUE) / length(x)
}


fftshift <- function(x) {
    stopifnot(is.double(x) || is.complex(x) || is.integer(x))
    if (!is.vector(x))
        stop("Argument 'x' must be a real or complex vector.")

    m <- length(x)
    p <- ceiling(m/2)
    idx <- c((p+1):m, 1:p)
    x[idx]
}


ifftshift <- function(x) {
    stopifnot(is.double(x) || is.complex(x) || is.integer(x))
    if (!is.vector(x))
        stop("Argument 'x' must be a real or complex vector.")

    m <- length(x)
    p <- floor(m/2)
    idx <- c((p+1):m, 1:p)
    x[idx]
}
