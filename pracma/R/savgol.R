##
##  s a v g o l . R  Savitzky-Golay Smoothing
##


savgol <- function(T, fl, forder = 4, dorder = 0) {
    stopifnot(is.numeric(T), is.numeric(fl))
    if (fl <= 1 || fl %% 2 == 0)
        stop("Argument 'fl' must be an odd integer greater than 1.")
    n <- length(T)

    # -- calculate filter coefficients --
    fc <- (fl-1)/2                          # index: window left and right
    X <- outer(-fc:fc, 0:forder, FUN="^")   # polynomial terms and coeffs
    Y <- pinv(X);                           # pseudoinverse

    # -- filter via convolution and take care of the end points --
    T2 <- convolve(T, rev(Y[(dorder+1),]), type="o")   # convolve(...)
    T2 <- T2[(fc+1):(length(T2)-fc)]

    Tsg <- (-1)^dorder * T2
    return( Tsg )
}
