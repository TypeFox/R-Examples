bits <- function(x, k = 54, pos_sign = FALSE, break0 = FALSE) {
    stopifnot(is.numeric(x), length(x) == 1)
    if (ceiling(k) != floor(k) || k <= 0)
        stop("Argument 'k' must be a (positive) integer.")
  
    if (x >= 0) {
        b <- if (pos_sign) "+" else ""
    } else {
        b <- "-"
        x <- -x
    }
      
    xn <- trunc(x)
    xf <- x - xn
      
    if (xn >= 1) {
        m2 <- nextpow2(xn)
        if (2^m2 > xn) m2 <- m2 - 1
        for (i in seq(m2, 0, by=-1)) {
            s <- 2^i
            if (xn >= s) {
                b <- paste(b, "1", sep="")
                xn <- xn - s
            } else
                b <- paste(b, "0", sep="")
        }
    } else {
          b <- paste(b, "0", sep="")
    }
      
    if (xf > 0) {
        b <- paste(b, ".", sep="")
        for (i in 1:k) {
            s <- 1/2^i
            if (xf >= s) {
                b <- paste(b, "1", sep="")
                xf <- xf - s
                if (break0 && xf == 0) break
            } else
                b <- paste(b, "0", sep="")
        }
    }
      
    return(b)
}

