##
##  l a m b e r t W . R  Lambert W Function
##


lambertWp <- function(x) {
    if (!is.numeric(x))
        stop("Argument 'x' must be a numeric (real) vector.")

    if (length(x) == 1) {
        if (x  < -1/exp(1)) return(NaN)
        if (x == -1/exp(1)) return(-1)

        # compute first iteration of $W_0$
	    if (x <= 1) {
            eta <- 2 + 2*exp(1)*x;
            f2 <- 3*sqrt(2) + 6 - (((2237+1457*sqrt(2))*exp(1) - 4108*sqrt(2) - 5764)*sqrt(eta)) /
                                  ((215+199*sqrt(2))*exp(1) - 430*sqrt(2)-796)
            f1 <- (1-1/sqrt(2))*(f2+sqrt(2));
            w0 <- -1 + sqrt(eta)/(1 + f1*sqrt(eta)/(f2 + sqrt(eta)));
        } else {
            w0 = log( 6*x/(5*log( 12/5*(x/log(1+12*x/5)) )) )
        }

        # w0 <- 1
        w1 <- w0 - (w0*exp(w0)-x)/((w0+1)*exp(w0)-(w0+2)*(w0*exp(w0)-x)/(2*w0+2))
        while(abs(w1-w0) > 1e-15) {
            w0 <- w1
            w1 <- w0 - (w0*exp(w0)-x)/((w0+1)*exp(w0)-(w0+2)*(w0*exp(w0)-x)/(2*w0+2))
        }
        return(w1)
    } else {
        sapply(x, lambertWp)
    }
}


lambertWn <- function(x) {
    if (!is.numeric(x))
        stop("Argument 'x' must be a numeric (real) scalar.")

    if (length(x) == 1) {
        if (x < -exp(-1) || x >= 0) return(NaN)
        if (x == exp(-1)) return(-1)

        # compute first approximation of $W_-1$
        m1 <-  0.3361; m2 <- -0.0042; m3 <- -0.0201
        sigma <- -1 - log(-x)
        w0 <- -1 - sigma - 2/m1*(1-1/(1 + (m1*sqrt(sigma/2)) / 
                                 (1 + m2*sigma*exp(m3*sqrt(sigma)))))
        r <- abs(x - w0*exp(w0))
        while (r > 1e-15) {
            w1 <- w0 - (w0*exp(w0)-x) /
                       (exp(w0)*(w0+1)-(w0+2)*(w0*exp(w0)-x)/(2*w0+2))
            r <- abs(x - w1*exp(w1) )
            w0 <- w1
        }
        return(w0)
    } else {
        sapply(x, lambertWn)
    }
}

