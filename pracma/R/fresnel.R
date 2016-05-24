##
##  f r e s n e l . R  Fresnel Integrals
##


fresnelS <- Vectorize(function(x) .fresnel(x)$S)

fresnelC <- Vectorize(function(x) .fresnel(x)$C)


.fresnel <- function(x) {
    if (!is.numeric(x) || length(x) != 1)
        stop("Argument 'x' must be a numeric scalar.")

    eps <- .Machine$double.eps
    xa <- abs(x)
    px <- pi * xa
    t <- 0.5 * px * xa
    t2 <- t^2

    fc <- fs <- 0
    if (xa == 0) {
        return(list(C = fc, S = fs))

    } else if (xa < 2.5) {
        r <- xa
        fc <- r
        for (k in 1:50) {
            r <- -.5 *r *(4.0*k-3.0)/k/(2.0*k-1.0)/(4.0*k+1.0)*t2
            fc <- fc + r
            if (abs(r) < abs(fc)*eps) break
        }
        fs <- xa*t/3.0
        r <- fs
        for (k in 1:50) {
            r <- -0.5*r*(4.0*k-1.0)/k/(2.0*k+1.0)/(4.0*k+3.0)*t2
            fs <- fs + r
            if (abs(r) < abs(fs)*eps) break
        }

    } else if (xa < 4.5) {
        m <- trunc(42.0+1.75*t)
        su <- 0.0
        f1 <- 0.0
        f0 <- 1.0e-100
        for (k in m:0) {
            f <- (2.0*k+3.0)*f0/t - f1
            if (k == trunc(k/2)*2) {
                fc <- fc + f
            } else {
                fs <- fs + f
            }
            su <- su+(2.0*k+1.0)*f*f
            f1 <- f0
            f0 <- f
        }
        q <- sqrt(su)
        fc <- fc*xa/q
        fs <- fs*xa/q

    } else {
        r=1.0
        f=1.0
        for (k in 1:20) {
            r <- -0.25*r*(4.0*k-1.0)*(4.0*k-3.0)/t2
            f <- f + r
        }
        k <- 20+1
        r <- 1.0/(px*xa)
        g <- r
        for (k in 1:12) {
            r <- -0.25*r*(4.0*k+1.0)*(4.0*k-1.0)/t2
            g <- g + r
        }
        k <- 12+1
        t0 <- t-trunc(t/(2.0*pi))*2.0*pi
        fc <- 0.5+(f*sin(t0)-g*cos(t0))/px
        fs <- 0.5-(f*cos(t0)+g*sin(t0))/px
        
    }

    if (x < 0) {
        fc <- -fc
        fs <- -fs
    }
    return(list(C = fc, S = fs))
}
