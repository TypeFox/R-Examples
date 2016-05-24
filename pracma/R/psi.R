##
##  p s i . R  Psi (Polygamma) Function
##


psi <- function(k, z) {
    if (missing(z)) {
        z <- k; k <- 0
    }
    stopifnot(is.numeric(z) || is.complex(z))
    if (length(k) > 1 || k < 0)
        stop("Argument 'k': Invalid Polygamma order, or 'k' not a scalar.")
    k <- floor(k)

    sz <- dim(z)
    zz <- z <- c(z)
    f <- 0.0 * z

    if (k == 0) {
        # reflection point
        p <- which(Real(z) < 0.5)
        if (length(p) > 0)
            z[p] <- 1 - z[p]
    
        # Lanczos approximation for the complex plane
        g <- 607/128  # best results when 4 <= g <= 5
    
        cc <-  c(0.99999999999999709182, 57.156235665862923517, -59.597960355475491248,
                14.136097974741747174, -0.49191381609762019978, 0.33994649984811888699e-4,
                0.46523628927048575665e-4, -0.98374475304879564677e-4, 0.15808870322491248884e-3,
               -0.21026444172410488319e-3, 0.21743961811521264320e-3, -0.16431810653676389022e-3,
                0.84418223983852743293e-4, -0.26190838401581408670e-4, 0.36899182659531622704e-5)
    
        n <- d <- 0
        for (j in length(cc):2) {
            dz <- 1/(z+j-2)
            dd <- cc[j] * dz
            d  <- d + dd
            n  <- n - dd * dz
        }
        d  <- d + cc[1]
        gg <- z + g - 0.5
    
        # log is accurate to about 13 digits...
        f <- log(gg) + (n/d - g/gg)
    
        if (length(p) > 0)
           f[p] <- f[p] - pi*cot(pi*zz[p])
    
        p <- which(round(zz) == zz && Real(zz) <= 0 && Imag(zz) == 0)
        if (length(p) > 0)
           f[p] <- Inf

    } else if (k > 0) {
        isneg <- which(Real(z) < 0)
        isok  <- which(Real(z) >= 0)
        n <- k

        negmethod <- 1
        if (length(isneg) > 0) {
           if (negmethod == 0) {
              zneg <- z[isneg]
              gneg <- psi(n, zneg+1)  # recurse if to far to the left...
              hneg <- -(-1)^n * gamma(n+1) * zneg^(-(n+1))
              fneg <- gneg + hneg
           } else {
              zneg <- z[isneg]  # shift by, say, 500, to speed things up
              m <- 500
              gneg <- psi(n, zneg+m)
              hneg <- 0
              for (k in (m-1):0)
                  hneg <- hneg + (zneg+k)^(-(n+1))
              hneg <- -(-1)^n * gamma(n+1) * hneg
              fneg <- gneg + hneg
           }
        }
        if (length(isok) > 0)
           z <- z[isok]

        # the zeros of the Lanczos PFE series when g=607/128 are:
        r <- c(-4.1614709798720630 - 0.14578107125196249*1i,
               -4.1614709798720630 + 0.14578107125196249*1i,
               -4.3851935502539474 - 0.19149326909941256*1i,
               -4.3851935502539474 + 0.19149326909941256*1i,
               -4.0914355423005926, -5.0205261882982271,
               -5.9957952053472399, -7.0024851819328395,
               -7.9981186370233868, -9.0013449037361806,
               -9.9992157162305535,-11.0003314815563886,
              -11.9999115102434217,-13.0000110489923175587)
        # the poles of the Lanczos PFE series are:
        # p <- c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13)
        e <- exp(1) 
        g <- 607/128  # best results when 4<=g<=5
        h <- 1/2

        s <- 0
        for (k in (length(r)-1):0)
            s <- s + (1/((z-r[k+1])^(n+1)) - 1/((z+k)^(n+1)))
        # what happens if n is not a positive integer?
        s <- (-1)^n * gamma(n+1) * s

        zgh <- z + (g-h)
        if (n == 0) {
        #  s=log(zgh)+(-g./zgh + s);
        #  use existing more accurate digamma function if n=0
        #  should never reach this code since we trapped it above 
           f <- psi(z)
        } else {
        # do derivs of front end stuff
           s <- (-1)^(n+1) * (gamma(n) *     zgh^(-n) + 
                              g*gamma(n+1) * zgh^-(n+1)) + s
        }

        if (length(isneg) > 0)
           f[isneg] <- fneg
        if (length(isok) > 0)
           f[isok] <- s

    } else
        stop("Argument 'k': Invalid Polygamma order.")

    if (is.numeric(z)) f <- Re(f)
    dim(f) <- sz

    return(f)
}
