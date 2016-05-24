beamDiskOverlap <-
function(s=0, r=1, theta.fwhm=1){
    len <- length(s)
    val <- numeric(len)
    alpha <- log(2)/theta.fwhm^2
    f <- function(x, s, alpha) exp(-alpha*x^2) * besselI(2*alpha*s*x, 0) * x
    for (i in 1:len) {
        val[i] <- 2*alpha*exp(-alpha*s[i]^2)/(1 - exp(-alpha*r^2))
        val[i] <- val[i] * integrate(f, 0, r, s=s[i], alpha=alpha)$value
    }
    val
}

