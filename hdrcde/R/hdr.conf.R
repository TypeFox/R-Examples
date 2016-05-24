hdrconf <- function(x,den,prob=95,conf=95)
{
    # Returns hdr with confidence limits
    # Assumes x is a sorted sample from the density den.

    alpha <- 1-conf/100
    info <- calc.falpha(x,den,alpha)
    hdr <- hdr.ends(den,info$falpha)$hdr
    nint <- length(hdr)
    delta <- range(x)/1000
    f2 <- approx(den$x,den$y,hdr+delta)$y
    fprime <- (f2-info$falpha)/delta
    g <- info$falpha*sum(abs(1/fprime))
    var.falpha <- alpha*(1-alpha)/(length(x)*g*g)
    z <- abs(qnorm(0.5-conf/200))
    falpha.ci <- z*sqrt(var.falpha)
    falpha.ci <- info$falpha + c(-falpha.ci,falpha.ci)
    if(falpha.ci[1] < 0)
        hdr1 <- c(NA,NA)
    else
        hdr1 <- hdr.ends(den,falpha.ci[1])$hdr
    if(falpha.ci[2] < 0)
        hdr2 <- c(NA,NA)
    else
        hdr2 <- hdr.ends(den,falpha.ci[2])$hdr
    return(structure(list(hdr=hdr, hdr.lo=hdr1, hdr.hi=hdr2,
            falpha=info$falpha, falpha.ci=falpha.ci),class="hdrconf"))
}
