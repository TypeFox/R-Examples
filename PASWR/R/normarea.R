normarea <-
function (lower = -Inf, upper = Inf, m, sig)
{
    Altblue <- "#CDCDED"
    Fontcol <-"#3333B3"
    par(mar=c(4,1,4,1))
    area <- pnorm(upper, m, sig) - pnorm(lower, m, sig)
    ra <- round(area,4)
    x <- seq(m - 4 * sig, m + 4 * sig, length = 1000)
    y <- dnorm(x, m, sig)
    par(pty = "m")
    plot(x, y, type = "n", xaxt = "n", yaxt = "n", xlab = "",
        ylab = "",main="")

mtext(substitute(paste("The area between ",lower," and ",upper," is ",ra)),
side=3,line=1,font=2,cex=1.15)
mtext(substitute(paste("X~Normal (" ,mu==m,", ",sigma==sig,")" ),
list(m=m,sig=sig)),side=1,line=3,col=Fontcol)


    if (lower == -Inf || lower < m - 4 * sig) {
        lower <- m - 4 * sig
    }
    if (upper == Inf || upper > m + 4 * sig) {
        upper <- m + 4 * sig
    }
    axis(1, at = c(m, lower, upper), labels = c(m, lower, upper))
    xaxis1 <- seq(lower, upper, length = 200)
    yaxis1 <- dnorm(xaxis1, m, sig)
    xaxis1 <- c(lower, xaxis1, upper)
    yaxis1 <- c(0, yaxis1, 0)
    polygon(xaxis1, yaxis1, density = -1, col = Altblue)
    lines(x, y, lwd = 2)
    lines(c(m - 4 * sig, m + 4 * sig), c(0, 0), lwd = 2)
    lines(c(lower, lower), c(0, dnorm(lower, m, sig)), lwd = 2)
    lines(c(upper, upper), c(0, dnorm(upper, m, sig)), lwd = 2)
    par(mar=c(5.1,4.1,4.1,2.1))
}

