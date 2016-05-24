spl.coef.conv = function (z)
{
    n <- length(z$x)
    h <- z$x[2:n] - z$x[1:(n - 1)]
    y0 <- z$y[1:(n - 1)]
    y1 <- z$y[2:n]
    b0 <- z$b[1:(n - 1)]
    b1 <- z$b[2:n]
    cc <- -(3 * (y0 - y1) + (2 * b0 + b1) * h)/h^2
    c1 <- (3 * (y0[n - 1] - y1[n - 1]) + (b0[n - 1] + 2 * b1[n -
          1]) * h[n - 1])/h[n - 1]^2
    dd <- (2 * (y0 - y1)/h + b0 + b1)/h^2
    d1 <- dd[n - 1]
    z$c <- c(cc, c1)
    z$d <- c(dd, d1)
    z
}
