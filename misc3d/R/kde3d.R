## Based on kde2d in MASS.
kde3d <- function (x, y, z, h, n = 20, lims = c(range(x), range(y), range(z))) 
{
    nx <- length(x)
    if (length(y) != nx || length(z) != nx) 
        stop("data vectors must be the same length")
    if (missing(h)) 
        h <- c(MASS::bandwidth.nrd(x),
               MASS::bandwidth.nrd(y),
               MASS::bandwidth.nrd(z)) / 6
    else if (length(h) != 3)
        h <- rep(h, length = 3)
    if (length(n) != 3)
        n <- rep(n, length = 3)
    if (length(lims) == 2)
        lims <- rep(lims, length = 6)
    gx <- seq(lims[1], lims[2], length = n[1])
    gy <- seq(lims[3], lims[4], length = n[2])
    gz <- seq(lims[5], lims[6], length = n[3])
    mx <- matrix(outer(gx, x, dnorm, h[1]), n[1], nx)
    my <- matrix(outer(gy, y, dnorm, h[2]), n[2], nx)
    mz <- matrix(outer(gz, z, dnorm, h[3]), n[3], nx)
    v <- array(0, n)
    tmy.nx <- t(my) / nx
    for (k in 1:n[3]) {
        tmy.nz.zk <- tmy.nx * mz[k,] # uses recycling to scale the rows
        v[,,k] <- mx %*%  tmy.nz.zk
    }
    return(list(x = gx, y = gy, z = gz, d = v))
}

