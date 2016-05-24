## segments.circle.R

segments.circle <- function(rho, shift, circle=360,
                            ncp=1000, gp=gpar()) {
    ## Creates a grob for segments of a circle

    ## Author: Rene Locher
    ## Version: 2009-03-16

    Dphi <- circle/nrow(rho)   ## size of segment
    dphi0 <- circle/ncp        ## delta wanted
    n <- floor(Dphi/dphi0)     ## max. number of pieces of lines per segment
    dphi <- circle/n/nrow(rho) ## delta used

    x <- y <- matrix(rep(NA,n*prod(dim(rho))),ncol=ncol(rho))

    phi <- seq(-Dphi/2+shift,length.out=n*nrow(rho),by=dphi)
    for (ii in 0:(nrow(rho)-1)) {
        x[(ii*n+1):((ii+1)*n),] <-
            outer(sin(2*pi*phi[(ii*n+1):((ii+1)*n)]/circle),rho[ii+1,])
        y[(ii*n+1):((ii+1)*n),] <-
            outer(cos(2*pi*phi[(ii*n+1):((ii+1)*n)]/circle),rho[ii+1,])
    }
    x <- as.vector(x)
    y <- as.vector(y)
    id <- rep(1:ncol(rho), rep(nrow(rho)*n,ncol(rho)))

    polygonGrob(name = "data",
                 x = x,
                 y = y,
                 id = id,
                 default.units = "native",
                 gp = gp)
} ## segments.circle
