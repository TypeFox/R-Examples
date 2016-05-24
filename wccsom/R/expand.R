## create a new grid, same topology, twice as large
## interpolate for new units

expand <- function(somnet, plotit=FALSE)
{
  if (plotit) cat("Creating new net...")
  oldxdim <- somnet$grid$xdim
  xd <- oldxdim*2
  oldydim <- somnet$grid$ydim
  yd <- oldydim*2
  topo <- somnet$grid$topo
  gr <- somgrid(xd, yd, topo)
  nhbrdist <- unit.distances(gr, somnet$toroidal)

  if (plotit) cat("redistributing code vectors...")
  ncodes <- gr$xdim * gr$ydim
  noldcodes <- oldxdim * oldydim
  colnrs <- (1:noldcodes) %% oldxdim
  colnrs[colnrs == 0] <- oldxdim ## start from 1
  rownrs <- floor(((1:noldcodes)-1) / oldxdim) ## start from 0
  newindices <- rownrs*2*xd + colnrs * 2 - 1
  acors <- rep(-1, ncodes)
  codes <- matrix(0, ncodes, ncol(somnet$codes))
  codes[newindices,] <- somnet$codes
  acors[newindices] <- somnet$acors
  mycolors <- rep(0, ncodes)
  mycolors[newindices] <- 1

  if (plotit) cat("interpolating unit codes...")
  nclose <- ifelse (topo == "hexagonal", 6, 4)
  for (i in 1:ncodes) {
    if (acors[i] < 0) {
      closeones <- order(nhbrdist[i, newindices])[1:nclose]
      closedists <- nhbrdist[i, newindices[closeones]]
      normfactor <- sum(1/closedists)
      codes[i,] <- colSums(codes[newindices[closeones],] * closedists) /
          normfactor
      ## codes[i,] <- colSums(sweep(codes[newindices[closeones],],
      ##                            2,
      ##                            closedists,
      ##                            FUN="*")) / normfactor
    }
  }
  
  if (plotit) cat("calculating autocovariances...")
  trwdth <- somnet$trwdth
  wghts <- 1 - (0:trwdth)/trwdth
  acors <- wacmat(codes, trwdth=trwdth, wghts=wghts)

  if (plotit) cat("ready.\n")

  somnet$codes <- codes
  somnet$acors <- acors
  somnet$changes <- NULL
  somnet$nhbrdist <- nhbrdist
  somnet$grid <- gr
  
  if (plotit)
    plot(somnet, type="property", property=mycolors)
  
  somnet
}


### Lower the dimensionality of the X matrix by a certain factor.
### Typically used to speed up the training of a network

bucket <- function(x, factor)
{
  if (is.vector(x)) {
    nx <- length(x)
    nout <- ceiling(nx/factor)
    
    nzeros <- factor * nout - nx
    xout <- c(x, rep(0, nzeros))
    xout <- matrix(xout, nrow=factor)
    
    return(colMeans(xout))
  }

  if (is.matrix(x)) {
    nx <- ncol(x)
    nout <- ceiling(nx/factor)
    
    nzeros <- factor * nout - nx
    zeromat <- matrix(0, nrow(x), nzeros)
    
    xout <- array(c(x, zeromat), dim=c(nrow(x), factor, nout))
    return(apply(xout, c(1,3), mean))
  }

  stop("x should be a vector or a matrix")
}


### Reverse procedure as bucket: now, linear interpolation is used.

debucket <- function(x, nout)
{
  t(apply(x, 1, function(y) approx(1:ncol(x), y, n=nout)$y))
}
