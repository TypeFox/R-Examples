
# ------------------------------------------------------------------------------
# R functions to support MKDE in 2D and 3D
# Author: Jeff A. Tracey, PhD: USGS-WERC, San Diego Field Station
# Created: 4 August 2011
# ------------------------------------------------------------------------------

#library(Rcpp)
#library(raster)
#library(sp)

# ---------------------------------------------------------------------
# Functions to initialize data structures
# ---------------------------------------------------------------------

# x, y, and z coordinates are CELL CENTERS
# arrays are filled in order of x, y, z (x=row, y=column, z=level)

# initialize 2D MKDE
initializeMKDE2D <- function(xLL, xCellSize, nX, yLL, yCellSize, nY) {
  xGrid <- xLL + xCellSize*(0:(nX - 1))
  yGrid <- yLL + yCellSize*(0:(nY - 1))
  nx <- length(xGrid)
  ny <- length(yGrid)
  z.min <- array(-Inf, c(nx, ny)) # include all z
  z.max <- array(Inf, c(nx, ny)) # include all z
  d <- array(NA, c(nx, ny, 1))
  res <- list(dimension=2, x=xGrid, y=yGrid, z=NA, z.min=z.min, z.max=z.max, nx=nx, ny=ny, nz=1, d=d)
  return(res)
}

# initialize 3D MKDE
initializeMKDE3D <- function(xLL, xCellSize, nX, yLL, yCellSize, nY, zLL,
                             zCellSize, nZ) {
  xGrid <- xLL + xCellSize*(0:(nX - 1))
  yGrid <- yLL + yCellSize*(0:(nY - 1))
  zGrid <- zLL + zCellSize*(0:(nZ - 1))
  nx <- length(xGrid)
  ny <- length(yGrid)
  nz <- length(zGrid)
  z.min <- array((min(zGrid) - 0.5*zCellSize), c(nx, ny)) # include all z
  z.max <- array((max(zGrid) + 0.5*zCellSize), c(nx, ny)) # include all z
  d <- array(NA, c(nx, ny, nz))
  res <- list(dimension=3, x=xGrid, y=yGrid, z=zGrid, z.min=z.min, z.max=z.max, nx=nx, ny=ny, nz=nz, d=d)
  return(res)
}

# Sets up the list with movement data
initializeMovementData <- function(t.obs, x.obs, y.obs, z.obs=NULL, 
                                   sig2obs=0.0, sig2obs.z=NA, t.max=max(diff(t.obs), na.rm=TRUE)) {
  # CHECK LENGTHS
  if (is.null(z.obs)) {
    dimension=2
  } else {
    dimension=3
  }
  n <- length(t.obs)
  a.obs <- rep(NA, n)
  if (length(sig2obs) == 1) {
    sig2obs.vec <- rep(sig2obs, n)
  } else if (length(sig2obs) == n) {
    sig2obs.vec <- sig2obs
  } else {
    stop("The length of sig2obs is not correct.")
  }
  if (is.na(sig2obs.z)) {
    if (length(sig2obs) == 1) {
      sig2obs.z.vec <- rep(sig2obs, n)
    } else if (length(sig2obs) == n) {
      sig2obs.z.vec <- sig2obs
    } else {
      stop("The length of sig2obs is not correct.")
    }
  } else {
    if (length(sig2obs) == 1) {
      sig2obs.z.vec <- rep(sig2obs.z, n)
    } else if (length(sig2obs) == n) {
      sig2obs.z.vec <- sig2obs.z
    } else {
      stop("The length of sig2obs is not correct.")
    }
  }
  too.much.time <- c((diff(t.obs) > t.max), TRUE)
  move.dat <- list(dimension=dimension, 
                   t.obs=t.obs, 
                   x.obs=x.obs,
                   y.obs=y.obs, 
                   z.obs=z.obs, 
                   a.obs=a.obs,
                   t.max=t.max,
                   sig2xy=rep(NA, n-1), 
                   sig2z=rep(NA, n-1),
                   sig2obs=sig2obs.vec, 
                   sig2obs.z=sig2obs.z.vec, 
                   n.excl.time=too.much.time,       # step-based (pre-computed)
                   n.excl.bound=rep(FALSE, n),      # location-based
                   n.excl.nomove=rep(FALSE, n),     # step-based
                   use.obs=(!too.much.time)         # overall indicator for each step
  )
  return(move.dat)
}

# Set z lower bound from raster
setMinimumZfromRaster <- function(mkde.obj, raster.layer) {
  xy <- expand.grid(x=mkde.obj$x, y=mkde.obj$y)
  z.tmp <- extract(raster.layer, xy)
  # now put in mkde.obj$z.min
  ij <- as.matrix(expand.grid(i=1:mkde.obj$nx, j=1:mkde.obj$ny))
  mkde.obj$z.min[ij] <- as.numeric(z.tmp)
  i <- which(is.na(mkde.obj$z.min))
  zCellSize <- mkde.obj$z[2] - mkde.obj$z[1]
  mkde.obj$z.min[i] <- (min(mkde.obj$z) - zCellSize)
  return(mkde.obj)
}

# Set z upper bound from raster
setMaximumZfromRaster <- function(mkde.obj, raster.layer) {
  xy <- expand.grid(x=mkde.obj$x, y=mkde.obj$y)
  z.tmp <- extract(raster.layer, xy)
  # now put in mkde.obj$z.max
  ij <- as.matrix(expand.grid(i=1:mkde.obj$nx, j=1:mkde.obj$ny))
  mkde.obj$z.max[ij] <- as.numeric(z.tmp)
  i <- which(is.na(mkde.obj$z.max))
  zCellSize <- mkde.obj$z[2] - mkde.obj$z[1]
  mkde.obj$z.max[i] <-(max(mkde.obj$z) + zCellSize)
  return(mkde.obj)
}

# Set z lower bound from constant
setMinimumZfromConstant <- function(mkde.obj, val) {
  mkde.obj$z.min[,] <- val
  return(mkde.obj)
}

# Set z upper bound from constant
setMaximumZfromConstant <- function(mkde.obj, val) {
  mkde.obj$z.max[,] <- val
  return(mkde.obj)
}

# ---------------------------------------------------------------------
# Functions to set up data and estimate variances from data
# ---------------------------------------------------------------------
estVarMKDE <- function(move.dat) {
  # set up data
  t.dat <- move.dat$t.obs[move.dat$use.obs]
  x.dat <- move.dat$x.obs[move.dat$use.obs]
  y.dat <- move.dat$y.obs[move.dat$use.obs]
  z.dat <- move.dat$z.obs[move.dat$use.obs]
  sig2 <- move.dat$sig2obs[move.dat$use.obs]
  sig2z <- move.dat$sig2obs.z[move.dat$use.obs]
  # nested objective functions
  objFunc <- function(x, vfac, vterm, d) {
    s2 <- x*vfac + vterm
    nm <- (1/(2*pi*s2))*exp(-d/(2*s2))
    nll <- -sum(log(nm), na.rm=T) # 
    return(nll)
  }
  
  # initial set up
  n <- length(move.dat$t.obs)
  i <- seq(2, (n-1), by=2) # indices of test points
  i.fm <- i-1
  i.to <- i+1
  dt.i <- t.dat[i.to] - t.dat[i.fm]
  j <- which(dt.i <= 2.0*move.dat$t.max) # x2 because using alt. locs
  i <- i[j]
  i.fm <- i-1
  i.to <- i+1
  dt.i <- t.dat[i.to] - t.dat[i.fm]
  if ((length(na.omit(dt.i)) > 0) & (all(dt.i > 0))) {
    alpha.i <- (t.dat[i] - t.dat[i.fm])/dt.i
    varfact <- dt.i*alpha.i*(1 - alpha.i)
    varterm <- sig2[i.fm]*(1 - alpha.i)^2 + sig2[i.to]*alpha.i^2
    # --- (X, Y) variance
    e.x <- x.dat[i.fm] + alpha.i*(x.dat[i.to] - x.dat[i.fm])
    e.y <- y.dat[i.fm] + alpha.i*(y.dat[i.to] - y.dat[i.fm])
    dsq.xy <- (x.dat[i] - e.x)^2 + (y.dat[i] - e.y)^2
    #
    mn.sig <- 1e-06
    mx.sig <- max((dsq.xy/(qnorm(0.95, 0, 1)^2) - varterm)/varfact, na.rm=TRUE)
    opt.xy <- optimize(objFunc, interval=c(mn.sig, mx.sig), varfact,
                       varterm, dsq.xy)
    move.dat$sig2xy <- rep(opt.xy$minimum, n-1)
    # --- Z variance
    move.dat$sig2z <- rep(NA, n-1)
    if (move.dat$dimension == 3) {
      varfact.z <- dt.i*alpha.i*(1 - alpha.i)
      varterm.z <- sig2z[i.fm]*(1 - alpha.i)^2 + sig2z[i.to]*alpha.i^2
      e.z <- z.dat[i.fm] + alpha.i*(z.dat[i.to] - z.dat[i.fm])
      dsq.z <- (z.dat[i] - e.z)^2
      mx.sig <- max((dsq.z/(qnorm(0.95, 0, 1)^2) - varterm.z)/varfact.z, na.rm=TRUE)
      opt.z <- optimize(objFunc, interval=c(mn.sig, mx.sig), varfact.z,
                        varterm.z, dsq.z)
      move.dat$sig2z <- rep(opt.z$minimum, n-1)
    }
  }
  return(move.dat)
}


# 
deselectLocationsOutsideBounds <- function(move.dat, mkde.obj) {
  # check to make sure move.dat$z.obs is present
  if (all(is.finite(mkde.obj$z.min)) & all(is.finite(mkde.obj$z.max)) & all(is.finite(move.dat$z.obs[move.dat$use.obs]))) {
    n <- length(move.dat$x.obs)
    nx <- length(mkde.obj$x)
    ny <- length(mkde.obj$y)
    xcs <- mkde.obj$x[2] - mkde.obj$x[1]
    ycs <- mkde.obj$y[2] - mkde.obj$y[1]
    move.dat$a.obs <- rep(NA, n)
    for (i in 1:n) {
      # get indices into mkde.obj$h, [ii and ij], make sure they are in bounds
      ii <- round((move.dat$x.obs[i] - mkde.obj$x[1] + 0.5*xcs)/xcs)
      ij <- round((move.dat$y.obs[i] - mkde.obj$y[1] + 0.5*xcs)/xcs)
      if (ii >= 1 & ii <= nx & ij >=1 & ij <= ny) {
        tmp.zmin <- mkde.obj$z.min[ii, ij]
        tmp.zmax <- mkde.obj$z.max[ii, ij]
         move.dat$a.obs[i] <- move.dat$z.obs[i] - tmp.zmin
        if (move.dat$z.obs[i] >= tmp.zmin & move.dat$z.obs[i] <= tmp.zmax) {
          move.dat$n.excl.bound[i] <- FALSE
        } else {
          move.dat$n.excl.bound[i] <- TRUE
        }
      } else { # outside density region
        move.dat$n.excl.bound[i] <- TRUE
      }
    }
    # update move.dat$use.obs
    for (i in 1:(n-1)) {
      if ((!move.dat$n.excl.bound[i]) & (!move.dat$n.excl.bound[i+1]) & 
            (!move.dat$n.excl.nomove[i]) & (!move.dat$n.excl.time[i])) {
        move.dat$use.obs[i] <- TRUE
      } else {
        move.dat$use.obs[i] <- FALSE
      }
    }
    move.dat$use.obs[n] <- FALSE
  } # else do nothing
  return(move.dat)
}

deselectNonMovementSteps <- function(move.dat, p=0.05) {
  # integrate normal distribution for location error in x, y, z
  coOccurProb <- function(x, m0, m1, sd0, sd1) {
    return(sqrt(dnorm(x, m0, sd0)*dnorm(x, m1, sd1)))
  }
  n <- length(move.dat$t.obs) - 1
  pp <- 1 - p/2
  for (i in 1:n) {
    p.x <- pnorm(abs(move.dat$x.obs[i+1]-move.dat$x.obs[i]), 0, sqrt(move.dat$sig2obs[i]))
    p.y <- pnorm(abs(move.dat$y.obs[i+1]-move.dat$y.obs[i]), 0, sqrt(move.dat$sig2obs[i]))
    if (move.dat$dimension == 3) {
      p.z <- pnorm(abs(move.dat$z.obs[i+1]-move.dat$z.obs[i]), 0, sqrt(move.dat$sig2obs.z[i]))
    } else {
      p.z <- 0.5
    }
    if (p.x > pp || p.y > pp || p.z > pp) {
      #
    } else {
      move.dat$n.excl.nomove[i] <- TRUE
      cat(paste("Move step", i, "considered a non-movement\n"))
    }
    move.dat$n.excl.nomove[n+1] <- TRUE # last obs is non-movement
  }
  # update move.dat$use.obs
  for (i in 1:(n-1)) {
    if ((!move.dat$n.excl.bound[i]) & (!move.dat$n.excl.bound[i+1]) & 
          (!move.dat$n.excl.nomove[i]) & (!move.dat$n.excl.time[i])) {
      move.dat$use.obs[i] <- TRUE
    } else {
      move.dat$use.obs[i] <- FALSE
    }
  }
  move.dat$use.obs[n] <- FALSE
  return(move.dat)
}

# ---------------------------------------------------------------------
# Functions to compute density (inc. wrappers for C++ functions)
# ---------------------------------------------------------------------

mkde2Dgrid <- function(mkde.obj, move.dat, t.step, d.thresh) {
  out <- .Call("mkde2dGrid02", 
               move.dat$t.obs, 
               move.dat$x.obs, 
               move.dat$y.obs, 
               as.integer(move.dat$use.obs), 
               mkde.obj$x, mkde.obj$y, 
               move.dat$sig2xy, move.dat$sig2obs, 
               t.step, d.thresh, 
               PACKAGE = "mkde")
  mkde.obj$d <- array(out, c(mkde.obj$nx, mkde.obj$ny, 1))
  return(mkde.obj)
}

mkde2Dinteraction <- function(mkde.obj, move.dat0, move.dat1, t.step, d.thresh) {
  # CHECK move.dat0 AND move.dat1 TO MAKE SURE THE TIMES MATCH
  if (all(is.na(move.dat0$sig2xy))) {
    move.dat0 <- estVarMKDE(move.dat0)
  }
  if (all(is.na(move.dat1$sig2xy))) {
    move.dat1 <- estVarMKDE(move.dat1)
  }
  out <- .Call("mkde2dGridv02interact", 
               move.dat0$t.obs, 
               move.dat0$x.obs, 
               move.dat0$y.obs, 
               move.dat1$x.obs, 
               move.dat1$y.obs,
               as.integer(move.dat0$use.obs & move.dat1$use.obs), 
               mkde.obj$x, mkde.obj$y, move.dat0$sig2xy, move.dat1$sig2xy, 
               move.dat0$sig2obs, move.dat1$sig2obs, 
               t.step, d.thresh, PACKAGE = "mkde")
  mkde.obj$d <- array(out, c(mkde.obj$nx, mkde.obj$ny, 1))
  return(list(mkde.obj=mkde.obj, move.dat0=move.dat0, move.dat1=move.dat1))
}

mkde3Dgrid <- function(mkde.obj, move.dat, t.step, d.thresh) {
  if (all(is.na(move.dat$sig2xy))) {
    move.dat <- estVarMKDE(move.dat)
  }
  out <- .Call("mkde3dGridv02", 
               move.dat$t.obs, 
               move.dat$x.obs, 
               move.dat$y.obs, 
               move.dat$z.obs,
               as.integer(move.dat$use.obs), 
               mkde.obj$x, mkde.obj$y, mkde.obj$z, 
               mkde.obj$z.min,  mkde.obj$z.max, 
               move.dat$sig2xy, move.dat$sig2z, 
               move.dat$sig2obs, move.dat$sig2obs.z, 
               t.step, d.thresh, 
               PACKAGE = "mkde")
  mkde.obj$d <- array(out, c(mkde.obj$nx, mkde.obj$ny, mkde.obj$nz))
  # may have to use aperm()
  return(mkde.obj)
}

mkde3Dinteraction <- function(mkde.obj, move.dat0, move.dat1, t.step, d.thresh) {
  # CHECK move.dat0 AND move.dat1 TO MAKE SURE THE TIMES MATCH
  if (all(is.na(move.dat0$sig2xy))) {
    move.dat0 <- estVarMKDE(move.dat0)
  }
  if (all(is.na(move.dat1$sig2xy))) {
    move.dat1 <- estVarMKDE(move.dat1)
  }
  res <- .Call("mkde3dGridv02interact", 
               move.dat0$t.obs, 
               move.dat0$x.obs, 
               move.dat0$y.obs, 
               move.dat0$z.obs, 
               move.dat1$x.obs, 
               move.dat1$y.obs, 
               move.dat1$z.obs, 
               as.integer(move.dat0$use.obs & move.dat1$use.obs), 
               mkde.obj$x, mkde.obj$y, mkde.obj$z, 
               mkde.obj$z.min,  mkde.obj$z.max, move.dat0$sig2xy, move.dat1$sig2xy, 
               move.dat0$sig2z, move.dat1$sig2z, move.dat0$sig2obs, move.dat1$sig2obs, 
               move.dat0$sig2obs.z, move.dat1$sig2obs.z, t.step, d.thresh, 
               PACKAGE = "mkde")
  mkde.obj$d <- array(res, c(mkde.obj$nx, mkde.obj$ny, mkde.obj$nz))
  return(list(mkde.obj=mkde.obj, move.dat0=move.dat0, move.dat1=move.dat1))
}

computeAreaRaster <- function(RelevMatrix, RcellSize) {
	out <- .Call( "computeCellSurfaceArea",RelevMatrix,
                     RcellSize, PACKAGE = "mkde" )
  out2 <- matrix(out, nrow=nrow(RelevMatrix), 
                 ncol=ncol(RelevMatrix), byrow=TRUE)
  return(out)
}


# initializeDensity(mkde.obj, move.dat)
# this function will check the dimension of mkde.obj and
# the existence/size of move.dat$z and then fill mkde.obj$d
# accordingly using mkde2Dgrid() or mkde3Dgrid()
initializeDensity <- function(mkde.obj, move.dat, integration.step=0.5, d.thresh=1e-25) {
  # calculate variances
  if (all(is.na(move.dat$sig2xy))) {
    move.dat <- estVarMKDE(move.dat)
  }
  if (mkde.obj$dimension == 2) {
    mkde.obj <- mkde2Dgrid(mkde.obj, move.dat, integration.step, d.thresh)
  } else if (mkde.obj$dimension == 2.5) {
    mkde.obj <- mkde2Dgrid(mkde.obj, move.dat, integration.step, d.thresh)
    if (!is.null(mkde.obj$area)) {
      mkde.obj$d[,,1] <- mkde.obj$d[,,1]*mkde.obj$area/sum(mkde.obj$d[,,1]*mkde.obj$area, na.rm=TRUE)
    }
  }else if (mkde.obj$dimension == 3) {
    mkde.obj <- mkde3Dgrid(mkde.obj, move.dat, integration.step, d.thresh)
  } else {
    stop("Dimension of mkde.obj is not 2, 2.5, or 3.")
  }
  return(list(mkde.obj=mkde.obj, move.dat=move.dat))
}

initializeAreaRaster <- function(mkde.obj) {
  if (mkde.obj$dimension == 2 | mkde.obj$dimension == 2.5) {
    cellSzX = mkde.obj$x[2] - mkde.obj$x[1]
    cellSzY = mkde.obj$y[2] - mkde.obj$y[1]
    if (abs(cellSzX - cellSzY) < 1e-8) {
      mkde.obj$area <- computeAreaRaster(mkde.obj$z.min, cellSzX)
      mkde.obj$dimension = 2.5
    } else {
      stop("Cannot create area raster if cell size in x and y dimensions unequal.")
    }
  }
  return(mkde.obj)
}

# ---------------------------------------------------------------------
# Functions to support analysis and visualization
# ---------------------------------------------------------------------

# mkde.obj - MKDE object (list)
# prob     - probabilities associated with quantiles
computeContourValues <- function(mkde.obj, prob) {
  a <- 1 - prob
  d2 <- sort(c(mkde.obj$d))
  d3 <- cumsum(d2)/sum(d2)
  nq <- length(a)
  thresh <- rep(NA, nq)
  for (i in 1:nq) {
    j <- na.omit(which(d3 <= a[i]))
    if (length(j) > 0) {
      thresh[i] <- d2[max(j)]
    }
  }
  out <- data.frame(prob=prob, threshold=thresh)
  return(out)
}

# Area or volume within some contour (incl. product of densities)
#   calculate volume or area of cells
#   compute number of cells >= each contour value
#   multipy each number by the volume/area
#   return the results
computeSizeMKDE <- function(mkde.obj, prob) {
  cntrs <- computeContourValues(mkde.obj, prob)
  nc <- nrow(cntrs)
  out <- rep(NA, nc)
  xSz <- mkde.obj$x[2] - mkde.obj$x[1]
  ySz <- mkde.obj$y[2] - mkde.obj$y[1]
  if (mkde.obj$dimension == 2) {
    av <- xSz*ySz
  } else if (mkde.obj$dimension == 3) {
    zSz <- mkde.obj$z[2] - mkde.obj$z[1]
    av <- xSz*ySz*zSz
  }
  for (i in 1:nc) {
    j <- which(mkde.obj$d >= cntrs$threshold[i])
    if (mkde.obj$dimension == 2.5) {
      out[i] <- sum(mkde.obj$area[j])
    } else {
      out[i] <- av*length(j)
    }
  }
  return(out)
}


# function to visualize densities
# 2D image plot
plotMKDE <- function(mkde.obj, z.index=1, probs=c(0.99, 0.95, 0.90, 0.75, 0.5, 0.0),
                     cmap=rev(rainbow(length(probs)-1)), add=FALSE, ...) {
  if (mkde.obj$dimension == 3) {
    dens.dat <- mkde.obj$d[,,z.index]
  } else {
    dens.dat <- mkde.obj$d[,,1]
  }
  cont.vals <- computeContourValues(mkde.obj, probs)
  image(mkde.obj$x, mkde.obj$y, dens.dat, breaks=cont.vals$threshold, col=cmap, add=add, ...)
}

# ---------------------------------------------------------------------
# Functions to write output in other formats
# ---------------------------------------------------------------------

# convert an mkde object to a RasterLayer or RasterStack (if 3D)
mkdeToRaster <- function(mkde.obj) {
  sx <- mkde.obj$x[2] - mkde.obj$x[1]
  sy <- mkde.obj$y[2] - mkde.obj$y[1]
  rst.res <- NA
  if (mkde.obj$dimension == 2 | mkde.obj$dimension == 2.5) {
    # make a raster layer
    rst.res <- raster(t(matrix(mkde.obj$d[,,1], nrow=mkde.obj$nx, ncol=mkde.obj$ny))[mkde.obj$ny:1,], 
                      xmn=(mkde.obj$x[1] - 0.5*sx), xmx=(mkde.obj$x[mkde.obj$nx] + 0.5*sx), 
                      ymn=(mkde.obj$y[1] - 0.5*sy), ymx=(mkde.obj$y[mkde.obj$ny] + 0.5*sy))
  } else if (mkde.obj$dimension == 3) {
    # make a raster stack
    rst.res <- stack()
    for (k in 1:mkde.obj$nz) {
      rst.tmp <- raster::raster(t(matrix(mkde.obj$d[,,k], nrow=mkde.obj$nx, ncol=mkde.obj$ny))[mkde.obj$ny:1,], 
                        xmn=(mkde.obj$x[1] - 0.5*sx), xmx=(mkde.obj$x[mkde.obj$nx] + 0.5*sx), 
                        ymn=(mkde.obj$y[1] - 0.5*sy), ymx=(mkde.obj$y[mkde.obj$ny] + 0.5*sy))
      rst.res <- addLayer(rst.res, rst.tmp)
    }
    names(rst.res) <- paste("z=", mkde.obj$z, sep="")
  }
  return(rst.res)
}

# Write 3D MKDE in VTK format
writeToVTK <- function(mkde.obj, fname, description="3D MKDE", cumprob=FALSE) {
  if (cumprob) {
    dtmp <- sort(c(mkde.obj$d), index.return=TRUE)
    dens <- rep(NA, length(dtmp))
    dens[dtmp$ix] <- cumsum(dtmp$x)/sum(dtmp$x, na.rm=TRUE)
  } else {
    dens <- c(mkde.obj$d)
  }
  if (mkde.obj$dimension == 3) {
    .Call( "writeMKDE3DtoVTK", mkde.obj$x, mkde.obj$y, mkde.obj$z,
          dens, fname, description, PACKAGE = "mkde" )
  } # don't do anything for other dimensions...YET
}

# write to GRASS GIS 3D ASCII raster
writeToGRASS <- function(mkde.obj, fname, nodat="NA", cumprob=FALSE) {
  if (mkde.obj$dimension == 3) {
    if (cumprob) {
      dtmp <- sort(c(mkde.obj$d), index.return=TRUE)
      dens <- rep(NA, length(dtmp))
      dens[dtmp$ix] <- cumsum(dtmp$x)/sum(dtmp$x, na.rm=TRUE)
    } else {
      dens <- c(mkde.obj$d)
    }
    .Call( "writeMKDE3DtoGRASS", mkde.obj$x, mkde.obj$y, mkde.obj$z,
          dens, fname, nodat, PACKAGE = "mkde" )
  } # don't do anything for other dimensions...YET
}

# write to XDMF files
writeToXDMF <- function(mkde.obj, fname, nodat="NA", cumprob=FALSE) {
  if (mkde.obj$dimension == 3) {
    fnXDMF <- paste(fname, ".xdmf", sep="")
    fnDAT <- paste(fname, ".dat", sep="")
    if (cumprob) {
      dtmp <- sort(c(mkde.obj$d), index.return=TRUE)
      dens <- rep(NA, length(dtmp))
      dens[dtmp$ix] <- cumsum(dtmp$x)/sum(dtmp$x, na.rm=TRUE)
    } else {
      dens <- c(mkde.obj$d)
    }
    .Call( "writeMKDE3DtoXDMF", mkde.obj$x, mkde.obj$y, mkde.obj$z,
           dens, fnXDMF, fnDAT, PACKAGE = "mkde" )
  } # don't do anything for other dimensions...YET
}

writeRasterToXDMF <- function(rast, fname, nodat="NA") {
  r.xy <- coordinates(rast)
  ext <- extent(rast)
  n.xy <- dim(rast)
  cell.sz <- res(rast)
  r.x <- ext@xmin + 0.5*cell.sz[1] + (0:(n.xy[2] - 1))*cell.sz[1]
  r.y <- ext@ymin + 0.5*cell.sz[2] + (0:(n.xy[1] - 1))*cell.sz[2]
  r.v <- values(rast, format="matrix")
  nrw <- nrow(r.v)
  ncl <- ncol(r.v)
  r.v <- t(r.v[nrw:1,]) # flip
  r.v <- as.vector(r.v, mode="numeric")
  fnXDMF <- paste(fname, ".xdmf", sep="")
  fnDAT <- paste(fname, ".dat", sep="")
  .Call("writeRasterToXDMF", r.x, r.y, r.v, fnXDMF, fnDAT, PACKAGE = "mkde")
}

writeRasterToVTK <- function(elev, r.rst, g.rst, b.rst, descr, fname) {
  # make sure all rasters are same dims, etc.
  r.xy <- coordinates(elev)
  ext <- extent(elev)
  n.xy <- dim(elev)
  cell.sz <- res(elev)
  r.x <- ext@xmin + 0.5*cell.sz[1] + (0:(n.xy[2] - 1))*cell.sz[1]
  r.y <- ext@ymin + 0.5*cell.sz[2] + (0:(n.xy[1] - 1))*cell.sz[2]
  r.v <- raster::values(elev, format="matrix")
  nrw <- nrow(r.v)
  ncl <- ncol(r.v)
  r.v <- t(r.v[nrw:1,]) # flip
  r.v <- as.vector(r.v, mode="numeric")
  #
  r.r <- values(r.rst, format="matrix")
  r.r <- t(r.r[nrw:1,]) # flip
  r.r <- as.vector(r.r, mode="numeric")
  r.g <- values(g.rst, format="matrix")
  r.g <- t(r.g[nrw:1,]) # flip
  r.g <- as.vector(r.g, mode="numeric")
  r.b <- raster::values(b.rst, format="matrix")
  r.b <- t(r.b[nrw:1,]) # flip
  r.b <- as.vector(r.b, mode="numeric")
  # (SEXP xgrid, SEXP ygrid, SEXP elev, SEXP rd, SEXP gr, SEXP bl, SEXP filenameVTK)
  .Call("writeRasterToVTK", r.x, r.y, r.v, r.r, r.g, r.b, descr, fname, PACKAGE="mkde")
}

writeObservedLocationVTK <- function(move.dat, mkde.obj, 
                                     description="Observed Locations", 
                                     filename="move_locations.vtk", 
                                     control=list(z.scale=1)) {
  # set up move data
  # ---------------------------------------------------------------------------
  n.loc <- length(move.dat$t.obs)
  n.use <- sum(move.dat$use.obs)
  
  # for future use
  if (mkde.obj$dimension == 3) {
    #
  } else { # 2 or 2.5; base height on DEM if not given
    #
  }
  
  # WRITE VTK FILE
  # ---------------------------------------------------------------------------
  # write header
  cat('# vtk DataFile Version 3.0\n', file=filename)
  cat(paste(description, '\n', sep=''), file=filename, append=TRUE)
  cat('ASCII\n', file=filename, append=TRUE)
  cat('DATASET POLYDATA\n',file=filename, append=TRUE)
  # Write points
  cat(paste('POINTS ', n.use, ' double\n', sep=''), file=filename, append=TRUE)
  for (i in 1:n.loc) {
    if(move.dat$use.obs[i]) {
      cat(paste(move.dat$x.obs[i], ' ', move.dat$y.obs[i], ' ', 
                control$z.scale*move.dat$z.obs[i], '\n', sep=''), file=filename, append=TRUE)
    }
  }
  # Write vertexes
  cat(paste('\nVERTICES ', n.use, ' ', 2*n.use, '\n', sep=''), file=filename, append=TRUE)
  j <- 0
  for (i in 1:n.loc) {
    if(move.dat$use.obs[i]) {
      cat(paste("1 ", j, '\n', sep=''), file=filename, append=TRUE)
      j <- j + 1
    }
  }
}

writeInterpolatedPathVTK <- function(move.dat, mkde.obj, 
                        description="Interpolated Movement Path", 
                        filename="move_path.vtk", 
                        control=list(method="linear", method.par=list(n=1), z.scale=1)) {
  # set up move data
  # ---------------------------------------------------------------------------
  n.mv <- length(move.dat$t.obs) - 1
  
  # anything to do?
  if (mkde.obj$dimension == 3) {
    #
  } else { # 2 or 2.5; interpolate along DEM
    #
  }
  
  # make line data
  # ---------------------------------------------------------------------------
  line.data <- list() # things to write are stored in a list of lists
  ii <- 1
  offset.current <- 0
  num.ind.total <- 0
  for (i in 1:n.mv) {
    if (move.dat$use.obs[i]) {
      x0 <- c(move.dat$x.obs[i], move.dat$x.obs[i + 1])
      y0 <- c(move.dat$y.obs[i], move.dat$y.obs[i + 1])
      z0 <- control$z.scale*c(move.dat$z.obs[i], move.dat$z.obs[i + 1])
      line.data[[ii]] <- list(ind.offset=offset.current, 
                              num.points=2, 
                              line.points=data.frame(x=x0, y=y0, z=z0), 
                              ind2points=c(1:2) + offset.current)
      offset.current <- line.data[[ii]]$ind.offset + 2
      ii <- ii + 1
      num.ind.total <- num.ind.total + 2
    }
  }
  
  # WRITE VTK FILE
  # ---------------------------------------------------------------------------
  n.lines <- length(line.data)
  # write header
  cat('# vtk DataFile Version 3.0\n', file=filename)
  cat(paste(description, '\n', sep=''), file=filename, append=TRUE)
  cat('ASCII\n', file=filename, append=TRUE)
  cat('DATASET POLYDATA\n',file=filename, append=TRUE)
  # Write points
  cat(paste('POINTS ', num.ind.total, ' double\n', sep=''), file=filename, append=TRUE)
  for (i in 1:n.lines) {
    n.pnts <- line.data[[i]]$num.points
    for (j in 1:n.pnts) {
      cat(paste(line.data[[i]]$line.points$x[j], ' ', line.data[[i]]$line.points$y[j], ' ', 
                line.data[[i]]$line.points$z[j], '\n', sep=''), file=filename, append=TRUE)
    }
  }
  # Write lines
  cat(paste('\nLINES ', n.lines, ' ', 3*n.lines, '\n', sep=''), file=filename, append=TRUE)
  for (i in 1:n.lines) {
    n.pnts <- line.data[[i]]$num.points
    cat(paste(n.pnts, ' ', sep=''), file=filename, append=TRUE)
    for (j in 1:n.pnts) {
      if (j < n.pnts) {
        cat(paste(line.data[[i]]$ind2points[j] - 1, ' ', sep=''), file=filename, append=TRUE)
      } else {
        cat(paste(line.data[[i]]$ind2points[j] - 1, '\n', sep=''), file=filename, append=TRUE)
      }
    }
  }
}







