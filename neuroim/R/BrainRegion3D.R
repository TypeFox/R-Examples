#' @import iterators

#' @include AllClass.R
{}
#' @include AllGeneric.R
{}


#' Create an instance of class ROIVolume
#' @param vspace the volume \code{BrainSpace}
#' @param coords matrix of voxel coordinates
#' @param data the data values
#' @return an instance of class \code{ROIVolume}
#' @rdname ROIVolume
#' @name ROIVolume
#' @export
ROIVolume <- function(vspace, coords, data=rep(length(indices),1)) {
  new("ROIVolume", space=vspace, coords=coords, data=data)
}
  

.makeSquareGrid <- function(bvol, centroid, surround, fixdim=3) {
  vspacing <- spacing(bvol)
  vdim <- dim(bvol)
  centroid <- as.integer(centroid)
  
  dimnums <- seq(1,3)[-fixdim]
  
  coords <- lapply(centroid, function(x) { round(seq(x-surround, x+surround)) })
  coords <- lapply(dimnums, function(i) {
    x <- coords[[i]]
    x[x > 0 & x <= vdim[i]]
  })
  
  if (all(sapply(coords, length) == 0)) {
    stop(paste("invalid cube for centroid", centroid, " with surround", surround, ": volume is zero"))
  }
  
  if (fixdim == 3) {
    grid <- as.matrix(expand.grid(x=coords[[1]],y=coords[[2]],z=centroid[3]))
  } else if (fixdim == 2) {
    grid <- as.matrix(expand.grid(x=coords[[1]],y=centroid[2],z=coords[[2]]))
  } else if (fixdim == 1) {
    grid <- as.matrix(expand.grid(x=centroid[1],y=coords[[1]],z=coords[[2]]))
  }
  
  grid
  
}
.makeCubicGrid <- function(bvol, centroid, surround) {
  vspacing <- spacing(bvol)
  vdim <- dim(bvol)
  centroid <- as.integer(centroid)
  
  coords <- lapply(centroid, function(x) { round(seq(x-surround, x+surround)) })
  coords <- lapply(1:3, function(i) {
    x <- coords[[i]]
    x[x > 0 & x <= vdim[i]]
  })
 
  if (all(sapply(coords, length) == 0)) {
    stop(paste("invalid cube for centroid", centroid, " with surround", surround, ": volume is zero"))
  }
  
  grid <- as.matrix(expand.grid(x=coords[[1]],y=coords[[2]],z=coords[[3]]))
}



#' Create a square region of interest where the z-dimension is fixed at one voxel coordinate.
#' @param bvol an \code{BrainVolume} or \code{BrainSpace} instance.
#' @param centroid the center of the cube in \emph{voxel} coordinates.
#' @param surround the number of voxels on either side of the central voxel.
#' @param fill optional value(s) to assign to data slot.
#' @param nonzero keep only nonzero elements from \code{bvol}. If \code{bvol} is A \code{BrainSpace} then this argument is ignored.
#' @param fixdim the fixed dimension is the third, or z, dimension.
#' @return an instance of class \code{ROIVolume}.
#' @examples
#'  sp1 <- BrainSpace(c(10,10,10), c(1,1,1))
#'  square <- RegionSquare(sp1, c(5,5,5), 1)
#'  vox <- coords(square)
#'  ## a 3 X 3 X 1 grid
#'  nrow(vox) == 9
#' @export
RegionSquare <- function(bvol, centroid, surround, fill=NULL, nonzero=FALSE, fixdim=3) {
  if (is.matrix(centroid)) {
    centroid <- drop(centroid)
  }
  
  if (length(centroid) != 3) {
    stop("RegionSquare: centroid must have length of 3 (x,y,z coordinates)")
  }
  
  if (surround < 0) {
    stop("'surround' argument cannot be negative")
  }
  
  if (is(bvol, "BrainSpace") && is.null(fill)) {
    fill = 1
  }
  
  grid <- .makeSquareGrid(bvol,centroid,surround,fixdim=fixdim)
  
  vals <- if (!is.null(fill)) {
    rep(fill, nrow(grid))
  } else {
    ## coercion to numeric shouldn't be necessary here.
    as.numeric(bvol[grid])
  }   
  
  keep <- if (nonzero) {
    which(vals != 0)    
  } else {
    seq_along(vals)
  }
  
  ### add central voxel
  new("ROIVolume", space = space(bvol), data = vals[keep], coords = grid[keep, ])
  
}

  
#' Create A Cuboid Region of Interest
#' @param bvol an \code{BrainVolume} or \code{BrainSpace} instance
#' @param centroid the center of the cube in \emph{voxel} coordinates
#' @param surround the number of voxels on either side of the central voxel. A \code{vector} of length 3.
#' @param fill optional value(s) to assign to data slot. 
#' @param nonzero keep only nonzero elements from \code{bvol}. If \code{bvol} is A \code{BrainSpace} then this argument is ignored.
#' @return an instance of class \code{ROIVolume}
#' @examples
#'  sp1 <- BrainSpace(c(10,10,10), c(1,1,1))
#'  cube <- RegionCube(sp1, c(5,5,5), 3)
#'  vox <- coords(cube)
#'  cube2 <- RegionCube(sp1, c(5,5,5), 3, fill=5)
#'  
#'  
#' @export
RegionCube <- function(bvol, centroid, surround, fill=NULL, nonzero=FALSE) {
  if (is.matrix(centroid)) {
    centroid <- drop(centroid)
  }
  
  if (length(centroid) != 3) {
    stop("RegionCube: centroid must have length of 3 (x,y,z coordinates)")
  }
  
  if (surround < 0) {
    stop("'surround' argument cannot be negative")
  }
  
  if (is(bvol, "BrainSpace") && is.null(fill)) {
    fill = 1
  }
  
  grid <- .makeCubicGrid(bvol,centroid,surround)
  
  vals <- if (!is.null(fill)) {
    rep(fill, nrow(grid))
  } else {
    ## coercion to numeric shouldn't be necessary here.
    as.numeric(bvol[grid])
  }   
  
  keep <- if (nonzero) {
    which(vals != 0)    
  } else {
    seq_along(vals)
  }
  
  ### add central voxel
  new("ROIVolume", space = space(bvol), data = vals[keep], coords = grid[keep, ])
  
}


.makeSphericalGrid <- function(bvol, centroid, radius) {
  vspacing <- spacing(bvol)
  vdim <- dim(bvol)
  centroid <- as.integer(centroid)
  mcentroid <- ((centroid-1) * vspacing + vspacing/2)
  cubedim <- ceiling(radius/vspacing)
  
  nsamples <- max(cubedim) * 2 + 1
  vmat <- apply(cbind(cubedim, centroid), 1, function(cdim) {
    round(seq(cdim[2] - cdim[1], cdim[2] + cdim[1], length.out=nsamples))
  })
  
  vlist <- lapply(1:NCOL(vmat), function(i) {
    v <- vmat[,i]
    unique(v[v >= 1 & v <= vdim[i]])
  })
  
  
  if (all(sapply(vlist, length) == 0)) {
    stop(paste("invalid sphere for centroid", paste(centroid, collapse=" "), " with radius",
               radius))
  }

  ## should be: need to add dvals
  ## as.matrix(expand.grid(x = vlist[[1]], y = vlist[[2]], z = vlist[[3]]))[which(dvals <= radius),]
  grid <- as.matrix(expand.grid(x = vlist[[1]], y = vlist[[2]], z = vlist[[3]]))
  
  dvals <- apply(grid, 1, function(gvals) {
    coord <- (gvals-1) * vspacing + vspacing/2
    sqrt(sum((coord - mcentroid)^2))
  })
  
  grid[which(dvals <= radius),]
  
}


#' Create A Spherical Region of Interest
#' @param bvol an \code{BrainVolume} or \code{BrainSpace} instance
#' @param centroid the center of the sphere in voxel space
#' @param radius the radius in real units (e.g. millimeters) of the spherical ROI
#' @param fill optional value(s) to assign to data slot
#' @param nonzero keep only nonzero elements from \code{bvol}
#' @return an instance of class \code{ROIVolume}
#' @examples
#'  sp1 <- BrainSpace(c(10,10,10), c(1,1,1))
#'  cube <- RegionSphere(sp1, c(5,5,5), 3.5)
#'  vox = coords(cube)
#' @export
RegionSphere <- function (bvol, centroid, radius, fill=NULL, nonzero=FALSE) {
  if (is.matrix(centroid)) {
    centroid <- drop(centroid)
  }
  if (length(centroid) != 3) {
    stop("RegionSphere: centroid must have length of 3 (x,y,z coordinates)")
  }
  
  if (is.null(fill) && is(bvol, "BrainSpace")) {
    fill = 1
  }
  
  bspace <- space(bvol)
  vspacing <- spacing(bvol)
  vdim <- dim(bvol)
  centroid <- as.integer(centroid)
  mcentroid <- ((centroid-1) * vspacing + vspacing/2)
 
  grid <- .makeSphericalGrid(bvol, centroid, radius)
   
  #dvals <- apply(grid, 1, function(gvals) {
  #  coord <- (gvals-1) * vspacing + vspacing/2
  #  sqrt(sum((coord - mcentroid)^2))
  #})
  
  #idx <- which(dvals <= radius)
   
  vals <- if (!is.null(fill)) {
    rep(fill, nrow(grid))
  } else {    
    ## coercion to numeric shouldn't be necessary here.
    as.numeric(bvol[grid])
  }   
  
  keep <- if (nonzero) {
    which(vals != 0)    
  } else {
    seq_along(vals)
  }
  
  new("ROIVolume", space = bspace, data = vals[keep], coords = grid[keep, ])
}

.resample <- function(x, ...) x[sample.int(length(x), ...)]

#' Create an Random Searchlight iterator
#' @param mask an image volume containing valid central voxels for roving searchlight
#' @param radius in mm of spherical searchlight
#' @export
RandomSearchlight <- function(mask, radius) {
  done <- array(FALSE, dim(mask))
  mask.idx <- which(mask != 0)
  grid <- indexToGrid(mask, mask.idx)
  


  prog <- function() { sum(done)/length(mask.idx) }
  
  nextEl <- function() {
    if (!all(done[mask.idx])) {
      center <- .resample(which(!done[mask.idx]), 1)
      done[center] <<- TRUE
      search <- RegionSphere(mask, grid[center,], radius, nonzero=TRUE) 
      vox <- coords(search)
      vox <- vox[!done[vox],,drop=FALSE]
      done[vox] <<- TRUE
      attr(vox, "center") <- grid[center,]
      vox
      
    } else {
      stop('StopIteration')
    }
  }
  obj <- list(nextElem=nextEl, progress=prog)
  class(obj) <- c("RandomSearchlight", 'abstractiter', 'iter')
  obj
}

#' Create a searchlight iterator that samples regions from within a mask.
#' Searchlight centers are sampled *without* replacement, but the same voxel can belong to multiple searchlight samples.
#' It is in the latter sense that this is a bootstrap resampling scheme.
#' @param mask an image volume containing valid central voxels for roving searchlight
#' @param radius in mm of spherical searchlight
#' @param iter the total number of searchlights to sample (default is 100)
#' @export
BootstrapSearchlight <- function(mask, radius, iter=100) {
  mask.idx <- which(mask != 0)
  grid <- indexToGrid(mask, mask.idx)
  index <- 0
  
  sample.idx <- sample(1:nrow(grid))
  
  prog <- function() { index/length(mask.idx) }
  
  nextEl <- function() {
    if (index <= iter & length(sample.idx) > 0) { 
      index <<- index + 1
      
      cenidx <- sample.idx[1]
      sample.idx <<- sample.idx[-1]
      
      search <- RegionSphere(mask, grid[cenidx,], radius, nonzero=TRUE) 
      vox <- search@coords
      attr(vox, "center") <- grid[index,]
      vox
    } else {
      stop('StopIteration')
    }
  }
  
  obj <- list(nextElem=nextEl, progress=prog)
  class(obj) <- c("BootstrapSearchlight", 'abstractiter', 'iter')
  obj
  
}


#' Create an exhaustive searchlight iterator
#' @param mask an image volume containing valid central voxels for roving searchlight
#' @param radius in mm of spherical searchlight
#' @export
Searchlight <- function(mask, radius) {
  mask.idx <- which(mask != 0)
	grid <- indexToGrid(mask, mask.idx)
	index <- 0
  
	prog <- function() { index/length(mask.idx) }
  
	nextEl <- function() {
		if (index < nrow(grid)) { 
			 index <<- index + 1
		 	 search <- RegionSphere(mask, grid[index,], radius, nonzero=TRUE) 
       vox <- search@coords
			 attr(vox, "center") <- grid[index,]
       vox
		} else {
			stop('StopIteration')
		}
	}
	
	obj <- list(nextElem=nextEl, progress=prog)
  class(obj) <- c("Searchlight", 'abstractiter', 'iter')
	obj
			
}


#' @name as
#' @rdname as-methods
setAs(from="ROIVolume", to="DenseBrainVolume", function(from) {
  dat <- array(0, dim(from@space))
  dat[from@coords] <- from@data
  ovol <- DenseBrainVolume(dat, from@space, from@source)
})


#' @rdname values-methods
#' @export 
setMethod("values", signature(x="ROIVolume"),
          function(x, ...) {
             x@data
          })



#' @rdname indices-methods
#' @export 
setMethod("indices", signature(x="ROIVolume"),
          function(x) {
			  gridToIndex(x@space, x@coords)
		  })
            


#' @export
#' @rdname coords-methods
setMethod(f="coords", signature=signature(x="ROIVolume"),
          function(x) {
            x@coords
          })


#' @export 
#' @rdname length-methods
#' @param x the object to get \code{length}
setMethod(f="length", signature=signature(x="ROIVolume"),
          function(x) {
            length(x@data)
		})

#' extract data from \code{ROIVolume}
#' @export
#' @param x the object
#' @param i first index
#' @param j second index
#' @param drop drop dimension
setMethod("[", signature=signature(x = "ROIVolume", i = "numeric", j = "missing", drop = "ANY"),
          function (x, i, j, drop) {
            x@data[i]
          })

#' show an \code{ROIVolime} 
#' @param object the object
#' @export
setMethod("show", signature=signature(object = "ROIVolume"),
		  function (object) {
			  cat("\n\n\tROIVolume", "\n")
			  cat("\tsize: ", length(object), "\n")
			  cat("\tparent dim:", dim(object), "\n")
			  cat("\tvoxel center of mass: ", colMeans(coords(object)), "\n")
		  })
  
      
.distance <- function(p1, p2) {
  diffs = (p1 - p2)
  sqrt(sum(diffs*diffs))
}


#' Create a Kernel object
#' @param kerndim the dimensions in voxels of the kernel
#' @param vdim the dimensions of the voxels in real units
#' @param FUN the kernel function taking as its first argument representing the distance from the center of the kernel
#' @param ... additional parameters to the kernel FUN
#' @importFrom stats dnorm
#' @export
Kernel <- function(kerndim, vdim, FUN=dnorm, ...) {
  if (length(kerndim) < 2) {
    stop("kernel dim length must be greater than 1")
  }
  
  #kern <- array(0, kerndim)
  
  ## the half-width for each dimensions
  hwidth <- sapply(kerndim, function(d) ceiling(d/2 -1))
  
  ## note, if a kernel dim is even, this will force it to be odd numbered
  grid.vec <- lapply(hwidth, function(sv) seq(-sv, sv))

  # compute relative voxel locations (i.e. centered at 0,0,0)
  voxel.ind <- as.matrix(do.call("expand.grid", grid.vec))
  
  # fractional voxel locations so that the location of a voxel coordinate is centered within the voxel
  cvoxel.ind <- t(apply(voxel.ind, 1, function(vals) sign(vals)* ifelse(vals == 0, 0, abs(vals)-.5)))
  
  ## the coordinates ofthe voxels (i.e. after multiplying by pixel dims)
  coords <- t(apply(cvoxel.ind, 1, function(v) (v * vdim)))
  
  ## distance of coordinate from kernel center
  coord.dist <- apply(coords, 1, .distance, c(0,0,0))
  
  wts <- FUN(coord.dist, ...)
  wts <- wts/sum(wts)

  
  kern.weights <- wts
  
  new("Kernel", width=kerndim, weights=kern.weights, voxels=voxel.ind, coords=coords)

}



#' @param centerVoxel the absolute location of the center of the voxel, default is (0,0,0)
#' @rdname voxels-methods
#' @export
setMethod(f="voxels", signature=signature(x="Kernel"),
          function(x, centerVoxel=NULL) {
            if (is.null(centerVoxel)) {
              x@voxels
            } else {
              sweep(x@voxels, 2, centerVoxel, "+")
            }
          })


  

  
