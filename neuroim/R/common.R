#' @import abind
NULL


#' matrixToVolumeList
#' converts a matrix to a list of BrainVolumes with values filled at grid coordinates determined by the \code{vox} argument.
#' @param voxmat an N by 3 matrix of voxel coordinates
#' @param mat an N by M matrix of values where M is the number of volumes to create (e.g. one volume per column in \code{mat})
#' @param mask a reference volume defining the geometry of the output volumes. This can either be of type \code{BrainSpace} or \code{BrainVolume}
#' @param default the value that will be used for voxels not contained within voxmat (defualt is \code{NA})
#' @return a \code{list} of \code{BrainVolume} instances, one for each column of \code{mat}
#' @export  
matrixToVolumeList <- function(voxmat, mat, mask, default=NA) {
  if (nrow(voxmat) != nrow(mat)) {
    stop("mismatching dimensions: nrow(voxmat) must equal nrow(mat)")
  }
  lapply(1:ncol(mat), function(i) {
    vol <- array(default, dim(mask))   
    vol[voxmat] <- mat[,i]
    
    if (is(mask, "BrainSpace")) {
      BrainVolume(vol, mask)
    } else {
      BrainVolume(vol, space(mask))
    }
  })
}  

matrixToVols <- matrixToVolumeList

#' @export
#' @rdname splitReduce-methods
setMethod(f="splitReduce", signature=signature(x = "matrix", fac="integer", FUN="function"),
          def=function(x, fac, FUN) {
            callGeneric(x,as.factor(fac), FUN)
          })

#' @export
#' @rdname splitReduce-methods
setMethod(f="splitReduce", signature=signature(x = "matrix", fac="integer", FUN="missing"),
          def=function(x, fac) {
            callGeneric(x,as.factor(fac))
          })

#' @export
#' @rdname splitReduce-methods
setMethod(f="splitReduce", signature=signature(x = "matrix", fac="factor", FUN="missing"),
          def=function(x, fac) {
            if (length(fac) != nrow(x)) {
              stop(paste("x must be same length as factor used for splitting rows"))
            }
            
            ind <- split(seq_along(fac), fac)
            out <- do.call(rbind, lapply(levels(fac), function(lev) {
              colMeans(x[ind[[lev]],])
            }))
            
            row.names(out) <- levels(fac)           
            out
          })

#' @export
#' @rdname splitReduce-methods
setMethod(f="splitReduce", signature=signature(x = "matrix", fac="factor", FUN="function"),
          def=function(x, fac, FUN) {
            if (length(fac) != nrow(x)) {
              stop(paste("x must be same length as split variable"))
            }        
            
            ind <- split(seq_along(fac), fac)
            out <- do.call(rbind, lapply(names(ind), function(lev) {
              apply(x[ind[[lev]],], 2, FUN)
            }))
            
            row.names(out) <- levels(fac)           
            out
          })

#' @export
#' @rdname splitReduce-methods
setMethod(f="splitReduce", signature=signature(x = "BrainVector", fac="factor", FUN="function"),
          def=function(x, fac, FUN) {
            if (length(fac) != prod(dim(x)[1:3])) {
              stop(paste("fac must have as many elements as the number of voxels"))
            }
            
            ind <- split(seq_along(fac), fac)
            out <- do.call(rbind, lapply(names(ind), function(lev) {
              #idx <- which(fac == lev)
              mat <- series(x, ind[[lev]])
              apply(mat, 1, FUN)
            }))
            
            row.names(out) <- levels(fac)           
            out
          })

#' @export
#' @rdname splitReduce-methods
setMethod(f="splitReduce", signature=signature(x = "BrainVector", fac="factor", FUN="missing"),
          def=function(x, fac, FUN) {
            if (length(fac) != prod(dim(x)[1:3])) {
              stop(paste("fac must have as many elements as the number of voxels"))
            }
            
            ind <- split(seq_along(fac), fac)
            out <- do.call(rbind, lapply(names(ind), function(lev) {
              rowMeans(series(x, ind[[lev]]))
            }))
            
            row.names(out) <- levels(fac)           
            out
          })



#' @export
#' @rdname splitScale-methods
setMethod(f="splitScale", signature=signature(x = "matrix", f="factor", center="logical", scale="logical"),
          def=function(x, f, center=TRUE, scale=TRUE) {
            if (length(f) != nrow(x)) {
              stop(paste("x must be same length as split variable"))
            }
            
            out <- matrix(0, nrow(x), ncol(x))
            ind <- split(seq_along(f), f)
            
            for (lev in names(ind)) {
              keep <- ind[[lev]]
              xs <- scale(x[keep,,drop=FALSE], center=center, scale=scale)
              out[keep,] <- xs         
            }
            
            out
          })


#' @export
#' @rdname splitScale-methods
setMethod(f="splitScale", signature=signature(x = "matrix", f="factor", center="missing", scale="missing"),
          def=function(x, f) {
            callGeneric(x,f, TRUE, TRUE)
          })


#' .isExtension
#' @rdname internal-methods
#' @keywords internal
.isExtension <- function(fname, extension) {
  last <- substr(fname, nchar(fname)+1 - nchar(extension), nchar(fname))
  return(last==extension)
}

#' .concat4D
#' @rdname internal-methods
#' @keywords internal
.concat4D <- function(x,y, ...) {
	rest <- list(...)
	
	D <- dim(x)[1:3]
	
	lapply(rest, function(z) {
			stopifnot(length(dim(z)) >= 3)
			stopifnot(identical(D, dim(z)[1:3]))
	})

	D4 <- function(vol) { if (length(dim(vol)) == 3) 1 else dim(vol)[4] }
		
	NVOLS <- D4(x) + D4(y)
	

	ndat <- abind(as.array(x), as.array(y), along=4)
	
	new.dim <- c(D, NVOLS)
	
	nspace <- BrainSpace(new.dim, origin=origin(x@space), spacing=spacing(x@space),
			axes=axes(x@space), trans=trans(x@space))
	
	ret <- DenseBrainVector(ndat, nspace)
	
  ## TODO fix me ridiculously slow
	if (length(rest) > 0) {
		for (i in seq_along(rest)) {
			ret <- concat(ret, rest[[i]])
		}
	}
	
	ret
	
}

#' .extract.array 
#' @rdname internal-methods
#' @keywords internal
.extract.array <- function(x, ..., drop=FALSE, indices=list(...)) {
  nindices <- length(indices)
  if (nindices == 0) {
    stop("Argument 'indices' is empty.")
  }
  dims <- names(indices)
  if (is.null(dims)) {
    dims <- seq(length = nindices)
  }
  
  ndim <- length(dim(x))
  if (any(dims < 1 | dims > ndim)) {
    stop("Argument 'dims' is out of bounds [1,", ndim, "]: ",
         paste(dims, collapse = ", "))
  }
  if (is.null(ndim)) {
    stop("Argument 'x' is not an array: ", class(x)[1])
  }

  args <- rep("", ndim)
  for (kk in seq(length = length(indices))) {
    dd <- dims[kk]
    ii <- sprintf("indices[[%d]]", kk)
    args[dd] <- ii
  }

  if (ndim > 1) {
    args <- c(args, "drop=drop")
  }

  args <- paste(args, collapse = ",")
  code <- paste("x[", args, "]", sep = "")
  expr <- parse(text = code)
  eval(expr)
}

#' .gridToIndex3D
#' @rdname internal-methods
#' @importFrom assertthat assert_that
#' @keywords internal
.gridToIndex3D <- function(dimensions, voxmat) {
	assert_that(length(dimensions) == 3)
  if (is.vector(voxmat)) {
    assert_that(length(voxmat) == 3)
    voxmat <- matrix(voxmat, 1,3)
  }
  
  assert_that(ncol(voxmat) == 3)
  gridToIndex3DCpp(dimensions, voxmat)

}

#' .gridToIndex
#' @rdname internal-methods
#' @keywords internal
.gridToIndex <- function(dimensions, vmat) {
	D <- Reduce("*", dimensions, accumulate=TRUE)
	apply(vmat, 1, function(vox) {
		sum(sapply(length(D):2, function(i) {
			D[i-1]*(vox[i]-1)
		})) + vox[1]
	})
	
}

#' .indexToGrid
#' @rdname internal-methods
#' @keywords internal
.indexToGrid <- function(idx, array.dim) {
  assert_that(all(idx > 0 & idx <= prod(array.dim)))
  assert_that(length(array.dim) <= 5)
  indexToGridCpp(idx, array.dim)
  
}



#' .getRStorage
#' @rdname internal-methods
#' @keywords internal
.getRStorage <- function(dataType) {
  if (any(toupper(dataType) == c("BINARY", "BYTE", "UBYTE", "SHORT", "INTEGER", "INT", "LONG"))) {
    "integer"
  } else if (any(dataType == c("FLOAT", "DOUBLE"))) {
    "double"
  } else {
	  stop(paste("unrecognized data type", dataType))
  }
}

# @nord
# .getMMapMode <- function(code) {	
#	if (code == "UNKNOWN") {
#		stop(paste(".getMMapMode: no memory map mode for UNKNOWN data type: ", code))
#	} else if (code == "BINARY") {
#		int8()
#	} else if (code == "UBYTE") {
#		uint8()
#	} else if(code == "SHORT") {
#		int16()
#	} else if(code == "INT") {
#		int32()
#	} else if (code == "FLOAT") {
#		real32()
#	} else if (code == "DOUBLE") {
#		real64()
#	} else {
#		stop(paste(".getMMapMode: unsupported data type: ", code))
#	}
#}
	

#' .getDataStorage
#' @rdname internal-methods
#' @keywords internal
.getDataStorage <- function(code) {
  if (code == 0) {
    return("UNKNOWN")
  } else if (code == 1) {
    return("BINARY")
  } else if (code == 2) {
    return("UBYTE")
  } else if(code == 4) {
    return("SHORT")
  } else if(code == 8) {
    return("INT")
  } else if (code == 16) {
    return("FLOAT")
  } else if (code == 32) {
    return("DOUBLE")
  } else {
    stop(paste("nifti(getDataStorage): unsupported data type: ", code))
  }
}

#' .getDataCode
#' @rdname internal-methods
#' @keywords internal
.getDataCode <- function(dataType) {
  if (dataType == "UNKNOWN") {
    return(0)
  }else if (dataType == "BINARY") {
    return(1)
  } else if (dataType == "UBYTE") {
    return(2)
  } else if(dataType == "SHORT") {
    return(4)
  } else if(dataType == "INT") {
    return(8)
  } else if (dataType == "FLOAT") {
    return(16)
  } else if (dataType == "DOUBLE") {
    return(32)
  } else {
    stop(paste("getDataCode: unsupported data type: ", dataType))
  }
}

#' .getDataSize
#' @rdname internal-methods
#' @keywords internal
.getDataSize <- function(dataType) {
  if (dataType == "BINARY") {
    return(1)
  } else if (dataType == "UBYTE") {
	  return(1)
  } else if (dataType == "UBYTE") {
    return(1)
  }
  else if (dataType == "SHORT") {
    return(2)
  }
  else if (dataType == "INTEGER") {
    return(4)
  }
  else if (dataType == "INT") {
    return(4)
  }
  else if (dataType == "FLOAT") {
    return(4)
  }
  else if (dataType == "DOUBLE") {
    return(8)
  }
  else if (dataType == "LONG") {
    return(8)
  }

  stop(paste("unrecognized data type: ", dataType))
}

#' .getEndian
#' @rdname internal-methods
#' @keywords internal
.getEndian <- function(conn) {
  #try little endian
  endian <- "little"

  hsize <- readBin(conn, integer(), 1, endian=endian)
  if (hsize != 348) {
    # might be bigendian
    endian <- "big"
    seek(conn, 0)
    hsize <- readBin(conn, integer(), 1, endian=endian)
    if (hsize != 348) {
      stop("nifti(getEndian): header size is not 348, invalid header.")
    }
  }

  return(endian)
}

#' @rdname internal-methods
#' @name .niftiExt
#' @keywords internal
.niftiExt <- function(filetype) {

  extensions <- list()

  if (filetype == "nifti-single") {
    extensions[["header"]]  <- "nii"
    extensions[["data"]] <- "nii"
  }
  else if (filetype == "nifti-pair") {
    extensions[["header"]]  <- "hdr"
    extensions[["data"]] <- "img"
  }
  else if (filetype == "nifti-gz") {
    extensions[["header"]]  <- "nii.gz"
    extensions[["data"]] <- "nii.gz"
  } else {
    stop(paste("unsupported filetype: ", filetype))
  }

  return(extensions)
}

#' @rdname internal-methods
#' @keywords internal
.matrixToQuatern <- function(mat) {
  xd <- sqrt(drop(crossprod(mat[1:3,1])))
  yd <- sqrt(drop(crossprod(mat[1:3,2])))
  zd <- sqrt(drop(crossprod(mat[1:3,3])))

  if (xd == 0) { mat[1,1] = 1; mat[2:3,1] = 0; xd = 1; }
  if (yd == 0) { mat[2,2] = 1; mat[c(1,3),2] = 0; yd = 1; }
  if (zd == 0) { mat[3,3] = 1; mat[1:2,3] = 0; zd = 1; }

  rmat = mat[1:3, 1:3]
  rmat[,1] = rmat[,1]/xd
  rmat[,2] = rmat[,2]/yd
  rmat[,3] = rmat[,3]/zd

  ####### skipping orthogonalization of columns

  #################################################

  zd = det(rmat)
  qfac = 1

  if (zd > 0) {
    qfac = 1
  } else {
    qfac = -1
    rmat[1:3,3] = -rmat[1:3,3]
  }

  # compute quaternion parameters

  a = rmat[1,1] + rmat[2,2] + rmat[3,3] + 1

  if (a > .5) {
    a = .5 * sqrt(a)
    b = 0.25 * (rmat[3,2]-rmat[2,3]) / a
    c = 0.25 * (rmat[1,3]-rmat[3,1]) / a
    d = 0.25 * (rmat[2,1]-rmat[1,2]) / a
   } else {
     xd = 1.0 + rmat[1,1] - (rmat[2,2]+rmat[3,3])
     yd = 1.0 + rmat[2,2] - (rmat[1,1]+rmat[3,3])
     zd = 1.0 + rmat[3,3] - (rmat[1,1]+rmat[2,2])
     if( xd > 1.0 ){
       b = 0.5 * sqrt(xd)
       c = 0.25* (rmat[1,2]+rmat[2,1]) / b
       d = 0.25* (rmat[1,3]+rmat[3,1]) / b
       a = 0.25* (rmat[3,2]-rmat[2,3]) / b
     } else if( yd > 1.0 ){
       c = 0.5 * sqrt(yd) ;
       b = 0.25* (rmat[1,2]+rmat[2,1]) / c
       d = 0.25* (rmat[2,3]+rmat[3,2]) / c
       a = 0.25* (rmat[1,3]-rmat[3,1]) / c
     } else {
       d = 0.5 * sqrt(zd) ;
       b = 0.25* (rmat[1,3]+rmat[3,1]) / d
       c = 0.25* (rmat[2,3]+rmat[3,2]) / d
       a = 0.25* (rmat[2,1]-rmat[1,2]) / d
     }
     if( a < 0.0 ){ b=-b ; c=-c ; d=-d; a=-a; }
   }

  return(list(quaternion=c(b,c,d), qfac=qfac))

}


.quaternToMatrix <- function(quat, origin, stepSize, qfac) {
  mat <- matrix(0, 4,4)
  mat[4,] <- c(0,0,0,1)

  a <- 1 - sum(quat^2)
  if (a < 1e-07) {
    a <- 1 /(sqrt(sum(quat^2)))
    quat <- quat*a
    a <- 0
  } else {
    a <- sqrt(a)
  }

  stepSize <- ifelse(stepSize > 0, stepSize, 1)
  xd <- stepSize[1]
  yd <- stepSize[2]
  zd <- stepSize[3]

  if (qfac < 0) {
    zd <- -zd
  }

  b <- quat[1]
  c <- quat[2]
  d <- quat[3]


  mat[1,1] <- (a*a+b*b-c*c-d*d) * xd
  mat[1,2] <- 2 * (b*c-a*d)     * yd
  mat[1,3] <- 2 * (b*d+a*c)     * zd
  mat[2,1] <- 2 * (b*c+a*d)     * xd
  mat[2,2] <- (a*a+c*c-b*b-d*d) * yd
  mat[2,3] <- 2 * (c*d-a*b)     * zd
  mat[3,1] <- 2 * (b*d-a*c)     * xd
  mat[3,2] <- 2 * (c*d+a*b)     * yd
  mat[3,3] <- (a*a+d*d-c*c-b*b) * zd

  mat[1:3,4] <- origin

  return(mat)
}

    
