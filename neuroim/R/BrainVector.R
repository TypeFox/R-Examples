#' @include AllClass.R
{}
#' @include AllGeneric.R
{}
#' @include common.R
{}
#' @include SparseBrainVector.R
{}



## TODO ought to be able to easily create BrainVector from list of vols


.BrainVectorFromMatrix <- function(data, space) {
	nvols <- dim(space)[4]
	nelements <-  prod(dim(space)[1:3])
	
	if ( (dim(data)[1] == nvols) && (dim(data)[2] == nelements) ) {
		#fourth dimension is rows
		DenseBrainVector(t(data), space)        
	} else if ((dim(data)[2] == nvols) && (dim(data)[1] == nelements )) {
		#fourth dimension is columns
		DenseBrainVector(data, space=space)
	} else {
		stop(paste("illegal matrix dimension ", dim(data)))
	}
}

#' makeVector
#' 
#' Construct a \code{\linkS4class{BrainVector}} instance, using default (dense) implementation
#' @param data a four-dimensional \code{array}
#' @param refdata an instance of class \code{\linkS4class{BrainVector}} or \code{\linkS4class{BrainVolume}} containing the reference space for the new vector.
#' @param label a \code{character} string
#' @param source an instance of class \code{\linkS4class{BrainSource}}
#' @return \code{\linkS4class{DenseBrainVector}} instance 
#' @export makeVector
makeVector <- function(data, refdata, source=NULL, label="") {	
	stopifnot(length(dim(refdata)) == 4)
	rspace <- if (ndim(space(refdata)) == 4) {
		dropDim(space(refdata))
	} else if (ndim(space(refdata)) == 3) {
		space(refdata)
	} else {
		stop("refdata must have 3 or 4 dimensions")
	}

	DenseBrainVector(data,addDim(rspace, dim(data)[4]),source, label)
	
}


#' BrainVector
#' 
#' constructor function for virtual class \code{\linkS4class{BrainVector}}
#' 
#' @param data the image data which can be a \code{matrix}, a 4d \code{array}, or a list of \code{BrainVolumes}. 
#'        If the latter, the geometric space of the data \code{BrainSpace} will be inferred from the constituent volumes, 
#'        which must all be identical.
#' @param space a \code{\linkS4class{BrainSpace}} object. Does not ned to be included if \code{data} argument is a list of \code{BrainVolumes}
#' @param mask an optional \code{array} of type \code{logical}
#' @param source an optional \code{\linkS4class{BrainSource}} object
#' @param label a label of type \code{character} 
#' @return a concrete instance of \code{\linkS4class{BrainVector}} class. 
#' If \code{mask} is provided then \code{\linkS4class{SparseBrainVector}}, otherwise \code{\linkS4class{DenseBrainVector}}
#' @export BrainVector
#' @rdname BrainVector-class
BrainVector <- function(data, space=NULL, mask=NULL, source=NULL, label="") {
  if (is.list(data)) {
    space <- space(data[[1]])
    space <- addDim(space, length(data))
    data <- do.call(cbind, lapply(data, function(x) as.vector(x)))
      
  }
    
	if (prod(dim(space)) != length(data)) {
		stop("dimensions of data argument do not match dimensions of space argument")
	}
	if (is.null(mask)) {
		DenseBrainVector(data,space, source, label)
	} else {
		SparseBrainVector(data,space,mask,source, label)
	}
	
}


#' DenseBrainVector
#' 
#' constructor function for class \code{\linkS4class{DenseBrainVector}}
#' 
#' @param data a 4-dimensional \code{array} or a 2-dimension \code{matrix} that is either nvoxels by ntime-points or ntime-points by nvoxels
#' @param space a \code{\linkS4class{BrainSpace}} object
#' @param source an optional \code{\linkS4class{BrainSource}} object
#' @param label a label of type \code{character} 
#' @return \code{\linkS4class{DenseBrainVector}} instance 
#' @export DenseBrainVector
#' @rdname DenseBrainVector-class
DenseBrainVector <- function(data, space, source=NULL, label="") {
	
	if (is.matrix(data)) {
		splen <- prod(dim(space)[1:3])
		data <- if (ncol(data) == splen) {
			t(data)
		} else if (nrow(data) == splen) {
			data
		}
    
    if (length(dim(space)) == 3) {
      ## add 4th dim to space arg
      space <- addDim(space, ncol(data))
    }

		dim(data) <- dim(space)
	}
  
  
	
	if (is.null(source)) {
		meta <- BrainMetaInfo(dim(data), spacing(space), origin(space), "FLOAT", label)
		source <- new("BrainSource", metaInfo=meta)	
	}
	
	new("DenseBrainVector", .Data=data, source=source, space=space)
	
}


#.makeMMap <- function(meta) {
#	mmap(meta@dataFile, mode=.getMMapMode(meta@dataType), off=meta@dataOffset)
	
#}


#' loadData
#' @return an instance of class \code{\linkS4class{BrainVector}} 
#' @param mmap use memory-mapped file
#' @rdname loadData-methods
setMethod(f="loadData", signature=c("BrainVectorSource"), 
		def=function(x, mmap=FALSE) {		
			
			meta <- x@metaInfo
      
      if (mmap) {
        stop("memory mapping not implemented")
      }
			if (mmap && (.Platform$endian != meta@endian)) {
				stop("cannot create memory mapped file when image endianness does not equal OS endianess")
			}
			
			if (mmap && .isExtension(meta@dataFile, ".gz")) {
				stop("cannot memory map to a gzipped file")		
			}
						
			stopifnot(length(meta@Dim) == 4)
						
			nels <- prod(meta@Dim[1:4]) 		
			ind <- x@indices
			
			#nels <- prod(meta@Dim[1:3]) 
      #datlist <- list()
			#for (i in 1:length(ind)) {
			#	offset <- prod(nels * (ind[i]-1)) * meta@bytesPerElement
			#	reader <- dataReader(meta, offset)		
			#	datlist[[i]] <- array(readElements(reader, nels), meta@Dim[1:3])
			#	close(reader)				
			#}
			
			
			#mappedData <- .makeMMap(meta)
			
			reader <- dataReader(meta, 0)	
			arr <- array(readElements(reader, nels), c(meta@Dim[1:4]))
      
			## bit of a hack to deal with scale factors
			if (.hasSlot(meta, "slope")) {
        
        if (meta@slope != 0) {		  
			    arr <- arr*meta@slope
        }
			}
      
			close(reader)
				
			#arr <- abind(datlist, along=4)			
			
      bspace <- BrainSpace(c(meta@Dim[1:3], length(ind)),meta@spacing, meta@origin, meta@spatialAxes, trans(meta))
			DenseBrainVector(arr[,,,ind,drop=FALSE], bspace, x)
			
		})



#' BrainVectorSource
#' 
#' Construct a \code{\linkS4class{BrainVectorSource}} object
#' 
#' @param fileName name of the 4-dimensional image file
#' @param indices the subset of integer volume indices to load -- if \code{NULL} then all volumes will be loaded
#' @param mask image volume indicating the subset of voxels that will be loaded. If provided, function returns \code{\linkS4class{SparseBrainVectorSource}}
#' @return a instance deriving from \code{\linkS4class{BrainVectorSource}}
#' 
#' @details If a \code{mask} is supplied then it should be a \code{\linkS4class{LogicalBrainVolume}} or \code{\linkS4class{BrainVolume}} instance. If the latter, then the mask will be defined by nonzero elements of the volume.
#'
#' @rdname BrainVectorSource
#' @importFrom assertthat assert_that
#' @export 
BrainVectorSource <- function(fileName, indices=NULL, mask=NULL) {
	assert_that(is.character(fileName))
	assert_that(file.exists(fileName))
	
	
	metaInfo <- readHeader(fileName)
	
	if (!is.null(indices) && max(indices) > 1) {
	  assert_that(length(dim(metaInfo)) == 4)
	  assert_that(max(indices) <= dim(metaInfo)[4])
	  assert_that(min(indices) > 0)
	}
	
  if (length(metaInfo@Dim) == 2) {
    stop(paste("cannot create BrainVector with only two dimensions: ", paste(metaInfo@Dim, collapse=" ")))  
  }
	
  if ( length(metaInfo@Dim) == 3) {
		indices <- 1
    metaInfo@Dim <- c(metaInfo@Dim,1)
	} else if (length(metaInfo@Dim) == 4 && is.null(indices)) {
		indices=seq(1, metaInfo@Dim[4])
	}
	

	if (is.null(mask)) {
		new("BrainVectorSource", metaInfo=metaInfo, indices=as.integer(indices))		
	} else {
		SparseBrainVectorSource(metaInfo, as.integer(indices), mask)		
	}
	
}

#' names
#' @rdname names-methods
#' @export
#' @param x the object to get \code{names} of
setMethod("names", signature=c("BrainBucketSource"),
		def=function(x) {
			x@metaInfo@label[x@indices]
		})


#' @rdname names-methods
#' @export
setMethod("names", signature=c("BrainBucket"),
		def=function(x) {
			x@labels
		})

#' Get length of \code{BrainVector}. This is the numbe rof volumes in the volume vector (e.g. the 4th image dimension)
#' @export
#' @rdname length-methods
setMethod("length", signature=c("BrainVector"),
		def=function(x) {
			dim(x)[4]
		})

#' Load data from a \code{\linkS4class{BrainBucketSource}}
#' @param key the name or index of the bucket to load
#' @return an instance of class \code{\linkS4class{BrainVolume}} 
#' @rdname loadData-methods
setMethod(f="loadData", signature=signature("BrainBucketSource"), 
		def=function(x, key) {

			if (is.numeric(key)) {
				labs <- names(x)
				if (any(key < 1) || any(key > length(labs))) {
					stop(paste("illegal index: ", key))
				}
				
				key <- labs[key]							
			}
      
      
			ret <- lapply(key, function(k) {				
				haskey <- exists(k, envir=x@cache, inherits=FALSE)
				if (!haskey) {
					idx <- which(names(x) == k) 
					vol <- loadData(x@sourceList[[idx]])
					assign(k, vol, envir=x@cache)
				} else {					
					vol <- get(k, envir=x@cache, inherits=FALSE)
				}		
				attr(vol, "label") <- k
				vol
			})
	
			if (length(ret) == 1) {
				ret[[1]]
			} else {
				ret
			}
			
		})
		

#' BrainBucketSource
#' 
#' Constructor function for \code{\linkS4class{BrainBucketSource}} class
#' 
#' @param fileName the name of the bucket file
#' @param pattern optional regular expression used to filter the sub-volumes using associated labels
#' @param indices optional set of sub-volume indices to load
#' @export BrainBucketSource
#' @rdname BrainBucketSource-class
BrainBucketSource <- function(fileName, pattern=NULL, indices=NULL) {
	stopifnot(is.character(fileName))
	stopifnot(file.exists(fileName))
	
	
	metaInfo <- readHeader(fileName)
	
	labels <- metaInfo@label
	nvols <- length(labels)	
	
	if (is.null(indices)) {
		indices <- seq_along(labels)
	} else {
		stopifnot(all(indices >0 & indices < nvols))
	}
	
	if (!is.null(pattern)) {
		idx <- grep(pattern, labels)
		
		if (length(idx) < 1) {
			stop(paste("pattern: ", pattern, "does not match any volume labels"))
		}
		
		indices <- intersect(idx, indices)
	}
			
	
	sourceList <- lapply(indices, function(i) { new("BrainVolumeSource", metaInfo=metaInfo, index=as.integer(i)) })	
	new("BrainBucketSource", metaInfo=metaInfo, indices=indices, sourceList=sourceList, cache=new.env(hash=TRUE))	
}

#' BrainBucket
#' 
#' Constructor function for \code{\linkS4class{BrainBucket}} class
#' 
#' @param volumeList a named list of \code{\linkS4class{BrainVolume}} instances
#' @export BrainBucket
#' @rdname BrainBucket-class
#' @examples 
#' vol1 <- BrainVolume(rnorm(24*24*24), BrainSpace(c(24,24,24), c(1,1,1)))
#' vol2 <- BrainVolume(rnorm(24*24*24), BrainSpace(c(24,24,24), c(1,1,1)))
#' vol3 <- BrainVolume(rnorm(24*24*24), BrainSpace(c(24,24,24), c(1,1,1)))
#' vlist <- list(vol1,vol2,vol3)
#' names(vlist) <- paste0("V", 1:3)
#' bucket <- BrainBucket(vlist)
#' all.equal(dim(bucket[[1]]), dim(vol1))
#' @return an instance of class \code{\linkS4class{BrainBucket}}
#' 
BrainBucket <- function(volumeList) {
  
  isvol <- sapply(volumeList, function(vol) inherits(vol, "BrainVolume"))
  
  if (length(volumeList) < 1 || any(!isvol)){
    stop("BrainBucket: 'volumeList' must be a nonempty list of instances of or extending 'BrainVolume'")
  }
  
  vnames <- names(volumeList)
  if (is.null(vnames)) {
    vnames <- paste0("V", 1:length(volumeList))
  } else if (any(vnames == "")) {
    whichEmpty <- which(vnames == "")
    vnames[whichEmpty] <- paste0("V", whichEmpty)
  }
  
  sp <- space(volumeList[[1]])
  D <- c(dim(sp), length(volumeList))
  
  meta <- BrainMetaInfo(D, spacing(sp), origin=origin(sp), dataType="FLOAT", label=vnames, spatialAxes=axes(sp))

  sourceList <- lapply(1:length(volumeList), function(i) { 
    volumeList[[i]]@source

  })
  
  bsource <- new("BrainBucketSource", metaInfo=meta, indices=1:length(volumeList), sourceList=sourceList, cache=new.env(hash=TRUE))
  
  new("BrainBucket", source=bsource,space=addDim(sp, length(volumeList)), labels=vnames, data=volumeList)
}
  
  



#' loadBucket
#' 
#' load a BrainBucket object from file
#' 
#' @param fileName the name of the file to load
#' @param pattern optional regular expression used to filter the sub-volumes using associated labels
#' @param indices optional set of sub-volume indices to load
#' @export loadBucket
loadBucket <- function(fileName, pattern=NULL, indices=NULL) {
	bsource <- BrainBucketSource(fileName, pattern, indices)
	
	meta <- bsource@metaInfo	
	labels <- meta@label
	idx <- bsource@indices
	
	D <- c(meta@Dim[1:3], length(idx))
	bspace <- BrainSpace(D, meta@spacing, meta@origin, meta@spatialAxes)
  
  ## TODO chceck: does not actually provide data ....
	buck <- new("BrainBucket", source=bsource, space=bspace, labels=labels[idx])
}

#' loadVolList
#' 
#' load a list of image volumes and return a \code{\linkS4class{BrainVector}} instance
#' 
#' @param fileNames a list of files to load
#' @param mask an optional mask indicating subset of voxels to load
#' @return an instance of class \code{\linkS4class{BrainVector}}
#' @export loadVolumeList
loadVolumeList <- function(fileNames, mask=NULL) {
	stopifnot(all(sapply(fileNames, file.exists)))
	metaInfo <- lapply(fileNames, readHeader)
	
	dims <- do.call(rbind, lapply(metaInfo, dim))
	if (!all(sapply(1:nrow(dims), function(i) all.equal(dims[1,], dims[i,])))) {
		stop("list of volumes must all have same dimensions")
	}
	
	if (!all(apply(dims, 1, length) == 3)) {
		stop("all volumes in list must have dim = 3")
	}
	
	nvols <- length(fileNames)	
	sourceList <- lapply(fileNames, function(fname) {
		BrainVolumeSource(fname, 1)
	})

	vols <- lapply(sourceList, loadData)
	if (is.null(mask)) {
		mat <- do.call(cbind, vols)
		dspace <- addDim(space(vols[[1]]), length(vols))	
		DenseBrainVector(mat, dspace, label=sapply(metaInfo, function(m) m@label))
	} else {
		mat <- do.call(cbind, vols)
		dspace <- addDim(space(vols[[1]]), length(vols))
		if (is.vector(mask)) {
			## mask supplied as index vector, convert to logical
			M <- array(logical(prod(dim(dspace)[1:3])), dim(dspace)[1:3])
			M[mask] <- TRUE
			mask <- M
		} else {
			mask <- as.logical(mask)
		}
		
		
		SparseBrainVector(mat[mask,], dspace, mask=mask, label=sapply(metaInfo, function(m) m@label))
		
	}
}

#' extract labeled volume from \code{BrainBucket}
#' @param x the object
#' @param i the first index
setMethod(f="[[", signature=signature(x="BrainBucket", i = "index", j = "missing"),
		def=function(x, i) {
		  x@data[[i]]
			#loadData(x@source, i)
		})


#' extract labeled volume from \code{BrainBucket}
#' @export
#' @param x the object
#' @param i first index
setMethod(f="[", signature=signature(x="BrainBucket", i = "index", j = "missing"),
          def=function(x, i) {
            x@data[i]
            #loadData(x@source, i)
          })




setAs("DenseBrainVector", "array", function(from) from@.Data)

setAs("BrainVector", "array", function(from) from[,,,])

#' show a \code{BrainVectorSource}
#' @param object the object
#' @export
setMethod(f="show",
		signature=signature(object="BrainVectorSource"),
		def=function(object) {
			cat("an instance of class",  class(object), "\n\n")
			cat("   indices: ", object@indices, "\n\n")
			cat("   metaInfo: \n")
			show(object@metaInfo)
			cat("\n\n")
			
		})




#' show a \code{BrainVector}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("BrainVector"),
          def=function(object) {
            sp <- space(object)
            cat(class(object), "\n")
            cat("  Type           :", class(object), "\n")
            cat("  Dimension      :", dim(object), "\n")
            cat("  Spacing        :", paste(paste(sp@spacing[1:(length(sp@spacing)-1)], " X ", collapse=" "), 
                                            sp@spacing[length(sp@spacing)], "\n"))
            cat("  Origin         :", paste(paste(sp@origin[1:(length(sp@origin)-1)], " X ", collapse=" "), 
                                            sp@origin[length(sp@origin)], "\n"))
            cat("  Axes           :", paste(sp@axes@i@axis, sp@axes@j@axis,sp@axes@k@axis), "\n")
            cat("  Coordinate Transform :", zapsmall(sp@trans), "\n")
            
          }
)



 
#' @rdname eachVolume-methods
#' @export
setMethod(f="eachVolume", signature=signature(x="BrainVector", FUN="function", withIndex="missing", mask="missing"),
		def=function(x, FUN, withIndex, mask, ...) {
			lapply(1:(dim(x)[4]), function(tt) FUN(x[,,,tt], ...))				
		})


#' @rdname eachVolume-methods
#' @export
setMethod(f="eachVolume", signature=signature(x="BrainVector", FUN="function", withIndex="missing", mask="BrainVolume"),
          def=function(x, FUN, withIndex, mask, ...) {
            mask.idx <- which(mask > 0)
            lapply(1:(dim(x)[4]), function(tt) {
              vals <- x[,,,tt]
              FUN(vals[mask.idx], ...)
            })
          })

 
#' @rdname eachVolume-methods
#' @export
setMethod(f="eachVolume", signature=signature(x="BrainVector", FUN="function", withIndex="missing", mask="missing"),
          def=function(x, FUN, withIndex, mask, ...) {   
            lapply(1:(dim(x)[4]), function(tt) {
              vals <- x[,,,tt]
              FUN(vals, ...)
            })
          })



#' @rdname eachVolume-methods
#' @export
setMethod(f="eachVolume", signature=signature(x="BrainBucket", FUN="function", withIndex="missing",mask="missing"),
		def=function(x, FUN, withIndex, ...) {
			lapply(1:(dim(x)[4]), function(tt) FUN(x[[tt]], ...))				
		})



#' @rdname eachVolume-methods
#' @export
setMethod("eachVolume", signature=signature(x="BrainBucket", FUN="function", withIndex="logical"),
		def=function(x, FUN, withIndex, ...) {
			lapply(1:(dim(x)[4]), function(tt) {					
						vol <- x[[tt]]
						if (withIndex) FUN(vol,tt,...) else FUN(vol,...)
					})
		})




#' @rdname eachVolume-methods
#' @export
setMethod("eachVolume", signature=signature(x="BrainVector", FUN="function", withIndex="logical"),
		def=function(x, FUN, withIndex, ...) {
			lapply(1:(dim(x)[4]), function(tt) {					
						vol <- x[,,,tt]
						if (withIndex) FUN(vol,tt,...) else FUN(vol,...)
					})
		})



#' @rdname subVector-methods
#' @export
setMethod(f="subVector", signature=signature(x="DenseBrainVector", i="numeric"),
          def=function(x, i) {
            xs <- space(x)
            dat <- x@.Data[,,,i]
            
            newdim <- c(dim(x)[1:3], length(i))
            bspace <- BrainSpace(newdim, spacing=spacing(xs), origin=origin(xs), axes(xs), trans(xs))
            DenseBrainVector(dat, bspace)
          })


#' @rdname takeVolume-methods
#' @param merge concatenate extracted volumes
#' @export
setMethod(f="takeVolume", signature=signature(x="BrainVector", i="numeric"),
		def=function(x, i, merge=FALSE) {
			## TODO this is VERY slow
			## TODO should be renamed "volSlice"
		  
		  
			xs <- space(x)
			bspace <- BrainSpace(dim(x)[1:3], spacing=spacing(xs), origin=origin(xs), axes(xs), trans(xs))
			
			makevol <- function(i) {				
				BrainVolume(x@.Data[,,,i], bspace)
			}
			
			res <- lapply(i, makevol)
			
			if (length(res) > 1 && merge) {
				res <- do.call("concat", res)				
			}
			
			if (length(res) == 1) {
			  ## TODO should be consistent, e.g. always return list
				res[[1]]
			} else {
				res
			}											
		})


#' @rdname eachSeries-methods
#' @export
setMethod(f="eachSeries", signature=signature(x="BrainVector", FUN="function", withIndex="missing"),
		def=function(x, FUN, withIndex=FALSE, ...) {
			
			NX <- dim(x)[1]
			NY <- dim(x)[2]
			NZ <- dim(x)[3]
			ret <- vector("list", prod(NX, NY, NZ))
			index <- 1
			for (i in 1:NZ) {
				for (j in 1:NY) {
					for (k in 1:NX) {
						ret[[index]] <- FUN(x[k,j,i,])
						index <- index+1
					}
				}
			}
			
			ret
			
		})




#' loadVector
#' 
#' load an image volume from a file
#' 
#' @param fileName the name of the file to load
#' @param indices the indices of the sub-volumes to load (e.g. if the file is 4-dimensional)
#' @param mask a mask defining the spatial elements to load 
#' @return an \code{\linkS4class{BrainVector}} object
#' @export 
loadVector  <- function(fileName, indices=NULL, mask=NULL) {
	src <- BrainVectorSource(fileName, indices, mask)
	loadData(src)
}




#setMethod("sliceMeans", signature(x="BrainVector"),
#          function(x) {
#             t(colMeans(x, dims=2))
#           })


#' @rdname concat-methods
#' @export
setMethod(f="concat", signature=signature(x="BrainVector", y="BrainVolume"),
		def=function(x,y, ...) {
			.concat4D(x,y,...)			
		})


#' @rdname concat-methods
#' @export
setMethod(f="concat", signature=signature(x="BrainVolume", y="BrainVector"),
  def=function(x,y, ...) {
    .concat4D(x,y,...)			
  })

#' @rdname scaleSeries-methods
#' @export
setMethod(f="scaleSeries", signature=signature(x="BrainVector", center="logical", scale="logical"),
          def=function(x, center, scale) {
            M <- as.matrix(x)
            Ms <- scale(t(M), center, scale)
            BrainVector(Ms, space(x))
             		
          })

#' @rdname scaleSeries-methods
#' @export
setMethod(f="scaleSeries", signature=signature(x="BrainVector", center="missing", scale="logical"),
          def=function(x, center, scale) {
            callGeneric(x, TRUE, scale)
          })


#' @rdname scaleSeries-methods
#' @export
setMethod(f="scaleSeries", signature=signature(x="BrainVector", center="missing", scale="missing"),
          def=function(x, center, scale) {
            callGeneric(x, TRUE, TRUE)
          })

#' @rdname scaleSeries-methods
#' @export
setMethod(f="scaleSeries", signature=signature(x="BrainVector", center="logical", scale="missing"),
          def=function(x, center, scale) {
            callGeneric(x, center, TRUE)
          })


#' @rdname concat-methods
#' @export
setMethod(f="concat", signature=signature(x="BrainVector", y="BrainVector"),
		def=function(x,y,...) {
			.concat4D(x,y,...)
		})


#' @rdname series-methods
#' @export
setMethod("series", signature(x="BrainVector", i="matrix"),
		def=function(x,i) {
      ## TODO not necessary, has to be matrix based on type of i arg...
			if (!is.matrix(i) && length(i) == 3) {
				i <- matrix(i, 1, 3)
			}
			
			stopifnot(ncol(i) == 3)
			apply(i, 1, function(i) x[i[1], i[2], i[3],])
	
		})



#' @rdname series-methods
#' @export
setMethod("series", signature(x="BrainVector", i="numeric"),
		def=function(x,i, j, k) {	
			if (missing(j) && missing(k)) {
				vdim <- dim(x)[1:3]
				mat <- arrayInd(i, vdim)
				apply(mat, 1, function(i) x[i[1], i[2], i[3],])			
			} else {
				x[i,j,k,]	
			}
		})



#' @describeIn seriesIter get a series iterator for a \code{\linkS4class{BrainVector}} instance
#' @export
setMethod(f="seriesIter", signature=signature(x="BrainVector"), 
		def=function(x) {
			len <- prod(dim(x)[1:3])
			vdim <- dim(x)[1:3]
			i <- 1
			nextEl <- function() {
				if (i <= len) {
					vox <- .indexToGrid(i, vdim)
					i <<- i + 1
					x[vox[1], vox[2], vox[3],]
					
				} else {
					stop("StopIteration") 
				}		
			}
			
			hasNx <- function() {
				i <= len
			}
			
			obj <- list(nextElem = nextEl, hasNext=hasNx) 
			class(obj) <- c("seriesIter", "abstractiter", "iter") 
			obj
			
			
		})


#' @export
setAs(from="DenseBrainVector", to="matrix",
		function(from) {
			data <- from@.Data
			dm <- dim(data)
			d123 <- prod(dm[1:3])
			d4 <- dm[4]
			
			dim(data) <- c(d123,d4)
			return(data)
			
		})


#' convert a \code{BrainVector} to \code{list} of volumes. 
#' @rdname as.list-methods
#' @param x the object
#' @export 
setMethod(f="as.list", signature=signature(x = "BrainVector"), def=function(x) {
  out = list()
  for (i in 1:dim(x)[4]) {
    out[[i]] <- takeVolume(x,i)
  }
  
  out
})


#' convert a \code{DenseBrainVector} to a matrix
#' 
#' @rdname as.matrix-methods
#' @param x the object
#' @export 
setMethod(f="as.matrix", signature=signature(x = "DenseBrainVector"), def=function(x) {
			as(x, "matrix")						
		})

 
#' @rdname as.sparse-methods
#' @export
setMethod(f="as.sparse", signature=signature(x="DenseBrainVector", mask="LogicalBrainVolume"),
          def=function(x, mask) {
            assert_that(all(dim(x)[1:3] == dim(mask)))
            assert_that(all(spacing(x) == spacing(mask)))
            
            vdim <- dim(x)[1:3]
            dat <- as.matrix(x)[mask == TRUE,]
            bvec <- SparseBrainVector(dat, space(x), mask)
            
          })

 
#' @rdname as.sparse-methods
setMethod(f="as.sparse", signature=signature(x="DenseBrainVector", mask="numeric"),
		def=function(x, mask) {
			vdim <- dim(x)[1:3]
			m <- array(0, vdim)
			m[mask] <- TRUE
			
			logivol <- LogicalBrainVolume(m, dropDim(space(x)))
			
			dat <- as.matrix(x)[mask,]
			
			bvec <- SparseBrainVector(dat, space(x), logivol)
			
		})

# @export gridToIndex
# @rdname gridToIndex-methods
# setMethod(f="gridToIndex", signature=signature(x="BrainVector", coords="matrix"),
#		def=function(x, coords) {
#			stopifnot(ncol(coords) == 4)
#			array.dim <- dim(x)
#			ind3d <- .gridToIndex(dim(x), coords[,1:3])
#		})

#setMethod(f="takeSeries", signature=signature(x="BrainVector", indices="numeric"),
#		def=function(x, indices) {
#			
#			D <- dim(x)[1:3]
#			if (all.equal(dim(indices), D)) {
#				indices <- which(indices>0)
#			}
			
#			vox <- t(sapply(indices, .indexToGrid, D)) 
			
#			apply(vox, 1, function(v) {
#				x[v[1], v[2], v[3],]
#			})
			
#		})

#setMethod(f="takeSeries", signature=signature(x="BrainVector", indices="BrainVolume"),
#		def=function(x, indices) {				
#			D <- dim(x)[1:3]
#			stopifnot(all.equal(dim(indices), D))
#			
#			idx <- which(indices > 0)
#			callGeneric(x, idx)			
#			
#		})          



#' @export
#' @rdname writeVector-methods
setMethod(f="writeVector",signature=signature(x="BrainVector", fileName="character", format="missing", dataType="missing"),
		def=function(x, fileName) {
			write.nifti.vector(x, fileName)           
		})


#' @export 
#' @rdname writeVector-methods
setMethod(f="writeVector",signature=signature(x="BrainVector", fileName="character", format="character", dataType="missing"),
		def=function(x, fileName, format) {
			if (toupper(format) == "NIFTI" || toupper(format) == "NIFTI1" || toupper(format) == "NIFTI-1") {
				callGeneric(x, fileName)
			} else {
				stop(paste("sorry, cannot write format: ", format))
			}      
		})


#' @export writeVector
#' @rdname writeVector-methods
#' @aliases writeVector,BrainVector,character,missing,character,ANY-method
setMethod(f="writeVector",signature=signature(x="BrainVector", fileName="character", format="missing", dataType="character"),
		def=function(x, fileName, dataType) {
			write.nifti.vector(x, fileName, dataType)   
			
		})




#loadSeries <- function(filenames, indices, volidx=NULL, reduce=T, demean=F, verbose=F, bulk.thresh=100 ) {
#stopifnot(all(sapply(filenames, .isNIFTI)))

#ret <- lapply(filenames, function(filename) {
#			bvec <- loadVector(filename, indices)
#			
#			if (demean) {
#				cmeans <- colMeans(retmat)
#				retmat <- sweep(retmat, 2, cmeans)
#			}
#			
#			if (reduce) {
#				retmat <- rowMeans(retmat)
#			}
#			
#			if (is.matrix(retmat) && NCOL(retmat) == 1) {
#				retmat <- retmat[,1]
#			}
#			
#			close(conn)
#			
#			retmat
#		})
#	
#	
#	if (length(ret) == 1) {
#		ret[[1]]
#	} else {
#		ret
#	}
#}



#.loadSparseVector <- function(filename, mask) {
#	if (!.isNIFTI(filename)) {
#		stop("only support NIFTI files at present")
#	}
#	
#	nfile <- NIFTIFile(filename)
#	header <- readHeader(nfile)
#	ddim <- dataDim(header)
#	
#	if (length(ddim) != 4) {
#		stop("Error: file does not have 4 dimensions, which is required for a BrainVector object")
#	}
#	
#	
#	if (inherits(mask, "array") || inherits(mask, "BrainData")) {
#		# should check that dimensions are equal
#		ind <- which(mask > 0)
#	} else {
#		stop("mask must be of type array or BrainData")
#	}
#	
#	mat <- readData(nfile, ind)
#	
#	space <- createSpace(header)
#	SparseBrainVector(mat, space, indices=ind)
#}

