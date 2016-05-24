#' @import hash
#' @importFrom assertthat assert_that
#' @importFrom Matrix sparseVector
#' @importFrom yaImpute ann
NULL

#' @include AllClass.R
{}
#' @include BrainVector.R
{}
#' @include common.R
{}
#' @include BrainMetaInfo.R
{}
#' @include NIFTI_IO.R
{}
#' @include Axis.R
{}

#' makeVolume
#' 
#' Construct a \code{\linkS4class{BrainVolume}} instance, using default (dense) implementation
#' @param data an optional one- or three-dimensional \code{vector} or \code{array}
#' @param refvol an instance of class \code{\linkS4class{BrainVolume}} containing the reference space for the new volume.
#' @param label an optional \code{character} string
#' @param source an optional instance of class \code{\linkS4class{BrainSource}}
#' @param indices an optional 1d vector of indices in to the 3d space
#' @return \code{\linkS4class{DenseBrainVolume}} instance 
#' @examples
#' bspace <- BrainSpace(c(64,64,64), spacing=c(1,1,1))
#' dat <- array(rnorm(64*64*64), c(64,64,64))
#' bvol <- BrainVolume(dat,bspace, label="test")
#' bvol2 <- makeVolume(dat, bvol)
#' all.equal(as.array(bvol),as.array(bvol2))
#' data <- 1:10
#' indices = seq(1,1000, length.out=10)
#' bvol3 <- makeVolume(data,bvol,indices=indices)
#' sum(bvol3) == sum(data)
#' @export makeVolume
makeVolume <- function(data=NULL, refvol, source=NULL, label="", indices=NULL) {
  if (is.null(data)) {
	  DenseBrainVolume(array(0, dim(refvol)),space(refvol),source, label,indices)	
  } else {
    DenseBrainVolume(data,space(refvol),source, label,indices)  
  }
}

#' BrainVolume
#' 
#' Construct a \code{\linkS4class{BrainVolume}} instance, using default (dense) implementation
#' @param data a three-dimensional \code{array}
#' @param space an instance of class \code{\linkS4class{BrainSpace}}
#' @param source an instance of class \code{\linkS4class{BrainSource}}
#' @param label a \code{character} string to identify volume
#' @param indices an 1D vector that gives the linear indices of the associated \code{data} vector
#' @return a \code{\linkS4class{DenseBrainVolume}} instance 
#' @examples
#' bspace <- BrainSpace(c(64,64,64), spacing=c(1,1,1))
#' dat <- array(rnorm(64*64*64), c(64,64,64))
#' bvol <- BrainVolume(dat,bspace, label="test")
#' print(bvol) 
#' @export BrainVolume
#' @rdname BrainVolume
BrainVolume <- function(data, space, source=NULL, label="", indices=NULL) {
	DenseBrainVolume(data,space,source=source, label=label, indices=indices)	
}

#' DenseBrainVolume
#' 
#' Construct a \code{\linkS4class{DenseBrainVolume}} instance
#' @param data a three-dimensional \code{array}
#' @param space an instance of class \code{\linkS4class{BrainSpace}}
#' @param source an instance of class \code{\linkS4class{BrainSource}}
#' @param label a \code{character} string
#' @param indices an optional 1-d index vector
#' @return \code{\linkS4class{DenseBrainVolume}} instance 
#' @export DenseBrainVolume
#' @rdname DenseBrainVolume-class
DenseBrainVolume <- function(data, space, source=NULL, label="", indices=NULL) {
	
	if (length(dim(space)) != 3) {
		stop("DenseBrainVolume: space argument must have three dimensions")
	} 
  
  if (is.matrix(data)) {
    if (nrow(data) == 1 || ncol(data) == 1) {
      data <- as.vector(data)
    }
  }
	
	if (length(data) == prod(dim(space)) && is.vector(data)) {
		dim(data) <- dim(space)
	}
	
	if (ndim(space) != 3) {
		stop("DenseBrainVolume: space argument must have three dimensions")
	} 
	
	if (!all(dim(space) == dim(data))) {
		stop("DenseBrainVolume: data and space argument must have equal dimensions")
	} 
	
	if (is.null(source)) {
		meta <- BrainMetaInfo(dim(space), spacing(space), origin(space), "FLOAT", label)
		source <- new("BrainSource", metaInfo=meta)
	}
	
	if (!is.null(indices)) {
		newdat <- array(0, dim(space))
		newdat[indices] <- data
		data <- newdat
	}
	
			
	#new("DenseBrainVolume", .Data=data, source=source, space=space)
	new("DenseBrainVolume", data, source=source, space=space)

}

#' ClusteredBrainVolume
#' 
#' Construct a \code{\linkS4class{ClusteredBrainVolume}} instance
#' @param mask an instance of class \code{\linkS4class{LogicalBrainVolume}}
#' @param clusters a vector of clusters ids with length equal to number of nonzero voxels in mask \code{mask}
#' @param labelMap an optional \code{list} that maps from cluster id to a cluster label, e.g. (1 -> "FFA", 2 -> "PPA")
#' @param source an optional instance of class \code{\linkS4class{BrainSource}}
#' @param label an optional \code{character} string used to label of the volume
#' @return \code{\linkS4class{ClusteredBrainVolume}} instance 
#' @export ClusteredBrainVolume
#' @examples
#' 
#' bspace <- BrainSpace(c(16,16,16), spacing=c(1,1,1))
#' grid <- indexToGrid(bspace, 1:(16*16*16))
#' kres <- kmeans(grid, centers=10)
#' mask <- BrainVolume(rep(1, 16^3),bspace)
#' clusvol <- ClusteredBrainVolume(mask, kres$cluster)
#' @rdname ClusteredBrainVolume-class
ClusteredBrainVolume <- function(mask, clusters, labelMap=NULL, source=NULL, label="") {
  mask <- as(mask, "LogicalBrainVolume")
  space <- space(mask)
  ids <- sort(unique(clusters)) 
  
  stopifnot(length(clusters) == sum(mask))
  
  if (length(ids) == 1) {
    warning("clustered volume only contains 1 partition")
  }
  
  if (is.null(labelMap)) {   
    labs <- paste("Clus_", ids, sep="")
    labelMap <- as.list(ids)
    names(labelMap) <- labs    
  } else {
    stopifnot(length(labelMap) == length(ids))
    stopifnot(all(unlist(labelMap) %in% ids))
  }
  
  if (is.null(source)) {
    meta <- BrainMetaInfo(dim(space), spacing(space), origin(space), "FLOAT", label)
    source <- new("BrainSource", metaInfo=meta)
  }
  

  clus.idx <- which(mask == TRUE)
  clus.split <- split(clus.idx, clusters)
  clus.names <- names(clus.split)
  clusterMap <- hash()
  
  for (i in 1:length(clus.split)) {
    clusterMap[[clus.names[[i]]]] <- clus.split[[clus.names[[i]]]]
  }
  
  new("ClusteredBrainVolume", mask=mask, clusters=as.integer(clusters), labelMap=labelMap, clusterMap=clusterMap, source=source, space=space)
}
  
  
  

#' SparseBrainVolume
#' 
#' Construct a \code{\linkS4class{SparseBrainVolume}} instance
#' @param data a numeric vector 
#' @param space an instance of class \code{\linkS4class{BrainSpace}}
#' @param indices a index vector indicating the 1-d coordinates of the data values
#' @param source an instance of class \code{\linkS4class{BrainSource}}
#' @param label a \code{character} string
#' @return \code{\linkS4class{SparseBrainVolume}} instance 
#' @export SparseBrainVolume
#' @details 
#' Image data is backed by \code{Matrix::sparseVector}. 
#' @examples
#' data <- 1:10
#' indices <- seq(1,1000, length.out=10)
#' bspace <- BrainSpace(c(64,64,64), spacing=c(1,1,1))
#' sparsevol <- SparseBrainVolume(data,bspace,indices=indices)
#' densevol <- BrainVolume(data,bspace,indices=indices)
#' sum(sparsevol) == sum(densevol)
#' 
#' 
#' 
#' @rdname SparseBrainVolume-class
SparseBrainVolume <- function(data, space, indices=NULL, source=NULL, label="") {
  if (length(indices) != length(data)) {
    stop(paste("length of 'data' must equal length of 'indices'"))
  }
  
  sv <- Matrix::sparseVector(x=data, i=indices, length=prod(dim(space)))
  
  if (is.null(source)) {
    meta <- BrainMetaInfo(dim(space), spacing(space), origin(space), "FLOAT", label)
    source <- new("BrainSource", metaInfo=meta)
  }
  
  new("SparseBrainVolume", data=sv, source=source, space=space)
}
  
  

#' LogicalBrainVolume
#' 
#' Construct a \code{\linkS4class{LogicalBrainVolume}} instance
#' @param data a three-dimensional \code{array}, a 1D vector with length equal to \code{prod(dim(space))}, or a set of \code{indices} where elements are \code{TRUE}
#' @param space an instance of class \code{\linkS4class{BrainSpace}}
#' @param source an instance of class \code{\linkS4class{BrainSource}}
#' @param label a \code{character} string
#' @param indices an optional 1-d index vector
#' @return \code{\linkS4class{LogicalBrainVolume}} instance 
#' @export LogicalBrainVolume
#' @rdname LogicalBrainVolume-class
LogicalBrainVolume <- function(data, space, source=NULL, label="", indices=NULL) {
	
	if (is.null(dim(data)) && length(data) == prod(dim(space))) {
		data <- array(data, dim(space))
	} else if (is.null(dim(data)) && !is.null(indices)) {
	  newdat <- array(FALSE, dim(space))
	  newdat[indices] <- data
	  data <- newdat
	}
	
	if (length(dim(data)) != 3) {
		stop("LogicalBrainVolume: data argument must have three dimensions")
	} 
	
	if (ndim(space) != 3) {
		stop("LogicalVolume: space argument must have three dimensions")
	} 
	
	
	if (!is.logical(data)) {
		D <- dim(data)
		data <- as.logical(data)
		dim(data) <- D
	}
	
	if (is.null(source)) {
		meta <- BrainMetaInfo(dim(data), spacing(space), origin(space), "BINARY", label)
		source <- new("BrainSource", metaInfo=meta)
	}
	
	new("LogicalBrainVolume", source=source, .Data=data, space=space)
}



#' conversion from DenseBrainVolume to array
#' @rdname as-methods
#' @name as
setAs(from="DenseBrainVolume", to="array", def=function(from) from@.Data)



#' conversion from \code{SparseBrainVolume} to \code{array}
#' @rdname as-methods
#' @name as
setAs(from="SparseBrainVolume", to="array", def=function(from) {
  vals <- as.numeric(from@data)
  array(vals, dim(from))
})


#' conversion from SparseBrainVolume to numeric
#' @rdname as-methods
#' @name as
setAs(from="SparseBrainVolume", to="numeric", def=function(from) {
  as.numeric(from@data)
})


#' Convert SparseBrainVolume to numeric
#' 
#' @rdname as.numeric-methods
#' @param x the object to convert
#' @export 
setMethod(f="as.numeric", signature=signature(x = "SparseBrainVolume"), def=function(x) {
  as(x, "numeric")			
})


#' conversion from BrainVolume to LogicalBrainVolume
#' @rdname as-methods
#' @name as
setAs(from="BrainVolume", to="LogicalBrainVolume", def=function(from) {
	LogicalBrainVolume(as.array(from), space(from), from@source)
})

#' conversion from DenseBrainVolume to LogicalBrainVolume
#' @name as
#' @rdname as-methods    
setAs(from="DenseBrainVolume", to="LogicalBrainVolume", def=function(from) {
	LogicalBrainVolume(as.array(from), space(from), from@source)
})

#' conversion from ClusteredBrainVolume to LogicalBrainVolume
#' @name as
#' @rdname as-methods
setAs(from="ClusteredBrainVolume", to="DenseBrainVolume", def=function(from) {
  data = from@clusters
  indices <- which(from@mask == TRUE)
  DenseBrainVolume(data, space(from), from@source, indices=indices)
})

#' conversion from BrainVolume to array
#' @rdname as-methods
#' @name as
setAs(from="BrainVolume", to="array", def=function(from) from[,,])

#' show a \code{BrainVolume}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("BrainVolume"),
          def=function(object) {
            sp <- space(object)
            cat("BrainVolume\n")
            cat("  Type           :", class(object), "\n")
            cat("  Dimension      :", dim(object), "\n")
            cat("  Spacing        :", paste(paste(signif(sp@spacing[1:(length(sp@spacing)-1)],2), " X ", collapse=" "), 
                                            sp@spacing[length(sp@spacing)], "\n"))
            cat("  Origin         :", paste(paste(signif(sp@origin[1:(length(sp@origin)-1)],2), " X ", collapse=" "), 
                                            sp@origin[length(sp@origin)], "\n"))
            cat("  Axes           :", paste(sp@axes@i@axis, sp@axes@j@axis,sp@axes@k@axis), "\n")
           
          }
)


#' load a BrainVolume
#' @export loadData
#' @rdname loadData-methods
## TODO reduce code duplication with loadData#BrainVectorSource
setMethod(f="loadData", signature=c("BrainVolumeSource"), 
		def=function(x) {
			
			meta <- x@metaInfo
			nels <- prod(meta@Dim[1:3]) 
			
			### for brain buckets, this offset needs to be precomputed ....
			offset <- (nels * (x@index-1)) * meta@bytesPerElement
			
			reader <- dataReader(meta, offset)
			dat <- readElements(reader, nels)
			## bit of a hack to deal with scale factors
      if (.hasSlot(meta, "slope")) {
        if (meta@slope != 0) {  	  
          dat <- dat*meta@slope
        }
      }
      
			close(reader)
			arr <- array(dat, meta@Dim[1:3])
			
			bspace <- BrainSpace(meta@Dim[1:3], meta@spacing, meta@origin, meta@spatialAxes, trans(meta))
			DenseBrainVolume(arr, bspace, x)
					
		})

#' Constructor for BrainVolumeSource
#' 
#' @param input the input file name
#' @param index the image subvolume index
#' @export 
#' @rdname BrainVolumeSource-class
BrainVolumeSource <- function(input, index=1) {
	stopifnot(index >= 1)
	stopifnot(is.character(input))
	
	if (!file.exists(input)) {
		candidates <- Sys.glob(paste(input, "*", sep=""))
		if (length(candidates) > 0) {
			input <- candidates[1]
		}
	}
	
	stopifnot(file.exists(input))
		
	metaInfo <- readHeader(input)
	
	if ( length(metaInfo@Dim) < 4 && index > 1) {
		stop("index cannot be greater than 1 for a image with dimensions < 4")
	}
	
	new("BrainVolumeSource", metaInfo=metaInfo, index=as.integer(index))								
	
}

#' Load an image volume from a file
#' 
#' @param fileName the name of the file to load
#' @param index the index of the volume (e.g. if the file is 4-dimensional)
#' @return an instance of the class \code{\linkS4class{DenseBrainVolume}}
#' 
#' @examples
#' fname <- system.file("extdata", "global_mask.nii", package="neuroim")
#' x <- loadVolume(fname)
#' print(dim(x))
#' space(x)
#' 
#' @export loadVolume
loadVolume  <- function(fileName, index=1) {
	src <- BrainVolumeSource(fileName, index)
	loadData(src)
}


#' @note dimensions of x and y must be equal
#' @export concat
#' @rdname concat-methods
setMethod(f="concat", signature=signature(x="DenseBrainVolume", y="DenseBrainVolume"),
		def=function(x,y,...) {
			.concat4D(x,y,...)			
		})


#' @export
#' @rdname fill-methods
setMethod(f="fill", signature=signature(x="BrainVolume", lookup="list"),
          def=function(x,lookup) {
            out <- lookup[unlist(x)]
            DenseBrainVolume(unlist(out), space(x))
          })

#' @export
#' @rdname fill-methods
setMethod(f="fill", signature=signature(x="BrainVolume", lookup="matrix"),
          def=function(x,lookup) {
            if (ncol(lookup) != 2) {
              stop("fill: lookup matrix must have two columns: column 1 is key, column 2 is value")
            } else if (nrow(lookup) < 1) {
              stop("fill: lookup matrix have at least one row")
            }
            
            out <- array(0, dim(x))
            
            for (i in 1:nrow(lookup)) {
              idx <- which(x == lookup[i,1])
              out[idx] <- as.vector(lookup[i,2])             
            }
            
            BrainVolume(out, space(x))
          })


#' @export splitFill
#' @rdname splitFill-methods
setMethod(f="splitFill", signature=signature(x="BrainVolume", fac="factor", FUN="function"),
		def=function(x,fac,FUN) {
			stopifnot(length(x) == length(fac))
			S <- split(1:length(x), fac, drop=TRUE)
			res <- sapply(S, function(ind) {
						X <- FUN(x[ind])
						if (length(X) == length(ind)) {
							X
						} else {
							rep(X, length(ind))
						}
					}, simplify=FALSE)
			
			ovol <- x
			ovol[1:length(x)] <- unsplit(res, fac)
			ovol
					
		})

 
#' @export 
#' @rdname slice-methods
setMethod(f="slice", signature=signature(x="BrainVolume", zlevel="numeric", along="numeric", orientation="character"),
          def=function(x, zlevel, along, orientation) {
            imslice <- switch(as.character(along),
                              "1"=x[zlevel,,],
                              "2"=x[,zlevel,],
                              "3"=x[,,zlevel])
            
            BrainSlice(imslice, dropDim(space(x), along))
            
            
          })


 
#' @export 
#' @rdname eachSlice-methods
setMethod(f="eachSlice", signature=signature(x="BrainVolume", FUN="function", withIndex="missing"),
		def=function(x, FUN) {
			lapply(1:(dim(x)[3]), function(z) FUN(x[,,z]))				
		})



#' @export 
#' @rdname eachSlice-methods
setMethod(f="eachSlice", signature=signature(x="BrainVolume", FUN="function", withIndex="logical"),
		def=function(x, FUN, withIndex) {
			lapply(1:(dim(x)[3]), function(z) {					
				slice <- x[,,z]
				if (withIndex) FUN(slice,z) else FUN(slice)
			})
		})


#' @export 
#' @rdname indexToCoord-methods
setMethod(f="indexToCoord", signature=signature(x="BrainVolume", idx="index"),
          def=function(x, idx) {
            callGeneric(space(x),idx)
          })


#' @export 
#' @rdname coordToIndex-methods
setMethod(f="coordToIndex", signature=signature(x="BrainVolume", coords="matrix"),
          def=function(x, coords) {
            assert_that(ncol(coords) == 3)
            callGeneric(space(x), coords)
          })


 
#' @export 
#' @rdname indexToGrid-methods
setMethod(f="indexToGrid", signature=signature(x="BrainVector", idx="index"),
		  def=function(x, idx) {
			  callGeneric(space(x), idx)
		  })

 
#' @export 
#' @rdname indexToGrid-methods
setMethod(f="indexToGrid", signature=signature(x="BrainVolume", idx="index"),
		  def=function(x, idx) {
			  callGeneric(space(x), idx)
		  })


#' @export 
#' @rdname gridToIndex-methods
setMethod(f="gridToIndex", signature=signature(x="BrainVolume", coords="matrix"),
          def=function(x, coords) {
            assert_that(ncol(coords) == 3)
            array.dim <- dim(x)
            .gridToIndex3D(dim(x), coords)
          })
  

#' @export 
#' @rdname gridToIndex-methods
  setMethod(f="gridToIndex", signature=signature(x="BrainVolume", coords="numeric"),
		  def=function(x, coords) {
		    assert_that(length(coords) == 3)
			  array.dim <- dim(x)
			  .gridToIndex3D(dim(x), matrix(coords, nrow=1, byrow=TRUE))
		  })

 
#' @export 
#' @rdname as.mask-methods
setMethod(f="as.mask", signature=signature(x="BrainVolume", indices="missing"),
          def=function(x) {
            LogicalBrainVolume(x > 0, space(x))
          })

#' @export 
#' @rdname as.mask-methods
setMethod(f="as.mask", signature=signature(x="BrainVolume", indices="numeric"),
          def=function(x, indices) {
            M <- array(0, dim(x))
            M[indices] <- 1
            LogicalBrainVolume(M, space(x))
          })

 
#' @export 
#' @rdname coordToGrid-methods
setMethod(f="coordToGrid", signature=signature(x="BrainVolume", coords="matrix"),
          def=function(x, coords) {
            callGeneric(space(x), coords)            
          })  


.pruneCoords <- function(coord.set,  vals,  mindist=10) {

	if (NROW(coord.set) == 1) {
		1
	}
	
	.prune <- function(keepIndices) {
		if (length(keepIndices) == 1) {
			keepIndices
		} else {
			ret <- yaImpute::ann(coord.set[keepIndices,], coord.set[keepIndices,], verbose=F,  k=2)$knn
			ind <- ret[, 2]
			ds <- sqrt(ret[, 4])
			v <- vals[keepIndices] 
			ovals <- v[ind]		
			pruneSet <- ifelse(ds < mindist & ovals > v,  TRUE, FALSE)
			
			if (any(pruneSet)) {
				Recall(keepIndices[!pruneSet])
			} else {
				keepIndices
			}
		}
		
		
	}
		
	.prune(1:NROW(coord.set))
	
	
	  
}

#' apply a kernel function to a \code{\linkS4class{BrainVolume}}
#' 
#' @rdname map-methods
#' @param mask restrict application of kernel to maksed area
#' @export
setMethod(f="map", signature=signature(x="BrainVolume", m="Kernel"),
          def=function(x, m, mask=NULL) {          
            ovol <- array(0, dim(x))
            hwidth <- sapply(m@width, function(d) ceiling(d/2 -1)) + 1
            xdim <- dim(x)[1]
            ydim <- dim(x)[2]
            zdim <- dim(x)[3]
            
            if (!is.null(mask)) {
              if (!all.equal(dim(mask), dim(ovol))) {
                stop(paste("mask must have same dimensions as input volume"))
              }
              # TODO check that mask is same shape as volume
              grid <- indexToGrid(mask, which(mask != 0))
            } else {
              grid <- as.matrix(expand.grid(i=hwidth[1]:(xdim - hwidth[1]), j=hwidth[2]:(ydim - hwidth[2]), k=hwidth[3]:(zdim - hwidth[3])))
            }
                      
            res <- apply(grid, 1, function(vox) {
              loc <- sweep(m@voxmat, 2, vox, "+")
              ivals <- x[loc]
              if (all(ivals == 0)) {
                0
              } else {
                sum(ivals * m@weights)
              }
            })
            ovol[grid] <- res                  
      
            BrainVolume(ovol, space(x))
          })

#' tesselate a LogicalBrainVolume into K spatial disjoint components
#' @param features use addiitonal feature set to tesselate volume
#' @param spatialWeight weight voxels according to distance
#' @importFrom stats kmeans
#' @importFrom stats sd
#' @rdname tesselate-methods
setMethod(f="tesselate", signature=signature(x="LogicalBrainVolume", K="numeric"), 
          def=function(x, K, features=NULL, spatialWeight=4) {
            voxgrid <- indexToGrid(x, which(x == TRUE))
            voxgrid <- sweep(voxgrid, 2, spacing(x), "*")
            
            if (!is.null(features)) {
              if (nrow(features) == length(x)) {
                features <- features[which(x == TRUE),]
              }
              
        
              features <- apply(features,2, scale)
              avg.sd <- sum(apply(voxgrid,2, sd))/ncol(features)
              sfeatures <- (features*avg.sd)/spatialWeight
              voxgrid <- cbind(voxgrid, sfeatures)
            }
    
            kgrid <- kmeans(voxgrid, centers=K, iter.max=100)
            ClusteredBrainVolume(x, kgrid$cluster, kgrid$centers)
          })


#' get number of clusters in a ClusteredBrainVolume
#' @rdname numClusters-methods
#' @export
setMethod(f="numClusters", signature=signature(x="ClusteredBrainVolume"), 
          def=function(x) {
            length(x@clusterMap)
          })


#' @rdname clusterCenters-methods
#' @import parallel
#' @export
setMethod(f="clusterCenters", signature=signature(x="ClusteredBrainVolume", features="matrix", FUN="missing"), 
          def=function(x, features) {
            cmap <- x@clusterMap
            res <- parallel::mclapply(sort(as.integer(names(cmap))), function(cnum) {
              idx <- cmap[[as.character(cnum)]]
              mat <- features[idx,]
              colMeans(mat)
            })
            
            
            mat <- t(do.call(cbind, res))
            row.names(mat) <- as.character(sort(as.integer(names(cmap))))
            mat
          })


#' merge partititons in a ClusteredBrainVolume
#' @rdname mergePartitions-methods
#' @importFrom stats kmeans
#' @export
setMethod(f="mergePartitions", signature=signature(x="ClusteredBrainVolume", K="numeric", features="matrix"), 
          def=function(x, K, features) {
            centers <- clusterCenters(x, features)
            kres <- kmeans(centers, centers=K)
            oclusters <- numeric(prod(dim(x)))
            for (cnum in 1:nrow(centers)) {
              idx <- x@clusterMap[[as.character(cnum)]]
              oclusters[idx] <- kres$cluster[cnum]
            }
            
            oclusters <- oclusters[x@mask == TRUE]
            ClusteredBrainVolume(x@mask, as.integer(oclusters))
            
          })
          

#' partition a \code{ClusteredBrainVolume} into K spatial disjoint components for every existing partition in the volume
#' @param method clustering method
#' @import parallel
#' @importFrom stats kmeans
#' @rdname partition-methods
#' @export
setMethod(f="partition", signature=signature(x="ClusteredBrainVolume", K="numeric", features="matrix"), 
          def=function(x, K, features, method="kmeans") {
            cmap <- x@clusterMap
            cnums <- sort(as.integer(names(cmap)))
            
            kres <- parallel::mclapply(cnums, function(cnum) {
              idx <- cmap[[as.character(cnum)]]
              fmat <- features[idx,]
              kmeans(fmat, centers=K)
            })
            
            oclusters <- numeric(prod(dim(x@mask)))
  
            for (id in 1:length(kres)) {
              idx <- cmap[[as.character(id)]]
              clus <- kres[[id]]$cluster
              new.clus <- (id-1)*(K) + clus
              oclusters[idx] <- new.clus
            }
            oclusters <- oclusters[x@mask == TRUE]
            perm.clusters <- sample(unique(oclusters[oclusters!= 0]))
            oclusters <- perm.clusters[oclusters]
            ClusteredBrainVolume(x@mask, as.integer(oclusters))
            
          })            
            
#' find connected components in BrainVolume
#' @export
#' @param threshold threshold defining lower intensity bound for image mask
#' @param clusterTable return clusterTable
#' @param localMaxima return table of local maxima
#' @param localMaximaDistance the distance used to define minum distance between local maxima 
#' @rdname connComp-methods
setMethod(f="connComp", signature=signature(x="BrainVolume"), 
	def=function(x, threshold=0, clusterTable=TRUE, localMaxima=TRUE, localMaximaDistance=15) {
		mask <- (x > threshold)
		stopifnot(any(mask))
	
		comps <- connComp3D(mask)
		
		grid <- as.data.frame(indexToGrid(mask, which(mask>0)))
		colnames(grid) <- c("x", "y", "z")
		locations <- split(grid, comps$index[comps$index>0])
		
		ret <- list(size=BrainVolume(comps$size, space(x)), index=BrainVolume(comps$index, space(x)), voxels=locations)
		
				
		
		if (clusterTable) {
			maxima <- do.call(rbind, lapply(locations, function(loc) {
				if (nrow(loc) == 1) {
					loc
				} else {
					vals <- x[as.matrix(loc)]
					loc[which.max(vals),]
				}			
			}))
			N <- comps$size[as.matrix(maxima)]
			Area <- N * prod(spacing(x))
			maxvals <- x[as.matrix(maxima)]
			ret$clusterTable <- data.frame(index=1:NROW(maxima), x=maxima[,1], y=maxima[,2], z=maxima[,3], N=N, Area=Area, value=maxvals)			
		}
		
		if (localMaxima) {	
			if (all(sapply(locations, NROW) == 1)) {
				
			}	
			coord.sets <- lapply(locations, function(loc) {
				sweep(as.matrix(loc), 2, spacing(x), "*")
			})
			
	
		    
			loc.max <- do.call(rbind, mapply(function(cset, i) {	
				idx <- .pruneCoords(as.matrix(cset), x[as.matrix(locations[[i]])], mindist=localMaximaDistance)
				maxvox <- as.matrix(locations[[i]])[idx,,drop=F]
				cbind(i, maxvox)
			}, coord.sets, 1:length(coord.sets), SIMPLIFY=FALSE))
			
			
			loc.max <- cbind(loc.max, x[loc.max[, 2:4, drop=F]])
			
			row.names(loc.max) <- 1:NROW(loc.max)
			colnames(loc.max) <- c("index", "x", "y", "z", "value")
			ret$localMaxima <- loc.max
		}
		
		ret
		
		
	})
	
	
	
    
### TODO when the source voulme is an AFNI BRIK the output file does not preserve orientation
### e.g. writeVolume(x, space(BRIK), ...)
### if BRIK was RAI, output NIFTI file reverts to LPI


#' @export
#' @rdname writeVolume-methods
setMethod(f="writeVolume",signature=signature(x="BrainVolume", fileName="character", format="missing", dataType="missing"),
		def=function(x, fileName) {
			write.nifti.volume(x, fileName)           
		})

 
#' @export
#' @rdname writeVolume-methods
setMethod(f="writeVolume",signature=signature(x="ClusteredBrainVolume", fileName="character", format="missing", dataType="missing"),
          def=function(x, fileName) {
            callGeneric(as(x, "DenseBrainVolume"), fileName)        
          })




#' @export
#' @rdname writeVolume-methods
setMethod(f="writeVolume",signature=signature(x="BrainVolume", fileName="character", format="character", dataType="missing"),
		def=function(x, fileName, format) {
			if (toupper(format) == "NIFTI" || toupper(format) == "NIFTI1" || toupper(format) == "NIFTI-1") {
				callGeneric(x, fileName)
			} else {
				stop(paste("sorry, cannot write format: ", format))
			}      
		})


#' @export writeVolume
#' @rdname writeVolume-methods
setMethod(f="writeVolume",signature=signature(x="BrainVolume", fileName="character", format="missing", dataType="character"),
		def=function(x, fileName, dataType) {
			write.nifti.volume(x, fileName, dataType)   
			
		})

#' as.logical
#' 
#' Convert BrainVolume to \code{linkS4class{LogicalBrainVolume}}
#' 
#' the image values will be converted to using R base function \code{as.logical} and wrapped in \code{LogicalBrainVolume}
#' 
#' @param x the object
#' @return an instance of \code{linkS4class{LogicalBrainVolume}}
#' @rdname as.logical-methods
#' @export 
setMethod(f="as.logical", signature=signature(x = "BrainVolume"), def=function(x) {
			vals <- as.logical(as.vector(x))	
			LogicalBrainVolume(vals, space(x))
})

#' @rdname as.sparse-methods
setMethod(f="as.sparse", signature=signature(x="DenseBrainVolume", mask="LogicalBrainVolume"),
          def=function(x, mask) {
            assert_that(all(dim(x) == dim(mask)))
            assert_that(all(spacing(x) == spacing(mask)))
            
            dat <- x[mask]
            bvec <- SparseBrainVector(dat, space(x))
            
})

#' @rdname as.sparse-methods
setMethod(f="as.sparse", signature=signature(x="DenseBrainVolume", mask="numeric"),
          def=function(x, mask) {
            m <- as.integer(mask)
            bvec <- SparseBrainVolume(x[m], space(x), indices=m)
            
          })



#' extractor
#' @export 
#' @param x the object
#' @param i first index 
#' @param j second index 
#' @param k third index 
#' @param ... additional args
#' @param drop dimension
setMethod(f="[", signature=signature(x = "SparseBrainVolume", i = "numeric", j = "numeric", drop="ANY"),
          def=function (x, i, j, k, ..., drop=TRUE) {  
            if (missing(k)) {
              k <- seq(1, dim(x)[3])
            }
            
            n <- length(i) * length(j) * length(k)
            mind <- cbind(rep(i, length.out=n), rep(j, each=length(i)), rep(k, each=length(i) * length(j)))
            
            ind <- gridToIndex(x, mind)
            x@data[ind]
            
        })


#' extractor
#' @export
#' @param x the object
#' @param i first index 
#' @param j second index 
#' @param k third index 
#' @param ... additional args
#' @param drop drop dimension
setMethod(f="[", signature=signature(x = "SparseBrainVolume", i = "numeric", j = "missing", drop="missing"),
          def=function (x, i, j, k, ..., drop) {  
            if (missing(k) && nargs() == 4) {
              x@data[i]
            } else {
              callGeneric(x, i=i,  seq(1,dim(x)[2]), k, drop)
            } 
         }
)

#' extractor
#' @export 
#' @param x the object
#' @param i first index 
#' @param j second index 
#' @param k third index 
#' @param ... additional args
#' @param drop dimension
setMethod(f="[", signature=signature(x = "SparseBrainVolume", i = "matrix", j="missing", drop="ANY"),
          def=function (x, i, j, k, ..., drop=TRUE) {  
            ind <- gridToIndex(x,i)
            x@data[ind]
          }
)

#' extractor
#' @export 
#' @param x the object
#' @param i first index 
#' @param j second index 
#' @param k third index 
#' @param ... additional args
#' @param drop dimension
setMethod(f="[", signature=signature(x = "SparseBrainVolume", i = "missing", j = "missing", drop="ANY"),
          def=function (x, i, j, k, ..., drop=TRUE) {  
            if (missing(k)) {
              x@data
            } else {
              callGeneric(x, seq(1, dim(x)[1]), seq(1, dim(x)[2]), k, ...)
            }
          }
)

#' extractor
#' @export 
#' @param x the object
#' @param i first index 
#' @param j second index 
#' @param k third index 
#' @param ... additional args
#' @param drop dimension
setMethod(f="[", signature=signature(x = "SparseBrainVolume", i = "missing", j = "numeric", drop="ANY"),
          def=function (x, i, j, k,  ..., drop=TRUE) {  
            callGeneric(x, seq(1, dim(x)[1]), j, k,...)
          }
)



            
            
