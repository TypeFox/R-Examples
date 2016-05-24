#' @include AllClass.R
{}
#' @include AllGeneric.R
{}

#TiledBrainVector <- function(fname, mask, ntiles=5, capacity=.5) {
#  stopifnot(capacity <= 1 && capacity > 0)
#  
#  indices <- which(mask > 0)
  
#  indexList <- split(indices, cut(1:length(indices),ntiles))
#  names(indexList) <- seq(1,length(indexList))
#  new("TiledBrainVector", filename=fname, indexList=indexList, mask=mask, capacity=capacity) 

#}

#setMethod("initialize", "TiledBrainVector", function(.Object, filename, mask, indexList, capacity) {
#  .Object@filename <- filename
#  .Object@cache <- vector("list", length(indexList))
#  .Object@indexList <- indexList
#  .Object@mask <- mask
#  .Object@capacity <- capacity
#  .Object
#})
  
#' SparseBrainVectorSource
#' 
#' constructs a SparseBrainVectorSource object
#' 
#' @param metaInfo an object of class \code{\linkS4class{BrainMetaInfo}}
#' @param indices a vector of 1D indices
#' @param mask a 3D \code{array} of type \code{logical} 
#' @export
#' @rdname SparseBrainVectorSource-class 	  
SparseBrainVectorSource <- function(metaInfo, indices, mask) {
  
	
	stopifnot(length(dim(metaInfo)) >= 3)
	stopifnot(all(indices >= 1 & indices <= dim(metaInfo)[4]))
	
	D <- dim(metaInfo)[1:3]
	
  
	if (is.vector(mask) && length(mask) < prod(D)) {
    ### this is a vector of indices
		m <- array(FALSE, D)
		m[mask] <- TRUE
		mask <- m
	} else if (identical(dim(mask), as.integer(D))) {
		mask <- as.array(mask)
	} else if (is.vector(mask) && length(mask) == prod(D)) {
		mask <- array(mask, D)
	} else {
		stop("illegal mask argument with dim: ", paste(dim(mask), collapse=", "))
	}
	
  if (!inherits(mask, "LogicalBrainVolume")) {
    mspace <- BrainSpace(dim(mask),  metaInfo@spacing, metaInfo@origin, metaInfo@spatialAxes)
    mask <- LogicalBrainVolume(mask, mspace)  	
  }
	
	stopifnot(all(dim(mask) == D))
	
	new("SparseBrainVectorSource", metaInfo=metaInfo, indices=indices, mask=mask)				
}


#' SparseBrainVector
#' 
#' constructs a SparseBrainVector object
#' 
#' @param data an array which can be a \code{matrix} or 4-D \code{array}
#' @param space a BrainSpace instance
#' @param mask a 3D \code{array} of type \code{logical} 
#' @param source the data source -- an instance of class \code{\linkS4class{BrainSource}}
#' @param label associated sub-image labels
#' @export 
#' @examples 
#' 
#' bspace <- BrainSpace(c(10,10,10,100), c(1,1,1))
#' mask <- array(rnorm(10*10*10) > .5, c(10,10,10))
#' mat <- matrix(rnorm(sum(mask)), 100, sum(mask))
#' svec <- SparseBrainVector(mat, bspace,mask)
#' length(indices(svec)) == sum(mask)
#' @rdname SparseBrainVector-class 	
SparseBrainVector <- function(data, space, mask, source=NULL, label="") {
	stopifnot(inherits(space, "BrainSpace"))
	
	
	if (!inherits(mask, "LogicalBrainVolume")) {
		mspace <- BrainSpace(dim(space)[1:3], spacing(space), origin(space), axes(space), trans(space))
		mask <- LogicalBrainVolume(as.logical(mask), mspace)
	}
	
	stopifnot(inherits(mask, "LogicalBrainVolume"))
	
	
	D4 <- if (is.matrix(data)) {
		Nind <- sum(mask == TRUE)
		if (nrow(data) == Nind) {
			data <- t(data)
			nrow(data)	
		} else if (ncol(data) == Nind) {
			nrow(data)
		} else {
			stop(paste("matrix with dim:", dim(data), " does not match mask cardinality: ", Nind))
		}
	} else if (length(dim(data)) == 4) {		
		mat <- apply(data, 4, function(vals) vals)
		data <- t(mat[mask==TRUE,])
		dim(data)[4]
	}
	
	if (ndim(space) == 3) {
		space <- addDim(space, nrow(data))
	}
				
  stopifnot(ndim(space) == 4)
	
	if (is.null(source)) {
		meta <- BrainMetaInfo(dim(space), spacing(space), origin(space), "FLOAT", label)
		source <- new("BrainSource", metaInfo=meta)	
	}
	 
	new("SparseBrainVector", source=source, space=space, mask=mask, data=data, map=IndexLookupVolume(space(mask), as.integer(which(mask))))
  
}




#' @export
#' @rdname loadData-methods 
setMethod(f="loadData", signature=c("SparseBrainVectorSource"), 
		def=function(x) {		
			
			### considerable code duplication with BrainVectorSource#loadData
			meta <- x@metaInfo
			#stopifnot(length(meta@Dim) == 4)
			
			meta <- x@metaInfo
			nels <- prod(meta@Dim[1:3]) 		
			
			ind <- x@indices
			M <- x@mask > 0
			
			reader <- dataReader(meta, offset=0)		
			dat4D <- readElements(reader, prod(meta@Dim[1:4]))
			datlist <- lapply(1:length(ind), function(i) {
				offset <- (nels * (ind[i]-1))
				dat4D[(offset+1):(offset + nels)][M]
			})
	
			close(reader)
				
			#datlist <- list()
			#for (i in 1:length(ind)) {
			#	offset <- (nels * (ind[i]-1)) * meta@bytesPerElement
			#	reader <- dataReader(meta, offset)		
			#	datlist[[i]] <- readElements(reader, nels)[M]
			#	close(reader)				
			#}
			
			arr <- do.call(rbind, datlist)		
			bspace <- BrainSpace(c(meta@Dim[1:3], length(ind)), meta@spacing, meta@origin, meta@spatialAxes)
			SparseBrainVector(arr, bspace, x@mask)
			
		})

#' @export
#' @rdname indices-methods
setMethod(f="indices", signature=signature(x="SparseBrainVector"),
          def=function(x) {
            indices(x@map)
          })
          
#setMethod(f="indices", signature=signature(x="TiledBrainVector"),
#          def=function(x) {
#           ret <- unlist(x@indexList)
#            names(ret) <- NULL
#            ret
#          })
  

#setMethod(f="coords", signature=signature(x="TiledBrainVector"),
#          def=function(x,i) {
#            if (missing(i)) {
#              return(indexToGrid(space(mask), indices(x)))
#            }
#
#            indexToGrid(space(mask), i)
#          })
            
            
#' @export 
#' @rdname coords-methods             
setMethod(f="coords", signature=signature(x="SparseBrainVector"),
          def=function(x,i) {           
            if (missing(i)) {
              return(coords(x@map, indices(x@map)))
            }
            coords(x@map, i)            
          })            

 
#' @export
#' @rdname eachVolume-methods
setMethod("eachVolume", signature=signature(x="SparseBrainVector", FUN="function", withIndex="logical", mask="missing"),
		def=function(x, FUN, withIndex=FALSE, mask, ...) {
			lapply(1:nrow(x@data), function(i) {
				if (withIndex) {
					FUN(takeVolume(x, i), i,...)
				} else {
					FUN(takeVolume(x, i), ...)
				}
			})
		})




 
#' @export
#' @rdname eachVolume-methods
setMethod("eachVolume", signature=signature(x="SparseBrainVector", FUN="function", withIndex="missing", mask="missing"),
		def=function(x, FUN, withIndex, ...) {
			lapply(1:nrow(x@data), function(i) FUN(takeVolume(x, i), ...))					
		})


#' @export
#' @rdname eachVolume-methods
setMethod("eachVolume", signature=signature(x="SparseBrainVector", FUN="function", withIndex="missing", mask="LogicalBrainVolume"),
          def=function(x, FUN, withIndex, mask, ...) {
            mask.idx <- which(mask > 0)
            lapply(1:nrow(x@data), function(i) FUN(takeVolume(x, i)[mask.idx], ...))					
          })

  

#' @export  
#' @details when \code{x} is a \code{SparseBrainVector} \code{eachSeries} only iterates over nonzero series.
#' @rdname eachSeries-methods 
setMethod(f="eachSeries", signature=signature(x="SparseBrainVector", FUN="function", withIndex="logical"),
          def=function(x, FUN, withIndex=FALSE, ...) {
            ## eachSeries only iterates over nonzero entries ...
            
            ret <- list()
            if (withIndex) {
              idx <- indices(x)
              for (i in 1:NCOL(x@data)) {
                ret[[i]] <- FUN(x@data[,i], idx[i])
              }
            } else {
              for (i in 1:NCOL(x@data)) {
                ret[[i]] <- FUN(x@data[,i])
              }
            }

            ret
          })


#' @export  
#' @describeIn seriesIter get a seriesIter for a \code{\linkS4class{SparseBrainVector}} instance
setMethod(f="seriesIter", signature=signature(x="SparseBrainVector"), 
	def=function(x) {
		len <- NCOL(x@data)
		i <- 0
		nextEl <- function() {
			i <<- i+1
			if (i <= len) {
				x@data[,i]
			} else {
				stop("StopIteration") 
			}		
		}

		hasNx <- function() {
			i < len
		}

		obj <- list(nextElem = nextEl, hasNext=hasNx) 
		class(obj) <- c("seriesIter", "abstractiter", "iter") 
		obj


	})


            


#setMethod(f="series", signature=signature(x="TiledBrainVector", i="numeric"),
#	def=function(x,i, j, k) {
#		stop()         
#      })
              

#setMethod(f="series", signature=signature(x="TiledBrainVector", i="matrix"),
#         def=function(x,i) {
#          idx <- gridToIndex(x@mask, i)
#           callGeneric(x,idx)
#         })
  
 
#' @export
#' @rdname series-methods 
setMethod(f="series", signature=signature(x="SparseBrainVector", i="matrix"),
         def=function(x,i) {
           idx <- gridToIndex(x@mask, i)
           callGeneric(x,idx)
         })
 

 #' @export
 #' @rdname series-methods 
 #' @param j index for 2nd dimension
 #' @param k index for 3rd dimension
 setMethod("series", signature(x="SparseBrainVector", i="numeric"),
		 def=function(x,i, j, k) {	
			 if (missing(j) && missing(k)) { 
				 idx <- lookup(x, as.integer(i))
				 idx.nz <- idx[idx!=0]
				 if (length(idx.nz) == 0) {
					 matrix(0, dim(x)[4], length(i))
				 } else {
					 mat <- matrix(0, dim(x)[4], length(i))
					 mat[, idx != 0] <- x@data[,idx.nz]     
					 mat
				 }
			 } else {
				 vdim <- dim(x)
				 slicedim <- vdim[1] * vdim[2]
				 idx <- slicedim*(k-1) + (j-1)*vdim[1] + i
				 callGeneric(x, idx)
			 }
				 
			
		 })
 

#' @export           
#' @rdname concat-methods 
setMethod(f="concat", signature=signature(x="SparseBrainVector", y="SparseBrainVector"),
          def=function(x,y,...) {
            if (!all(indices(x) == indices(y))) {
              stop("cannot concatenate arguments with different index maps")
            }

            ndat <- rbind(x@data, y@data)
            d1 <- dim(x)
            d2 <- dim(y)
            
            ndim <- c(d1[1:3], d1[4] + d2[4])
            nspace <- BrainSpace(ndim, spacing(x@space),  origin(x@space), axes(x@space), trans(x@space))
  
            
            ret <- SparseBrainVector(ndat, nspace, mask=x@mask)
            rest <- list(...)
            
            if (length(rest) >= 1) {
              ret <- Reduce("concat", c(ret, rest))
            }

            return(ret)
            
          })
          

#' @export          
#' @rdname lookup-methods          
setMethod(f="lookup", signature=signature(x="SparseBrainVector", i="numeric"),
         def=function(x,i) {
            lookup(x@map, i)
          })


#' extractor
#' @export 
#' @param x the object
#' @param i first index 
#' @param j second index 
#' @param k third index 
#' @param m the fourth index
#' @param ... additional args
#' @param drop dimension
setMethod(f="[", signature=signature(x = "SparseBrainVector", i = "numeric", j = "missing"),
		  def=function (x, i, j, k, m, ..., drop=TRUE) {  
			  callGeneric(x, i, 1:(dim(x)[2]))
		  }
  )


#' extractor
#' @export 
#' @param x the object
#' @param i first index 
#' @param j second index 
#' @param k third index 
#' @param m the fourth index
#' @param ... additional args
#' @param drop dimension
setMethod(f="[", signature=signature(x = "SparseBrainVector", i = "missing", j = "missing"),
		  def=function (x, i, j, k, m, ..., drop=TRUE) {  
			  callGeneric(x, 1:(dim(x)[1]), 1:(dim(x)[2]))
		  }
  )
  

#' extractor
#' @export 
#' @param x the object
#' @param i first index 
#' @param j second index 
#' @param k third index 
#' @param m the fourth index
#' @param ... additional args
#' @param drop dimension
setMethod(f="[", signature=signature(x = "SparseBrainVector", i = "missing", j = "numeric"),
		  def=function (x, i, j, k, m, ..., drop=TRUE) {  
			  callGeneric(x, i:(dim(x)[1]), j)
		  }
  )


#' extractor
#' @export 
#' @param x the object
#' @param i first index 
#' @param j second index 
#' @param k third index 
#' @param m the fourth index
#' @param ... additional args
#' @param drop dimension
setMethod(f="[", signature=signature(x = "SparseBrainVector", i = "numeric", j = "numeric"),
          def = function (x, i, j, k, m, ..., drop = TRUE) {
            if (missing(k))
              k = 1:(dim(x)[3])
             
            vmat <- as.matrix(expand.grid(i,j,k,m))
            ind <- .gridToIndex3D(dim(x)[1:3], vmat[,1:3,drop = FALSE])
            
            
            mapped <- cbind(lookup(x, ind), m)
            
            
            
            vals <- unlist(apply(mapped, 1, function(i) {
              if (i[1] == 0) {
                0
              } else {
                x@data[i[2], i[1]]
                #x@data[i[1],i[2]]
              }
            }))
            
            dim(vals) <- c(length(i),length(j),length(k),length(m))
            
            if (drop) {
              drop(vals)
            } else {
              vals
            }
			           
})

#' @export
#' @rdname subVector-methods
setMethod(f="subVector", signature=signature(x="SparseBrainVector", i="numeric"),
          def=function(x, i) {
            idx <- which(x@mask > 0)      
            bspace <- dropDim(space(x))
            
            res <- lapply(i, function(i) x@data[i,])
            res <- do.call("cbind", res)			
            SparseBrainVector(res, bspace, x@mask)
          })


 
#' @export
#' @rdname takeVolume-methods
setMethod(f="takeVolume", signature=signature(x="SparseBrainVector", i="numeric"),
          def=function(x, i, merge=FALSE) {
            idx <- which(x@mask > 0)      
            bspace <- dropDim(space(x))
             
            res <- lapply(i, function(i) x@data[i,])
       
            if (length(res) > 1 && merge) {
              res <- do.call("cbind", res)			
              SparseBrainVector(res, bspace, x@mask)
            } else {
              if (length(res) == 1) {
                BrainVolume(res[[1]], bspace, indices=idx)
              } else {
                lapply(res, function(x) BrainVolume(x, bspace, indices=idx))
              }
            }
          })

#' @export
setAs(from="SparseBrainVector", to="matrix",
		  function(from) {
		    ## TODO this should return a dense matrix
		    ind <- indices(from)
		    out <- matrix(dim(from)[4], length(ind))
		    out[, ind] <- from@data
		    out
			  #from@data			  
		  })

#' as.matrix
#' 
#' convert SparseBrainVector to matrix
#' @rdname as.matrix-methods
#' @export 
setMethod(f="as.matrix", signature=signature(x = "SparseBrainVector"), def=function(x) {
			  as(x, "matrix")						
		  })

#' as.list
#' 
#' convert SparseBrainVector to list of \code{\linkS4class{DenseBrainVolume}}
#' @rdname as.list-methods
#' @export
setMethod(f="as.list", signature=signature(x = "SparseBrainVector"), def=function(x) {
			D4 <- dim(x)[4]
			lapply(1:D4, function(i) takeVolume(x,i))
			
})

#' show a \code{SparseBrainVector}
#' @param object the object
#' @export
setMethod("show",
          signature=signature(object="SparseBrainVector"),
          def=function(object) {
            cat("an instance of class",  class(object), "\n\n")
            cat("   dimensions: ", dim(object), "\n")
            cat("   voxel spacing: ", spacing(object), "\n")
            cat("   cardinality: ", length(object@map@indices))
            cat("\n\n")
            
          })


#		setAs(from="BrainVector", to="matrix",
#		      function(from) {
#		        data <- from@.Data
#		        dm <- dim(data)
#		        d123 <- prod(dm[1:3])
#		        d4 <- dm[4]
#
#		        dim(data) <- c(d123,d4)
#		        data
#
#		      })
            

          
          
