#' @include AllClass.R
roxygen()
#' @include AllGeneric.R
NULL
#' @include BrainVector.R
NULL
#' @include BrainVolume.R
NULL

#' @importFrom assertthat assert_that
checkDim <- function(e1,e2) {
  assert_that(all(dim(e1) == dim(e2)))
  assert_that(all(spacing(e1) == spacing(e2)))
 
}

setMethod(f="Arith", signature=signature(e1="SparseBrainVolume", e2="SparseBrainVolume"),
          def=function(e1, e2) {
            checkDim(e1,e2)          
            res <- callGeneric(e1@data,e2@data)   
            new("SparseBrainVolume", data=res, source=e1@source, space=space(e1))
             
          })


setMethod(f="Arith", signature=signature(e1="ROIVolume", e2="ROIVolume"),
          def=function(e1, e2) {
            checkDim(e1,e2)
            
            idx1 <- gridToIndex(e1@space, e1@coords)
            idx2 <- gridToIndex(e2@space, e2@coords)
            
            indices <- sort(union(idx1, idx2))   
            v1 <- numeric(length(indices))
            v2 <- numeric(length(indices))
            v1[indices %in% idx1] <- e1@data
            v2[indices %in% idx2] <- e2@data
            res <- callGeneric(v1,v2)   
          
            new("ROIVolume", space=space(e1), data=res, coords = indexToGrid(space(e1), indices))
            
          })





setMethod(f="Arith", signature=signature(e1="BrainVolume", e2="BrainVolume"),
          def=function(e1, e2) {
            checkDim(e1,e2)  
            
            ret <- callGeneric(e1@.Data,e2@.Data)
            bv <- DenseBrainVolume(ret, space(e1))
     
          })


setMethod(f="Arith", signature=signature(e1="BrainVector", e2="BrainVector"),
          def=function(e1, e2) {
            
      checkDim(e1,e2)  
            
			if (inherits(e1, "DenseBrainVector") && inherits(e2, "DenseBrainVector")) {
            	ret <- callGeneric(e1@.Data,e2@.Data)
            	DenseBrainVector(ret, space(e1))
			} else {
				D4 <- dim(e1)[4]		  
				vols <- list()
				for (i in 1:D4) {
					vols[[i]] <- callGeneric(takeVolume(e1,i), takeVolume(e2,i))
				}
				
				mat <- do.call(cbind, vols)
				dspace <- addDim(space(vols[[1]]), length(vols))	
				DenseBrainVector(mat, dspace)
				
			}
     
})
 

 setMethod(f="Arith", signature=signature(e1="BrainVector", e2="BrainVolume"),
		  def=function(e1, e2) {
			  if (!all(dim(e1)[1:3] == dim(e2))) {
				  stop("cannot perform operation on argument with different dimensions")
			  }
			  
			  D4 <- dim(e1)[4]	
			  vols <- list()
			  for (i in 1:D4) {
				  vols[[i]] <- callGeneric(takeVolume(e1,i), e2)
			  }
			  
			  mat <- do.call(cbind, vols)
			  dspace <- addDim(space(vols[[1]]), length(vols))	
			  DenseBrainVector(mat, dspace)
			  
		  
		  })
  

setMethod(f="Arith", signature=signature(e1="BrainVolume", e2="BrainVector"),
		def=function(e1, e2) {
			callGeneric(e2,e1)
		})


setMethod(f="Summary", signature=signature(x="SparseBrainVector"),
		def=function(x) {
			callGeneric(x@data)
		})


setMethod(f="Summary", signature=signature(x="SparseBrainVolume", na.rm="ANY"),
    def=function(x, ..., na.rm) {
      callGeneric(x@data)
    })
  

#setMethod("sum", signature()
          
