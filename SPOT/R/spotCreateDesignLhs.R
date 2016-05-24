####################################################################################
#' spotCreateDesignLhs
#' 
#' Uses the improvedLHS function from the lhs package to create a Latin Hypercube
#' Design. improvedLHS optimizes the euclidean distance between points in the sample.
#'
#' @param spotConfig list of spotConfiguration
#' @param noDesPoints number of design points, default is NaN
#' @param repeats number of repeats, default is NaN
#'
#' @return matrix \code{M} \cr
#' - \code{M} has \code{dimension} columns and \code{noDesPoints} rows
#' with entries corresponding to the region of interest.
#' @export
#' @seealso \code{\link{spotCreateDesignBasicDoe}}, \code{\link{spotCreateDesignFrF2}}, 
#' \code{\link{spotCreateDesignLhd}}, \code{\link{spotCreateDesignLhsOpt}}
####################################################################################
spotCreateDesignLhs <- function(spotConfig, noDesPoints = NaN, repeats=NaN){	
	spotWriteLines(spotConfig$io.verbosity,2,"  Entering spotCreateDesignLhs.R::spotCreateDesignLhs()");
	spotInstAndLoadPackages("lhs")	
	
	## use roi or aroi:
	if(spotConfig$spot.fileMode){
		if(file.exists(spotConfig$io.aroiFileName))
			roiConfig <- spotReadRoi(spotConfig$io.aroiFileName,spotConfig$io.columnSep,spotConfig$io.verbosity)
		else
			roiConfig <- spotReadRoi(spotConfig$io.roiFileName,spotConfig$io.columnSep,spotConfig$io.verbosity)
	}else{
		roiConfig <- spotConfig$alg.aroi
		if(is.null(roiConfig)) roiConfig <- spotConfig$alg.roi
	}
	
	
	pNames <- row.names(roiConfig);
	lowerBound <-  roiConfig[ ,"lower"];
	upperBound <-  roiConfig[ ,"upper"];
		
	A <- t(rbind(t(lowerBound), t(upperBound)))
			
  M<-as.matrix(lhs::improvedLHS(noDesPoints, length(pNames), dup=2))
	## M has entries in the range from 0 to 1, so we 	
	for (i in 1:nrow(M)){## transform these values into the regions of interest:
		for (j in 1: ncol(M)){
			M[i,j] <- A[j,1] + M[i,j] * (A[j,2] - A[j,1])
		}
	}
	M <- as.data.frame(M)
	colnames(M) <- pNames	
	spotWriteLines(spotConfig$io.verbosity,2,"  Leaving spotCreateDesignLhs.R::spotCreateDesignLhs")
	M
}
