
####################################################################################
#' spotCreateDesignBasicDoe
#'
#' Create a full factorial design based on the number of dimensions and the number of design points.
#' Uses the function gen.factorial from the AlgDesign package. Two levels are used for each variable.
#'
#' The dimension is driven from the number of rows of the .roi - file 
#' (each row in the roi file defines a variable, 
#' In case of a missing number of design points a value is calculated from the dimension 
#'
#' @param spotConfig list of spotConfiguration, where only the table $alg.roi (.roi-file) is used,
#' @param noDesPoints optional parameter gives the number of design points. If noDesPoints
#' is larger than the size of full factorial design, an error is reported.
#' @param repeats is obsolete for this type of design, it has no influence.
#'
#' @return Matrix \code{M} \cr
#' - \code{M} has \code{dimension} columns and \code{noDesPoints} rows
#' @export
#' @seealso \code{\link{spotCreateDesignLhd}}, \code{\link{spotCreateDesignFrF2}}, 
#' \code{\link{spotCreateDesignLhs}}, \code{\link{spotCreateDesignLhsOpt}}
####################################################################################
spotCreateDesignBasicDoe <- function(spotConfig, noDesPoints = 100, repeats=1){
	spotWriteLines(spotConfig$io.verbosity,2,"  Entering spotCreateDesignDoe.R::spotCreateDesignDoe()");
	#spotInstAndLoadPackages("AlgDesign")
	
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
	
	
	pNames <- row.names(spotConfig$alg.roi);
	lowerBound <-  roiConfig[ ,"lower"];
	upperBound <-  roiConfig[ ,"upper"];
	
	A <- t(rbind(t(lowerBound), t(upperBound)))
	
	M <-  gen.factorial(2, length(pNames))
	if (nrow(M)<noDesPoints){# JZ test must be made here, because now lines are added to M
		stop("Number of initial design points (init.design.size) too large for this design")		
	}
	for (i in 1:nrow(M)){
		for (j in 1: ncol(M)){
			M[i,j] <- A[j,1] + (M[i,j] +1)/2 * (A[j,2] - A[j,1])
		}
	}
	
	colnames(M) <- pNames
	x <- 1:nrow(M)
	
	M <- M[ sample(x,replace=FALSE)[1:noDesPoints], ]	
	spotWriteLines(spotConfig$io.verbosity,2,"  Leaving spotCreateDesignDoe.R::spotCreateDesignDoe")
	M		
}