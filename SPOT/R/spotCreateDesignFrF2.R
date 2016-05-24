####################################################################################
#' spotCreateDesignFrF2
#' 
#' Uses the function FrF2 from the FrF2-package to generate a (fractional) factorial design.
#' The factorial design is augmented with the augment.ccd function from DoE.wrapper package.
#' The dimension is determined from the number of rows of the .roi - file 
#' (each row in the roi file defines a variable). 
#'
#' @param spotConfig list of spot settings
#' @param noDesPoints  is obsolete for this type of design, it has no influence. The number of points is fixed.
#' @param repeats is obsolete for this type of design, it has no influence.
#'
#' @return matrix \code{M} \cr
#' - \code{M} has \code{dimension} columns and \code{noDesPoints} rows
#' @export
#' @seealso \code{\link{spotCreateDesignBasicDoe}}, \code{\link{spotCreateDesignLhd}}, 
#' \code{\link{spotCreateDesignLhs}}, \code{\link{spotCreateDesignLhsOpt}}
####################################################################################
spotCreateDesignFrF2 <- function(spotConfig, noDesPoints = NaN, repeats=NaN){	
	spotWriteLines(spotConfig$io.verbosity,2,"  Entering spotCreateDesignFrF2.R::spotCreateDesignFrF2()");
	spotInstAndLoadPackages(c('FrF2',  'DoE.wrapper'))
	#require(DoE.wrapper)#allready required in spotInstAndLoadPackages function
	
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
	
	# design: range from -1 to +1
	M<-FrF2::FrF2(nfactors=length(pNames), resolution=3)
	
	M <- DoE.wrapper::ccd.augment(M, ncenter = 1, columns="all",	alpha = sqrt(2)/2, randomize=TRUE)
	## remove block information
	M <- M[,-1]
	##delete replicates
	M <- unique(M)
	###############
	
	M<-as.matrix(M)	
		
	for (i in 1:nrow(M)){
		for (j in 1: ncol(M)){
			M[i,j] <- A[j,1] + (M[i,j] +1)/2 * (A[j,2] - A[j,1])
		}
	}		
	colnames(M) <- pNames	
	rownames(M) <- NULL
	spotWriteLines(spotConfig$io.verbosity,2,"  Leaving spotCreateDesignFrF2.R::spotCreateDesignFrF2")
	M <- as.data.frame(M)
	M		
}