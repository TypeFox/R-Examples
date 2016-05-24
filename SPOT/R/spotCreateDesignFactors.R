####################################################################################
#' spotCreateDesignFactors
#' 
#' Design function for designs consisting of factors only. Factorial designs for al factors.
#' Assumse that factors are encoded in the ROI, with integer notation.
#'
#' @param spotConfig list of spot settings
#' @param noDesPoints  is obsolete for this type of design, it has no influence. The number of points is fixed by the number of factors and their levels.
#' @param repeats is obsolete for this type of design, it has no influence.
#'
#' @return matrix \code{M} \cr
#' - \code{M} has \code{dimension} columns and \code{noDesPoints} rows
#' @export
#' @seealso \code{\link{spotCreateDesignBasicDoe}}, \code{\link{spotCreateDesignLhd}}, 
#' \code{\link{spotCreateDesignLhs}}, \code{\link{spotCreateDesignLhsOpt}}
####################################################################################
spotCreateDesignFactors <- function(spotConfig, noDesPoints = NaN, repeats=NaN){	
	spotWriteLines(spotConfig$io.verbosity,2,"  Entering spotCreateDesignFrF2");
	#spotInstAndLoadPackages("AlgDesign") #for gen.factorial

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
	loLvl <-  roiConfig[ ,"lower"];
	upLvl <-  roiConfig[ ,"upper"];
	
	nfactors <- length(pNames)
	nlevels <- upLvl-loLvl + 1	
	#create design
	M<- gen.factorial(nlevels,nfactors,center=FALSE,varNames=pNames)	
	## fix level
	for(i in 1:nfactors){
		M[,i]<-M[,i]+loLvl[i] -1
	}	
	#
	rownames(M) <- NULL
	spotWriteLines(spotConfig$io.verbosity,2,"  Leaving spotCreateDesignFactors")
	M <- as.data.frame(M)
	M		
}