convData <- function(linesMat, maps){
	
	linesDes <- vapply(linesMat, function(lin){
		goodLin <- lin
		for(mp in maps){
			goodLin <- mp(goodLin)
		}
		return(goodLin)
		
	}, "e")
	
	return(linesDes)
}

#' Make the maps for converting slice notation
#'
#' @inheritParams makeDataMap
#' @return A function that takes in a string and converts all the given slice
#'   notation
#' @details Slice notation for matrices are tricky because they can easily be
#'   confused with the requirements for conversion are the bounds given by both
#'   left and right symbols or the Matlab class. The Matlab class allows for the
#'   conversion of structures but is really just a dictionary for the different 
#'   bounds.
#' @examples
#'   sliceMap <- makeSliceMap("{", "}", "list")
#'   sliceMap("junk <- importData{300}")
#'   # "junk <- importData[[300]]"
#'   
#'   sliceMap <- makeSliceMap(matClass = "structure", rClass = "list")
#'   sliceMap("junk <- students.AP.GPA")
#'   # junk <- students[['AP']][['GPA']]
#' @export
makeSliceMap <- function(leftSym, rightSym, rClass, matClass = ""){
	
	if(!isClassName(rClass)){
		stop(sprintf("'%s' is not a valid R class", rClass))
	}
	
	if(!nzchar(matClass)){
		if(missing(leftSym) || missing(rightSym)){
			stop("Please provide either the bounds of the data or the Matlab class")
		}
	} else {
		syms <- getMatLabClassBounds(matClass)
		leftSym <- syms[1]
		rightSym <- syms[2]
	}
	
	if(matClass == "matrix"){
		stop(paste(
			"Matrix slicing will not be converted because it has the same",
			"notation as function use. Converting these would be destructive."))
	}
	
	rBounds <- switch(rClass,
		vector = c("[", "]"),
		data.frame = ,
		list = c("[[", "]]")
	)
	
	return(function(lin){
		
		goodLin <- lin
		
		guts <- shExtractData(lin, leftSym, rightSym, type = "slice")
		while(nzchar(guts)){
			bef <- goodLin
			
			goodLin <- if(matClass == "structure"){
				rout <- sprintf("%s'%s'%s", rBounds[1], guts, rBounds[2])
				bas <- getBetween(goodLin, leftSym, rightSym, insertChar = rout)
				sub("[.]", "", bas)
			} else {
				rout <- paste0(rBounds[1], guts, rBounds[2])
				getBetween(goodLin, leftSym, rightSym, insertChar = rout, shInclude = TRUE)
			}
			guts <- getBetween(removeStrings(goodLin), leftSym, rightSym)
		}
		return(goodLin)
		
	})
}

shExtractData <- function(lin, leftSym, rightSym, type = c("inst", "slice")[1]){
	preLin <- lin
	preLin <- removeStrings(preLin)
	guts <- getBetween(preLin, leftSym, rightSym)
	if(type == "inst") guts <- putBackStrings(guts, lin)
	
	equ <- regexpr("[<]-|=", preLin)
	st <- regexpr(paste0("\\", leftSym), preLin)
	
	varNamePos <- regexpr(paste0("\\w+\\", leftSym), preLin)
	
	
	if(type == "inst"){
		if(equ > st || varNamePos > 0) guts <- "" 
	} else {
		if(varNamePos < 0) guts <- ""
	}
	
	return(guts)
}

getMatLabClassBounds <- function(matClass){
	
	matDict <- list(
		string = c("[", "]"),
		structure = c(".", "W|$"),
		cell = c("{", "}"),
		matrix = c("[", "]")
		)
	
	bounds <- matDict[[matClass]]
	
	if(is.null(bounds)){
		stop(paste(
			"The class '", matClass,
			"is not supported. supported classes are",
			paste(names(matDict), sep = " | ")))
	}
	return(bounds)
}

#' Make the maps for the data
#'
#' @param leftSym The left symbol that contains the Matlab data
#' @param rightSym the right symbol that contains the Matlab data
#' @param rClass The formal r class name that defines what the R data is
#'   outputted as
#' @param matClass The name of the Matlab class that should be converted
#' @return A function that takes in a Matlab lines and changes the data into R
#'   data lines
#' @details The requirements for conversion are the bounds given by both left
#'   and right symbols or the MatLab class. The Matlab class allows for the
#'   conversion of structures but is really just a dictionary for the different
#'   bounds.
#' @examples
#' 	 dataMap <- makeDataMap("[", "]", "matrix")
#' 	 dataMap("thing <- [23,2, 3.2; 7, 6, 8]")
#' 	 # "thing <- matrix(c(23, 2, 3.2, 7, 6, 8), nrow = 2, ncol = 3)"
#' 	 
#' 	 dataMap <- makeDataMap(rClass = "list", matClass = "cell")
#' 	 dataMap("otherThing <- {23,2, '3.2'; NaN, 6, 8}")
#' 	 # "otherThing <- list(list(23, 2, '3.2'), list(NaN, 6, 8))"
#' @export
makeDataMap <- function(leftSym, rightSym, rClass, matClass = ""){

	if(!isClassName(rClass)){
		stop(sprintf("'%s' is not a valid R class", rClass))
	}
	
	if(!nzchar(matClass)){
		if(missing(leftSym) || missing(rightSym)){
			stop("Please provide either the bounds of the data or the matlab class")
		}
	} else {
		syms <- getMatLabClassBounds(matClass)
		leftSym <- syms[1]
		rightSym <- syms[2]
	}




	return(function(lin){
		guts <- shExtractData(lin, leftSym, rightSym, type = "inst")
		stringFlag <- !(removeStrings(lin) == lin)
		if(!stringFlag){
			if(matClass == "string") guts <- ""
		} else {
			if(matClass == "matrix") guts <- ""
		}
		
		if(!nzchar(guts)){
			return(lin)
		} else {
			
			rout <- switch(rClass,
				vector = sprintf("c(%s)",
					paste(splitMatVec(guts, stringFlag), collapse = ", ")),
				data.frame = sprintf("as.data.frame(%s)", matrixify(guts)),
				list = listify(guts),
				matrix = matrixify(guts)
			)
			return(
				getBetween(lin,leftSym, rightSym, insertChar = rout, shInclude = TRUE))
		}


	})
}

matrixify <- function(lin){
	noNums <- gsub("\\d+(\\.\\d+)?|NA|NaN", "|", lin)
	
	numVec <- splitMatVec(lin)
	numVec <- as.numeric(numVec)

	refMat <- lapply(strsplit(noNums, "[|]"), getRowColFromData)

	rMat <- simplify2array(refMat)
	maxes <- vapply(1:dim(rMat)[3], function(x){
		c(max(rMat[,1,x]), max(rMat[,2,x]))
	}, rep(1,2))

	outMat <- sprintf(
		"matrix(c(%s), nrow = %d, ncol = %d)",
		paste(numVec, collapse = ", "),
		maxes[1,],
		maxes[2,])

	return(outMat)

}

getRowColFromData <- function(vin){
	rows <- c(1, grep("[;]", vin))
	cols <- c(1, grep("\\s|[,]", vin))
	
	rowInd <- findInterval(1:length(vin), rows, rightmost.closed = FALSE)
	tem <- diff(c(rows, length(vin)+1))
	colInd <- c(lapply(tem, function(x){1:x}), recursive = TRUE)
	if(length(rowInd) - length(colInd) != 0) stop("non equal in matrixfy")
	return(cbind(rowInd, colInd))
}


listify <- function(lin){
	bef <- lin
	lin <- removeStrings(lin)
	
	coorVec <- regmatches(lin, gregexpr("\\,|\\;",lin))
	coorVec <- lapply(coorVec, function(x) c("",x))
	
	refMat <- lapply(coorVec, getRowColFromData)
	
	ele <- splitMatVec(lin)
	ele <- putBackStrings(ele, bef)
	
	out <- vapply(refMat, function(RCMat){
		innerLists <- c()
		for(rind in unique(RCMat[, "rowInd"])){
			thisRowSet <- (RCMat[, "rowInd"] == rind)
			innerLists[rind] <- sprintf("list(%s)",
				paste(ele[thisRowSet], collapse = ", "))
		}
		
		return(paste(innerLists, collapse = ", "))
		
	}, "e")
	return(sprintf("list(%s)", out))
	
}

splitMatVec <- function(sin, hasStrings = FALSE){
	bef <- sin
	sin <- removeStrings(sin)

	allNums <- gsub(",", " ", gsub(";", " ", sin))
	thingVec <- trimWhite(c(strsplit(allNums, " "), recursive = TRUE))
	
	
	thingVec <- thingVec[nzchar(thingVec)]
	
	return(if(hasStrings){
		sprintf("paste0(%s)", paste(putBackStrings(thingVec, bef), collapse = ", "))
	} else {
		thingVec
	})

}
