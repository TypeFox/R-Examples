convEasySyntax  <- function(linesMat){
	lcode <- getInlineComments(linesMat)
	linesOut <- lcode$sansInLineCom


	linesOut <- convSymbols(linesOut)
	linesOut <- convSemiColon(linesOut)
	linesOut <- convIfElse(linesOut)
	for (loopType in c("for", "while")){
		linesOut <- convLoops(linesOut,loopType)
	}

	linesOut <- convEqualsArrow(linesOut)

	linesOut <- paste0(linesOut, lcode$inLineCom)
	return(linesOut)
}

convSemiColon <- function(linesMat){
	return(gsub(";$", "", linesMat))
}

convEqualsArrow <- function(linesMat){

	linesOut <- gsub(" = ", " <- ", linesMat)


	#filter all the other cases matlab uses =
	logStrs <- c("==", "<=", ">=", "!=")
	logSet <- asRightMatrix(vapply(logStrs, function(x){ grepl(x, linesMat) },
		rep(TRUE, length(linesMat))))
	logSet <- as.logical(rowSums(logSet))
	linesOut[!logSet] <- gsub("=", " <- ", linesOut[!logSet])
	return(linesOut)
}

convLoops <- function(linesMat, loopStr){
	loopSet <- grepl(loopStr, linesMat)
	loopLines <- paraAroundThings(linesMat[loopSet], loopStr)
	loopLines <- gsub(" = ", " in ", loopLines)
	loopLines <- gsub("=", " in ", loopLines)

	linesMat[loopSet] <- loopLines
	return(linesMat)

}

convIfElse  <- function(linesMat){

	elseSet <- grepl("else$", linesMat)
	ifelseSet <- grepl("elseif", linesMat)
	ifSet <- grepl("if", linesMat)

	linesMat[ifelseSet] <- gsub("elseif", "} else if", linesMat[ifelseSet])
	linesMat[elseSet] <- gsub("else", "} else {", linesMat[elseSet])

	linesMat[ifSet] <- paraAroundThings(linesMat[ifSet], "if")

	return(linesMat)
}

convSymbols <- function(linesMat){
	
	#Change "end"s in slicing 
	sliceEndSet <- grepl("\\:*end\\s*(\\}|\\))", linesMat)
	
	rightLoc <- gregexpr("\\}|\\)", linesMat[sliceEndSet])
	leftLoc <- gregexpr("\\(|\\{", linesMat[sliceEndSet])
	endLoc <- gregexpr("end", linesMat[sliceEndSet])
	
	linesMat[sliceEndSet] <- vapply(linesMat[sliceEndSet], function(lin){
		ends <- regexpr("end", lin)
		endInd <-  ends + attr(ends, "match.length")
		rightChar <- substr(lin, endInd, endInd)
		leftChar <- switch(rightChar, `}` = "{", `)` = "(")
		rghts <- gregexpr(paste0("\\", rightChar), lin)[[1]]
		lfts <- gregexpr(paste0("\\", leftChar), lin)[[1]]
		
		hier <- findInterval(ends, rghts)
		endCap <- rev(lfts)[hier + 1]
		sliceCut <- substr(lin, 1, endCap - 1)
		varName <- regmatches(sliceCut, regexpr("\\w+$", sliceCut))
		gsub("end", sprintf("length(%s)", varName), lin)
		
	}, "e", USE.NAMES = FALSE)
	
	
	
	linesMat <- gsub("end", "}", linesMat)
	return(gsub("~", "!", linesMat))
}

getInlineComments <- function(linesMat){
	lout <- list()
	lout$inLineCom <- rep('', length(linesMat))
	lout$sansInLineCom <- linesMat

	endInd <- regexpr("\\;", linesMat)
	percentInd <- regexpr("\\%", linesMat)
	inLineSet <- (endInd < percentInd)


	lout$inLineCom[inLineSet] <- substr(linesMat[inLineSet],
																			percentInd[inLineSet],
																			nchar(linesMat[inLineSet]))
	lout$inLineCom <- gsub("\\%", "#", lout$inLineCom)

	lout$sansInLineCom[inLineSet] <-
		trimWhite(substr(linesMat[inLineSet],
										 1,
										 percentInd[inLineSet] - 1),
							"end")

	return(lout)
}

paraAroundThings <- function(sin, type){
	startChar <- regexpr(type, sin) + 2
	lastChar <- nchar(sin)
	aftThing <- substr(sin, startChar, lastChar)
	firstLet <- regexpr("\\w", aftThing)
	firstPara <- regexpr("\\(", aftThing)

	badSet <- (firstLet < firstPara | firstPara < 0)
	sin[badSet] <- paste0(
		gsub(paste0(type, ' '), paste0(type, " ("), sin[badSet]),
		"){")
	sin[!badSet] <- gsub(paste0(type, "\\("), paste(type, '('), sin[!badSet])
	sin[!badSet] <- paste0(sin[!badSet], '{')
	return(sin)

}
