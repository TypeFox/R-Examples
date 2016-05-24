getBetween <- function(sin, left, right,
	insertChar = NULL,
	whatIsEmpty = c("oneChar", "last", "first")[1],
	shInclude = FALSE){

	if(!nzchar(left)){
		rightPos <- regexpr(paste0("\\", right), sin)
		leftPos <- dealEmpty(rightPos, type = whatIsEmpty, fun = `-`, lin = sin)
	} else if(!nzchar(right)){
		leftPos <- regexpr(paste0("\\", left), sin)
		rightPos <- dealEmpty(leftPos, type = whatIsEmpty, fun = `+`, lin = sin)
	} else {
		leftPos <- regexpr(paste0("\\", left), sin)
		restSin <- substr(sin, leftPos + 1, nchar(sin))
		
		potRightPos <- gregexpr(paste0("\\", right), restSin)
		extLeftPos <- gregexpr(paste0("\\", left), restSin)
		if(length(potRightPos) > 0){
			rightPos <- mapply(function(potR, extL){
				ind <- 1
				while(ind <= length(potR) && any(extL < potR[ind] & extL > 0)){
					extL[ind] <- 99999
					potR[ind] <- 99999
					ind <- ind + 1
				}
				if(ind > length(potR)) return(-4)
				else return(potR[ind])
			}, potRightPos, extLeftPos)
			
			goodSet <- rightPos > 0
			rightPos[goodSet] <- rightPos[goodSet] + leftPos[goodSet]
			attr(rightPos, "match.length") <- attr(potRightPos[[1]], "match.length")[1]
			attr(rightPos, "useBytes") <- TRUE
		} else {
			rightPos <- leftPos
		}
		
		
		if(shInclude){
			leftPos <- leftPos - 1
			attr(leftPos, "match.length") <- 1

			rightPos <- rightPos + attr(rightPos, "match.length")
			attr(rightPos, "match.length") <- 1
		}
	}
	
	ind <- length(leftPos)
	while(ind > 0){
		if(leftPos[ind] < 0) rightPos[ind] <- leftPos[ind]
		if(rightPos[ind] < 0) leftPos[ind] <- rightPos[ind]
		ind <- (ind - 1)
	}


	if(is.null(insertChar)){
		cap <- substr(sin,
			leftPos + attr(leftPos, "match.length"),
			rightPos - 1
		)
		return(trimWhite(cap))
	} else {
		newStr <- paste0(
			substr(sin, 1, leftPos + attr(leftPos, "match.length") - 1 ),
			ifelse(leftPos >= 0, insertChar, ""),
			substr(sin, rightPos, nchar(sin))
		)
		return(newStr)
	}
}

dealEmpty <- function(pos, type, fun = NULL, lin = ""){
	out <- switch(type,
		oneChar = defaultOneChar(pos, fun),
		first = 0,
		last = nchar(lin))
	attr(out, "match.length") <- 1
	return(out)
}

trimWhite <- function(sin, where = "both"){
	return(switch(where,
                beg = gsub("^\\s+", "", sin),
                end = gsub("\\s+$", "", sin),
                both = gsub("^\\s+|\\s+$", "", sin)
	))
}

asRightMatrix <- function(vin){
	if(!is.matrix(vin)){
		t(as.matrix(vin))
	} else {
		vin
	}
}
defaultOneChar <- function(oppsMatch, func){
	defMatch <- func(oppsMatch, 2)
	attr(defMatch, "match.length") <- 1
	return(defMatch)
}

#' @importFrom methods new
isClassName <- function(sin){
	out <- tryCatch(new(sin),
		error = function(cond){
			!grepl("is not a defined class", cond)
		})
	return(!is.logical(out) || out)
}

removeStrings <- function(linesMat){
	
	linesDes <- vapply(linesMat, function(lin){
		bef <- noStringLin <- lin
		ind <- 1
		check <- TRUE
		rightQuote <- getMainQuote(lin)
		
		while(check){
			ins <- sprintf("#%s#", ind)
			noStringLin <- getBetween(noStringLin, rightQuote, rightQuote, 
				insertChar = ins, shInclude = TRUE)
			check <- !all(bef == noStringLin)
			bef <- noStringLin
			ind <- ind + 1
		}
		return(noStringLin)
	}, "sdf", USE.NAMES = FALSE)
	return(linesDes)
}

getMainQuote <- function(lin){
	rightQuote <- rep("'", length(lin))
	doubleSet <- (regexpr("'", lin) > regexpr('"', lin) & regexpr('"', lin) > 0)
	rightQuote[doubleSet] <- '"'
	return(rightQuote)
}

putBackStrings <- function(argVec, lin){
	bef <- lin
	if(length(argVec) == 0) return(lin)
	
	rightQuote <- getMainQuote(lin)
	
	ins <- regmatches(argVec, gregexpr("[#]\\d+[#]", argVec))
	needRep <- vapply(ins, function(x){ length(x) > 0 }, TRUE)
	for(ind in which(needRep)){
		for(pat in ins[[ind]]){
			realStr <- getBetween(lin, rightQuote, rightQuote, shInclude = TRUE)
			argVec[ind] <- gsub(pat, realStr, argVec[ind])
			
			lin <- getBetween(lin, rightQuote, rightQuote,
				insertChar = "", shInclude = TRUE)
		}
	}
	return(argVec)
}
