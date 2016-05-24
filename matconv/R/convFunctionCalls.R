convFunctionsCalls <- function(linesMat, maps){
	linesDes <- linesMat
	assignInd <- regexpr("[<]\\-", linesMat)
	leftParList <- gregexpr("\\(", linesMat)
	leftParInd <- vapply(leftParList, function(x){ rev(x)[1] }, 1)
	potSet <- (assignInd < leftParInd)
	
	noStringLin <- removeStrings(linesDes)
	mapNames <- paste0(names(maps), "\\(")
	
	linMapMat <- as.matrix(vapply(mapNames, function(pat){
		grepl(pat, noStringLin)
	},rep(TRUE, length(linesDes)), USE.NAMES = FALSE))
	
	linMapVec <- which(linMapMat, arr.ind = TRUE, useNames = FALSE)
	
	convSeq <- if(nrow(linMapVec) > 0){ 1:nrow(linMapVec) } else { NULL }
	for(convInd in convSeq){
		lin <- linesDes[linMapVec[convInd, 1]]
		map <- maps[[linMapVec[convInd, 2]]]
		pat <- mapNames[linMapVec[convInd, 2]]
		funcStart <- regexpr(pat, lin)
		restLin <- substr(lin, funcStart, nchar(lin)) 
		guts <- getBetween(restLin, "(", ")")
		matArgs <- strsplit(removeStrings(guts), ",")[[1]]
		matArgs <- putBackStrings(matArgs, guts)
		if(length(matArgs) ==  0) next
		
		if(!(is.null(map$flags$spaceSepMatArgs))){
			matArgs <- strsplit(lin, " ")[[1]]
		}
		matArgs <- trimWhite(matArgs)
		
		if(length(map$argMap) == 1){
			useMapInd <- 1
		} else {
			#Multiple dictionaries per matlab function
			#use fun switcher
			useMapInd <- map$flags$multSwitch(matArgs)
		}
		
		rargs <- map$argMap[[useMapInd]](matArgs)$rargs
		
		#Use other flags
		if(!is.null(map$flags$varOut)){
			sliceAdd <- ifelse(grepl("\\[", map$flags$varOut[1]), "", "$")
			reqVars <- strsplit(getBetween(lin, "[", "]"), " ")[[1]]
			addCalls <- paste(
				paste0(reqVars, " <- lout", sliceAdd, map$flags$varOut),
				collapse = "; ")
			out <- sprintf("lout <- %s); %s",
				rargs,
				addCalls)
			
		} else {
			out <- getBetween(restLin, '', ')',
				insertChar = rargs, whatIsEmpty = "first")
			out <- paste0(
				substr(lin, 1, funcStart - 1),
				out)
		}
		
		linesDes[linMapVec[convInd, 1]]	<- as.character(out)
	}


	#deal with space sep ones
	spaceArgSet <- vapply(maps, function(x){ !is.null(x$flags$spaceSepMatArgs) }, FALSE)
	potSpace <- strsplit(linesDes, "(?<=\\w)\\s", perl = TRUE )
	spaceLineSet <- vapply(potSpace, function(x){
		!is.na(match(x[1], names(maps[spaceArgSet])))
	}, TRUE)
	spaceArgs <- strsplit(linesMat[spaceLineSet], " ")
	if(any(spaceArgSet) && any(spaceLineSet)){
		linesDes[spaceLineSet] <- mapply(function(marg, mp){
			rout <- mp$argMap[[1]](marg[-1])$rargs

			out <- paste0(paste(rout, collapse = ", "), ")")
			return(out)
		}, spaceArgs, maps[spaceArgSet])
	}


	return(linesDes)
}

#' Turn dictionary lines into functions that map matlab to R function calls
#'
#' @param addDict An optional character vector with manufactored lines
#' @param pathDict The path to a text file with the dictionary lines written to it
#'
#' @return a list of functions to convert the arguments of a matlab function. It
#'   comes with the names of matlab functions.
#' 
#' @details The output of the individual maps consits of the actual map for the
#'   given matlab arguments as a vector and a list of flags included in the
#'   dictionary. The argMap itself is a list of potential functions that could
#'   be used if a some flags are detected in the dictionary line. A more
#'   expansive look at the different dictionaries that could be used can be seen
#'   in the base dictionary at "extdata/HiebelerDict.txt" or in the vignette
#'   "vignettes/functionCalls.rmd". It returns a list with the R version of the
#'   arguments with a left parentheisis.
#' @examples
#' 
#' funcMap <- makeFuncMaps("trace: sum, diag(%1)")
#' funcMap[['trace']]$argMap[[1]]("matThing")
#' #$rargs
#' # "sum(diag(matThing)"
#' 
#' funcMap <- makeFuncMaps("mod: , 1 %% 2")
#' funcMap[['mod']]$argMap[[1]](c(4, 2))
#' #$rargs
#' # "(4, %%, 2"
#' 
#' test1 <- "mat"
#' test2 <- c("mat", "2")
#' 
#' funcMap <- makeFuncMaps(c("size--if 1:dim, 1", "size--if 2: ,dim(%1)[%2]"))
#' rightConv <- funcMap$size$flags$multSwitch(test1)
#' funcMap$size$argMap[[rightConv]](test1)
#' #$rargs
#'  "dim(mat"
#'  
#' rightConv <- funcMap$size$flags$multSwitch(test2)
#' funcMap$size$argMap[[rightConv]](test2)
#' #$rargs
#'  "dim(mat)[2]"
#' @export
makeFuncMaps <- function(addDict = NULL, pathDict = ''){
	dictLines <- addDict
	if(nzchar(pathDict)){
		if(!file.exists(pathDict)){
			stop(paste(pathDict, "does not exist, please supply a dictionary file"))
		}

		dictFile <- readLines(pathDict)
		if (length(dictFile) == 0){
			stop(paste(pathDict, "is empty, please fill with matLab functions"))
		}
		dictLines <- c(dictLines, dictFile)
	}

	if(length(dictLines) == 0){
		stop(paste("No dictionaries supplied",
			"either feed in a character vector",
			"or a file with the dictionaries", sep = ", "))
	}



	lout <- parseFlags(dictLines)
	dictLines <- lout$strSansFlags

	keyVal <- strsplit(dictLines, ":")
	allFunNames <- vapply(keyVal, function(x){ x[1] }, "e")
	allDictArgs <- vapply(keyVal, function(x){ x[2] }, "e")
	finFunNames <- unique(allFunNames)

	maps <- lapply(1:length(finFunNames), function(x){
		list(argMap = list(), flags = list()) })
	names(maps) <- finFunNames

	argFuns <- lapply(allDictArgs, function(x){ parseArgs(x) })

	dupsMat <- (duplicated(allFunNames) | duplicated(allFunNames, fromLast = TRUE))

	anum <- 1
	while(anum <= length(argFuns)){
		nam <- allFunNames[anum]
		wantVec <- anum

		if(dupsMat[anum]){
			lastDup <- which(!dupsMat[anum:length(argFuns)])[1] - 2 + anum
			if(is.na(lastDup)){
				#All dups
				lastDup <- length(dupsMat)
			}
			wantVec <- anum:lastDup
			anum <- lastDup
		}

		maps[[nam]]$argMap <- argFuns[wantVec]
		anum <- anum + 1
	}

	for(nm in finFunNames){
		maps[[nm]]$flags <- lout$flags[[nm]]
	}

	return(maps)

}

`%isKey%` <- function(vals, ldict){
	return(is.element(names(ldict), vals))
}

parseArgs <- function(dictArg){
	sargs <- strsplit(dictArg, ',')
	sargs <- trimWhite(sargs[[1]])

	rname <- sargs[1]
	sargs <- sargs[-1]

	swiSet <- grepl("^[0-9]+$", sargs)
	literalNumSet <- grepl("^[0-9]+L$", sargs)
	strInsertSet <- grepl("\\%[0-9]", sargs)
	stringSet <- !literalNumSet & !swiSet & !strInsertSet

	return(function(matArg){
		rargs <- NULL
		rargs[swiSet] <- matArg[as.integer(sargs[swiSet])]
		rargs[literalNumSet] <- as.numeric(gsub("L", "", sargs[literalNumSet]))
		for(iar in which(strInsertSet)){
			arg <- sargs[iar]
			test <- TRUE
			while(test){
				ind <- as.numeric(getBetween(arg, '%', ''))
				arg <- sub("\\%[0-9]", matArg[ind], arg)
				test <- grepl("\\%[0-9]", arg)
			}
			rargs[iar] <- arg
		}

		rargs[stringSet] <- sargs[stringSet]

		return(list(
			rargs = paste0(rname, '(', paste(rargs, collapse = ", "))
			))
	})
}

parseFlags <- function(dictLines){


	flagStr <- lapply(1:length(dictLines), function(x){ list() })
	strSansFlags <- dictLines

	#separate flags
	stFlag <- gregexpr("\\-\\-", dictLines)
	stDiv <- regexpr("[:]", dictLines)
	flagSet <- vapply(stFlag, function(x){ x[1] > 0 }, TRUE)

	for(ind in which(flagSet)){
		left <- stFlag[[ind]] + 2
		right <- ifelse(stFlag[[ind]] > stDiv[[ind]],
			nchar(dictLines[ind]),
			stDiv[[ind]] - 1
		)
		for(flagInd in 1:length(left)){
			flagStr[[ind]] <- substr(dictLines[ind], left[flagInd], right[flagInd])
			strSansFlags[ind] <- paste0(
				substr(strSansFlags[ind], 1, left[flagInd] - 3),
				substr(strSansFlags[ind], right[flagInd] + 1, nchar(strSansFlags[ind]))
			)
		}
	}

	#make flags and funcSwitchers
	matName <- vapply(strsplit(strSansFlags, ":"), function(x){ x[1] },"e")

	uniMatName <- unique(matName)
	dupsSet <- vapply(uniMatName, function(x){
		sum(grepl(x, matName)) > 1
	}, TRUE)

	matNameswFlags <- unique(matName[flagSet])
	uniFlagNums <- match(matNameswFlags, uniMatName)

	flags <- lapply(1:length(uniMatName), function(x){ list() })
	names(flags) <- uniMatName

	for(unind in uniFlagNums){

		wantVec <- which(!is.na(match(matName, uniMatName[unind])))
		if(dupsSet[unind]){

			flags[[unind]] <- lapply(wantVec, function(x){
				makeFlag(flagStr[[x]], makeSwitch = FALSE)
			})
			flags[[unind]]$multSwitch <- makeFunSwitcher(flagStr[wantVec])

		} else {
			flags[[unind]] <- makeFlag(flagStr[[wantVec]])
		}
	}

	return(mget(c("strSansFlags", "flags")))
}

makeFlag <- function(vin, makeSwitch = TRUE){
	flag <- list()
	possFlags <- c("if", "out", "space-sep", "not-req")




	for(si in vin){
		para <- strsplit(si, " ")[[1]]
		flagName <- para[1]

		if(flagName == "if"){
			if(makeSwitch) flag$multSwitch <- makeFunSwitcher(list(si))
		} else if (flagName == "out"){
			flag$varOut <- para[-1]
		} else if (flagName == "space-sep"){
			flag$spaceSepMatArgs <- TRUE
		} else {
			stop(paste("The flag:", si, "is indecipherable", sep = "\n"))
		}
	}

	return(flag)
}

makeFunSwitcher <- function(lFlags){

	finallyInd <- NULL
	lengthOutVec <- lengthVec <- rep(NA, length(lFlags))
	matMap <- lapply(1:length(lFlags), function(x){
		list(arg = NULL, val = NULL)
	})

	for(dictNum in 1:length(lFlags)){
		para <- strsplit(lFlags[[dictNum]][1], ' ')[[1]][-1]
		if(length(para) == 1){
			if(para[1] == "finally"){
				finallyInd <- dictNum
			} else {
				lengthVec[dictNum] <- as.integer(para[1])
			}
		} else {
			if(para[1] == "length(out)"){
				lengthOutVec[dictNum] <- as.integer(gsub("L", "", para[3]))
			} else {
				matMap[[dictNum]]$arg <- para[1]
				matMap[[dictNum]]$val <- gsub("L", "", para[3])
			}
		}
	}

	return(function(matArgs, numOut = 1){
		useInd <- NULL
		if(numOut > 1){
			useInd <- which(lengthOutVec == numOut)
		}

		useInd <- c(useInd, which(lengthVec == length(matArgs)))

		test <- vapply(matMap, function(mp){
			check <- matArgs[as.integer(mp$arg)] == mp$val
			if(length(check) == 0) check <- FALSE
			return(check)
		}, TRUE)
		useInd <- c(useInd, which(test))

		if(length(useInd) == 0){
			if(!is.null(finallyInd)){
				useInd <- finallyInd
			} else {
				stop(paste("Do not have rule that supports:" , matArgs))
			}
		}

		return(useInd[1])
	})
}
