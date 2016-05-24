setClass("RENAMECOL",
	representation = representation(
						strEqcCommand	=	"character",
						colInRename		=	"character",
						colOutRename	=	"character"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						colInRename			=	"",
						colOutRename		=	""
						)
	#contains = c("EcfReader")
)

setGeneric("setRENAMECOL", function(object) standardGeneric("setRENAMECOL"))
setMethod("setRENAMECOL", signature = (object = "RENAMECOL"), function(object) {
	
	aEqcSlotNamesIn = c("colInRename", "colOutRename")
	
	objEqcReader <- EqcReader(object@strEqcCommand,aEqcSlotNamesIn)
	
	if(length(objEqcReader@lsEqcSlotsOut) > 0) {
		for(i in 1:length(objEqcReader@lsEqcSlotsOut)) {
			tmpSlot <- names(objEqcReader@lsEqcSlotsOut)[i]
			tmpSlotVal <- objEqcReader@lsEqcSlotsOut[[i]]
			
			#if(all(!is.na(tmpSlotVal))) slot(object, tmpSlot) <- tmpSlotVal
			if(any(!is.na(tmpSlotVal))) slot(object, tmpSlot) <- tmpSlotVal
		}
	}
	return(object)
})

#############################################################################################################################
validRENAMECOL <- function(objRNC) {

	if(objRNC@colInRename == "")
		stop(paste(" EASY ERROR:RENAMECOL\n No colInRename defined. Please set colInRename.", sep=""))
	
	if(objRNC@colOutRename == "") 
		stop(paste(" EASY ERROR:RENAMECOL\n No colOutRename defined. Please set colOutRename.", sep=""))
		
	return(TRUE)
}

#############################################################################################################################
RENAMECOL.run <- function(objRNC, objGWA, objREPORT) {
	
	iMatchIn = match(objRNC@colInRename,objGWA@aHeader)
	iMatchOut = match(objRNC@colOutRename,objGWA@aHeader)
	
	isAnyRenamed = FALSE
	
	if(!is.na(iMatchIn) & is.na(iMatchOut)) {
		### Column must be renamed, colOut is not in table!
		names(objGWA@tblGWA)[iMatchIn] <-  objGWA@aHeader[iMatchIn] <- objRNC@colOutRename
		isAnyRenamed = TRUE
	} else if(!is.na(iMatchIn) & !is.na(iMatchOut)) {
		### Column colInRename AND colOutRename are already in table
		
		names(objGWA@tblGWA)[iMatchOut] <-  objGWA@aHeader[iMatchOut] <- paste(objGWA@aHeader[iMatchOut],".old",sep="")
		names(objGWA@tblGWA)[iMatchIn] <-  objGWA@aHeader[iMatchIn] <- objRNC@colOutRename
		isAnyRenamed = TRUE
		
	} else {
		## No action
	}
	
	
	if(!("numColRenamed" %in% names(objREPORT@tblReport))) objREPORT <- REPORT.addval(objREPORT, "numColRenamed", 0)
	else if(REPORT.getval(objREPORT, "numColRenamed") == "NA") objREPORT <- REPORT.addval(objREPORT, "numColRenamed", 0)
	
	if(isAnyRenamed) {
		## Update Report numColRenamed = numColRenamed + 1

		strNumColRenamed = REPORT.getval(objREPORT, "numColRenamed")
		numColRenamedNew = as.numeric(strNumColRenamed) + 1
		objREPORT <- REPORT.setval(objREPORT, "numColRenamed", numColRenamedNew)

	}
	
	return(list(objGWA,objREPORT))
}

RENAMECOL <- function(strEqcCommand){ 
	## Wrapper for class definition
	RENAMECOLout <- setRENAMECOL(new("RENAMECOL", strEqcCommand = strEqcCommand))
	validRENAMECOL(RENAMECOLout)
	#RENAMECOLout.valid <- validRENAMECOL(RENAMECOLout)
	return(RENAMECOLout)
	#validECF(ECFout)
	#return(ECFout)
	
	## Identical:
	# ECFin <- new("ECF5", fileECF = fileECFIn) 
	# ECFout <- setECF5(ECFin)
	# return(ECFout)
}

# setValidity("RENAMECOL", function(object){
	# print("RENAMECOL-CHECK")
	
	
	
	# print(TRUE)
	# return(TRUE)
# })

