setClass("STRSPLITCOL",
	representation = representation(
						strEqcCommand		=	"character",
						colSplit			=	"character",
						strSplit			=	"character",
						numSplitIdx			=	"numeric",
						colOut				=	"character"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						colSplit			=	"",
						strSplit			=	"",
						numSplitIdx			=	-1,
						colOut				=	"NewCol"
						)
	#contains = c("EcfReader")
)

setGeneric("setSTRSPLITCOL", function(object) standardGeneric("setSTRSPLITCOL"))
setMethod("setSTRSPLITCOL", signature = (object = "STRSPLITCOL"), function(object) {
	
	aEqcSlotNamesIn = c("colSplit","strSplit", "numSplitIdx", "colOut")
	#aEcfSlotNamesIn = c("arcdAddCol", "acolOut")

	objEqcReader <- EqcReader(object@strEqcCommand,aEqcSlotNamesIn)
	
	if(length(objEqcReader@lsEqcSlotsOut) > 0) {
		for(i in 1:length(objEqcReader@lsEqcSlotsOut)) {
			tmpSlot <- names(objEqcReader@lsEqcSlotsOut)[i]
			tmpSlotVal <- objEqcReader@lsEqcSlotsOut[[i]]
			
			if(all(!is.na(tmpSlotVal))) slot(object, tmpSlot) <- tmpSlotVal
		}
	}
	return(object)
})

#############################################################################################################################
validSTRSPLITCOL <- function(objSTRSPLITCOL){ 
	
	if(objSTRSPLITCOL@numSplitIdx == -1)
		stop(paste("EASY ERROR:STRSPLITCOL\n Please set --numSplitIdx to define which part of each output array should be used. \n !!!", sep=""))
	
	return(TRUE)
}
#############################################################################################################################
STRSPLITCOL.GWADATA.valid <- function(objSTRSPLITCOL,objGWA){ 
	
	idxSplitCol  = match(objSTRSPLITCOL@colSplit, objGWA@aHeader)
	
	if(is.na(idxSplitCol))
		stop(paste("EASY ERROR:STRSPLITCOL\n Column colSplit\n",objSTRSPLITCOL@colSplit," does not exist in file\n",objGWA@fileInShortName,"\n !!!", sep=""))
	
	if(objGWA@aClasses[idxSplitCol] != "character")
		stop(paste("EASY ERROR:STRSPLITCOL\n Column colSplit\n",objSTRSPLITCOL@colSplit," must be of type 'character' in file\n",objGWA@fileInShortName,"\n !!!", sep=""))
	
	return(TRUE)
}

#############################################################################################################################
STRSPLITCOL.run <- function(objSTRSPLITCOL, objGWA) {

	colSplit <- objSTRSPLITCOL@colSplit
	strSplit <- objSTRSPLITCOL@strSplit
	numSplitIdx <- objSTRSPLITCOL@numSplitIdx
	colOut <- objSTRSPLITCOL@colOut
	
	idxSplitCol  = match(objSTRSPLITCOL@colSplit, objGWA@aHeader)
	
	out <- as.character(unlist(lapply(strsplit(objGWA@tblGWA[,idxSplitCol],strSplit),function(x) x[numSplitIdx])))
	
	blnColExist = any(names(objGWA@tblGWA) == colOut)
	 
	objGWA <- GWADATA.cbind(objGWA, out, colOut, TRUE)
		
	return(objGWA)
}

STRSPLITCOL <- function(strEqcCommand){ 
	## Wrapper for class definition
	STRSPLITCOLout <- setSTRSPLITCOL(new("STRSPLITCOL", strEqcCommand = strEqcCommand))
	validSTRSPLITCOL(STRSPLITCOLout)
	#STRSPLITCOLout.valid <- validSTRSPLITCOL(STRSPLITCOLout)
	return(STRSPLITCOLout)
	#validECF(ECFout)
	#return(ECFout)
	
	## Identical:
	# ECFin <- new("ECF5", fileECF = fileECFIn) 
	# ECFout <- setECF5(ECFin)
	# return(ECFout)
}

# setValidity("STRSPLITCOL", function(object){
	# print("STRSPLITCOL-CHECK")
	
	
	
	# print(TRUE)
	# return(TRUE)
# })

