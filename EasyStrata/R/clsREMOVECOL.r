setClass("REMOVECOL",
	representation = representation(
						strEqcCommand		=	"character",
						colRemove			=	"character"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						colRemove			=	""
						)
	#contains = c("EcfReader")
)

setGeneric("setREMOVECOL", function(object) standardGeneric("setREMOVECOL"))
setMethod("setREMOVECOL", signature = (object = "REMOVECOL"), function(object) {
	
	aEqcSlotNamesIn = c("colRemove")
	
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
validREMOVECOL <- function(objRC) {
	
	if(objRC@colRemove[1] == "")
		stop(paste(" EASY ERROR:REMOVECOL\n No colRemove defined.\n Please set colRemove or remove REMOVECOL function.", sep=""))
	
	return(TRUE)
}

#############################################################################################################################
REMOVECOL.run <- function(objRC, objGWA) {
	
	objGWA <- GWADATA.removecols(objGWA, objRC@colRemove)
	
	return(objGWA)
}

REMOVECOL <- function(strEqcCommand){ 
	## Wrapper for class definition
	REMOVECOLout <- setREMOVECOL(new("REMOVECOL", strEqcCommand = strEqcCommand))
	validREMOVECOL(REMOVECOLout)
	return(REMOVECOLout)
}
