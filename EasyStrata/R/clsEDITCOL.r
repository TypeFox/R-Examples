setClass("EDITCOL",
	representation = representation(
						strEqcCommand		=	"character",
						rcdEditCol			=	"character",
						colEdit				=	"character"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						rcdEditCol			=	"",
						colEdit				=	""
						)
	#contains = c("EcfReader")
)

setGeneric("setEDITCOL", function(object) standardGeneric("setEDITCOL"))
setMethod("setEDITCOL", signature = (object = "EDITCOL"), function(object) {
	
	aEqcSlotNamesIn = c("rcdEditCol", "colEdit")
	#aEcfSlotNamesIn = c("arcdEditCol", "acolEdit")

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
validEDITCOL <- function(objEDITCOL) {

	return(TRUE)
}

EDITCOL.GWADATA.valid <- function(objEDITCOL, objGWA){ 

	if(!(objEDITCOL@colEdit %in% objGWA@aHeader))
		stop(paste("EASY ERROR:EDITCOL\n Column colEdit\n",objEDITCOL@colEdit," does not exist in file\n",objGWA@fileInShortName,"\n !!!", sep=""))
	
}

#############################################################################################################################
EDITCOL.run <- function(objEDITCOL, objGWA) {

	rcdEditCol 		<- objEDITCOL@rcdEditCol
	colEdit 	<- objEDITCOL@colEdit
	
	objRCD 	<- RCD(rcdEditCol)
	out 	<- RCD.eval(objRCD, objGWA)
	
	objGWA <- GWADATA.cbind(objGWA, out, colEdit, blnOverwrite=TRUE)
		
	return(objGWA)
}

EDITCOL <- function(strEqcCommand){ 
	## Wrapper for class definition
	EDITCOLout <- setEDITCOL(new("EDITCOL", strEqcCommand = strEqcCommand))
	validEDITCOL(EDITCOLout)
	return(EDITCOLout)	
}
