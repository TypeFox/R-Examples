setClass("GETNUM",
	representation = representation(
						strEqcCommand				=	"character",
						rcdGetNum					=	"character",
						strGetNumName				=	"character"
						),
	prototype = prototype(
						strEqcCommand				=	"",
						rcdGetNum					=	"",
						strGetNumName				=	""
						)
	#contains = c("EcfReader")
)

setGeneric("setGETNUM", function(object) standardGeneric("setGETNUM"))
setMethod("setGETNUM", signature = (object = "GETNUM"), function(object) {
	
	aEqcSlotNamesIn = c("rcdGetNum", "strGetNumName")

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
validGETNUM <- function(objGETNUM) {
	
	if(objGETNUM@strGetNumName == "") 
		stop(paste(" EASY ERROR:GETNUM\n No strGetNumName defined for rcdGetNum \n ",objGETNUM@rcdGetNum," \n Please set strGetNumName.", sep=""))
	
	if(objGETNUM@rcdGetNum == "") 
		stop(paste(" EASY ERROR:GETNUM\n No rcdGetNum defined.\n Please set rcdGetNum or remove GETNUM function.", sep=""))
	
	return(TRUE)
}

#############################################################################################################################
GETNUM.run <- function(objGN, objGWA, objREPORT) {

	rcdGetNum 		<- objGN@rcdGetNum
	strGetNumName 	<- objGN@strGetNumName
	
	objRCD 	<- RCD(rcdGetNum)
	out 	<- RCD.eval(objRCD, objGWA)
	
	numOut = length(which(out))
	
	objREPORT <- REPORT.addval(objREPORT,strGetNumName,numOut)
	
	return(objREPORT)
}

GETNUM <- function(strEqcCommand){ 
	## Wrapper for class definition
	GETNUMout <- setGETNUM(new("GETNUM", strEqcCommand = strEqcCommand))
	validGETNUM(GETNUMout)
	#GETNUMout.valid <- validGETNUM(GETNUMout)
	return(GETNUMout)
}

