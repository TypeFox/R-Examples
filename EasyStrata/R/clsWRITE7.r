setClass("WRITE",
	representation = representation(
						strEqcCommand	=	"character",
						strMode			=	"character",
						strPrefix		=	"character",
						strSuffix		=	"character",
						strSep			=	"character",
						strMissing		=	"character",
						strTabixParam	=	"character"
						),
	prototype = prototype(
						strEqcCommand	=	"",
						strMode			=	"txt",
						strPrefix		=	"",
						strSuffix		=	"",
						strSep			=	"TAB",
						strMissing		=	"NA",
						strTabixParam	=	""
						)
	#contains = c("EcfReader")
)

setGeneric("setWRITE", function(object) standardGeneric("setWRITE"))
setMethod("setWRITE", signature = (object = "WRITE"), function(object) {
	
	aEqcSlotNamesIn = c( "strMode", "strPrefix","strSuffix", "strSep", "strMissing","strTabixParam")
	
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
validWRITE <- function(objWRITE) {
	
	if(all(objWRITE@strSep != c("TAB", "SPACE", "COMMA")))
		stop("EASY ERROR:\nWRITE\nWrong File separator defined\nPlease use TAB, WHITESPACE, SPACE or COMMA!")
	
	## Reset separator
	if(objWRITE@strSep == "TAB") 			objWRITE@strSep <- "\t"
	if(objWRITE@strSep == "SPACE") 		objWRITE@strSep <- " "
	if(objWRITE@strSep == "COMMA") 		objWRITE@strSep <- ","
	
	if(all(objWRITE@strMode != c("txt", "gz", "bgz")))
		stop("EASY ERROR:\nWRITE\nWrong Writing Mode strMode defined\nPlease use txt, gz or bgz!")
	
	return(objWRITE)
}

#############################################################################################################################
WRITE.run <- function(objWRITE, objGWA, objREPORT) {
	
	GWADATA.write(objGWA, 
					strMode = objWRITE@strMode,
					strPrefix = objWRITE@strPrefix,
					strSuffix = objWRITE@strSuffix,
					strSep = objWRITE@strSep,
					strMissing = objWRITE@strMissing,
					strTabixParam = objWRITE@strTabixParam
					) 
					
	objREPORT <- REPORT.setval(objREPORT,"numSNPsOut",nrow(objGWA@tblGWA))
	
	return(objREPORT)
}

WRITE <- function(strEqcCommand){ 
	## Wrapper for class definition
	WRITEout <- setWRITE(new("WRITE", strEqcCommand = strEqcCommand))
	#validWRITE(WRITEout)
	WRITEout.valid <- validWRITE(WRITEout)
	return(WRITEout.valid)
	#validECF(ECFout)
	#return(ECFout)
	
	## Identical:
	# ECFin <- new("ECF5", fileECF = fileECFIn) 
	# ECFout <- setECF5(ECFin)
	# return(ECFout)
}

# setValidity("WRITE", function(object){
	# print("WRITE-CHECK")
	
	
	
	# print(TRUE)
	# return(TRUE)
# })

