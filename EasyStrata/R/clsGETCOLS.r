setClass("GETCOLS",
	representation = representation(
						strEqcCommand		=	"character",
						acolOut				=	"character"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						acolOut				=	""
						)
	#contains = c("EcfReader")
)

setGeneric("setGETCOLS", function(object) standardGeneric("setGETCOLS"))
setMethod("setGETCOLS", signature = (object = "GETCOLS"), function(object) {
	
	aEqcSlotNamesIn = c("acolOut")
	
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
validGETCOLS <- function(objGC) {
	
	if(objGC@acolOut[1] == "")
		stop(paste(" EASY ERROR:GETCOLS\n No acolOut defined.\n Please set at least one column for acolOut or remove GETCOLS function.", sep=""))
	
	return(TRUE)
}

#############################################################################################################################
GETCOLS.run <- function(objGC, objGWA) {
		
	objGWA <- GWADATA.getcols(objGWA, objGC@acolOut, blnSuppressError=TRUE)
	
	return(objGWA)
}

GETCOLS <- function(strEqcCommand){ 
	## Wrapper for class definition
	GETCOLSout <- setGETCOLS(new("GETCOLS", strEqcCommand = strEqcCommand))
	validGETCOLS(GETCOLSout)
	#GETCOLSout.valid <- validGETCOLS(GETCOLSout)
	return(GETCOLSout)
	#validECF(ECFout)
	#return(ECFout)
	
	## Identical:
	# ECFin <- new("ECF5", fileECF = fileECFIn) 
	# ECFout <- setECF5(ECFin)
	# return(ECFout)
}

# setValidity("GETCOLS", function(object){
	# print("GETCOLS-CHECK")
	
	
	
	# print(TRUE)
	# return(TRUE)
# })

