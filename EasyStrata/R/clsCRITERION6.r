setClass("CRITERION",
	representation = representation(
						strEqcCommand			=	"character",
						rcdCrit					=	"character",
						strCritName				=	"character",
						blnWriteEmpty			=	"logical"
						),
	prototype = prototype(
						strEqcCommand			=	"",
						rcdCrit					=	"",
						strCritName				=	"",
						blnWriteEmpty			=	FALSE
						)
	#contains = c("EcfReader")
)

setGeneric("setCRITERION", function(object) standardGeneric("setCRITERION"))
setMethod("setCRITERION", signature = (object = "CRITERION"), function(object) {
	
	aEqcSlotNamesIn = c("rcdCrit", "strCritName","blnWriteEmpty")

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
validCRITERION <- function(objCRITERION) {
	
	if(objCRITERION@strCritName == "") 
		stop(paste(" EASY ERROR:CRITERION\n No strCritName defined for rcdCrit \n ",objCRITERION@rcdCrit," \n Please set strCritName.", sep=""))
	
	if(objCRITERION@rcdCrit == "") 
		stop(paste(" EASY ERROR:CRITERION\n No rcdCrit defined.\n Please set rcdCrit or remove CRITERION function.", sep=""))
	
	return(TRUE)
}

#############################################################################################################################
CRITERION.run <- function(objCRIT, objGWA, objREPORT) {
	
		rcdCrit 		<- objCRIT@rcdCrit
		strCritName 	<- objCRIT@strCritName
		
		objRCD 	<- RCD(rcdCrit)
		out 	<- RCD.eval(objRCD, objGWA)
		
		numOut = length(which(out))
		
		objGWA.iCrit <- GWADATA.getrows(objGWA, which(out))
		
		#if(!blnSupressOutput) GWADATA.write(objGWA.iCrit, strSuffix = ".criterion")
		
		objREPORT <- REPORT.addval(objREPORT,strCritName,numOut)
	
	return(list(objGWA.iCrit,objREPORT))
}

CRITERION <- function(strEqcCommand){ 
	## Wrapper for class definition
	CRITERIONout <- setCRITERION(new("CRITERION", strEqcCommand = strEqcCommand))
	validCRITERION(CRITERIONout)
	#CRITERIONout.valid <- validCRITERION(CRITERIONout)
	return(CRITERIONout)
}

