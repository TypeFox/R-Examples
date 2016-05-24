setClass("CLEAN",
	representation = representation(
						strEqcCommand				=	"character",
						rcdClean					=	"character",
						strCleanName				=	"character",
						blnWriteCleaned				=	"logical"
						),
	prototype = prototype(
						strEqcCommand				=	"",
						rcdClean					=	"",
						strCleanName				=	"",
						blnWriteCleaned				=	FALSE
						)
	#contains = c("EcfReader")
)

setGeneric("setCLEAN", function(object) standardGeneric("setCLEAN"))
setMethod("setCLEAN", signature = (object = "CLEAN"), function(object) {
	
	aEqcSlotNamesIn = c("rcdClean", "strCleanName","blnWriteCleaned")

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
validCLEAN <- function(objCLEAN) {
	
	if(objCLEAN@strCleanName == "") 
		stop(paste(" EASY ERROR:CLEAN\n No strCleanName defined for rcdClean \n ",objCLEAN@rcdClean," \n Please set strCleanName.", sep=""))
	
	if(objCLEAN@rcdClean == "") 
		stop(paste(" EASY ERROR:CLEAN\n No rcdClean defined.\n Please set rcdClean or remove CLEAN function.", sep=""))
	
	return(TRUE)
}

#############################################################################################################################
CLEAN.run <- function(objCLEAN, objGWA, objREPORT, blnSuppressCleaning) {
	
	
	#for(iFilter in 1:length(objCLEAN@arcdClean)) {
		
		#rcdClean		<- objCLEAN@arcdClean[iFilter]
		#strCleanName 	<- objCLEAN@astrCleanName[iFilter]
		
		rcdClean <- objCLEAN@rcdClean
		strCleanName <- objCLEAN@strCleanName
		
		objRCD 	<- RCD(rcdClean)
		out 	<- RCD.eval(objRCD, objGWA)
		
		numOutDrop = length(which(out))
		
		# objGWA.iFilter <- GWADATA.getrows(objGWA, which(out))
		# GWADATA.write(objGWA.Duplicates, "filter", "2")
		
		
		#if(any(out)&!is.na(any(out))) objGWA <- GWADATA.getrows(objGWA, -which(out))
		#objGWA <- GWADATA.getrows(objGWA, -which(out))
		
		if(objCLEAN@blnWriteCleaned) objGWA.Cleaned <- GWADATA.getrows(objGWA, which(out))
		else objGWA.Cleaned <- GWADATA.copy(objGWA)
		
		if(!blnSuppressCleaning) objGWA <- GWADATA.removerows(objGWA, -which(out))
		
		objREPORT <- REPORT.addval(objREPORT,strCleanName,numOutDrop)
	#}
	
	
	lsOut <- list(objGWA, objREPORT, objGWA.Cleaned)
	
	return(lsOut)
}

CLEAN <- function(strEqcCommand){ 
	## Wrapper for class definition
	CLEANout <- setCLEAN(new("CLEAN", strEqcCommand = strEqcCommand))
	validCLEAN(CLEANout)
	#CLEANout.valid <- validCLEAN(CLEANout)
	return(CLEANout)
}

