setClass("CALCULATE",
	representation = representation(
						strEqcCommand			=	"character",
						rcdCalc					=	"character",
						strCalcName				=	"character"
						),
	prototype = prototype(
						strEqcCommand			=	"",
						rcdCalc					=	"",
						strCalcName				=	""
						)
	#contains = c("EcfReader")
)

setGeneric("setCALCULATE", function(object) standardGeneric("setCALCULATE"))
setMethod("setCALCULATE", signature = (object = "CALCULATE"), function(object) {
	
	aEqcSlotNamesIn = c("rcdCalc", "strCalcName")

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
validCALCULATE <- function(objCALCULATE) {
	
	# if(length(objADDCOL@arcdAddCol)>length(objADDCOL@astrAddColNames)) 
		# stop(paste("EASY ERROR:ADDCOL\n No new Column name defined in astrAddColNames for RCD \n",objADDCOL@arcdAddCol[length(objADDCOL@astrAddColNames)+1]," !!!", sep=""))
	# if(length(objADDCOL@arcdAddCol)<length(objADDCOL@astrAddColNames)) 
		# stop(paste("EASY ERROR:ADDCOL\n No rcd defined in arcdAddCol for New Column \n",objADDCOL@astrAddColNames[length(objADDCOL@arcdAddCol)+1]," !!!", sep=""))
	
	if(objCALCULATE@strCalcName == "") 
		stop(paste(" EASY ERROR:CALCULATE\n No strCalcName defined for rcdCalc \n ",objCALCULATE@rcdCalc," \n Please set strCalcName.", sep=""))
	
	if(objCALCULATE@rcdCalc == "") 
		stop(paste(" EASY ERROR:CALCULATE\n No rcdCalc defined.\n Please set rcdCalc or remove CALCULATE function.", sep=""))
	
	
	return(TRUE)
}

#############################################################################################################################
CALCULATE.run <- function(objCALCULATE, objGWA, objREPORT) {
			
	rcdCalc 		<- objCALCULATE@rcdCalc
	strCalcName 	<- objCALCULATE@strCalcName
	
	objRCD 	<- RCD(rcdCalc)
	out 	<- RCD.eval(objRCD, objGWA)
	
	objREPORT <- REPORT.addval(objREPORT,strCalcName,as.numeric(out))

	return(objREPORT)
}

CALCULATE <- function(strEqcCommand){ 
	## Wrapper for class definition
	CALCULATEout <- setCALCULATE(new("CALCULATE", strEqcCommand = strEqcCommand))
	validCALCULATE(CALCULATEout)
	#CALCULATEout.valid <- validCALCULATE(CALCULATEout)
	return(CALCULATEout)
}

