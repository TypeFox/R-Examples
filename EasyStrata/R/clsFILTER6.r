setClass("FILTER",
	representation = representation(
						strEqcCommand				=	"character",
						rcdFilter					=	"character",
						strFilterName				=	"character"
						),
	prototype = prototype(
						strEqcCommand				=	"",
						rcdFilter					=	"",
						strFilterName				=	""
						)
	#contains = c("EcfReader")
)

setGeneric("setFILTER", function(object) standardGeneric("setFILTER"))
setMethod("setFILTER", signature = (object = "FILTER"), function(object) {
	
	aEqcSlotNamesIn = c("rcdFilter", "strFilterName")

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
validFILTER <- function(objFILTER) {
	
	if(objFILTER@strFilterName == "") 
		stop(paste(" EASY ERROR:FILTER\n No strFilterName defined for rcdFilter \n ",objFILTER@rcdFilter," \n Please set strFilterName.", sep=""))
	
	if(objFILTER@rcdFilter == "")
		stop(paste(" EASY ERROR:FILTER\n No rcdFilter defined.\n Please set rcdFilter or remove FILTER function.", sep=""))
	
	return(TRUE)
}

#############################################################################################################################
FILTER.run <- function(objFILTER, objGWA, objREPORT) {
	

	rcdFilter		<- objFILTER@rcdFilter
	strFilterName 	<- objFILTER@strFilterName
	
	objRCD 	<- RCD(rcdFilter)
	out 	<- RCD.eval(objRCD, objGWA)
	
	numOut = length(which(out))
	
	# objGWA.iFilter <- GWADATA.getrows(objGWA, which(out))
	# GWADATA.write(objGWA.Duplicates, "filter", "2")
	
	if(any(out)) objGWA <- GWADATA.getrows(objGWA, which(out))
	
	objREPORT <- REPORT.addval(objREPORT,strFilterName,numOut)

	lsOut <- list(objGWA, objREPORT)
	
	return(lsOut)
}

FILTER <- function(strEqcCommand){ 
	## Wrapper for class definition
	FILTERout <- setFILTER(new("FILTER", strEqcCommand = strEqcCommand))
	validFILTER(FILTERout)
	#FILTERout.valid <- validFILTER(FILTERout)
	return(FILTERout)
}

