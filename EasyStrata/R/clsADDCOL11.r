setClass("ADDCOL",
	representation = representation(
						strEqcCommand		=	"character",
						rcdAddCol			=	"character",
						colOut				=	"character",
						blnOverwrite		=	"logical"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						rcdAddCol			=	"",
						colOut				=	"NewCol",
						blnOverwrite		=	TRUE
						)
	#contains = c("EcfReader")
)

setGeneric("setADDCOL", function(object) standardGeneric("setADDCOL"))
setMethod("setADDCOL", signature = (object = "ADDCOL"), function(object) {
	
	aEqcSlotNamesIn = c("rcdAddCol", "colOut", "blnOverwrite")
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
validADDCOL <- function(objADDCOL) {
	
	# if(length(objADDCOL@arcdAddCol)>length(objADDCOL@acolOut)) 
		# stop(paste("EASY ERROR:ADDCOL\n No new Column name defined in acolOut for RCD \n",objADDCOL@arcdAddCol[length(objADDCOL@acolOut)+1]," !!!", sep=""))
	# if(length(objADDCOL@arcdAddCol)<length(objADDCOL@acolOut)) 
		# stop(paste("EASY ERROR:ADDCOL\n No rcd defined in arcdAddCol for New Column \n",objADDCOL@acolOut[length(objADDCOL@arcdAddCol)+1]," !!!", sep=""))
	
	if(objADDCOL@colOut == "") 
		stop(paste(" EASY ERROR:ADDCOL\n No new column name colOut defined for rcdAddCol \n ",objADDCOL@rcdAddCol," \n Please set colOut.", sep=""))
	
	if(objADDCOL@rcdAddCol == "") 
		warning(paste(" EASY WARNING:ADDCOL\n rcdAddCol is NA.\n Added column ",objADDCOL@colOut, " only contains NA values.", sep=""))
		#stop(paste(" EASY ERROR:ADDCOL\n No rcdAddCol defined.\n Please set rcdAddCol or remove ADDCOL function.", sep=""))
	
	return(TRUE)
}

#############################################################################################################################
ADDCOL.run <- function(objADDCOL, objGWA) {

	rcdAddCol 		<- objADDCOL@rcdAddCol
	colOut 	<- objADDCOL@colOut
	blnOverwrite 	<- objADDCOL@blnOverwrite
	
	objRCD 	<- RCD(rcdAddCol)
	out 	<- RCD.eval(objRCD, objGWA)
	
	#if(length(out) == 1) out = rep(out, dim(objGWA)[1])
	if(length(out) == 1) out = rep(out, dim(objGWA@tblGWA)[1])
	else if(is.null(out)) out = rep(NA, dim(objGWA@tblGWA)[1])
	
	blnColExist = any(names(objGWA@tblGWA) == colOut)
	
	if((blnColExist & blnOverwrite) | !blnColExist) 
		objGWA <- GWADATA.cbind(objGWA, out, colOut, blnOverwrite)
	## else No Action if column exists and should not be overwritten!
		
	return(objGWA)
}

ADDCOL <- function(strEqcCommand){ 
	## Wrapper for class definition
	ADDCOLout <- setADDCOL(new("ADDCOL", strEqcCommand = strEqcCommand))
	validADDCOL(ADDCOLout)
	#ADDCOLout.valid <- validADDCOL(ADDCOLout)
	return(ADDCOLout)
	#validECF(ECFout)
	#return(ECFout)
	
	## Identical:
	# ECFin <- new("ECF5", fileECF = fileECFIn) 
	# ECFout <- setECF5(ECFin)
	# return(ECFout)
}

# setValidity("ADDCOL", function(object){
	# print("ADDCOL-CHECK")
	
	
	
	# print(TRUE)
	# return(TRUE)
# })

