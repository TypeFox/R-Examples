setClass("FDR",
	representation = representation(
						strEqcCommand			=	"character",
						colPval					=	"character",
						strFdrMethod			=	"character",
						colOut					=	"character"
						),
	prototype = prototype(
						strEqcCommand			=	"",
						colPval					=	"",
						strFdrMethod			=	"BH",
						colOut					=	""
						)
)


setGeneric("setFDR", function(object) standardGeneric("setFDR"))
setMethod("setFDR", signature = (object = "FDR"), function(object) {
	
	#aEqcSlotNamesIn = c("colBeta1", "colSe1", "colBeta2" , "colSe2", "blnUseGoncalo","bln2sided","rcdTestDirection","strPdiffName")
	aEqcSlotNamesIn = c("colPval", "strFdrMethod","colOut")
	
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
validFDR <- function(objFDR) {
	
	
	if(objFDR@colPval == "")
		stop(paste(" EASY ERROR:FDR\n No P-Value column defined.\n Please set colPval.", sep=""))
	
	if(!objFDR@strFdrMethod %in% c("BH","BY"))
		stop(paste(" EASY ERROR:FDR\n strFdrMethod must either be 'BH' or 'BY'.\n Please correct strFdrMethod.", sep=""))
	
	return(TRUE)
}

#############################################################################################################################
FDR.run <- function(objFDR, objGWA) {
	
	colPval			 <- objFDR@colPval
	strFdrMethod	 <- objFDR@strFdrMethod
	colOut			 <- objFDR@colOut
	
	iP = match(colPval, objGWA@aHeader)
	if(is.na(iP)) stop(paste("EASY ERROR:FDR\nColumn \n",colPval,"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
	aP 	<- objGWA@tblGWA[,iP]
	if(class(aP) != "numeric" & class(aP) != "double" ) stop(paste("EASY ERROR:FDR\nColumn \n",colPval,"\nis not a numeric variable in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
	## remove NA????
	
	pfdr <- p.adjust(aP,strFdrMethod)
	
	strNewColName <- ifelse(colOut == "", paste(colPval,".",strFdrMethod,sep=""), colOut)
	
	objGWA <- GWADATA.cbind(objGWA, pfdr, strNewColName)
	
	return(objGWA)
}

FDR <- function(strEqcCommand){ 
	## Wrapper for class definition
	FDRout <- setFDR(new("FDR", strEqcCommand = strEqcCommand))
	validFDR(FDRout)
	return(FDRout)

}
