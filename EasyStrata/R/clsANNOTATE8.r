setClass("ANNOTATE",
	representation = representation(
						strEqcCommand	=	"character",
						colInChr		=	"character",
						colInPos		=	"character",
						colInPval		=	"character",
						numPvalLim		= 	"numeric",
						fileAnnot		=	"character",
						colRefChr		=	"character",
						colRefPos		=	"character",
						strAnnotTag		=	"character",
						colOutAnnot		=	"character",
						numPosLim		=	"numeric"
						),
	prototype = prototype(
						strEqcCommand	=	"",
						colInChr		=	"",
						colInPos		=	"",
						colInPval		=	"",
						numPvalLim		= 	1,
						fileAnnot		=	"",
						colRefChr		=	"Chr",
						colRefPos		=	"Pos",
						strAnnotTag		=	"",
						colOutAnnot		=	"AnnotTag",
						numPosLim		=	500000
						)
	#contains = c("EcfReader")
)

setGeneric("setANNOTATE", function(object) standardGeneric("setANNOTATE"))
setMethod("setANNOTATE", signature = (object = "ANNOTATE"), function(object) {
	
	#aEqcSlotNamesIn = c("colInChr","colInPos","colInPval","numPvalLim","fileAnnot","colRefChr","colRefPos","strAnnotTag","colOutAnnot","numPosLim")
	aEqcSlotNamesIn = c("colInChr","colInPos","colInPval","numPvalLim","fileAnnot","strAnnotTag","colOutAnnot","numPosLim")

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
validANNOTATE <- function(objANNOTATE) {
	
	if(objANNOTATE@colInChr == "") 
		stop(paste(" EASY ERROR:ANNOTATE\n No column colInChr defined. Please set colInChr.", sep=""))
	if(objANNOTATE@colInPos == "") 
		stop(paste(" EASY ERROR:ANNOTATE\n No column colInPos defined. Please set colInPos.", sep=""))
	if(objANNOTATE@fileAnnot == "") 
		stop(paste(" EASY ERROR:ANNOTATE\n No reference file defined. Please set fileAnnot.", sep=""))
	
	### Valid with GWADATA?
	
	if(objANNOTATE@fileAnnot != "") {
		if(!file.exists(objANNOTATE@fileAnnot))
			stop(paste("EASY ERROR:ANNOTATE\n File fileAnnot\n ",objANNOTATE@fileAnnot,"\n does not exist.", sep=""))
		### Cols exist?
		
		tblAnnot<-read.table(objANNOTATE@fileAnnot,header=T, sep="",  stringsAsFactors=FALSE)
		
			isAv <- objANNOTATE@colRefChr %in% names(tblAnnot)
			if(!isAv)
				stop(paste(" EASY ERROR:ANNOTATE\n Defined column colRefChr \n",objANNOTATE@colRefChr, "\n is not available in fileAnnot. PLease specify correct column name.", sep=""))
			isAv <- objANNOTATE@colRefPos %in% names(tblAnnot)	
			if(!isAv)
				stop(paste(" EASY ERROR:ANNOTATE\n Defined column colRefPos \n",objANNOTATE@colRefPos, "\n is not available in fileAnnot. PLease specify correct column name.", sep=""))
			#isAv <- objANNOTATE@colRefTag %in% names(tblAnnot)
			#if(!isAv)
			#	stop(paste(" EASY ERROR:ANNOTATE\n Defined column colRefTag \n",objANNOTATE@colRefTag, "\n is not available in fileAnnot. PLease specify correct column name.", sep=""))
		
	}
	
	return(TRUE)
}
ANNOTATE.GWADATA.valid <- function(objANNOTATE, objGWA) {
	
	isNotAv <- !(objANNOTATE@colInChr %in% objGWA@aHeader)
	if(isNotAv)
		stop(paste(" EASY ERROR:ANNOTATE\n Defined column colInChr \n",objANNOTATE@colInChr, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	isNotAv <- !(objANNOTATE@colInPos %in% objGWA@aHeader)
	if(isNotAv)
		stop(paste(" EASY ERROR:ANNOTATE\n Defined column colInPos \n",objANNOTATE@colInPos, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	iPos = match(objANNOTATE@colInPos, objGWA@aHeader)
	isPosNumeric <- objGWA@aClasses[iPos] == "numeric" | objGWA@aClasses[iPos] == "integer" | objGWA@aClasses[iPos] == "double"

	if(!isPosNumeric)
		stop(paste(" EASY ERROR:ANNOTATE\n Defined column colInPos \n",objANNOTATE@colInPos, "\n is not numeric for GWA data-set \n",objGWA@fileIn,"\n . Please cast colInPos to numeric, integer or double.", sep=""))
	
	
	
	isDefinedNotAv <- objANNOTATE@colInPval != "" & !(objANNOTATE@colInPval %in% objGWA@aHeader)
	if(isDefinedNotAv)
		stop(paste(" EASY ERROR:ANNOTATE\n Defined column colInPval \n",objANNOTATE@colInPval, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	if(objANNOTATE@numPvalLim < 1) {
		if(!(objANNOTATE@colInPval %in% objGWA@aHeader))
			stop(paste(" EASY ERROR:ANNOTATE\n To use the Pvalue threshold for the annotation of loci, you have to define colInPval as well.","\n PLease specify colInPval with the ANNOTATE command !!!", sep=""))
		if(objANNOTATE@numPvalLim <=0)
			stop(paste(" EASY ERROR:ANNOTATE\n PLease specify numPvalLim within the range ]0,1[ within the ANNOTATE command !!!", sep=""))
	}
	
}
#############################################################################################################################
ANNOTATE.run <- function(objANNOTATE, objGWA) {
	
	colInChr		=	objANNOTATE@colInChr
	colInPos		=	objANNOTATE@colInPos
	colInPval		=	objANNOTATE@colInPval
	numPvalLim		=	objANNOTATE@numPvalLim
	colRefChr		=	objANNOTATE@colRefChr
	colRefPos		=	objANNOTATE@colRefPos
	strAnnotTag		=	objANNOTATE@strAnnotTag
	colOutAnnot		= 	objANNOTATE@colOutAnnot
	numPosLim		=	objANNOTATE@numPosLim
	
	tblLoci<-read.table(objANNOTATE@fileAnnot, header=T, sep="\t", stringsAsFactors=FALSE)
	numLoci = dim(tblLoci)[1]
	
	aInChr = GWADATA.getcol(objGWA, colInChr)
	aInPos = GWADATA.getcol(objGWA, colInPos)
	
	if(numPvalLim < 1) aInPval = GWADATA.getcol(objGWA, colInPval)
	else aInPval = rep(0, dim(objGWA@tblGWA)[1])
	
	aRefChr = tblLoci[, which(names(tblLoci) == colRefChr)]
	aRefPos = tblLoci[, which(names(tblLoci) == colRefPos)]
	#aRefTag = tblLoci[, which(names(tblLoci) == colRefTag)]
	
	
	if(colOutAnnot%in%objGWA@aHeader) 
		aTagOut = GWADATA.getcol(objGWA, colOutAnnot) # is already annotated
	else	
		aTagOut = rep(NA, dim(objGWA@tblGWA)[1]) # not used before
	
	for(i in 1:numLoci) {
		chrTmp = aRefChr[i]
		posTmp = aRefPos[i]
		#tagTmp = aRefTag[i]
		tagTmp = strAnnotTag
		
		#isCurrentLocus = aInChr == chrTmp & abs(aInPos - posTmp) <= numPosLim & aInPval <= numPvalLim
		isCurrentLocus = aInChr == chrTmp & abs(aInPos - posTmp) <= numPosLim
		if(any(isCurrentLocus)) {
			isAnyPvalLow = any(aInPval[isCurrentLocus] <= numPvalLim)
			if(isAnyPvalLow)
				aTagOut[isCurrentLocus] = ifelse(is.na(aTagOut[isCurrentLocus]), tagTmp, paste(aTagOut[isCurrentLocus],tagTmp, sep=";"))
		}
	}
	
	objGWA <- GWADATA.cbind(objGWA, aTagOut, colOutAnnot, blnOverwrite=TRUE)
	
	return(objGWA)
}

ANNOTATE <- function(strEqcCommand){ 
	## Wrapper for class definition
	ANNOTATEout <- setANNOTATE(new("ANNOTATE", strEqcCommand = strEqcCommand))
	validANNOTATE(ANNOTATEout)
	#ANNOTATEout.valid <- validANNOTATE(ANNOTATEout)
	return(ANNOTATEout)
}

