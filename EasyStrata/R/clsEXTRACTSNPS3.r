setClass("EXTRACTSNPS",
	representation = representation(
						strEqcCommand	=	"character",
						colInMarker		=	"character",
						fileRef			=	"character",
						colRefMarker	=	"character", 
						strTag			=	"character"
						),
	prototype = prototype(
						strEqcCommand	=	"",
						colInMarker		=	"",
						fileRef			=	"",
						colRefMarker	=	"SNP",
						strTag			=	""
						)
	#contains = c("EcfReader")
)

setGeneric("setEXTRACTSNPS", function(object) standardGeneric("setEXTRACTSNPS"))
setMethod("setEXTRACTSNPS", signature = (object = "EXTRACTSNPS"), function(object) {
	
	aEqcSlotNamesIn = c("colInMarker","fileRef","colRefMarker","strTag")

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
validEXTRACTSNPS <- function(objEXTRACTSNPS) {
	
	if(objEXTRACTSNPS@colInMarker == "") 
		stop(paste(" EASY ERROR:EXTRACTSNPS\n No column colInMarker defined. Please set colInMarker.", sep=""))
	if(objEXTRACTSNPS@fileRef == "") 
		stop(paste(" EASY ERROR:EXTRACTSNPS\n No reference file defined. Please set fileRef.", sep=""))
	
	if(!file.exists(objEXTRACTSNPS@fileRef))
			stop(paste("EASY ERROR:EXTRACTSNPS\n File fileRef\n ",objEXTRACTSNPS@fileRef,"\n does not exist.", sep=""))
		### Cols exist?
		
	tblRef<-read.table(objEXTRACTSNPS@fileRef,header=T, sep="",  nrows=1, stringsAsFactors=FALSE)
		
	isAv <- objEXTRACTSNPS@colRefMarker %in% names(tblRef)
	if(!isAv)
		stop(paste(" EASY ERROR:EXTRACTSNPS\n Defined column colRefMarker \n",objEXTRACTSNPS@colRefMarker, "\n is not available in fileRef. PLease specify correct column name.", sep=""))
	
	return(TRUE)
}
EXTRACTSNPS.GWADATA.valid <- function(objEXTRACTSNPS, objGWA) {
	
	isNotAv <- !(objEXTRACTSNPS@colInMarker %in% objGWA@aHeader)
	if(isNotAv)
		stop(paste(" EASY ERROR:EXTRACTSNPS\n Defined column colInMarker \n",objEXTRACTSNPS@colInMarker, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	return(TRUE)
}
#############################################################################################################################
EXTRACTSNPS.run <- function(objEXTRACTSNPS, objGWA, objREPORT) {
	
	colInMarker		=	objEXTRACTSNPS@colInMarker
	fileRef 		=	objEXTRACTSNPS@fileRef
	colRefMarker	=	objEXTRACTSNPS@colRefMarker
	strTag			=	objEXTRACTSNPS@strTag
	
	tblRef <- read.table(objEXTRACTSNPS@fileRef, header=T, sep="\t", stringsAsFactors=FALSE)
	
	iMarker = match(colRefMarker, names(tblRef))
	aRefMarker = tblRef[,iMarker]
	aInMarker = GWADATA.getcol(objGWA, colInMarker)
	numMiss = length(which(!aRefMarker%in%aInMarker))
	
	if(nchar(strTag)>0) strTag = paste(strTag,".",sep="") 
	objREPORT <- REPORT.addval(objREPORT,paste(objEXTRACTSNPS@strTag,"numExtractMissing",sep=""),numMiss)
	
	tblOut <- merge(tblRef, objGWA@tblGWA, by.x = colRefMarker, by.y = colInMarker, all.x = TRUE, all.y = FALSE)
	
	return(list(tblOut, objREPORT))

}

EXTRACTSNPS <- function(strEqcCommand){ 
	## Wrapper for class definition
	EXTRACTSNPSout <- setEXTRACTSNPS(new("EXTRACTSNPS", strEqcCommand = strEqcCommand))
	validEXTRACTSNPS(EXTRACTSNPSout)
	#EXTRACTSNPSout.valid <- validEXTRACTSNPS(EXTRACTSNPSout)
	return(EXTRACTSNPSout)
}

