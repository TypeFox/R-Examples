setClass("MERGE",
	representation = representation(
						strEqcCommand				=	"character",
						colInMarker					=	"character",
						strInSuffix					=	"character",
						blnInAll					=	"logical",
						fileRef						=	"character",
						colRefMarker				=	"character",
						strRefSuffix				=	"character",
						blnRefAll					=	"logical",
						fileRefTag					=	"character",
						blnWriteNotInRef			=	"logical",
						blnWriteNotInIn				=	"logical",
						strTag						=	"character"
						),
	prototype = prototype(
						strEqcCommand				=	"",
						colInMarker					=	"",
						strInSuffix					=	"",
						blnInAll					=	TRUE,
						fileRef						=	"",
						colRefMarker				=	"",
						strRefSuffix				=	"",
						blnRefAll					=	FALSE,
						fileRefTag					=	"1",
						blnWriteNotInRef			=	FALSE,
						blnWriteNotInIn				=	FALSE,
						strTag						=	""
						),
	contains = c("GWADATA")
)

setGeneric("setMERGE", function(object) standardGeneric("setMERGE"))
setMethod("setMERGE", signature = (object = "MERGE"), function(object) {
	
	aEqcSlotNamesIn = c("colInMarker", "strInSuffix", "blnInAll", "fileRef", "colRefMarker", "strRefSuffix", "blnRefAll","fileRefTag","blnWriteNotInRef","blnWriteNotInIn","strTag",
						"strMissing", "strSeparator", "acolIn", "acolInClasses","acolNewName") 
	
	### Last 4 are inherited from class GWADATA and can be used with MERGE for reference file!
						

	objEqcReader <- EqcReader(object@strEqcCommand,aEqcSlotNamesIn)
	
	if(length(objEqcReader@lsEqcSlotsOut) > 0) {
		for(i in 1:length(objEqcReader@lsEqcSlotsOut)) {
			tmpSlot <- names(objEqcReader@lsEqcSlotsOut)[i]
			tmpSlotVal <- objEqcReader@lsEqcSlotsOut[[i]]
			
			if(all(!is.na(tmpSlotVal))) slot(object, tmpSlot) <- tmpSlotVal
		}
	}
	
	#object@pathOut <- getwd()
	object@pathOut <- "NA"
	
	return(object)
})

#############################################################################################################################
MERGE.init <- function(objMERGE) {
	
	
	
	if(objMERGE@fileRef == "")
		stop(paste(" EASY ERROR:MERGE\n No reference file fileRef defined. \n Please set fileRef that will be merged to alll GWA data-sets.", sep=""))
	
	if(!file.exists(objMERGE@fileRef))
		stop(paste("EASY ERROR:MERGE\n Reference file \n ",objMERGE@fileRef,"\n that should be merged to input, does not exist!!!\n", sep=""))
	
	if(objMERGE@colRefMarker == "")
		stop(paste(" EASY ERROR:MERGE\n No reference Marker column defined for fileRef. \n Please define colRefMarker that will be used for merging the data-set.", sep=""))
		
	if(objMERGE@colInMarker == "")
		stop(paste(" EASY ERROR:MERGE\n No input Marker column defined. \n Please define colInMarker that will be used for merging the data-set.", sep=""))
	
	objMERGE@fileIn <- objMERGE@fileRef
	
	objMERGE@pathOut = "NA"
	
	# if(objMERGE@strTag == "") 
		# stop(paste(" EASY ERROR:MERGE\n Requested parameter 'strTag' missing. \n Please ensure that unique 'strTag' values are set for all your MERGE,ADJUSTALLELES and AFCHECK commands.", sep=""))
	
	objMERGE.init <- GWADATA.init(objMERGE)
	
	if(!(objMERGE.init@colRefMarker%in%objMERGE.init@aHeader)) 
		stop(paste(" EASY ERROR:MERGE\n Reference Marker column colRefMarker \n",objMERGE@colRefMarker," \n cannot be found in fileRef \n", objMERGE@fileRef," \n Please correct colRefMarker.", sep=""))
		
	return(objMERGE.init)
}

MERGE.GWADATA.valid <- function(objMERGE, objGWA) {
	
	isAv <- objMERGE@colInMarker %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:MERGE\n Defined column colInMarker \n",objMERGE@colInMarker, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))

	### Check duplicate colnames in merged data set

						
	aHeadIn <- objGWA@aHeader
	iMarkerIn = match(objMERGE@colInMarker, aHeadIn)
	aHeadIn[-iMarkerIn] <- paste(aHeadIn[-iMarkerIn],objMERGE@strInSuffix,sep="")
	
	aHeadRef <- objMERGE@aHeader
	iMarkerRef = match(objMERGE@colRefMarker, aHeadRef)
	aHeadRef[-iMarkerRef] <- paste(aHeadRef[-iMarkerRef],objMERGE@strRefSuffix,sep="")
	
	aHeadIntersect = intersect(aHeadIn[-iMarkerIn],aHeadRef[-iMarkerRef])
	if(length(aHeadIntersect)>0)
		stop(paste(" EASY ERROR:MERGE\n Column \n",paste(aHeadIntersect,collapse="\n"), 
					"\n will be duplicated after merging data to \n",objGWA@fileIn,
					"\n PLease specify correct suffixes strInSuffix and strRefSuffix or rename column in one file.", sep="")
			)
}


#############################################################################################################################
# MERGE.read.ref <- function(objMERGE) {

	# objMERGE <- GWADATA.read(objMERGE)

	# return(objMERGE)
# }

# MERGE.read.10rows <- function(objMERGE) {

	# objMERGE <- GWADATA.read.10rows(objMERGE)

	# return(objMERGE)
# }
#############################################################################################################################
MERGE.run <- function(objMERGE, objGWA, objREPORT, isValidScript) {
	
	
	colInMarker 	<- objMERGE@colInMarker
	colRefMarker 	<- objMERGE@colRefMarker
	
	strPrefix <- ifelse(objMERGE@strTag!="", paste(objMERGE@strTag,".",sep=""), "")
	
	objGWA.NotInRef <- GWADATA.copy(objGWA)
	objMERGE.NotInIn <- GWADATA.copy(objMERGE)
	objMERGE.NotInIn@pathOut <- objGWA@pathOut
	objMERGE.NotInIn@fileInShortName <- rev(strsplit(objMERGE.NotInIn@fileRef,"/",fixed=T)[[1]])[1]
			
	#aMarkerIn <- GWADATA.getcol(objGWA,objMERGE@colInMarker)
	#aMarkerRef <- GWADATA.getcol(objMERGE,objMERGE@colRefMarker)
	
	isNotInRef 	<- !objGWA@tblGWA[,colInMarker]%in%objMERGE@tblGWA[,colRefMarker]
	isNotInIn 	<- !objMERGE@tblGWA[,colRefMarker]%in%objGWA@tblGWA[,colInMarker]
	
	objREPORT <- REPORT.addval(objREPORT,paste(strPrefix,"NotInRef",sep=""),length(which(isNotInRef)))
	objREPORT <- REPORT.addval(objREPORT,paste(strPrefix,"NotInIn",sep=""),length(which(isNotInIn)))
	
	if(isValidScript) {
		if(objMERGE@blnWriteNotInRef & any(isNotInRef)) {
			objGWA.NotInRef <- GWADATA.getrows(objGWA,which(isNotInRef))
			GWADATA.write(objGWA.NotInRef, strSuffix = paste(".",strPrefix,"notinref",sep=""))
			rm(objGWA.NotInRef)
		}
		if(objMERGE@blnWriteNotInIn & any(isNotInIn)) {
			objMERGE.NotInIn <- GWADATA.getrows(objMERGE,which(isNotInIn))		
			GWADATA.write(objMERGE.NotInIn, strSuffix = paste(".",strPrefix,"notinin",sep=""))
			rm(objMERGE.NotInIn)
		}
	}
	objGWA <- GWADATA.merge(objGWA, objMERGE, 
				strSuffix.In = objMERGE@strInSuffix, 
				strSuffix.Add = objMERGE@strRefSuffix, 
				blnAll.In = objMERGE@blnInAll, 
				blnAll.Add = objMERGE@blnRefAll, 
				strBy.In = objMERGE@colInMarker, 
				strBy.Add = objMERGE@colRefMarker) 
	
	return(list(objGWA,objREPORT))
}
#############################################################################################################################
MERGE <- function(strEqcCommand){ 
	## Wrapper for class definition
	MERGEout <- setMERGE(new("MERGE", strEqcCommand = strEqcCommand))
	MERGEout.init <- MERGE.init(MERGEout) 
	return(MERGEout.init)

}
