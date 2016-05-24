setClass("GC",
	representation = representation(
						strEqcCommand				=	"character",
						colPval						=	"character",
						colSE						=	"character",
						numLambda					=	"numeric",
						fileGcSnps					=	"character",
						strTag						=	"character",
						colGcSnpsMarker				=	"character",
						colInMarker					=	"character",
						blnSuppressCorrection		=	"logical"
						),
	prototype = prototype(
						strEqcCommand				=	"",
						colPval						=	"",
						colSE						=	"",
						numLambda					=	-1,
						fileGcSnps					=	"",
						strTag						=	"",
						colGcSnpsMarker				=	"",
						colInMarker					=	"",
						blnSuppressCorrection		=	FALSE
						)
	#contains = c("EcfReader")
)

setGeneric("setGC", function(object) standardGeneric("setGC"))
setMethod("setGC", signature = (object = "GC"), function(object) {
	
	aEqcSlotNamesIn = c("colPval", "colSE", "numLambda" , "fileGcSnps", "colGcSnpsMarker", "colInMarker", "blnSuppressCorrection","strTag")
	#aEcfSlotNamesIn = c("arcdAddCol", "astrAddColNames")

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
validGC <- function(objGC) {
	## Paste validity checks specific for GC
	## Pval <-> se
	## fileGcSnps exists
	## colMarkerGcSnps exists
	
	if(objGC@colPval == "") 
		stop(paste(" EASY ERROR:GC\n No P-Value column colPval defined for GC.\n Please set colPval.", sep=""))
	
	if(objGC@fileGcSnps != "") {
		if(!file.exists(objGC@fileGcSnps))
			stop(paste("EASY ERROR:GC\n File fileGcSnps\n ",objGC@fileGcSnps,"\n does not exist.\n Please check path or remove --fileGcSnps.", sep=""))
		if(objGC@colGcSnpsMarker == "")
			stop(paste(" EASY ERROR:GC\n No SNP file Marker column defined. \n Please define colGcSnpsMarker that will be used for merging the SNP file to the data-set.", sep=""))
		if(objGC@colInMarker == "")
			stop(paste(" EASY ERROR:GC\n No input Marker column defined. \n Please define colInMarker that will be used for merging the data-set to the SNP file.", sep=""))
	}
	
	return(TRUE)
}

GC.GWADATA.valid <- function(objGC, objGWA) {
	
	if(objGC@fileGcSnps != "") {
		isAv <- objGC@colInMarker %in% objGWA@aHeader
		if(!isAv)
			stop(paste(" EASY ERROR:GC\n Defined column colInMarker \n",objGC@colInMarker, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	}
	
}

#############################################################################################################################
GC.run <- function(objGC, objGWA, objREPORT) {
	
	strTag <- objGC@strTag
	
	tblReport <- data.frame()

	colPval <- objGC@colPval
	colSE 	<- objGC@colSE
	
	iPval	= match(colPval,objGWA@aHeader)
	if(is.na(iPval)) 
		stop(paste("EASY ERROR:GC\nColumn \n",colPval,"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
	
	iSe		= match(colSE,objGWA@aHeader)

	aPvalNoGc <- aPvalNoGcCalc <- objGWA@tblGWA[,iPval]
	aChisqNoGc<- aChisqNoGcCalc<- qchisq(aPvalNoGc,df=1,lower.tail=F)
	
	if(objGC@numLambda == -1) {
		## Calculate Lambda
		
		if(objGC@fileGcSnps != "") {
			# If file is defined, use SNPs in file for calculation of Lambda
			icolMarker = match(objGC@colInMarker, objGWA@aHeader)
			aMarker <- objGWA@tblGWA[,icolMarker]
			
			tblGcSnps<-read.table(objGC@fileGcSnps,header=T, sep="",  stringsAsFactors=FALSE)
			iMarkerGcSnps = match(objGC@colGcSnpsMarker,names(tblGcSnps))[1]
			
			if(is.na(iMarkerGcSnps))
				stop(paste(" EASY ERROR:GC\n Column colGcSnpsMarker \n ",objGC@colGcSnpsMarker," does not exist in fileGcSnps.\n Please check the Marker column name.", sep=""))
			
			aMarkerGcSnps = tblGcSnps[,iMarkerGcSnps]
			
			# iGcCalc1 = match(aMarker,aMarkerGcSnps)
			# iGcCalc2 = which(!is.na(iGcCalc1))
			
			isUsedForGc = aMarker%in%aMarkerGcSnps
			
			aPvalNoGcCalc = aPvalNoGc[which(isUsedForGc)]
			aChisqNoGcCalc= aChisqNoGcCalc[which(isUsedForGc)]
		}
	
		# strAlgLambdaCalc == "median"
		# use median (like metal does)
		denom<-qchisq(0.5, df=1)
		objGC@numLambda<-median(aChisqNoGcCalc,na.rm=T)/denom
		
	}
	
	strNewColNameReal <- strNewColName <- paste(colPval,".GC",sep="")
	
	
	if(!objGC@blnSuppressCorrection) {
		
		if(objGC@numLambda > 1 & !is.na(objGC@numLambda)) {
			## Correct, if Lambda>1
			aChisq1Gc<-aChisqNoGc/objGC@numLambda
			aPval1Gc<-signif(pchisq(aChisq1Gc,df=1,lower.tail=F),4)
		} else {
			aPval1Gc <- signif(pchisq(aChisqNoGc,df=1,lower.tail=F),4)
		}
		
		
		objGWA <- GWADATA.cbind(objGWA, aPval1Gc, strNewColName)
		strNewColNameReal <- objGWA@aHeader[length(objGWA@aHeader)] ## Maybe change to Pvalue.GC.1
		
		if(!is.na(iSe)) {
			aSeNoGc<-objGWA@tblGWA[,iSe]
			if(objGC@numLambda > 1 & !is.na(objGC@numLambda)) {
				aChisq1Gc<-aChisqNoGc/objGC@numLambda
				aSe1Gc<-signif(aSeNoGc*sqrt(objGC@numLambda),4)
			} else {
				aSe1Gc <- aSeNoGc
			}
			strAddLikePval = sub(colPval,"",strNewColNameReal)
			strNewColName <- paste(colSE,strAddLikePval,sep="")
			objGWA <- GWADATA.cbind(objGWA, aSe1Gc, strNewColName)
			
			# tblOut<-cbind(tblOut,aSe1Gc)
			# names(tblOut)[ncol(tblOut)] = paste(colSe,strColAdd,sep="")
		}
	}
#	tblReport = data.frame(objGC@numLambda)		
#	names(tblReport)[1] = paste(objGC@strTag,"Lambda.",strNewColNameReal,sep="")
#	objREPORT <- REPORT.addval(objREPORT,tblReport)
	if(nchar(strTag)>0) strTag = paste(strTag,".",sep="")
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"Lambda.",strNewColNameReal,sep=""), objGC@numLambda)

	lsOut <- list(objGWA, objREPORT)
	
	return(lsOut)
}

GC <- function(strEqcCommand){ 
	## Wrapper for class definition
	GCout <- setGC(new("GC", strEqcCommand = strEqcCommand))
	validGC(GCout)
	#ADDCOLout.valid <- validADDCOL(ADDCOLout)
	return(GCout)
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

