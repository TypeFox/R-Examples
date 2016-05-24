setClass("EVALSTAT",
	representation = representation(
						strEqcCommand	=	"character",
						colStat			=	"character",
						strTag			=	"character"
						),
	prototype = prototype(
						strEqcCommand	=	"",
						colStat			=	"",
						strTag			=	""
						),
	#contains = c("EcfReader")
)

setGeneric("setEVALSTAT", function(object) standardGeneric("setEVALSTAT"))
setMethod("setEVALSTAT", signature = (object = "EVALSTAT"), function(object) {
	
	aEqcSlotNamesIn = c("colStat","strTag")

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


EVALSTAT.run <- function(objEVALSTAT, objGWA, objREPORT) {
	### ... 
	tblStat <- data.frame()
	
	strTag = objEVALSTAT@strTag

	colStat = objEVALSTAT@colStat
	iCol=match(colStat,names(objGWA@tblGWA))
	
	if(!is.na(iCol)) {
		
		aTmp<- as.numeric(as.character(objGWA@tblGWA[,iCol]))		
		
		if(all(is.na(aTmp))) {
			tblStat[1,1]	=	NA
			tblStat[1,2]	=	NA
			tblStat[1,3]	=	NA
			tblStat[1,4]	=	NA
			tblStat[1,5]	=	NA
			tblStat[1,6]	=	NA
			tblStat[1,7]	=	NA
			tblStat[1,8]	=	NA
			tblStat[1,9]	=	NA
		
		} else {
			# vals[1]	=	length(aTmp)
			# vals[2]	=	length(which(is.na(aTmp)))
			
			tblStat[1,1]	=	length(aTmp)
			tblStat[1,2]	=	length(which(is.na(aTmp)))
			tblStat[1,3]	=	min(aTmp,na.rm=TRUE)
			tblStat[1,4]	=	max(aTmp,na.rm=TRUE)
			tblStat[1,5]	=	median(aTmp,na.rm=TRUE)
			tblStat[1,6]	=	quantile(aTmp,na.rm=TRUE)[[2]]
			tblStat[1,7]	=	quantile(aTmp,na.rm=TRUE)[[4]]
			tblStat[1,8]	=	mean(aTmp,na.rm=TRUE)
			tblStat[1,9]	=	sd(aTmp,na.rm=TRUE)
		}
		
		if(nchar(strTag)>0) strTag = paste(strTag,".",sep="") 
		names(tblStat)[1] = paste(strTag,colStat,"_num",sep="")
		names(tblStat)[2] = paste(strTag,colStat,"_NA",sep="")
		names(tblStat)[3] = paste(strTag,colStat,"_min",sep="")
		names(tblStat)[4] = paste(strTag,colStat,"_max",sep="")
		names(tblStat)[5] = paste(strTag,colStat,"_median",sep="")
		names(tblStat)[6] = paste(strTag,colStat,"_p25",sep="")
		names(tblStat)[7] = paste(strTag,colStat,"_p75",sep="")
		names(tblStat)[8] = paste(strTag,colStat,"_mean",sep="")
		names(tblStat)[9] = paste(strTag,colStat,"_sd",sep="")
		
		objREPORT <- REPORT.addval(objREPORT,tblStat)
		
	} else {
		stop(paste(" EASY ERROR: EVALSTAT\n Column \n",colStat,"\n is not available in file \n",objGWA@fileInShortName,"\n !!!" ,sep="" ))			
	}

	return(objREPORT)

}

EVALSTAT <- function(strEqcCommand){ 
	## Wrapper for class definition
	EVALSTATout <- setEVALSTAT(new("EVALSTAT", strEqcCommand = strEqcCommand))
	
	return(EVALSTATout)
	#validECF(ECFout)
	#return(ECFout)
	
	## Identical:
	# ECFin <- new("ECF5", fileECF = fileECFIn) 
	# ECFout <- setECF5(ECFin)
	# return(ECFout)
}

# setValidity("EVALSTAT", function(object){
	# print("EVALSTAT-CHECK")
	
	
	
	# print(TRUE)
	# return(TRUE)
# })

