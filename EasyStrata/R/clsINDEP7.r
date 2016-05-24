setClass("INDEP",
	representation = representation(
						strEqcCommand		=	"character",
						rcdCriterion		=	"character",
						colIndep			=	"character",
						blnIndepMin			=	"logical",
						colInChr			=	"character",
						colInPos			=	"character",
						numPosLim			=	"numeric",
						blnAddIndepInfo 	= 	"logical",
						strTag				=	"character",
						blnStepDown		 	= 	"logical"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						rcdCriterion		=	"",
						colIndep			=	"",
						blnIndepMin			=	TRUE,
						colInChr			=	"",
						colInPos			=	"",
						numPosLim			=	500000,
						blnAddIndepInfo 	= 	FALSE,
						strTag				=	"",
						blnStepDown		 	= 	FALSE
						)
	#contains = c("EcfReader")
)

setGeneric("setINDEP", function(object) standardGeneric("setINDEP"))
setMethod("setINDEP", signature = (object = "INDEP"), function(object) {
	
	aEqcSlotNamesIn = c("rcdCriterion","colIndep","blnIndepMin","colInChr","colInPos","numPosLim","blnAddIndepInfo","strTag","blnStepDown")

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
validINDEP <- function(objINDEP) {
	
	### Valid with GWADATA?
	
	if(objINDEP@rcdCriterion == "") 
		cat(paste(" EASY WARNING:INDEP\n No criterion rcdCriterion defined. All data will be used for independentisation.", sep=""))
	if(objINDEP@colIndep == "") 
		stop(paste(" EASY ERROR:INDEP\n No column colIndep defined. Please set colIndep.", sep=""))
	if(objINDEP@colInChr == "") 
		stop(paste(" EASY ERROR:INDEP\n No column colInChr defined. Please set colInChr.", sep=""))
	if(objINDEP@colInPos == "") 
		stop(paste(" EASY ERROR:INDEP\n No column colInPos defined. Please set colInPos.", sep=""))
	
	return(TRUE)
}
INDEP.GWADATA.valid <- function(objINDEP, objGWA) {
	
	isAv <- objINDEP@colIndep %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:INDEP\n Defined column colIndep \n",objINDEP@colIndep, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	isAv <- objINDEP@colInChr %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:INDEP\n Defined column colInChr \n",objINDEP@colInChr, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	isAv <- objINDEP@colInPos %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:INDEP\n Defined column colInPos \n",objINDEP@colInPos, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))

}

#############################################################################################################################
INDEP.run <- function(objINDEP, objGWA, objREPORT) {

		rcdCriterion 	<- objINDEP@rcdCriterion
		colIndep		<- objINDEP@colIndep
		blnIndepMin		<- objINDEP@blnIndepMin
		colInChr		<- objINDEP@colInChr
		colInPos		<- objINDEP@colInPos
		numPosLim		<- objINDEP@numPosLim
		blnAddIndepInfo	<- objINDEP@blnAddIndepInfo
		strTag			<- objINDEP@strTag
		blnStepDown		<- objINDEP@blnStepDown
		
		if(rcdCriterion != "") {
			objRCD 	<- RCD(rcdCriterion)
			out 	<- RCD.eval(objRCD, objGWA)
			numIndepCrit = length(which(out))
			objGWA.indep <- GWADATA.getrows(objGWA, which(out))
		} else {
			numIndepCrit <- nrow(objGWA@tblGWA)
			objGWA.indep <- objGWA
		}
		
		is_na_indep = is.na(GWADATA.getcol(objGWA.indep, colIndep)) | is.na(GWADATA.getcol(objGWA.indep, colInChr)) | is.na(GWADATA.getcol(objGWA.indep, colInPos))
	
		objGWA.indep@tblGWA <- 	objGWA.indep@tblGWA[!is_na_indep,]
	
		numIndepNA = length(which(is_na_indep))
		
		aTopHit <- aLociTag <- aNumLocusSNPs <- rep(NA, nrow(objGWA.indep@tblGWA))
		
		hitcount=1
		
		aIndepVal 	= GWADATA.getcol(objGWA.indep, colIndep)
		aChr 		= GWADATA.getcol(objGWA.indep, colInChr)
		aPos 		= GWADATA.getcol(objGWA.indep, colInPos)
		
		aIndepValBak <- aIndepVal
		
		while(any(is.na(aLociTag))) {
			
			if(blnIndepMin) iTmpExtr = which(aIndepVal == min(aIndepVal))[1]
			else iTmpExtr = which(aIndepVal == max(aIndepVal))[1]
			
			chrTmp =  aChr[iTmpExtr]
			posTmp =  aPos[iTmpExtr]
			
			isCurrentLocus = aChr == chrTmp & abs(aPos - posTmp) <= numPosLim
			isCurrentLocusAlreadySet = any(!is.na(aLociTag[which(isCurrentLocus)]))
			
			if(blnStepDown & isCurrentLocusAlreadySet) {
				### step down
				## get set Locus IDs
				aLocusIDsSet = unique(aLociTag[which(isCurrentLocus)])
				aLocusIDsSet = aLocusIDsSet[!is.na(aLocusIDsSet)]
				## usually this should be a single locus
				## however in theory it is possible that loci overlap to the left and to the right
				## get lociNum with lowest P and reset all SNPs in the locus to that 
				## this is always the lower locus number because the data set is sorted
				isCurrentLocusWide = isCurrentLocus | aLociTag%in%aLocusIDsSet
				locusIdNew = min(aLocusIDsSet)
				
				aLociTag[which(isCurrentLocusWide)] = locusIdNew
				aTopHit[which(isCurrentLocusWide & aTopHit!=locusIdNew)] = NA
				
			} else {
				aTopHit[iTmpExtr] = hitcount
				aLociTag[which(isCurrentLocus & is.na(aLociTag))] = hitcount
				hitcount = hitcount + 1
			}
			
			if(blnIndepMin) aIndepVal[which(isCurrentLocus)] = Inf
			else aIndepVal[which(isCurrentLocus)] = -Inf
			
		}
	
	objGWA.indep = GWADATA.cbind(objGWA.indep, aTopHit, "aTopHit")
	objGWA.indep = GWADATA.cbind(objGWA.indep, aLociTag, "aLociTag")
	### add number of locus SNPs
	for(locusNum in unique(na.omit(aTopHit))) {
		aNumLocusSNPs[aLociTag == locusNum] <- length(which(aLociTag == locusNum))
	}
	objGWA.indep = GWADATA.cbind(objGWA.indep, aNumLocusSNPs, "aNumLocusSNPs")
	
	objGWA.indep.x <- GWADATA.getrows(objGWA.indep, which(!is.na(aTopHit)))
	
	if(nchar(strTag)>0) strTag = paste(strTag,".",sep="") 
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numSNPIndepCrit",sep=""),numIndepCrit)
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numSNPIndepNA",sep=""),numIndepNA)
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numIndepLoci",sep=""),hitcount-1)
	
	if(blnAddIndepInfo) {
		### merges by intersect(names)
		objGWA.indep.tmp <- objGWA.indep
		objGWA.indep.tmp <- GWADATA.renamecol(objGWA.indep.tmp, "aLociTag", paste(strTag,"aLociTag",sep=""))
		objGWA.indep.tmp <- GWADATA.renamecol(objGWA.indep.tmp, "aTopHit", paste(strTag,"aTopHit",sep=""))
		objGWA.indep.tmp <- GWADATA.renamecol(objGWA.indep.tmp, "aNumLocusSNPs", paste(strTag,"aNumLocusSNPs",sep=""))
		
		if(paste(strTag,"aLociTag",sep="") %in% objGWA@aHeader) 
			stop("EASY ERROR:INDEP\nAdding column aLociTag to large data-set failed because column\n",paste(strTag,"aLociTag",sep=""),"\n is already present GWADATA. PLease use --strTag to influence the column name.")
		if(paste(strTag,"aTopHit",sep="") %in% objGWA@aHeader) 
			stop("EASY ERROR:INDEP\nAdding column aTopHit to large data-set failed because column\n",paste(strTag,"aTopHit",sep=""),"\n is already present GWADATA. PLease use --strTag to influence the column name.")
		if(paste(strTag,"aNumLocusSNPs",sep="") %in% objGWA@aHeader) 
			stop("EASY ERROR:INDEP\nAdding column aNumLocusSNPs to large data-set failed because column\n",paste(strTag,"aNumLocusSNPs",sep=""),"\n is already present GWADATA. PLease use --strTag to influence the column name.")
			
		### merges by intersect(names)		
		objGWA <- GWADATA.merge(objGWA,objGWA.indep.tmp, strSuffix.In = "", strSuffix.Add = "", blnAll.In = TRUE, blnAll.Add = TRUE, strBy.In = NA, strBy.Add = NA)
	}
	
	return(list(objGWA,objGWA.indep,objGWA.indep.x,objREPORT))
}

INDEP <- function(strEqcCommand){ 
	## Wrapper for class definition
	INDEPout <- setINDEP(new("INDEP", strEqcCommand = strEqcCommand))
	validINDEP(INDEPout)
	#INDEPout.valid <- validINDEP(INDEPout)
	return(INDEPout)
}

