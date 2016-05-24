setClass("METAANALYSIS",
	representation = representation(
						strEqcCommand				=	"character",
						acolBETAs					=	"character",
						acolSEs						=	"character",
						colOutBeta					=	"character",
						colOutSe					=	"character",
						acolZSCOREs					=	"character",
						acolNs						=	"character",
						colOutZscore				=	"character",
						colOutN						=	"character",
						colOutP						=	"character",
						acolA1s						=	"character",
						acolA2s						=	"character",
						strMethod					=	"character"
						),
	prototype = prototype(
						strEqcCommand				=	"",
						acolBETAs					=	"",
						acolSEs						=	"",
						colOutBeta					=	"bmeta",
						colOutSe					=	"semeta",
						acolZSCOREs					=	"",
						acolNs						=	"",
						colOutZscore				=	"zmeta",
						colOutN						=	"nmeta",
						colOutP						=	"pmeta",
						acolA1s						=	"",
						acolA2s						=	"",
						strMethod					=	""
						)
)
setGeneric("setMETAANALYSIS", function(object) standardGeneric("setMETAANALYSIS"))
setMethod("setMETAANALYSIS", signature = (object = "METAANALYSIS"), function(object) {
	
	#aEqcSlotNamesIn = c("colBeta1", "colSe1", "colBeta2" , "colSe2", "blnUseGoncalo","bln2sided","rcdTestDirection","strPdiffName")
	aEqcSlotNamesIn = c("acolBETAs", "acolSEs", "colOutBeta", "colOutSe", "acolZSCOREs", "acolNs", "colOutZscore", "colOutN", "colOutP","acolA1s","acolA2s")
	
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
validMETAANALYSIS <- function(objMA) {
	
	if(any(objMA@acolBETAs != "") | any(objMA@acolSEs != "")) strMethod = "FIXEDEFFECT"
	else if(any(objMA@acolZSCOREs != "")) strMethod = "ZSCOREBASED"
	else stop(paste(" EASY ERROR:METAANALYSIS\n Please either set (--acolBETAs AND --acolSEs) OR (--acolZSCOREs AND --acolNs) !!!", sep=""))

	if(strMethod == "FIXEDEFFECT") {
		if(length(objMA@acolBETAs)<2) 
			stop(paste(" EASY ERROR:METAANALYSIS\n Please set at least 2 beta columns for acolBETAs (e.g. --acolBETAs beta1;beta2) !!!", sep=""))
		if(length(objMA@acolSEs)<2) 
			stop(paste(" EASY ERROR:METAANALYSIS\n Please set at least 2 se columns for acolSEs (e.g. --acolSEs se1;se2) !!!", sep=""))
		if(length(objMA@acolBETAs) != length(objMA@acolSEs))
			stop(paste(" EASY ERROR:METAANALYSIS\n Number of defined Beta columns acolBETAs does not match number of defined Se columns acolSEs. \n PLease check!!!", sep=""))
		if(any(objMA@acolBETAs == ""))
			stop(paste(" EASY ERROR:METAANALYSIS\n No Effect estimate columns acolBETAs defined.\n Please set acolBETAs.", sep=""))
		if(any(objMA@acolSEs == ""))
			stop(paste(" EASY ERROR:METAANALYSIS\n No Effect estimate columns acolSEs defined.\n Please set acolSEs.", sep=""))
		if(any(objMA@acolA1s!="")) {
			if(length(objMA@acolA1s) != length(objMA@acolBETAs))
				stop(paste(" EASY ERROR:METAANALYSIS\n Number of defined allele A1 columns acolA1s does not match number of defined Beta columns acolBETAs. \n PLease check!!!", sep=""))
		}
		if(any(objMA@acolA2s!="")) {
			if(length(objMA@acolA2s) != length(objMA@acolBETAs))
				stop(paste(" EASY ERROR:METAANALYSIS\n Number of defined allele A2 columns acolA1s does not match number of defined Beta columns acolBETAs. \n PLease check!!!", sep=""))
		}
	} else {
		if(length(objMA@acolZSCOREs)<2) 
			stop(paste(" EASY ERROR:METAANALYSIS\n Please set at least 2 zscore columns for acolZSCOREs (e.g. --acolZSCOREs z1;z2) !!!", sep=""))
		if(length(objMA@acolNs)<2) 
			stop(paste(" EASY ERROR:METAANALYSIS\n Please set at least 2 n columns for acolNs (e.g. --acolNs n1;n2) !!!", sep=""))
		if(length(objMA@acolZSCOREs) != length(objMA@acolNs))
			stop(paste(" EASY ERROR:METAANALYSIS\n Number of defined Zscore columns acolZSCOREs does not match number of defined N columns acolNs. \n PLease check!!!", sep=""))
		if(any(objMA@acolZSCOREs == ""))
			stop(paste(" EASY ERROR:METAANALYSIS\n No Zscore columns acolZSCOREs defined.\n Please set acolZSCOREs.", sep=""))
		if(any(objMA@acolNs == ""))
			stop(paste(" EASY ERROR:METAANALYSIS\n No N columns acolNs defined.\n Please set acolNs.", sep=""))
		if(any(objMA@acolA1s!="")) {
			if(length(objMA@acolA1s) != length(objMA@acolZSCOREs))
				stop(paste(" EASY ERROR:METAANALYSIS\n Number of defined allele A1 columns acolA1s does not match number of defined Zscore columns acolZSCOREs. \n PLease check!!!", sep=""))
		}
		if(any(objMA@acolA2s!="")) {
			if(length(objMA@acolA2s) != length(objMA@acolZSCOREs))
				stop(paste(" EASY ERROR:METAANALYSIS\n Number of defined allele A2 columns acolA1s does not match number of defined Zscore columns acolZSCOREs. \n PLease check!!!", sep=""))
		}
	}
	
	objMA@strMethod <- strMethod
	
	return(objMA)
}

#############################################################################################################################
METAANALYSIS.run <- function(objMA, objGWA) {
	
	
	if(objMA@strMethod == "FIXEDEFFECT") {
		nDf = length(objMA@acolBETAs)
		
		#b: ((Effect.MEN/(StdErr.MEN^2))+(Effect.WOMEN/(StdErr.WOMEN^2)))/(1/StdErr.MEN^2+1/StdErr.WOMEN^2)
		#se: sqrt(1/(1/StdErr.MEN^2+1/StdErr.WOMEN^2))
		#p: pchisq((Effect.OVERALL/StdErr.OVERALL)^2,df=1,lower.tail=FALSE)
		

		a1ref <- a2ref <- c()
		
		isUseAlleles <- isUseN <- FALSE 
		
		if(length(objMA@acolA1s)==length(objMA@acolBETAs)&length(objMA@acolA1s)>1)
			isUseAlleles = TRUE
		
		if(all(objMA@acolNs!="")) 
			isUseN <- TRUE
				
		for(i in 1:nDf) {
			
			iBeta = match(objMA@acolBETAs[i], objGWA@aHeader)
			iSe = match(objMA@acolSEs[i], objGWA@aHeader)
			if(isUseN) {
				iN = match(objMA@acolNs[i], objGWA@aHeader)
				if(is.na(iN))  stop(paste("EASY ERROR:METAANALYSIS\nColumn \n",objMA@acolNs[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
				aN 	<- objGWA@tblGWA[,iN]
				if(class(aN) != "numeric" & class(aN) != "double" ) stop(paste("EASY ERROR:METAANALYSIS\nColumn \n",objMA@acolNs[i],"\nis not a numeric variable in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
				if(i == 1) nmeta <- rep(0,numSNPs)
			}
			
			if(is.na(iBeta)) stop(paste("EASY ERROR:METAANALYSIS\nColumn \n",objMA@acolBETAs[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			if(is.na(iSe))  stop(paste("EASY ERROR:METAANALYSIS\nColumn \n",objMA@acolSEs[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			
			aBeta 	<- objGWA@tblGWA[,iBeta]
			aSe 	<- objGWA@tblGWA[,iSe]

			if(i == 1) {
				numSNPs = length(aBeta)
				bmeta.enum <- bmeta.denom <- nmetadf <- rep(0,numSNPs)
			}
			
			isOk = !is.na(aBeta)
			# rep(TRUE,numSNPs)
			
			if(class(aBeta) != "numeric" & class(aBeta) != "double" ) stop(paste("EASY ERROR:METAANALYSIS\nColumn \n",objMA@acolBETAs[i],"\nis not a numeric variable in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			if(class(aSe) != "numeric" & class(aSe) != "double" ) stop(paste("EASY ERROR:METAANALYSIS\nColumn \n",objMA@acolSEs[i],"\nis not a numeric variable in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			
			if(isUseAlleles) {
				iA1 = match(objMA@acolA1s[i], objGWA@aHeader)
				iA2 = match(objMA@acolA2s[i], objGWA@aHeader)
				if(is.na(iA1)) stop(paste("EASY ERROR:METAANALYSIS\nColumn \n",objMA@acolA1s[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
				if(is.na(iA2))  stop(paste("EASY ERROR:METAANALYSIS\nColumn \n",objMA@acolA2s[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
				
				if(i == 1) {
					a1ref <- objGWA@tblGWA[,iA1]
					a2ref <- objGWA@tblGWA[,iA2]
				} else {
					a1in <- objGWA@tblGWA[,iA1]
					a2in <- objGWA@tblGWA[,iA2]
					isSame = a1in==a1ref & a2in==a2ref
					isChange = a1in==a2ref & a2in==a1ref
					isChange[is.na(isChange)]<-FALSE
					aBeta[isChange] <- -aBeta[isChange]
					isOk = isSame|isChange
					isOk[is.na(isOk)] <- FALSE
				}
			}
			
			bmeta.enum[isOk] = bmeta.enum[isOk] + (aBeta[isOk]/aSe[isOk]^2)
			bmeta.denom[isOk] = bmeta.denom[isOk] + (1/aSe[isOk]^2)
			nmetadf[isOk] = nmetadf[isOk] + 1
			if(isUseN) nmeta[isOk] = nmeta[isOk] + aN[isOk]
		}
		
		bmeta = bmeta.enum / bmeta.denom
		semeta = sqrt(1/bmeta.denom)
		pmeta = pchisq((bmeta/semeta)^2, df = 1, lower.tail = FALSE)
		
		objGWA <- GWADATA.cbind(objGWA, bmeta, objMA@colOutBeta)
		objGWA <- GWADATA.cbind(objGWA, semeta, objMA@colOutSe)
		objGWA <- GWADATA.cbind(objGWA, pmeta, objMA@colOutP)
		objGWA <- GWADATA.cbind(objGWA, nmetadf, "nmeta.df")
		if(isUseN) objGWA <- GWADATA.cbind(objGWA, nmeta, objMA@colOutN)
	} else {
		## z score based
		nDf = length(objMA@acolZSCOREs)
		
		## Zmeta = sum(Zi*wi)/sqrt(sum(wi^2)); wi = sqrt(Ni)
		## P = 2*phi(abs(Z))

		a1ref <- a2ref <- c()
		
		isUseAlleles = FALSE 
		
		if(length(objMA@acolA1s)==length(objMA@acolBETAs)&length(objMA@acolA1s)>1)
			isUseAlleles = TRUE
		
		for(i in 1:nDf) {
			
			iZ = match(objMA@acolZSCOREs[i], objGWA@aHeader)
			iN = match(objMA@acolNs[i], objGWA@aHeader)
			
			if(is.na(iZ)) stop(paste("EASY ERROR:METAANALYSIS\nColumn \n",objMA@acolZSCOREs[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			if(is.na(iN)) stop(paste("EASY ERROR:METAANALYSIS\nColumn \n",objMA@acolNs[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			
			aZ 	<- objGWA@tblGWA[,iZ]
			aN 	<- objGWA@tblGWA[,iN]
			
			if(i == 1) {
				numSNPs = length(aZ)
				enum <- denom <- nmeta <- nmetadf <- rep(0,numSNPs)
			}
			
			isOk = !is.na(aZ)
			# rep(TRUE,numSNPs)
			
			if(class(aZ) != "numeric" & class(aZ) != "double" ) stop(paste("EASY ERROR:METAANALYSIS\nColumn \n",objMA@acolZSCOREs[i],"\nis not a numeric variable in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			if(class(aN) != "numeric" & class(aN) != "double" & class(aN) != "integer") stop(paste("EASY ERROR:METAANALYSIS\nColumn \n",objMA@acolNs[i],"\nis not a numeric variable in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			
			if(isUseAlleles) {
				iA1 = match(objMA@acolA1s[i], objGWA@aHeader)
				iA2 = match(objMA@acolA2s[i], objGWA@aHeader)
				if(is.na(iA1)) stop(paste("EASY ERROR:METAANALYSIS\nColumn \n",objMA@acolA1s[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
				if(is.na(iA2))  stop(paste("EASY ERROR:METAANALYSIS\nColumn \n",objMA@acolA2s[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
				
				if(i == 1) {
					a1ref <- objGWA@tblGWA[,iA1]
					a2ref <- objGWA@tblGWA[,iA2]
				} else {
					a1in <- objGWA@tblGWA[,iA1]
					a2in <- objGWA@tblGWA[,iA2]
					isSame = a1in==a1ref & a2in==a2ref
					isChange = a1in==a2ref & a2in==a1ref
					isChange[is.na(isChange)]<-FALSE
					aZ[isChange] <- -aZ[isChange]
					isOk = isSame|isChange
					isOk[is.na(isOk)] <- FALSE
				}
			}
			
			enum[isOk] = enum[isOk] + aZ[isOk]*sqrt(aN[isOk])
			denom[isOk] = denom[isOk] + sqrt(aN[isOk])^2
			nmetadf[isOk] = nmetadf[isOk] + 1
			nmeta[isOk] = nmeta[isOk] + aN[isOk]
		}
		
		zmeta = enum / sqrt(denom)
		# semeta = sqrt(1/bmeta.denom)
		pmeta = 2*pnorm(abs(zmeta), lower.tail = FALSE)
		
		objGWA <- GWADATA.cbind(objGWA, zmeta, objMA@colOutZscore)
		objGWA <- GWADATA.cbind(objGWA, nmeta, objMA@colOutN)
		objGWA <- GWADATA.cbind(objGWA, pmeta, objMA@colOutP)
		objGWA <- GWADATA.cbind(objGWA, nmetadf, "nmeta.df")
	}
	return(objGWA)
}

METAANALYSIS <- function(strEqcCommand){ 
	## Wrapper for class definition
	METAANALYSISout <- setMETAANALYSIS(new("METAANALYSIS", strEqcCommand = strEqcCommand))
	METAANALYSISout <- validMETAANALYSIS(METAANALYSISout)
	return(METAANALYSISout)

}
