setClass("CALCPHET",
	representation = representation(
						strEqcCommand				=	"character",
						acolBETAs					=	"character",
						acolSEs						=	"character",
						acolZSCOREs					=	"character",
						acolNs						=	"character",
						colOutPhet					=	"character",
						strMethod					=	"character"
						),
	prototype = prototype(
						strEqcCommand				=	"",
						acolBETAs					=	"",
						acolSEs						=	"",
						acolZSCOREs					=	"",
						acolNs						=	"",
						colOutPhet					=	"phet",
						strMethod					=	""
						)
)
setGeneric("setCALCPHET", function(object) standardGeneric("setCALCPHET"))
setMethod("setCALCPHET", signature = (object = "CALCPHET"), function(object) {
	
	#aEqcSlotNamesIn = c("colBeta1", "colSe1", "colBeta2" , "colSe2", "blnUseGoncalo","bln2sided","rcdTestDirection","strPdiffName")
	aEqcSlotNamesIn = c("acolBETAs", "acolSEs", "colOutPhet", "acolZSCOREs", "acolNs")
	
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
validCALCPHET <- function(objCP) {
	
	if(any(objCP@acolBETAs != "") | any(objCP@acolSEs != "")) strMethod = "BETA"
	else if(any(objCP@acolZSCOREs != "") | any(objCP@acolNs != "")) strMethod = "ZSCORE"
	else stop(paste(" EASY ERROR:CALCPHET\n Please either set (--acolBETAs AND --acolSEs) OR (--acolZSCOREs AND --acolNs) !!!", sep=""))
	
	if(strMethod == "BETA") {
		if(length(objCP@acolBETAs)<2) 
			stop(paste(" EASY ERROR:CALCPHET\n Please set at least 2 beta columns for acolBETAs (e.g. --acolBETAs beta1;beta2) !!!", sep=""))
		if(length(objCP@acolSEs)<2) 
			stop(paste(" EASY ERROR:CALCPHET\n Please set at least 2 se columns for acolSEs (e.g. --acolSEs se1;se2) !!!", sep=""))
		if(length(objCP@acolBETAs) != length(objCP@acolSEs))
			stop(paste(" EASY ERROR:CALCPHET\n Number of defined Beta columns acolBETAs does not match number of defined Se columns acolSEs. \n PLease check!!!", sep=""))
		if(any(objCP@acolBETAs == ""))
			stop(paste(" EASY ERROR:CALCPHET\n No Effect estimate columns acolBETAs defined.\n Please set acolBETAs.", sep=""))
		if(any(objCP@acolSEs == ""))
			stop(paste(" EASY ERROR:CALCPHET\n No standard error columns acolSEs defined.\n Please set acolSEs.", sep=""))
	} else {
		if(length(objCP@acolZSCOREs)<2) 
			stop(paste(" EASY ERROR:CALCPHET\n Please set at least 2 Z score columns for acolZSCOREs (e.g. --acolZSCOREs Z1;Z2) !!!", sep=""))
		if(length(objCP@acolNs)<2) 
			stop(paste(" EASY ERROR:CALCPHET\n Please set at least 2 sample size columns for acolNs (e.g. --acolNs N1;N2) !!!", sep=""))
		if(length(objCP@acolZSCOREs) != length(objCP@acolNs))
			stop(paste(" EASY ERROR:CALCPHET\n Number of defined Z score columns acolZSCOREs does not match number of defined sample size columns acolNs. \n PLease check!!!", sep=""))
		if(any(objCP@acolZSCOREs == ""))
			stop(paste(" EASY ERROR:CALCPHET\n No Z score columns acolZSCOREs defined.\n Please set acolZSCOREs.", sep=""))
		if(any(objCP@acolNs == ""))
			stop(paste(" EASY ERROR:CALCPHET\n No sample size columns acolNs defined.\n Please set acolNs.", sep=""))		
	}
	
	objCP@strMethod <- strMethod
		
	return(objCP)
}

#############################################################################################################################
CALCPHET.run <- function(objCP, objGWA) {
	
	if(objCP@strMethod == "BETA") {
		nDf = length(objCP@acolBETAs)
		
		### Meta
		bmeta.enum <- bmeta.denom <- 0 
		
		for(i in 1:nDf) {
			### Meta
			iBeta = match(objCP@acolBETAs[i], objGWA@aHeader)
			iSe = match(objCP@acolSEs[i], objGWA@aHeader)
			
			if(is.na(iBeta)) stop(paste("EASY ERROR:CALCPHET\nColumn \n",objCP@acolBETAs[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			if(is.na(iSe))  stop(paste("EASY ERROR:CALCPHET\nColumn \n",objCP@acolSEs[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			
			aBeta 	<- objGWA@tblGWA[,iBeta]
			aSe 	<- objGWA@tblGWA[,iSe]
			
			if(class(aBeta) != "numeric" & class(aBeta) != "double" ) stop(paste("EASY ERROR:CALCPHET\nColumn \n",objCP@acolBETAs[i],"\nis not a numeric variable in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			if(class(aSe) != "numeric" & class(aSe) != "double" ) stop(paste("EASY ERROR:CALCPHET\nColumn \n",objCP@acolSEs[i],"\nis not a numeric variable in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			
			bmeta.enum = bmeta.enum + (aBeta/aSe^2)
			bmeta.denom = bmeta.denom + (1/aSe^2)
		}
		
		bmeta = bmeta.enum / bmeta.denom
		
		### Het
		
		qhet = 0
		
		for(i in 1:nDf) {
			### Het
			iBeta = match(objCP@acolBETAs[i], objGWA@aHeader)
			iSe = match(objCP@acolSEs[i], objGWA@aHeader)
			
			aBeta 	<- objGWA@tblGWA[,iBeta]
			aSe 	<- objGWA@tblGWA[,iSe]
		
			qhet = qhet + ((aBeta - bmeta)^2)/(aSe^2)
		}
		
		phet = pchisq(qhet,(nDf-1),lower.tail=FALSE)

		objGWA <- GWADATA.cbind(objGWA, phet, objCP@colOutPhet)
	} else {
		
		nDf = length(objCP@acolZSCOREs)
		
		### Meta
		enum <- denom <- nmeta <- 0 
		
		
		for(i in 1:nDf) {
			### Meta
			
			iZ = match(objCP@acolZSCOREs[i], objGWA@aHeader)
			iN = match(objCP@acolNs[i], objGWA@aHeader)

			if(is.na(iZ)) stop(paste("EASY ERROR:CALCPHET\nColumnZ \n",objCP@acolZSCOREs[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			if(is.na(iN))  stop(paste("EASY ERROR:CALCPHET\nColumnN \n",objCP@acolNs[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			
			aZ 	<- objGWA@tblGWA[,iZ]
			aN 	<- objGWA@tblGWA[,iN]
			
			if(class(aZ) != "numeric" & class(aZ) != "double" ) stop(paste("EASY ERROR:CALCPHET\nColumn \n",objCP@acolZSCOREs[i],"\nis not a numeric variable in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			if(class(aN) != "numeric" & class(aN) != "double" & class(aN) != "integer" ) stop(paste("EASY ERROR:CALCPHET\nColumn \n",objCP@acolNs[i],"\nis not a numeric variable in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			
			enum = enum + aZ*sqrt(aN)
			denom = denom + sqrt(aN)^2
			nmeta = nmeta + aN
		}
		
		zmeta = enum / sqrt(denom)
		
		### Het
		
		qhet = 0
		
		for(i in 1:nDf) {
			### Het
			iZ = match(objCP@acolZSCOREs[i], objGWA@aHeader)
			iN = match(objCP@acolNs[i], objGWA@aHeader)
			
			aZ 	<- objGWA@tblGWA[,iZ]
			aN 	<- objGWA@tblGWA[,iN]
		
			qhet = qhet + ((aZ/sqrt(aN) - zmeta/sqrt(nmeta))^2)*aN
		}
		
		phet = pchisq(qhet,(nDf-1),lower.tail=FALSE)

		objGWA <- GWADATA.cbind(objGWA, phet, objCP@colOutPhet)
	
	
	}
	
	return(objGWA)
}

CALCPHET <- function(strEqcCommand){ 
	## Wrapper for class definition
	CALCPHETout <- setCALCPHET(new("CALCPHET", strEqcCommand = strEqcCommand))
	CALCPHETout<-validCALCPHET(CALCPHETout)
	return(CALCPHETout)

}
