setClass("JOINTTEST",
	representation = representation(
						strEqcCommand				=	"character",
						acolBETAs					=	"character",
						acolSEs						=	"character",
						acolZSCOREs					=	"character",
						colOutPjoint				=	"character",
						strMethod					=	"character"
						),
	prototype = prototype(
						strEqcCommand				=	"",
						acolBETAs					=	"",
						acolSEs						=	"",
						acolZSCOREs					=	"",
						colOutPjoint				=	"pjoint",
						strMethod					=	""
						)
)


setGeneric("setJOINTTEST", function(object) standardGeneric("setJOINTTEST"))
setMethod("setJOINTTEST", signature = (object = "JOINTTEST"), function(object) {
	
	#aEqcSlotNamesIn = c("colBeta1", "colSe1", "colBeta2" , "colSe2", "blnUseGoncalo","bln2sided","rcdTestDirection","strPdiffName")
	aEqcSlotNamesIn = c("acolBETAs", "acolSEs", "acolZSCOREs" , "colOutPjoint")
	
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
validJOINTTEST <- function(objJT) {
	
	if(any(objJT@acolBETAs != "") | any(objJT@acolSEs != "")) strMethod = "BETA"
	else if(any(objJT@acolZSCOREs != "")) strMethod = "ZSCORE"
	else stop(paste(" EASY ERROR:JOINTTEST\n Please either set (--acolBETAs AND --acolSEs) OR (--acolZSCOREs AND --acolNs) !!!", sep=""))
	
	if(strMethod == "BETA") {
		if(length(objJT@acolBETAs)<2)
			stop(paste(" EASY ERROR:JOINTTEST\n Please set at least 2 beta columns for acolBETAs (e.g. --acolBETAs beta1;beta2) !!!", sep=""))
		if(length(objJT@acolSEs)<2) 
			stop(paste(" EASY ERROR:JOINTTEST\n Please set at least 2 se columns for acolSEs (e.g. --acolSEs se1;se2) !!!", sep=""))
		if(length(objJT@acolBETAs) != length(objJT@acolSEs))
			stop(paste(" EASY ERROR:JOINTTEST\n Number of defined Beta columns acolBETAs does not match number of defined Se columns acolSEs. \n PLease check!!!", sep=""))
		if(any(objJT@acolBETAs == ""))
			stop(paste(" EASY ERROR:JOINTTEST\n No Effect estimate columns acolBETAs defined.\n Please set acolBETAs.", sep=""))
		if(any(objJT@acolSEs == ""))
			stop(paste(" EASY ERROR:JOINTTEST\n No Effect estimate columns acolSEs defined.\n Please set acolSEs.", sep=""))
	} else {
		if(length(objJT@acolZSCOREs)<2) 
			stop(paste(" EASY ERROR:JOINTTEST\n Please set at least 2 Z score columns for acolZSCOREs (e.g. --acolZSCOREs z1;z2) !!!", sep=""))
		if(any(objJT@acolZSCOREs == ""))
			stop(paste(" EASY ERROR:JOINTTEST\n No Z score columns acolZSCOREs defined.\n Please set acolZSCOREs.", sep=""))
	}
	
	objJT@strMethod <- strMethod
	
	return(objJT)
}

#############################################################################################################################
JOINTTEST.run <- function(objJT, objGWA) {
	
	if(objJT@strMethod == "BETA") {
		nDf = length(objJT@acolBETAs)
		
		joint.chisq = 0 
		
		for(i in 1:nDf) {
			
			iBeta = match(objJT@acolBETAs[i], objGWA@aHeader)
			iSe = match(objJT@acolSEs[i], objGWA@aHeader)
			
			if(is.na(iBeta)) stop(paste("EASY ERROR:JOINTTEST\nColumn \n",objJT@acolBETAs[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			if(is.na(iSe))  stop(paste("EASY ERROR:JOINTTEST\nColumn \n",objJT@acolSEs[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			
			aBeta 	<- objGWA@tblGWA[,iBeta]
			aSe 	<- objGWA@tblGWA[,iSe]
			
			if(class(aBeta) != "numeric" & class(aBeta) != "double" ) stop(paste("EASY ERROR:JOINTTEST\nColumn \n",objJT@acolBETAs[i],"\nis not a numeric variable in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			if(class(aSe) != "numeric" & class(aSe) != "double" ) stop(paste("EASY ERROR:JOINTTEST\nColumn \n",objJT@acolSEs[i],"\nis not a numeric variable in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			
			joint.chisq = joint.chisq + ((aBeta/aSe)^2)
		}
		
		pjoint = pchisq(joint.chisq, df = nDf, lower.tail = FALSE)
		
		objGWA <- GWADATA.cbind(objGWA, pjoint, objJT@colOutPjoint)
	} else {
		nDf = length(objJT@acolZSCOREs)
		
		joint.chisq = 0 
		
		for(i in 1:nDf) {
			
			iZ = match(objJT@acolZSCOREs[i], objGWA@aHeader)
			
			if(is.na(iZ)) stop(paste("EASY ERROR:JOINTTEST\nColumn \n",objJT@acolZSCOREs[i],"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			
			aZ 	<- objGWA@tblGWA[,iZ]
			
			if(class(aZ) != "numeric" & class(aZ) != "double" ) stop(paste("EASY ERROR:JOINTTEST\nColumn \n",objJT@aZ[i],"\nis not a numeric variable in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
			
			joint.chisq = joint.chisq + aZ^2
		}
		
		pjoint = pchisq(joint.chisq, df = nDf, lower.tail = FALSE)
		
		objGWA <- GWADATA.cbind(objGWA, pjoint, objJT@colOutPjoint)
		
	}
	
	return(objGWA)
}

JOINTTEST <- function(strEqcCommand){ 
	## Wrapper for class definition
	JOINTTESTout <- setJOINTTEST(new("JOINTTEST", strEqcCommand = strEqcCommand))
	JOINTTESTout<-validJOINTTEST(JOINTTESTout)
	return(JOINTTESTout)

}
