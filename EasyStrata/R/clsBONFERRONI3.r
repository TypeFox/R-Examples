setClass("BONFERRONI",
	representation = representation(
						strEqcCommand			=	"character",
						colPval					=	"character",
						numIndepTest			=	"numeric",
						blnUseLengthPval		=	"logical",
						colOut					=	"character",
						colaTopHit				=	"character",
						blnEstimateDim			=	"logical",
						numEstimateDimPosLim	=	"numeric",
						numEstimateDimR2Lim		=	"numeric",
						colInMarker				=	"character",
						filePLINK				=	"character",
						fileBfile				=	"character"
						),
	prototype = prototype(
						strEqcCommand			=	"",
						colPval					=	"",
						numIndepTest			=	1000000,
						blnUseLengthPval		=	FALSE,
						colOut					=	"",
						colaTopHit					=	"",
						blnEstimateDim			=	FALSE,
						numEstimateDimPosLim	=	500000,
						numEstimateDimR2Lim		=	0.5,
						colInMarker				=	"",
						filePLINK				=	"",
						fileBfile				=	""
						)
)


setGeneric("setBONFERRONI", function(object) standardGeneric("setBONFERRONI"))
setMethod("setBONFERRONI", signature = (object = "BONFERRONI"), function(object) {
	
	#aEqcSlotNamesIn = c("colBeta1", "colSe1", "colBeta2" , "colSe2", "blnUseGoncalo","bln2sided","rcdTestDirection","strPdiffName")
	aEqcSlotNamesIn = c("colPval", "numIndepTest", "blnUseLengthPval","colOut", "colaTopHit", "blnEstimateDim","numEstimateDimPosLim", "numEstimateDimR2Lim", "colInMarker", "filePLINK", "fileBfile")
	
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
validBONFERRONI <- function(objBONF) {
	
	
	if(objBONF@colPval == "")
		stop(paste(" EASY ERROR:BONFERRONI\n No P-Value column defined.\n Please set colPval.", sep=""))
	
	if(objBONF@blnEstimateDim) {
		if(objBONF@colInMarker == "") 
			stop(paste(" EASY ERROR:BONFERRONI\n No column colInMarker defined. Please set colInMarker.", sep=""))
		if(objBONF@filePLINK == "") 
			stop(paste(" EASY ERROR:BONFERRONI\n No path to PLINK software --filePLINK defined. Please set filePLINK.", sep=""))
		if(objBONF@fileBfile == "") 
			stop(paste(" EASY ERROR:BONFERRONI\n No path to reference file --fileBfile defined. Please set fileBfile.", sep=""))		
		if(!file.exists(objBONF@filePLINK))
			stop(paste("EASY ERROR:BONFERRONI\n File filePLINK\n ",objBONF@filePLINK,"\n does not exist.", sep=""))
		if(!file.exists(paste(objBONF@fileBfile,".bed",sep="")))
			stop(paste("EASY ERROR:BONFERRONI\n File fileBfile\n ",objBONF@fileBfile,"\n does not exist.", sep=""))
	}
	
	return(TRUE)
}

#############################################################################################################################
BONFERRONI.run <- function(objBONF, objGWA, objREPORT,isValidScript) {
	
	colPval				<- objBONF@colPval
	numIndepTest		<- objBONF@numIndepTest
	blnUseLengthPval	<- objBONF@blnUseLengthPval
	blnEstimateDim		<- objBONF@blnEstimateDim
	colOut				<- objBONF@colOut
	colaTopHit			<- objBONF@colaTopHit
	
	iP = match(colPval, objGWA@aHeader)
	if(is.na(iP)) 
		stop(paste("EASY ERROR:BONFERRONI\nColumn \n",colPval,"\nis not available in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
	
	aP 	<- objGWA@tblGWA[,iP]
	if(class(aP) != "numeric" & class(aP) != "double" ) 
		stop(paste("EASY ERROR:BONFERRONI\nColumn \n",colPval,"\nis not a numeric variable in file \n",objGWA@fileIn,"\n!!!" ,sep="" ))
	## remove NA????
	
	if(blnUseLengthPval) {
	
		nTests <- length(which(!is.na(aP)))
		
	} else if(blnEstimateDim) {
		if(isValidScript) {
			### start plink to estimate number of independent tests
			colInMarker 	<- objBONF@colInMarker
			numPosLim		<- objBONF@numEstimateDimPosLim
			numR2Lim 		<- objBONF@numEstimateDimR2Lim
			filePLINK 		<- objBONF@filePLINK
			fileBfile 		<- objBONF@fileBfile
			
			objGWA.clumpin <- objGWA
			is_na = is.na(GWADATA.getcol(objGWA.clumpin, objBONF@colPval))
# print(length(is_na))
# print(length(which(is_na)))
# print(objGWA.clumpin@tblGWA[1:10,])
			objGWA.clumpin@tblGWA <- objGWA.clumpin@tblGWA[!is_na,]
			objGWA.clumpin <- GWADATA.renamecol(objGWA.clumpin, colInMarker, "SNP")
			
			fileBasePlink = paste(objGWA.clumpin@pathOut,"/",objGWA@fileInShortName,sep="")
			fileBasePlink = sub(".txt","",fileBasePlink)
			fileBasePlink = sub(".gz","",fileBasePlink)
			# if(nchar(strTag)>0) fileBasePlink = paste(fileBasePlink,".",strTag,sep="")
							
			fileInPlink	<- paste(fileBasePlink ,".tmp_plink_in",sep="")
			write.table(objGWA.clumpin@tblGWA, fileInPlink, row.names=F, quote=F, sep="\t")		
			fileOutPlink<- paste(fileBasePlink ,".tmp_plink_out",sep="")

			#### PLINK clump
			## bp to kb
			numPosLim <- round(numPosLim/1000)
			
			plink_Call = paste(filePLINK,"--noweb","--bfile",fileBfile,"--clump",fileInPlink,"--out",fileOutPlink,"--clump-kb",numPosLim, "--clump-r2",numR2Lim,"--clump-p1 1 --clump-p2 1 --clump-field",colPval,"--clump-verbose")
			system(plink_Call)
			
			## creates three files fileOutPlink.clumped, fileOutPlink.log, fileOutPlink.nosex

			vDir<-scan(file = paste(fileOutPlink,".clumped",sep=""), what=character(0), n = -1, sep = "\n",quiet=TRUE,blank.lines.skip=FALSE)
			
			strSplitClump = "------------------------------------------------------------------"

			iClumpChange = grep(strSplitClump, vDir)

			numClumps = length(iClumpChange)
			
			nTests <- numClumps
			
			if(file.exists(fileInPlink)) 
				file.remove(fileInPlink)
			if(file.exists(paste(fileOutPlink,".clumped",sep=""))) 
				file.remove(paste(fileOutPlink,".clumped",sep=""))
			if(file.exists(paste(fileOutPlink,".log",sep=""))) 
				file.remove(paste(fileOutPlink,".log",sep=""))
			if(file.exists(paste(fileOutPlink,".nosex",sep=""))) 
				file.remove(paste(fileOutPlink,".nosex",sep=""))
		
		} else {
			## just for validity check
			nTests <- numIndepTest
		}	
	} else if(colaTopHit != "") {
		
		nTests <- length(which(!is.na(objGWA@tblGWA[,colaTopHit])))
		
	} else {
		nTests <- numIndepTest
	}
	
	#pbonf <- p.adjust(aP,objBONF@strFdrMethod)
	pbonf <- aP*nTests
	pbonf[pbonf>1]<-1
	
	strNewColName <- ifelse(colOut == "", paste(colPval,".bonf",sep=""), colOut)
	objREPORT <- REPORT.addval(objREPORT,paste(strNewColName,".numDim",sep=""),nTests)
	
	objGWA <- GWADATA.cbind(objGWA, pbonf, strNewColName)
	
	return(list(objGWA,objREPORT))
}

BONFERRONI <- function(strEqcCommand){ 
	## Wrapper for class definition
	BONFERRONIout <- setBONFERRONI(new("BONFERRONI", strEqcCommand = strEqcCommand))
	validBONFERRONI(BONFERRONIout)
	return(BONFERRONIout)

}
