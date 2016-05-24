setClass("MIAMIPLOT",
	representation = representation(
						strEqcCommand			=	"character",
						colInChr				=	"character",
						colInPos				=	"character",
						colMIAMIPlotUp			=	"character",
						colMIAMIPlotDown		=	"character",
# global
						astrDefaultColourChr	=	"character",
						arcdColourCrit			=	"character",
						astrColour				=	"character",
						numDefaultSymbol		=	"numeric",
						numDefaultCex			=	"numeric",
						arcdSymbolCrit			=	"character",
						anumSymbol				=	"numeric",
						arcdCexCrit				=	"character",
						anumCex					=	"numeric",
#plot params (do not change btw up/down)			
						fileAnnot				=	"character",
						numAnnotPosLim			=	"numeric",
						numAnnotPvalLim			=	"numeric",
						arcdAdd2Plot			=	"character",
						anumAddPvalLine			=	"numeric",
						astrAddPvalLineCol		=	"character",
						anumAddPvalLineLty		=	"numeric",
						blnYAxisBreak			=	"logical",
						numYAxisBreak			=	"numeric",
						numPvalOffset			=	"numeric",
						strPlotName				= 	"character",
						strFormat				=	"character",
						numCexAxis				=	"numeric",
						numCexLab				=	"numeric",
						numWidth				=	"numeric",
						numHeight				=	"numeric",
						anumParMar				=	"numeric",
						anumParMgp				=	"numeric",
						strParBty				=	"character",
						blnLogPval				=	"logical",
#Up						
						astrDefaultColourChrUp	=	"character",
						arcdColourCritUp		=	"character",
						astrColourUp			=	"character",
						numDefaultSymbolUp		=	"numeric",
						numDefaultCexUp			=	"numeric",
						arcdSymbolCritUp		=	"character",
						anumSymbolUp			=	"numeric",
						arcdCexCritUp			=	"character",
						anumCexUp				=	"numeric",
#Down
						astrDefaultColourChrDown=	"character",
						arcdColourCritDown		=	"character",
						astrColourDown			=	"character",
						numDefaultSymbolDown	=	"numeric",
						numDefaultCexDown		=	"numeric",
						arcdSymbolCritDown		=	"character",
						anumSymbolDown			=	"numeric",
						arcdCexCritDown			=	"character",
						anumCexDown				=	"numeric"
						),
	prototype = prototype(
						strEqcCommand			=	"",
						colInChr				=	"",
						colInPos				=	"",
						colMIAMIPlotUp			=	"",
						colMIAMIPlotDown		=	"",
						astrDefaultColourChr	= 	c("gray51","gray66"),
						arcdColourCrit			=	"",
						astrColour				=	"",
						numDefaultSymbol		=	19,
						numDefaultCex			=	0.4,
						arcdSymbolCrit			=	"",
						anumSymbol				=	-1,
						arcdCexCrit				=	"",
						anumCex					=	-1,
						fileAnnot				=	"",
						numAnnotPosLim			=	500000,
						numAnnotPvalLim			=	1,
						arcdAdd2Plot			=	"",
						anumAddPvalLine			=	5e-8,
						astrAddPvalLineCol		=	"red",
						anumAddPvalLineLty		=	6,
						blnYAxisBreak			=	FALSE,
						numYAxisBreak			=	22,
						numPvalOffset			=	1,
						strPlotName				= 	"miami",
						strFormat				=	"png",
						numCexAxis				=	1.5,
						numCexLab				=	2,
						numWidth				=	1600,	# pixel px # needs change for pdf!
						numHeight				=	800,	# pixel px # needs change for pdf!
						anumParMar				=	c(7, 6, 4, 10),
						anumParMgp				=	c(3, 1, 0),
						strParBty				=	"o",
						blnLogPval				=	FALSE,
						astrDefaultColourChrUp	= 	"",
						arcdColourCritUp		=	"",
						astrColourUp			=	"",
						numDefaultSymbolUp		=	-1,
						numDefaultCexUp			=	-1,
						arcdSymbolCritUp		=	"",
						anumSymbolUp			=	-1,
						arcdCexCritUp			=	"",
						anumCexUp				=	-1,
						astrDefaultColourChrDown= 	"",
						arcdColourCritDown		=	"",
						astrColourDown			=	"",
						numDefaultSymbolDown	=	-1,
						numDefaultCexDown		=	-1,
						arcdSymbolCritDown		=	"",
						anumSymbolDown			=	-1,
						arcdCexCritDown			=	"",
						anumCexDown				=	-1
						)
	#contains = c("EcfReader")
)
setGeneric("setMIAMIPLOT", function(object) standardGeneric("setMIAMIPLOT"))
setMethod("setMIAMIPLOT", signature = (object = "MIAMIPLOT"), function(object) {
	
	aEqcSlotNamesIn = c("colInChr","colInPos",
						"colMIAMIPlotUp","colMIAMIPlotDown",
						"astrDefaultColourChr","arcdColourCrit","astrColour","numDefaultSymbol","numDefaultCex","arcdSymbolCrit","anumSymbol","arcdCexCrit","anumCex",
						"fileAnnot","numAnnotPosLim","numAnnotPvalLim",
						"numPvalOffset",
						"arcdAdd2Plot","anumAddPvalLine","astrAddPvalLineCol","anumAddPvalLineLty",
						"blnYAxisBreak","numYAxisBreak","strPlotName",
						"strFormat","numCexAxis","numCexLab","numWidth","numHeight","anumParMar","anumParMgp","strParBty",
						"blnLogPval",
						"astrDefaultColourChrUp","arcdColourCritUp","astrColourUp","numDefaultSymbolUp","numDefaultCexUp","arcdSymbolCritUp","anumSymbolUp","arcdCexCritUp","anumCexUp",
						"astrDefaultColourChrDown","arcdColourCritDown","astrColourDown","numDefaultSymbolDown","numDefaultCexDown","arcdSymbolCritDown","anumSymbolDown","arcdCexCritDown","anumCexDown")
	
	objEqcReader <- EqcReader(object@strEqcCommand,aEqcSlotNamesIn)
	
	if(length(objEqcReader@lsEqcSlotsOut) > 0) {
		for(i in 1:length(objEqcReader@lsEqcSlotsOut)) {
			tmpSlot <- names(objEqcReader@lsEqcSlotsOut)[i]
			tmpSlotVal <- objEqcReader@lsEqcSlotsOut[[i]]
			
			if(all(!is.na(tmpSlotVal))) slot(object, tmpSlot) <- tmpSlotVal
		}
	}
	isDefaultPngSize = object@numWidth == 1600 & object@numHeight == 1200
	if(object@strFormat == "pdf" & isDefaultPngSize) {
		object@numWidth <- 12
		object@numHeight <- 10
	}
	
	### Recompile Up/Down params
	if(all(object@astrDefaultColourChrDown == "")) 	object@astrDefaultColourChrDown <- object@astrDefaultColourChr
	if(all(object@arcdColourCritDown == "")) 		object@arcdColourCritDown 		<- object@arcdColourCrit
	if(all(object@astrColourDown == "")) 			object@astrColourDown 			<- object@astrColour
	if(object@numDefaultSymbolDown==-1) 			object@numDefaultSymbolDown 	<- object@numDefaultSymbol
	if(all(object@arcdSymbolCritDown == "")) 		object@arcdSymbolCritDown 		<- object@arcdSymbolCrit
	if(all(object@anumSymbolDown == -1)) 			object@anumSymbolDown 			<- object@anumSymbol	
	if(object@numDefaultCexDown==-1) 				object@numDefaultCexDown 		<- object@numDefaultCex
	if(all(object@arcdCexCritDown == "")) 			object@arcdCexCritDown 			<- object@arcdCexCrit
	if(all(object@anumCexDown == -1)) 				object@anumCexDown 				<- object@anumCex	
	
	if(all(object@astrDefaultColourChrUp == "")) 	object@astrDefaultColourChrUp 	<- object@astrDefaultColourChr
	if(all(object@arcdColourCritUp == "")) 			object@arcdColourCritUp 		<- object@arcdColourCrit
	if(all(object@astrColourUp == "")) 				object@astrColourUp 			<- object@astrColour
	if(object@numDefaultSymbolUp==-1) 				object@numDefaultSymbolUp 		<- object@numDefaultSymbol
	if(all(object@arcdSymbolCritUp == "")) 			object@arcdSymbolCritUp 		<- object@arcdSymbolCrit
	if(all(object@anumSymbolUp == -1)) 				object@anumSymbolUp 			<- object@anumSymbol	
	if(object@numDefaultCexUp==-1) 					object@numDefaultCexUp 			<- object@numDefaultCex
	if(all(object@arcdCexCritUp == "")) 			object@arcdCexCritUp 			<- object@arcdCexCrit
	if(all(object@anumCexUp == -1)) 				object@anumCexUp 				<- object@anumCex	
	
	return(object)
})
#############################################################################################################################
MIAMIPLOT.GWADATA.valid <- function(objMIAMIPLOT, objGWA) {
	
	isMatch = objMIAMIPLOT@colMIAMIPlotUp %in% objGWA@aHeader
	if(!isMatch)
		stop(paste("EASY ERROR:MIAMIPLOT\n Column \n",objMIAMIPLOT@colMIAMIPlotUp," does not exist in file\n",objGWA@fileInShortName,"\n !!!" ,sep=""))
	
	isMatch = objMIAMIPLOT@colMIAMIPlotDown %in% objGWA@aHeader
	if(!isMatch)
		stop(paste("EASY ERROR:MIAMIPLOT\n Column \n",objMIAMIPLOT@colMIAMIPlotDown," does not exist in file\n",objGWA@fileInShortName,"\n !!!" ,sep=""))
	
	
	isMatch = objMIAMIPLOT@colInChr %in% objGWA@aHeader
	if(!isMatch)
		stop(paste("EASY ERROR:MIAMIPLOT\n Column \n",objMIAMIPLOT@colInChr," does not exist in file\n",objGWA@fileInShortName,"\n !!!" ,sep=""))
	
	isMatch = objMIAMIPLOT@colInPos %in% objGWA@aHeader
	if(!isMatch)
		stop(paste("EASY ERROR:MIAMIPLOT\n Column \n",objMIAMIPLOT@colInPos," does not exist in file\n",objGWA@fileInShortName,"\n !!!" ,sep=""))
	
	if(objMIAMIPLOT@numPvalOffset<=0 | objMIAMIPLOT@numPvalOffset>1)
		stop(paste("EASY ERROR:MIAMIPLOT\n numPvalOffset must be within ]0;1] \n !!!" ,sep=""))
	
	if(any(objMIAMIPLOT@anumAddPvalLine<=0 | objMIAMIPLOT@anumAddPvalLine>=1)) 
		stop(paste("EASY ERROR:MIAMIPLOT\n Each value in anumAddPvalLine must be within ]0;1[ \n !!!" ,sep=""))
	
	if(length(objMIAMIPLOT@anumAddPvalLine) != length(objMIAMIPLOT@astrAddPvalLineCol)) 
		objMIAMIPLOT@astrAddPvalLineCol = rep("red",length(objMIAMIPLOT@anumAddPvalLine))
		
	if(length(objMIAMIPLOT@anumAddPvalLine) != length(objMIAMIPLOT@anumAddPvalLineLty)) 
		objMIAMIPLOT@anumAddPvalLineLty = rep(6,length(objMIAMIPLOT@anumAddPvalLine))

	### Default symbols:
	if(!(objMIAMIPLOT@numDefaultSymbol%in%c(0:25)))
		stop(paste("EASY ERROR:MIAMIPLOT\n Wrong numDefaultSymbol defined. PLease use integer value between 0 and 25.", sep=""))
	if(!(objMIAMIPLOT@numDefaultSymbolUp%in%c(0:25)))
		stop(paste("EASY ERROR:MIAMIPLOT\n Wrong numDefaultSymbolUp defined. PLease use integer value between 0 and 25.", sep=""))
	if(!(objMIAMIPLOT@numDefaultSymbolDown%in%c(0:25)))
		stop(paste("EASY ERROR:MIAMIPLOT\n Wrong numDefaultSymbolDown defined. PLease use integer value between 0 and 25.", sep=""))
	
	### Default symbol sizes:
	if(objMIAMIPLOT@numDefaultCex <= 0)
		stop(paste("EASY ERROR:MIAMIPLOT\n Wrong numDefaultCex defined. PLease use positive numeric value.", sep=""))
	if(objMIAMIPLOT@numDefaultCexUp <= 0)
		stop(paste("EASY ERROR:MIAMIPLOT\n Wrong numDefaultCexUp defined. PLease use positive numeric value.", sep=""))
	if(objMIAMIPLOT@numDefaultCexDown <= 0)
		stop(paste("EASY ERROR:MIAMIPLOT\n Wrong numDefaultCexDown defined. PLease use positive numeric value.", sep=""))
	
	### Default colours:
	if(length(objMIAMIPLOT@astrDefaultColourChr)!=2) 
		stop(paste("EASY ERROR:MIAMIPLOT\n Length of astrDefaultColourChr must be 2. PLease use two colours separated by ';'.", sep=""))
	for(i in 1:2) {
		strColour = objMIAMIPLOT@astrDefaultColourChr[i]
		isOkColors =  strColour %in% colors()
		isOkColorHex =  nchar(strColour)==7 & substring(strColour,1,1)=="#"
		sTemp = substring(strColour,2,7) 
		isOkColorHex = isOkColorHex & all(strsplit(sTemp,"",fixed=TRUE)[[1]]%in%c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"))
		if(!isOkColors & !isOkColorHex)
			stop(paste("EASY ERROR:MIAMIPLOT\n Wrong astrDefaultColourChr defined. PLease use color from R colors() or hexadecimal nomenclature, eg #FFFFFF.", sep=""))
	}
	if(length(objMIAMIPLOT@astrDefaultColourChrUp)!=2) 
		stop(paste("EASY ERROR:MIAMIPLOT\n Length of astrDefaultColourChrUp must be 2. PLease use two colours separated by ';'.", sep=""))
	for(i in 1:2) {
		strColour = objMIAMIPLOT@astrDefaultColourChrUp[i]
		isOkColors =  strColour %in% colors()
		isOkColorHex =  nchar(strColour)==7 & substring(strColour,1,1)=="#"
		sTemp = substring(strColour,2,7) 
		isOkColorHex = isOkColorHex & all(strsplit(sTemp,"",fixed=TRUE)[[1]]%in%c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"))
		if(!isOkColors & !isOkColorHex)
			stop(paste("EASY ERROR:MIAMIPLOT\n Wrong astrDefaultColourChrUp defined. PLease use color from R colors() or hexadecimal nomenclature, eg #FFFFFF.", sep=""))
	}
	if(length(objMIAMIPLOT@astrDefaultColourChrDown)!=2) 
		stop(paste("EASY ERROR:MIAMIPLOT\n Length of astrDefaultColourChrDown must be 2. PLease use two colours separated by ';'.", sep=""))
	for(i in 1:2) {
		strColour = objMIAMIPLOT@astrDefaultColourChrDown[i]
		isOkColors =  strColour %in% colors()
		isOkColorHex =  nchar(strColour)==7 & substring(strColour,1,1)=="#"
		sTemp = substring(strColour,2,7) 
		isOkColorHex = isOkColorHex & all(strsplit(sTemp,"",fixed=TRUE)[[1]]%in%c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"))
		if(!isOkColors & !isOkColorHex)
			stop(paste("EASY ERROR:MIAMIPLOT\n Wrong astrDefaultColourChrDown defined. PLease use color from R colors() or hexadecimal nomenclature, eg #FFFFFF.", sep=""))
	}
	
	### Colours:
	if(objMIAMIPLOT@astrColour[1] != "") {
		for(i in 1:length(objMIAMIPLOT@astrColour)) {
			strColour = objMIAMIPLOT@astrColour[i]
			isOkColors = strColour %in% colors()
			isOkColorHex =  nchar(strColour)==7 & substring(strColour,1,1)=="#"
			sTemp = substring(strColour,2,7) 
			isOkColorHex = isOkColorHex & all(strsplit(sTemp,"",fixed=TRUE)[[1]]%in%c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"))
			if(!isOkColors & !isOkColorHex)
				stop(paste("EASY ERROR:MIAMIPLOT\n Wrong astrColour defined. PLease use color from R colors() or hexadecimal nomenclature, eg #FFFFFF.", sep=""))
		}
	}
	if(objMIAMIPLOT@astrColourUp[1] != "") {
		for(i in 1:length(objMIAMIPLOT@astrColourUp)) {
			strColour = objMIAMIPLOT@astrColourUp[i]
			isOkColors = strColour %in% colors()
			isOkColorHex =  nchar(strColour)==7 & substring(strColour,1,1)=="#"
			sTemp = substring(strColour,2,7) 
			isOkColorHex = isOkColorHex & all(strsplit(sTemp,"",fixed=TRUE)[[1]]%in%c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"))
			if(!isOkColors & !isOkColorHex)
				stop(paste("EASY ERROR:MIAMIPLOT\n Wrong astrColourUp defined. PLease use color from R colors() or hexadecimal nomenclature, eg #FFFFFF.", sep=""))
		}
	}
	if(objMIAMIPLOT@astrColourDown[1] != "") {
		for(i in 1:length(objMIAMIPLOT@astrColourDown)) {
			strColour = objMIAMIPLOT@astrColourDown[i]
			isOkColors = strColour %in% colors()
			isOkColorHex =  nchar(strColour)==7 & substring(strColour,1,1)=="#"
			sTemp = substring(strColour,2,7) 
			isOkColorHex = isOkColorHex & all(strsplit(sTemp,"",fixed=TRUE)[[1]]%in%c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"))
			if(!isOkColors & !isOkColorHex)
				stop(paste("EASY ERROR:MIAMIPLOT\n Wrong astrColourDown defined. PLease use color from R colors() or hexadecimal nomenclature, eg #FFFFFF.", sep=""))
		}
	}
	
	### Symbols:
	if(objMIAMIPLOT@anumSymbol[1] != -1) {
		for(i in 1:length(objMIAMIPLOT@anumSymbol)) {
			numSymbol = objMIAMIPLOT@anumSymbol[i]
			if(!(numSymbol%in%c(0:25)))
				stop(paste("EASY ERROR:MIAMIPLOT\n Wrong anumSymbol defined. PLease use integer value between 0 and 25.", sep=""))
		}
	}
	if(objMIAMIPLOT@anumSymbolUp[1] != -1) {
		for(i in 1:length(objMIAMIPLOT@anumSymbolUp)) {
			numSymbol = objMIAMIPLOT@anumSymbolUp[i]
			if(!(numSymbol%in%c(0:25)))
				stop(paste("EASY ERROR:MIAMIPLOT\n Wrong anumSymbolUp defined. PLease use integer value between 0 and 25.", sep=""))
		}
	}
	if(objMIAMIPLOT@anumSymbolDown[1] != -1) {
		for(i in 1:length(objMIAMIPLOT@anumSymbolDown)) {
			numSymbol = objMIAMIPLOT@anumSymbolDown[i]
			if(!(numSymbol%in%c(0:25)))
				stop(paste("EASY ERROR:MIAMIPLOT\n Wrong anumSymbolDown defined. PLease use integer value between 0 and 25.", sep=""))
		}
	}
	### Symbol sizes:
	if(objMIAMIPLOT@anumCex[1] != -1) {
		for(i in 1:length(objMIAMIPLOT@anumCex)) {
			numCex = objMIAMIPLOT@anumCex[i]
			if(numCex <= 0)
				stop(paste("EASY ERROR:MIAMIPLOT\n Wrong anumCex defined. PLease use positive numeric value.", sep=""))			
		}
	}
	if(objMIAMIPLOT@anumCexUp[1] != -1) {
		for(i in 1:length(objMIAMIPLOT@anumCexUp)) {
			numCex = objMIAMIPLOT@anumCexUp[i]
			if(numCex <= 0)
				stop(paste("EASY ERROR:MIAMIPLOT\n Wrong anumCexUp defined. PLease use positive numeric value.", sep=""))			
		}
	}
	if(objMIAMIPLOT@anumCexDown[1] != -1) {
		for(i in 1:length(objMIAMIPLOT@anumCexDown)) {
			numCex = objMIAMIPLOT@anumCexDown[i]
			if(numCex <= 0)
				stop(paste("EASY ERROR:MIAMIPLOT\n Wrong anumCexDown defined. PLease use positive numeric value.", sep=""))			
		}
	}
	
	if(!(objMIAMIPLOT@strFormat == "png" | objMIAMIPLOT@strFormat == "pdf"))
		stop(paste("EASY ERROR:MIAMIPLOT\n Wrong strFormat defined, PLease use 'png' or 'pdf' \n !!!", sep=""))
	
	if(objMIAMIPLOT@numCexAxis <= 0) 
		stop(paste("EASY ERROR:MIAMIPLOT\n Wrong numCexAxis defined, PLease use numCexAxis > 0 \n !!!", sep=""))
		
	if(objMIAMIPLOT@numCexLab <= 0) 
		stop(paste("EASY ERROR:MIAMIPLOT\n Wrong numCexLab defined, PLease use numCexLab > 0 \n !!!", sep=""))				
	
	if(length(objMIAMIPLOT@anumParMar) != 4) 
		stop(paste("EASY ERROR:MIAMIPLOT\n Wrong anumParMar defined, Please use anumParMar of length 4 (Default=5.1;4.1;4.1;2.1) \n !!!", sep=""))				
	
	if(length(objMIAMIPLOT@anumParMgp) != 3) 
		stop(paste("EASY ERROR:MIAMIPLOT\n Wrong anumParMgp defined, Please use anumParMgp of length 3 (Default=3;1;0) \n !!!", sep=""))				
	
	if(!(objMIAMIPLOT@strParBty%in%c("o","l","7","c","u","]")))
		stop(paste("EASY ERROR:MIAMIPLOT\n Wrong strParBty defined, Please use 'o','l','7','c','u' or ']' \n !!!", sep=""))				
		
	if(objMIAMIPLOT@fileAnnot != "") {
		if(!file.exists(objMIAMIPLOT@fileAnnot))
			stop(paste("EASY ERROR:MIAMIPLOT\n File fileAnnot\n ",objMIAMIPLOT@fileAnnot,"\n does not exist.", sep=""))
		### Cols exist?
			tblAnnot10<-read.table(objMIAMIPLOT@fileAnnot,header=T, sep="",  stringsAsFactors=FALSE, nrows = 10)
			isAv <- "Chr" %in% names(tblAnnot10)
			if(!isAv)
				stop(paste(" EASY ERROR:MIAMIPLOT\n Column 'Chr' is not available in fileAnnot. PLease correct file header.", sep=""))
			isAv <- "Pos" %in% names(tblAnnot10)
			if(!isAv)
				stop(paste(" EASY ERROR:MIAMIPLOT\n Column 'Pos' is not available in fileAnnot. PLease correct file header.", sep=""))
			isAv <- "Colour" %in% names(tblAnnot10)
			if(!isAv)
				stop(paste(" EASY ERROR:MIAMIPLOT\n Column 'Colour' is not available in fileAnnot. PLease correct file header.", sep=""))	
	}
	
}
MIAMIPLOT.run <- function(objMIAMIPLOT, objGWA) {

	colInChr 				<- 	objMIAMIPLOT@colInChr
	colInPos 				<- 	objMIAMIPLOT@colInPos
	colMIAMIPlotUp 			<- 	objMIAMIPLOT@colMIAMIPlotUp
	colMIAMIPlotDown 		<- 	objMIAMIPLOT@colMIAMIPlotDown
	astrDefaultColourChrUp	<- 	objMIAMIPLOT@astrDefaultColourChrUp
	astrDefaultColourChrDown<- 	objMIAMIPLOT@astrDefaultColourChrDown
	numPvalOffset 			<-	objMIAMIPLOT@numPvalOffset
	arcdAdd2Plot			<-	objMIAMIPLOT@arcdAdd2Plot
	## 130829
	anumAddPvalLine			<-	objMIAMIPLOT@anumAddPvalLine
	astrAddPvalLineCol		<-	objMIAMIPLOT@astrAddPvalLineCol
	anumAddPvalLineLty		<-	objMIAMIPLOT@anumAddPvalLineLty
	blnYAxisBreak			<-	objMIAMIPLOT@blnYAxisBreak
	numYAxisBreak			<-	objMIAMIPLOT@numYAxisBreak
	numCexAxis				<-	objMIAMIPLOT@numCexAxis
	numCexLab				<-	objMIAMIPLOT@numCexLab
	blnLogPval				<-	objMIAMIPLOT@blnLogPval	
	## 130927
	arcdColourCritUp		<-	objMIAMIPLOT@arcdColourCritUp
	astrColourUp			<-	objMIAMIPLOT@astrColourUp
	numDefaultSymbolUp		<-	objMIAMIPLOT@numDefaultSymbolUp
	numDefaultCexUp			<-	objMIAMIPLOT@numDefaultCexUp
	arcdSymbolCritUp		<-	objMIAMIPLOT@arcdSymbolCritUp
	anumSymbolUp			<-	objMIAMIPLOT@anumSymbolUp
	arcdCexCritUp			<-	objMIAMIPLOT@arcdCexCritUp
	anumCexUp				<-	objMIAMIPLOT@anumCexUp
	
	arcdColourCritDown		<-	objMIAMIPLOT@arcdColourCritDown
	astrColourDown			<-	objMIAMIPLOT@astrColourDown
	numDefaultSymbolDown	<-	objMIAMIPLOT@numDefaultSymbolDown
	numDefaultCexDown		<-	objMIAMIPLOT@numDefaultCexDown
	arcdSymbolCritDown		<-	objMIAMIPLOT@arcdSymbolCritDown
	anumSymbolDown			<-	objMIAMIPLOT@anumSymbolDown
	arcdCexCritDown			<-	objMIAMIPLOT@arcdCexCritDown
	anumCexDown				<-	objMIAMIPLOT@anumCexDown
	
	## 140211
	fileAnnot 				=	objMIAMIPLOT@fileAnnot
	numAnnotPosLim			=	objMIAMIPLOT@numAnnotPosLim
	numAnnotPvalLim			=	objMIAMIPLOT@numAnnotPvalLim
	############################
	
	iValUp = match(colMIAMIPlotUp, names(objGWA@tblGWA))
	objGWA@tblGWA[,iValUp] = as.numeric(objGWA@tblGWA[,iValUp])
	iValDown = match(colMIAMIPlotDown, names(objGWA@tblGWA))
	objGWA@tblGWA[,iValDown] = as.numeric(objGWA@tblGWA[,iValDown])
	iChr = match(colInChr, names(objGWA@tblGWA))
	objGWA@tblGWA[,iChr] = as.numeric(objGWA@tblGWA[,iChr])
	iPos = match(colInPos, names(objGWA@tblGWA))
	objGWA@tblGWA[,iPos] = as.numeric(objGWA@tblGWA[,iPos])
	
	#isRemove = (is.na(objGWA@tblGWA[,iValUp])&is.na(objGWA@tblGWA[,iValDown])) | is.na(objGWA@tblGWA[,iChr]) | is.na(objGWA@tblGWA[,iPos])
	isRemove = is.na(objGWA@tblGWA[,iChr]) | is.na(objGWA@tblGWA[,iPos])
	if(any(isRemove)) objGWA <- GWADATA.removerows(objGWA, which(isRemove))
	
	############################
	#### Compile initital colour/symbol/cex array 
	
	colplotUp <- plotCritKeyUp <- rep(NA, dim(objGWA@tblGWA)[1])
	symplotUp <- rep(numDefaultSymbolUp, dim(objGWA@tblGWA)[1])
	cexplotUp <- rep(numDefaultCexUp, dim(objGWA@tblGWA)[1])
	
	if(arcdColourCritUp[1] != "") {
		for(iColourCrit in 1:length(arcdColourCritUp)) {
			rcdColourCrit 		<- arcdColourCritUp[iColourCrit]
			strColour 			<- astrColourUp[iColourCrit]
			objRCD 	<- RCD(rcdColourCrit)
			out 	<- RCD.eval(objRCD, objGWA)
			colplotUp[out] <- strColour
			plotCritKeyUp[out]= iColourCrit
		}
	} 
	if(arcdSymbolCritUp[1] != "") {
		for(iSymbolCrit in 1:length(arcdSymbolCritUp)) {
			rcdSymbolCrit 		<- arcdSymbolCritUp[iSymbolCrit]
			numSymbol 			<- anumSymbolUp[iSymbolCrit]
			objRCD 	<- RCD(rcdSymbolCrit)
			out 	<- RCD.eval(objRCD, objGWA)
			symplotUp[out] <- numSymbol
		}
	}
	if(arcdCexCritUp[1] != "") {
		for(iCexCrit in 1:length(arcdCexCritUp)) {
			rcdCexCrit 		<- arcdCexCritUp[iCexCrit]
			numCex 			<- anumCexUp[iCexCrit]
			objRCD 	<- RCD(rcdCexCrit)
			out 	<- RCD.eval(objRCD, objGWA)
			cexplotUp[out] <- numCex
		}
	}
	
	colplotDown <- plotCritKeyDown <- rep(NA, dim(objGWA@tblGWA)[1])
	symplotDown <- rep(numDefaultSymbolDown, dim(objGWA@tblGWA)[1])
	cexplotDown <- rep(numDefaultCexDown, dim(objGWA@tblGWA)[1])
	
	if(arcdColourCritDown[1] != "") {
		for(iColourCrit in 1:length(arcdColourCritDown)) {
			rcdColourCrit 		<- arcdColourCritDown[iColourCrit]
			strColour 			<- astrColourDown[iColourCrit]
			objRCD 	<- RCD(rcdColourCrit)
			out 	<- RCD.eval(objRCD, objGWA)
			colplotDown[out] <- strColour
			plotCritKeyDown[out]= iColourCrit
		}
	} 
	if(arcdSymbolCritDown[1] != "") {
		for(iSymbolCrit in 1:length(arcdSymbolCritDown)) {
			rcdSymbolCrit 		<- arcdSymbolCritDown[iSymbolCrit]
			numSymbol 			<- anumSymbolDown[iSymbolCrit]
			objRCD 	<- RCD(rcdSymbolCrit)
			out 	<- RCD.eval(objRCD, objGWA)
			symplotDown[out] <- numSymbol
		}
	}
	if(arcdCexCritDown[1] != "") {
		for(iCexCrit in 1:length(arcdCexCritDown)) {
			rcdCexCrit 		<- arcdCexCritDown[iCexCrit]
			numCex 			<- anumCexDown[iCexCrit]
			objRCD 	<- RCD(rcdCexCrit)
			out 	<- RCD.eval(objRCD, objGWA)
			cexplotDown[out] <- numCex
		}
	}
	
	############################
	#### Log on P-Values
	if(blnLogPval) {
		isInf = objGWA@tblGWA[,iValUp] == Inf
		if(any(isInf)) 
			warning(paste("EASY WARNING:\n There are ",length(which(isInf)), " SNPs with -log10_P=Inf being excluded from the MH plot for file \n ", objGWA@fileIn ,sep="" ))
		isOutRange = objGWA@tblGWA[,iValUp] < 0
		if(any(isOutRange))
			warning(paste("EASY WARNING:\n There are ",length(which(isOutRange)), " SNPs with -log10_P<0 being excluded from the MH plot for file \n ", objGWA@fileIn ,sep="" ))
		isBelowOffset = objGWA@tblGWA[,iValUp] < -log10(numPvalOffset)
		isRemove = isBelowOffset | isOutRange | isInf
		if(any(isRemove)) {
			objGWA@tblGWA[which(isRemove),iValUp] <- NA
		}
		##
		isInf = objGWA@tblGWA[,iValDown] == Inf
		if(any(isInf)) 
			warning(paste("EASY WARNING:\n There are ",length(which(isInf)), " SNPs with -log10_P=Inf being excluded from the MH plot for file \n ", objGWA@fileIn ,sep="" ))
		isOutRange = objGWA@tblGWA[,iValDown] < 0
		if(any(isOutRange))
			warning(paste("EASY WARNING:\n There are ",length(which(isOutRange)), " SNPs with -log10_P<0 being excluded from the MH plot for file \n ", objGWA@fileIn ,sep="" ))
		isBelowOffset = objGWA@tblGWA[,iValDown] < -log10(numPvalOffset)
		isRemove = isBelowOffset | isOutRange | isInf
		if(any(isRemove)) {
			objGWA@tblGWA[which(isRemove),iValDown] <- NA
		}
	} else {
		isZero = objGWA@tblGWA[,iValUp] == 0
		isZero[is.na(isZero)]<-FALSE

		if(any(isZero)) 
			warning(paste("EASY WARNING:\n There are ",length(which(isZero)), " SNPs with P=0 being excluded from the MH plot for file \n ", objGWA@fileIn ,sep="" ))
			
		isOutRange = objGWA@tblGWA[,iValUp] < 0 | objGWA@tblGWA[,iValUp] > 1
		isOutRange[is.na(isOutRange)]<-FALSE
		if(any(isOutRange))
			warning(paste("EASY WARNING:\n There are ",length(which(isOutRange)), " SNPs with P<0 or P>1 being excluded from the MH plot for file \n ", objGWA@fileIn ,sep="" ))
			
		isBelowOffset = objGWA@tblGWA[,iValUp] > numPvalOffset
		isBelowOffset[is.na(isBelowOffset)]<-FALSE
		
		isRemove = isBelowOffset | isOutRange | isZero
		if(any(isRemove)) 
			objGWA@tblGWA[which(isRemove),iValUp] <- NA
		
		objGWA@tblGWA[,iValUp] 		= -log10(objGWA@tblGWA[,iValUp])
		##
		isZero = objGWA@tblGWA[,iValDown] == 0
		isZero[is.na(isZero)]<-FALSE
		
		if(any(isZero)) 
			warning(paste("EASY WARNING:\n There are ",length(which(isZero)), " SNPs with P=0 being excluded from the MH plot for file \n ", objGWA@fileIn ,sep="" ))
		isOutRange = objGWA@tblGWA[,iValDown] < 0 | objGWA@tblGWA[,iValDown] > 1
		isOutRange[is.na(isOutRange)]<-FALSE
		
		if(any(isOutRange))
			warning(paste("EASY WARNING:\n There are ",length(which(isOutRange)), " SNPs with P<0 or P>1 being excluded from the MH plot for file \n ", objGWA@fileIn ,sep="" ))
		isBelowOffset = objGWA@tblGWA[,iValDown] > numPvalOffset
		isBelowOffset[is.na(isBelowOffset)]<-FALSE
		
		isRemove = isBelowOffset | isOutRange | isZero
		if(any(isRemove)) {
			objGWA@tblGWA[which(isRemove),iValDown] <- NA
		}
		objGWA@tblGWA[,iValDown] 		= -log10(objGWA@tblGWA[,iValDown])
	}
	
	objGWA <- GWADATA.getcols(objGWA, c(colMIAMIPlotUp,colMIAMIPlotDown, colInChr, colInPos))
	tblPlot <- objGWA@tblGWA
	names(tblPlot)[1] = "yUp"
	names(tblPlot)[2] = "yDown"
	names(tblPlot)[3] = "chr"
	names(tblPlot)[4] = "pos"

	tblPlot = cbind(tblPlot,colplotUp,colplotDown,plotCritKeyUp,plotCritKeyDown,symplotUp,symplotDown,cexplotUp,cexplotDown, stringsAsFactors = FALSE)

	### Annot
	if(fileAnnot != "") {
		tblAnnot<-read.table(fileAnnot,header=T, sep="",  stringsAsFactors=FALSE)
		isSignifUp <- tblPlot$yUp > -log10(numAnnotPvalLim)
		isSignifUp[is.na(isSignifUp)] <- FALSE
		isSignifDown <- tblPlot$yDown > -log10(numAnnotPvalLim)
		isSignifDown[is.na(isSignifDown)] <- FALSE
		for(iLocus in 1:nrow(tblAnnot)) {
			isLocus = (tblPlot$chr == tblAnnot$Chr[iLocus]) & (abs(tblPlot$pos-tblAnnot$Pos[iLocus])<numAnnotPosLim)
			#isColourUp = isLocus & is.na(tblPlot$colplotUp) & isSignifUp
			isColourUp = isLocus & isSignifUp
			if(any(isColourUp)) {
				# [isLocus] so that all SNPs in thelocus get coloured
				tblPlot$colplotUp[isLocus] <- tblAnnot$Colour[iLocus]
				tblPlot$plotCritKeyUp[isLocus] <- Inf
			}
			#isColourDown = isLocus & is.na(tblPlot$colplotDown) & isSignifDown
			isColourDown = isLocus & isSignifDown
			if(any(isColourDown)) {
				tblPlot$colplotDown[isLocus] <- tblAnnot$Colour[iLocus]
				tblPlot$plotCritKeyDown[isLocus] <- Inf
			}
		}
	}
	
	############################
	#### Compile x-axes 
	
	iOrderChrPos = order(tblPlot$chr,tblPlot$pos)
	tblPlot<-tblPlot[iOrderChrPos,] 
	
	##### Calc. position - vector
	arr_pos<-tblPlot$pos
	
	arr_pos_min<-tapply(tblPlot$pos,tblPlot$chr,min)
	arr_pos_max<-tapply(tblPlot$pos,tblPlot$chr,max)
	
	## tapply sortiert Ausgabe nach chr !!!
	
	arr_chr_uni<-sort(unique(tblPlot$chr))
	
	arr_dpos<-arr_pos_max-arr_pos_min+1
	num_chr1<-min(arr_chr_uni)
	
	# Set first appearing chr to minimum position
	arr_pos[tblPlot$chr==num_chr1]<-arr_pos[tblPlot$chr==num_chr1]-arr_pos_min[1]+1
	
	arr_chr2end<-tblPlot$chr[tblPlot$chr>num_chr1]
	arr_chr2end<-unique(arr_chr2end)
	arr_chr2end<-sort(arr_chr2end)

	count=2
	
	pos_lim = arr_labpos = rep(NA,length(arr_chr2end)+1)
	pos_lim[1] = max(arr_pos[tblPlot$chr==num_chr1])
	arr_labpos[1]<- min(arr_pos[tblPlot$chr==num_chr1]) + (max(arr_pos[tblPlot$chr==num_chr1]) - min(arr_pos[tblPlot$chr==num_chr1]))/2

	for(i in arr_chr2end)  #23?
	{
		is_chr = tblPlot$chr==i
		arr_pos[is_chr]<-sum(arr_dpos[1:(count-1)])+arr_pos[is_chr]-arr_pos_min[count]+1
	
		pos_lim[count] = max(arr_pos[is_chr])
		# arr_labpos[count]<-tapply(arr_pos,tblPlot$chr,mean)
		# arr_labname<-paste("chr",arr_chr,sep="")
		arr_labpos[count]<- min(arr_pos[is_chr]) + (max(arr_pos[is_chr]) - min(arr_pos[is_chr]))/2
	
		count=count+1
	}
	
	tblPlot = cbind(tblPlot, "x" = arr_pos, stringsAsFactors = FALSE)
	
	############################
	##### ReCheck Colouring and symbols:
	# create array, with unique, sorted, appearing chromosomes:
	arr_chr<-c(num_chr1,arr_chr2end)
	arr_labname<-paste("chr",arr_chr,sep="")
	arr_idx2=seq(1,length(arr_chr),by=2)
	# e.g. arr_idx2 = 1,3,5, ...
	#### Set default colours:
	isSetDefaultColour = is.na(tblPlot$colplotUp)
	tblPlot$colplotUp[isSetDefaultColour]	<-	ifelse(tblPlot$chr[isSetDefaultColour]%in%arr_chr[arr_idx2],astrDefaultColourChrUp[1],astrDefaultColourChrUp[2])
	isSetDefaultColour = is.na(tblPlot$colplotDown)
	tblPlot$colplotDown[isSetDefaultColour]	<-	ifelse(tblPlot$chr[isSetDefaultColour]%in%arr_chr[arr_idx2],astrDefaultColourChrDown[1],astrDefaultColourChrDown[2])	
	plot.title = objGWA@fileInShortName
	strXlab = "Chromosome"
	if(!blnLogPval) {
		strYlabUp = paste("-log10 (",colMIAMIPlotUp,")",sep="")
		strYlabDown = paste("+log10 (",colMIAMIPlotDown,")",sep="")
	}
	else {
		strYlabUp = colMIAMIPlotUp
		strYlabDown = colMIAMIPlotDown
	}

	############################

	yAxisTickPos <- NULL
	yAxisTickLab <- TRUE
	anumAddPvalLinePos = -log10(anumAddPvalLine)
	
	#### miami axis lim and pos
	yMax = max(pmax(tblPlot$yUp, tblPlot$yDown, na.rm=TRUE), na.rm=TRUE)
	if(blnYAxisBreak & numYAxisBreak<yMax) {
		## transform yUp
			tblPlot = tblPlot[order(tblPlot$yUp),]
			isUsed = !is.na(tblPlot$yUp)
			yUp = tblPlot$yUp[isUsed]
			yu = yUp[yUp<=numYAxisBreak] # [0;22] -> [0;80]
			yo = yUp[yUp>numYAxisBreak] # [22;yMax] -> [80;100]
			yus = yu*80/numYAxisBreak
			yos = 80 - 20*(numYAxisBreak-yo)/(yMax-numYAxisBreak)
			tblPlot$yUp[isUsed] = c(yus,yos)
		## transform yDown
			tblPlot = tblPlot[order(tblPlot$yDown),]
			isUsed = !is.na(tblPlot$yDown)
			yDown = tblPlot$yDown[isUsed]
			yu = yDown[yDown<=numYAxisBreak] # [0;22] -> [0;80]
			yo = yDown[yDown>numYAxisBreak] # [22;yMax] -> [80;100]
			yus = yu*80/numYAxisBreak
			yos = 80 - 20*(numYAxisBreak-yo)/(yMax-numYAxisBreak)
			tblPlot$yDown[isUsed] = c(yus,yos)
		## transform yTicks
		yAxisTickLabU = pretty(c(0,numYAxisBreak),n=4,min.n=0)
		yAxisTickLabU = yAxisTickLabU[-which(yAxisTickLabU>=numYAxisBreak)]
		yAxisTickPosU = yAxisTickLabU*80/numYAxisBreak
		yAxisTickLabO = pretty(c(numYAxisBreak,yMax),n=4,min.n=0)
		yAxisTickLabO = yAxisTickLabO[-which(yAxisTickLabO<=numYAxisBreak)]
		yAxisTickPosO = 80 - 20*(numYAxisBreak-yAxisTickLabO)/(yMax-numYAxisBreak)
		yAxisTickLab = c(yAxisTickLabU,yAxisTickLabO)
		yAxisTickPos = c(yAxisTickPosU,yAxisTickPosO)
		
		yAxisTickLab = c(-rev(yAxisTickLab),yAxisTickLab)
		yAxisTickPos = c(-rev(yAxisTickPos),yAxisTickPos)
		
		isLineU = anumAddPvalLinePos<numYAxisBreak
		anumAddPvalLinePos[isLineU] = anumAddPvalLinePos[isLineU]*80/numYAxisBreak
		anumAddPvalLinePos[!isLineU] = 80 - 20*(numYAxisBreak-anumAddPvalLinePos[!isLineU])/(yMax-numYAxisBreak)
	}
	
	## 
	yMaxScaled = max(pmax(tblPlot$yUp, tblPlot$yDown, na.rm=TRUE), na.rm=TRUE)
	arr_ylim = c(-yMaxScaled,yMaxScaled)
	arr_labpos_y = c((arr_ylim[1]/4)*3,arr_ylim[2]/2)
	############################	############################
	####### Start plot
	
	tblPlot = tblPlot[order(tblPlot$plotCritKeyUp,na.last=FALSE),]
	
	#plot( tblPlot$x, tblPlot$y, col=tblPlot$colplot, pch=tblPlot$symplot, cex=tblPlot$cexplot, cex.lab = numCexLab, xaxt="n",yaxt="n", ylab=strYlab, xlab="",ylim=c(0,max(tblPlot$y)))
	plot( tblPlot$x, tblPlot$yUp, col=tblPlot$colplotUp, pch=tblPlot$symplotUp, cex=tblPlot$cexplotUp, xaxt="n",yaxt="n", ylab="", xlab="",ylim=arr_ylim)
	#tblPlot = tblPlot[order(tblPlot$plotCritKeyDown,na.last=FALSE),]
	points( tblPlot$x, -tblPlot$yDown, col=tblPlot$colplotDown, pch=tblPlot$symplotDown, cex=tblPlot$cexplotDown)
	
	### yAxis:
	mtext(strYlabUp, side=2, at=arr_labpos_y[2], line=2.5, cex=numCexLab)
	mtext(strYlabDown, side=2, at=arr_labpos_y[1], line=2.5, cex=numCexLab)			
	axis(2, at =yAxisTickPos,  labels = yAxisTickLab, cex.axis = numCexAxis)
	### xAxisLab:
	mtext("Chromosome", side=1, line=2.5, cex=numCexLab, padj=2)
	axis(1,at=arr_labpos,labels=arr_labname,las=3,cex.axis = numCexAxis)
	
	abline(a=0,b=0,col="gray3",lty=6)
	for(k in 1:length(pos_lim)) abline(v=pos_lim[k],col="gray3",lty=3)
	grid(nx = NA, ny = NA, col = "gray", lty = "dotted")
	
	for(k in 1:length(anumAddPvalLine)) {
		if(is.na(anumAddPvalLineLty[k])) anumAddPvalLineLty[k] <- 6
		if(is.na(astrAddPvalLineCol[k])) astrAddPvalLineCol[k] <- "red"
	
		abline(anumAddPvalLinePos[k],0,col = astrAddPvalLineCol[k],lty = anumAddPvalLineLty[k])	
		axis(4,at=anumAddPvalLinePos[k],labels=anumAddPvalLine[k],las=3,cex.axis = numCexAxis)
		
		abline(-anumAddPvalLinePos[k],0,col = astrAddPvalLineCol[k],lty = anumAddPvalLineLty[k])	
		axis(4,at=-anumAddPvalLinePos[k],labels=anumAddPvalLine[k],las=3,cex.axis = numCexAxis)
		
	}
	#abline(yGwsLine,0,col = "red",lty = 6)	
	#axis(4,at=yGwsLine,labels="5e-8",las=3,cex.axis = 1.5)
	
	if(blnYAxisBreak & numYAxisBreak<yMax) {
		axis.break(2, 80, style = "zigzag")
		abline(80,0,col = "grey",lty = 1)	
		
		axis.break(2, -80, style = "zigzag")
		abline(-80,0,col = "grey",lty = 1)	
	}
	
	if(all(arcdAdd2Plot != "")) {
		for(rcdTmp in arcdAdd2Plot) {
			# for(colTmp in names(tblPlot)) {
				# if(length(grep(colTmp,rcdTmp)) > 0) {
					# icolTmp=which(names(tblPlot)==colTmp)
					# assign(colTmp,as.numeric(as.character(tblPlot[,icolTmp])))
					# # create arrays, that are used within criterion and name them accordingly
					# # tblIn$pmen -> create array pmen if criterion e.g. pmen>0.3
				# }
			# }
			eval(parse(text=rcdTmp))
		}
	}
	
}
MIAMIPLOT <- function(strEqcCommand){ 
	## Wrapper for class definition
	MIAMIPLOTout <- setMIAMIPLOT(new("MIAMIPLOT", strEqcCommand = strEqcCommand))
	return(MIAMIPLOTout)
}

