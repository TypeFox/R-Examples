setClass("MHPLOT",
	representation = representation(
						strEqcCommand			=	"character",
						colMHPlot				=	"character",
						colInChr				=	"character",
						colInPos				=	"character",
						astrDefaultColourChr	=	"character",
						arcdColourCrit			=	"character",
						astrColour				=	"character",
						numPvalOffset			=	"numeric",
						arcdAdd2Plot			=	"character",
						anumAddPvalLine			=	"numeric",
						astrAddPvalLineCol		=	"character",
						anumAddPvalLineLty		=	"numeric",
						blnYAxisBreak			=	"logical",
						numYAxisBreak			=	"numeric",
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
						numDefaultSymbol		=	"numeric",
						numDefaultCex			=	"numeric",
						arcdSymbolCrit			=	"character",
						anumSymbol				=	"numeric",
						arcdCexCrit				=	"character",
						anumCex					=	"numeric",
						fileAnnot				=	"character",
						numAnnotPosLim			=	"numeric",
						numAnnotPvalLim			=	"numeric"
						),
	prototype = prototype(
						strEqcCommand			=	"",
						colMHPlot				=	"",
						colInChr				=	"",
						colInPos				=	"",
						astrDefaultColourChr	= 	c("gray51","gray66"),
						arcdColourCrit			=	"",
						astrColour				=	"",
						numPvalOffset			=	1,
						arcdAdd2Plot			=	"",
						anumAddPvalLine			=	5e-8,
						astrAddPvalLineCol		=	"red",
						anumAddPvalLineLty		=	6,
						blnYAxisBreak			=	FALSE,
						numYAxisBreak			=	22,
						strPlotName				= 	"mh",
						strFormat				=	"png",
						numCexAxis				=	1.5,
						numCexLab				=	2,
						numWidth				=	1600,	# pixel px # needs change for pdf!
						numHeight				=	600,	# pixel px # needs change for pdf!
						anumParMar				=	c(7, 6, 4, 10),
						anumParMgp				=	c(3, 1, 0),
						strParBty				=	"o",
						blnLogPval				=	FALSE,
						numDefaultSymbol		=	19,
						numDefaultCex			=	0.4,
						arcdSymbolCrit			=	"",
						anumSymbol				=	-1,
						arcdCexCrit				=	"",
						anumCex					=	-1,
						fileAnnot				=	"",
						numAnnotPosLim			=	500000,
						numAnnotPvalLim			=	1
						)
	#contains = c("EcfReader")
)
setGeneric("setMHPLOT", function(object) standardGeneric("setMHPLOT"))
setMethod("setMHPLOT", signature = (object = "MHPLOT"), function(object) {
	
	aEqcSlotNamesIn = c("colMHPlot","colInChr","colInPos",
						"astrDefaultColourChr",
						"arcdColourCrit","astrColour",
						"numPvalOffset",
						"arcdAdd2Plot","anumAddPvalLine","astrAddPvalLineCol","anumAddPvalLineLty",
						"blnYAxisBreak","numYAxisBreak","strPlotName",
						"strFormat","numCexAxis","numCexLab","numWidth","numHeight","anumParMar","anumParMgp","strParBty",
						"blnLogPval",
						"numDefaultSymbol","numDefaultCex","arcdSymbolCrit","anumSymbol","arcdCexCrit","anumCex",
						"fileAnnot","numAnnotPosLim","numAnnotPvalLim")
	
	objEqcReader <- EqcReader(object@strEqcCommand,aEqcSlotNamesIn)
	
	if(length(objEqcReader@lsEqcSlotsOut) > 0) {
		for(i in 1:length(objEqcReader@lsEqcSlotsOut)) {
			tmpSlot <- names(objEqcReader@lsEqcSlotsOut)[i]
			tmpSlotVal <- objEqcReader@lsEqcSlotsOut[[i]]
			
			if(all(!is.na(tmpSlotVal))) slot(object, tmpSlot) <- tmpSlotVal
		}
	}
	isDefaultPngSize = object@numWidth == 1600 & object@numHeight == 600
	if(object@strFormat == "pdf" & isDefaultPngSize) {
		object@numWidth <- 12
		object@numHeight <- 5
	}
	
	return(object)
})
#############################################################################################################################
MHPLOT.GWADATA.valid <- function(objMHPLOT, objGWA) {
	
	isMatch = objMHPLOT@colMHPlot %in% objGWA@aHeader
	if(!isMatch)
		stop(paste("EASY ERROR:MHPLOT\n Column \n",objMHPLOT@colMHPlot," does not exist in file\n",objGWA@fileInShortName,"\n !!!" ,sep=""))
	
	isMatch = objMHPLOT@colInChr %in% objGWA@aHeader
	if(!isMatch)
		stop(paste("EASY ERROR:MHPLOT\n Column \n",objMHPLOT@colInChr," does not exist in file\n",objGWA@fileInShortName,"\n !!!" ,sep=""))
	
	isMatch = objMHPLOT@colInPos %in% objGWA@aHeader
	if(!isMatch)
		stop(paste("EASY ERROR:MHPLOT\n Column \n",objMHPLOT@colInPos," does not exist in file\n",objGWA@fileInShortName,"\n !!!" ,sep=""))
	
	if(objMHPLOT@numPvalOffset<=0 | objMHPLOT@numPvalOffset>1)
		stop(paste("EASY ERROR:MHPLOT\n numPvalOffset must be within ]0;1] \n !!!" ,sep=""))
	
	if(any(objMHPLOT@anumAddPvalLine<=0 | objMHPLOT@anumAddPvalLine>=1)) 
		stop(paste("EASY ERROR:MHPLOT\n Each value in anumAddPvalLine must be within ]0;1[ \n !!!" ,sep=""))
	
	if(length(objMHPLOT@anumAddPvalLine) != length(objMHPLOT@astrAddPvalLineCol)) 
		objMHPLOT@astrAddPvalLineCol = rep("red",length(objMHPLOT@anumAddPvalLine))
		
	if(length(objMHPLOT@anumAddPvalLine) != length(objMHPLOT@anumAddPvalLineLty)) 
		objMHPLOT@anumAddPvalLineLty = rep(6,length(objMHPLOT@anumAddPvalLine))

	#####
	
	if(!(objMHPLOT@numDefaultSymbol%in%c(0:25)))
		stop(paste("EASY ERROR:MHPLOT\n Wrong numDefaultSymbol defined. PLease use integer value between 0 and 25.", sep=""))

	if(objMHPLOT@numDefaultCex <= 0)
		stop(paste("EASY ERROR:MHPLOT\n Wrong numDefaultCex defined. PLease use positive numeric value.", sep=""))
	
	astrDefaultColourChr = objMHPLOT@astrDefaultColourChr
	
	if(length(astrDefaultColourChr)!=2) 
		stop(paste("EASY ERROR:MHPLOT\n Length of astrDefaultColourChr must be 2. PLease use two colours separated by ';'.", sep=""))
		
	for(i in 1:2) {
		strColour = astrDefaultColourChr[i]
		isOkColors =  strColour %in% colors()
		
		isOkColorHex =  nchar(strColour)==7 & substring(strColour,1,1)=="#"
		sTemp = substring(strColour,2,7) 
		isOkColorHex = isOkColorHex & all(strsplit(sTemp,"",fixed=TRUE)[[1]]%in%c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"))
		
		if(!isOkColors & !isOkColorHex)
			stop(paste("EASY ERROR:MHPLOT\n Wrong astrDefaultColourChr defined. PLease use color from R colors() or hexadecimal nomenclature, eg #FFFFFF.", sep=""))
	}
	
	if(objMHPLOT@astrColour[1] != "") {
		for(i in 1:length(objMHPLOT@astrColour)) {
			
			strColour = objMHPLOT@astrColour[i]
			isOkColors = strColour %in% colors()
			
			isOkColorHex =  nchar(strColour)==7 & substring(strColour,1,1)=="#"
			sTemp = substring(strColour,2,7) 
			isOkColorHex = isOkColorHex & all(strsplit(sTemp,"",fixed=TRUE)[[1]]%in%c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"))
			
			if(!isOkColors & !isOkColorHex)
				stop(paste("EASY ERROR:MHPLOT\n Wrong astrColour defined. PLease use color from R colors() or hexadecimal nomenclature, eg #FFFFFF.", sep=""))
		}
	}
	if(objMHPLOT@anumSymbol[1] != -1) {
		for(i in 1:length(objMHPLOT@anumSymbol)) {
			
			numSymbol = objMHPLOT@anumSymbol[i]
			if(!(numSymbol%in%c(0:25)))
				stop(paste("EASY ERROR:MHPLOT\n Wrong anumSymbol defined. PLease use integer value between 0 and 25.", sep=""))
		}
	}
	if(objMHPLOT@anumCex[1] != -1) {
		for(i in 1:length(objMHPLOT@anumCex)) {
			numCex = objMHPLOT@anumCex[i]
			if(numCex <= 0)
				stop(paste("EASY ERROR:MHPLOT\n Wrong anumCex defined. PLease use positive numeric value.", sep=""))			
		}
	}
	
	if(!(objMHPLOT@strFormat == "png" | objMHPLOT@strFormat == "pdf"))
		stop(paste("EASY ERROR:MHPLOT\n Wrong strFormat defined, PLease use 'png' or 'pdf' \n !!!", sep=""))
	
	if(objMHPLOT@numCexAxis <= 0) 
		stop(paste("EASY ERROR:MHPLOT\n Wrong numCexAxis defined, PLease use numCexAxis > 0 \n !!!", sep=""))
		
	if(objMHPLOT@numCexLab <= 0) 
		stop(paste("EASY ERROR:MHPLOT\n Wrong numCexLab defined, PLease use numCexLab > 0 \n !!!", sep=""))				
	
	if(length(objMHPLOT@anumParMar) != 4) 
		stop(paste("EASY ERROR:MHPLOT\n Wrong anumParMar defined, Please use anumParMar of length 4 (Default=5.1;4.1;4.1;2.1) \n !!!", sep=""))				
	
	if(length(objMHPLOT@anumParMgp) != 3) 
		stop(paste("EASY ERROR:MHPLOT\n Wrong anumParMgp defined, Please use anumParMgp of length 3 (Default=3;1;0) \n !!!", sep=""))				
	
	if(!(objMHPLOT@strParBty%in%c("o","l","7","c","u","]")))
		stop(paste("EASY ERROR:MHPLOT\n Wrong strParBty defined, Please use 'o','l','7','c','u' or ']' \n !!!", sep=""))				
		
	if(objMHPLOT@fileAnnot != "") {
		if(!file.exists(objMHPLOT@fileAnnot))
			stop(paste("EASY ERROR:MHPLOT\n File fileAnnot\n ",objMHPLOT@fileAnnot,"\n does not exist.", sep=""))
		### Cols exist?
			tblAnnot10<-read.table(objMHPLOT@fileAnnot,header=T, sep="",  stringsAsFactors=FALSE, nrows = 10)
		
			isAv <- "Chr" %in% names(tblAnnot10)
			if(!isAv)
				stop(paste(" EASY ERROR:MHPLOT\n Column 'Chr' is not available in fileAnnot. PLease correct file header.", sep=""))
			isAv <- "Pos" %in% names(tblAnnot10)
			if(!isAv)
				stop(paste(" EASY ERROR:MHPLOT\n Column 'Pos' is not available in fileAnnot. PLease correct file header.", sep=""))
			isAv <- "Colour" %in% names(tblAnnot10)
			if(!isAv)
				stop(paste(" EASY ERROR:MHPLOT\n Column 'Colour' is not available in fileAnnot. PLease correct file header.", sep=""))	
	}
	
}
MHPLOT.run <- function(objMHPLOT, objGWA) {

	colMHPlot 				<- 	objMHPLOT@colMHPlot
	colInChr 				<- 	objMHPLOT@colInChr
	colInPos 				<- 	objMHPLOT@colInPos
	astrDefaultColourChr	<- 	objMHPLOT@astrDefaultColourChr
	arcdColourCrit			<-	objMHPLOT@arcdColourCrit
	astrColour				<-	objMHPLOT@astrColour
	numPvalOffset 			<-	objMHPLOT@numPvalOffset
	arcdAdd2Plot			<-	objMHPLOT@arcdAdd2Plot
	## 130829
	anumAddPvalLine		<-	objMHPLOT@anumAddPvalLine
	astrAddPvalLineCol	<-	objMHPLOT@astrAddPvalLineCol
	anumAddPvalLineLty	<-	objMHPLOT@anumAddPvalLineLty
	blnYAxisBreak			<-	objMHPLOT@blnYAxisBreak
	numYAxisBreak			<-	objMHPLOT@numYAxisBreak
	numCexAxis				<-	objMHPLOT@numCexAxis
	numCexLab				<-	objMHPLOT@numCexLab
	blnLogPval				<-	objMHPLOT@blnLogPval	
	## 130927
	numDefaultSymbol		=	objMHPLOT@numDefaultSymbol
	numDefaultCex			=	objMHPLOT@numDefaultCex
	arcdSymbolCrit			=	objMHPLOT@arcdSymbolCrit
	anumSymbol				=	objMHPLOT@anumSymbol
	arcdCexCrit				=	objMHPLOT@arcdCexCrit
	anumCex					=	objMHPLOT@anumCex
	fileAnnot 				=	objMHPLOT@fileAnnot
	numAnnotPosLim			=	objMHPLOT@numAnnotPosLim
	numAnnotPvalLim			=	objMHPLOT@numAnnotPvalLim
	############################
	
	iVal = match(colMHPlot, names(objGWA@tblGWA))
	objGWA@tblGWA[,iVal] = as.numeric(objGWA@tblGWA[,iVal])
	iChr = match(colInChr, names(objGWA@tblGWA))
	objGWA@tblGWA[,iChr] = as.numeric(objGWA@tblGWA[,iChr])
	iPos = match(colInPos, names(objGWA@tblGWA))
	objGWA@tblGWA[,iPos] = as.numeric(objGWA@tblGWA[,iPos])
	
	isRemove = is.na(objGWA@tblGWA[,iVal]) | is.na(objGWA@tblGWA[,iChr]) | is.na(objGWA@tblGWA[,iPos])
	if(any(isRemove)) objGWA <- GWADATA.removerows(objGWA, which(isRemove))
	
	############################
	#### Compile initital colour/symbol/cex array 
	
	colplot <- plotCritKey <- rep(NA, dim(objGWA@tblGWA)[1])
	symplot <- rep(numDefaultSymbol, dim(objGWA@tblGWA)[1])
	cexplot <- rep(numDefaultCex, dim(objGWA@tblGWA)[1])
	
	if(arcdColourCrit[1] != "") {
		for(iColourCrit in 1:length(arcdColourCrit)) {
			rcdColourCrit 		<- arcdColourCrit[iColourCrit]
			strColour 			<- astrColour[iColourCrit]
			objRCD 	<- RCD(rcdColourCrit)
			out 	<- RCD.eval(objRCD, objGWA)
			colplot[out] <- strColour
			plotCritKey[out]= iColourCrit
		}
	} 
	if(arcdSymbolCrit[1] != "") {
		for(iSymbolCrit in 1:length(arcdSymbolCrit)) {
			rcdSymbolCrit 		<- arcdSymbolCrit[iSymbolCrit]
			numSymbol 			<- anumSymbol[iSymbolCrit]
			objRCD 	<- RCD(rcdSymbolCrit)
			out 	<- RCD.eval(objRCD, objGWA)
			symplot[out] <- numSymbol
		}
	}
	if(arcdCexCrit[1] != "") {
		for(iCexCrit in 1:length(arcdCexCrit)) {
			rcdCexCrit 		<- arcdCexCrit[iCexCrit]
			numCex 			<- anumCex[iCexCrit]
			objRCD 	<- RCD(rcdCexCrit)
			out 	<- RCD.eval(objRCD, objGWA)
			cexplot[out] <- numCex
		}
	}
	############################
	#### Log on P-Values
	if(blnLogPval) {
		isInf = objGWA@tblGWA[,iVal] == Inf
		if(any(isInf)) 
			warning(paste("EASY WARNING:\n There are ",length(which(isInf)), " SNPs with -log10_P=Inf being excluded from the MH plot for file \n ", objGWA@fileIn ,sep="" ))
		isOutRange = objGWA@tblGWA[,iVal] < 0
		if(any(isOutRange))
			warning(paste("EASY WARNING:\n There are ",length(which(isOutRange)), " SNPs with -log10_P<0 being excluded from the MH plot for file \n ", objGWA@fileIn ,sep="" ))
		isBelowOffset = objGWA@tblGWA[,iVal] < -log10(numPvalOffset)
		isRemove = isBelowOffset | isOutRange | isInf
		if(any(isRemove)) {
			objGWA <- GWADATA.removerows(objGWA, which(isRemove))
			colplot <- colplot[-which(isRemove)]
			symplot <- symplot[-which(isRemove)]
			cexplot <- cexplot[-which(isRemove)]
			plotCritKey <- plotCritKey[-which(isRemove)]
		}
	} else {
		isZero = objGWA@tblGWA[,iVal] == 0
		if(any(isZero)) 
			warning(paste("EASY WARNING:\n There are ",length(which(isZero)), " SNPs with P=0 being excluded from the MH plot for file \n ", objGWA@fileIn ,sep="" ))
		isOutRange = objGWA@tblGWA[,iVal] < 0 | objGWA@tblGWA[,iVal] > 1
		if(any(isOutRange))
			warning(paste("EASY WARNING:\n There are ",length(which(isOutRange)), " SNPs with P<0 or P>1 being excluded from the MH plot for file \n ", objGWA@fileIn ,sep="" ))
		isBelowOffset = objGWA@tblGWA[,iVal] > numPvalOffset
		isRemove = isBelowOffset | isOutRange | isZero
		if(any(isRemove)) {
			objGWA <- GWADATA.removerows(objGWA, which(isRemove))
			colplot <- colplot[-which(isRemove)]
			symplot <- symplot[-which(isRemove)]
			cexplot <- cexplot[-which(isRemove)]
			plotCritKey <- plotCritKey[-which(isRemove)]
		}
		objGWA@tblGWA[,iVal] 		= -log10(objGWA@tblGWA[,iVal])	
	}
	
	objGWA <- GWADATA.getcols(objGWA, c(colMHPlot, colInChr, colInPos))
	tblPlot <- objGWA@tblGWA

	names(tblPlot)[1] = "y"
	names(tblPlot)[2] = "chr"
	names(tblPlot)[3] = "pos"
	tblPlot = cbind(tblPlot, colplot, plotCritKey, symplot, cexplot, stringsAsFactors = FALSE)
	
	### Annot
	if(fileAnnot != "") {
		tblAnnot<-read.table(fileAnnot,header=T, sep="",  stringsAsFactors=FALSE)
		isSignif <- tblPlot$y > -log10(numAnnotPvalLim)
		isSignif[is.na(isSignif)] <- FALSE
		for(iLocus in 1:nrow(tblAnnot)) {
			#isLocus = is.na(tblPlot$colplot) & tblPlot$chr == tblAnnot$Chr[iLocus] & abs(tblPlot$pos-tblAnnot$Pos[iLocus])<numAnnotPosLim			
			isLocus = (tblPlot$chr == tblAnnot$Chr[iLocus]) & (abs(tblPlot$pos-tblAnnot$Pos[iLocus])<numAnnotPosLim)
			isColour = isLocus & isSignif
			if(any(isColour)) {
				tblPlot$colplot[isLocus] <- tblAnnot$Colour[iLocus]
				tblPlot$plotCritKey[isLocus] <- Inf
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
	isSetDefaultColour = is.na(tblPlot$colplot)

	tblPlot$colplot[isSetDefaultColour]	<-	ifelse(tblPlot$chr[isSetDefaultColour]%in%arr_chr[arr_idx2],astrDefaultColourChr[1],astrDefaultColourChr[2])
	
	plot.title = objGWA@fileInShortName
	strXlab = "Chromosome"
	if(!blnLogPval) strYlab = paste("-log10 (",colMHPlot,")",sep="")
	else strYlab = colMHPlot

	############################

	yAxisTickPos = NULL
	yAxisTickLab = TRUE
	anumAddPvalLinePos = -log10(anumAddPvalLine)
	
	############################	
	####### 130829 - Rescale y axis
	yMax = max(tblPlot$y)
	
	if(blnYAxisBreak & numYAxisBreak<yMax) {
	
		tblPlot = tblPlot[order(tblPlot$y),]
		y = tblPlot$y
		
		yu = y[y<=numYAxisBreak] # [0;22] -> [0;80]
		yo = y[y>numYAxisBreak] # [22;yMax] -> [80;100]
		
		yus = yu*80/numYAxisBreak
		yos = 80 - 20*(numYAxisBreak-yo)/(yMax-numYAxisBreak)
		
		tblPlot$y = c(yus,yos)
		
		yAxisTickLabU = pretty(c(0,numYAxisBreak),n=4,min.n=0)
		yAxisTickLabU = yAxisTickLabU[-which(yAxisTickLabU>=numYAxisBreak)]
		yAxisTickPosU = yAxisTickLabU*80/numYAxisBreak
		
		yAxisTickLabO = pretty(c(numYAxisBreak,yMax),n=4,min.n=0)
		yAxisTickLabO = yAxisTickLabO[-which(yAxisTickLabO<=numYAxisBreak)]
		yAxisTickPosO = 80 - 20*(numYAxisBreak-yAxisTickLabO)/(yMax-numYAxisBreak)
		
		yAxisTickLab = c(yAxisTickLabU,yAxisTickLabO)
		yAxisTickPos = c(yAxisTickPosU,yAxisTickPosO)
			
		isLineU = anumAddPvalLinePos<numYAxisBreak
		anumAddPvalLinePos[isLineU] = anumAddPvalLinePos[isLineU]*80/numYAxisBreak
		anumAddPvalLinePos[!isLineU] = 80 - 20*(numYAxisBreak-anumAddPvalLinePos[!isLineU])/(yMax-numYAxisBreak)
		
	}
	############################
	####### Start plot
	
	tblPlot = tblPlot[order(tblPlot$plotCritKey,na.last=FALSE),]
	
	plot( tblPlot$x, tblPlot$y, col=tblPlot$colplot, pch=tblPlot$symplot, cex=tblPlot$cexplot, cex.lab = numCexLab, xaxt="n",yaxt="n", ylab=strYlab, xlab="",ylim=c(0,max(tblPlot$y)))
	
	# Achsenlinien:
	axis(2, at =yAxisTickPos,  labels = yAxisTickLab, cex.axis = numCexAxis)
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
	}
	#abline(yGwsLine,0,col = "red",lty = 6)	
	#axis(4,at=yGwsLine,labels="5e-8",las=3,cex.axis = 1.5)
	
	if(blnYAxisBreak & numYAxisBreak<yMax) {
		axis.break(2, 80, style = "zigzag")
		abline(80,0,col = "grey",lty = 1)	
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
MHPLOT <- function(strEqcCommand){ 
	## Wrapper for class definition
	MHPLOTout <- setMHPLOT(new("MHPLOT", strEqcCommand = strEqcCommand))
	return(MHPLOTout)
}

