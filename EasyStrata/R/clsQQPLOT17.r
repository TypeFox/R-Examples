setClass("QQPLOT",
	representation = representation(
						strEqcCommand		=	"character",
						acolQQPlot			=	"character",
						astrColour			=	"character",
						anumSymbol			=	"numeric",
						anumCex				=	"numeric",
						blnCombined			=	"logical",
						blnLegend			=	"logical",
						astrLegendText		=	"character",
						arcdExclude			=	"character",
						numPvalOffset		=	"numeric",
						blnYAxisBreak		=	"logical",
						numYAxisBreak		=	"numeric",
						blnLogPval			=	"logical",
						blnPlotCI			=	"logical",
						fileRemove			=	"character",
						numRemovePosLim		=	"numeric",
						strRemovedColour	=	"character",
						numRemovedSymbol	=	"numeric",
						numRemovedCex		=	"numeric",
						colInChr			=	"character",
						colInPos			=	"character",
						strAxes				=	"character", ## added 130610
						strXlab				=	"character",
						strYlab				=	"character",
						strTitle			=	"character",
						arcdAdd2Plot		=	"character",
						strMode				=	"character",
						strFormat			=	"character",
						numCexAxis			=	"numeric",
						numCexLab			=	"numeric",
						numWidth			=	"numeric",
						numHeight			=	"numeric",
						anumParMar			=	"numeric",
						anumParMgp			=	"numeric",
						strParBty			=	"character",
						strPlotName			= 	"character",
						numParLas			=	"numeric"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						acolQQPlot			=	"",
						astrColour			=	"black",
						anumSymbol			=	20,
						anumCex				=	1,
						blnCombined			=	FALSE,
						blnLegend			=	FALSE,
						astrLegendText		=	"",
						arcdExclude			=	"",
						numPvalOffset		=	1,
						blnYAxisBreak		=	FALSE,
						numYAxisBreak		=	22,
						blnLogPval			=	FALSE,
						blnPlotCI			=	FALSE,
						fileRemove			=	"",
						numRemovePosLim		=	500000,
						strRemovedColour	=	"green",
						numRemovedSymbol	=	20,
						numRemovedCex		=	1,
						colInChr			=	"",
						colInPos			=	"",
						strAxes				=	"",
						strXlab				=	"",
						strYlab				=	"",
						strTitle			=	"",
						arcdAdd2Plot		=	"",
						strMode				=	"singleplot",
						strFormat			=	"png",
						numCexAxis			=	1,
						numCexLab			=	1,
						numWidth			=	640,	# pixel px # needs change for pdf!
						numHeight			=	640,	# pixel px # needs change for pdf!
						anumParMar			=	c(5.1, 4.1, 4.1, 2.1),
						anumParMgp			=	c(3, 1, 0),
						strParBty			=	"o",
						strPlotName			= 	"qq",
						numParLas			=	0
						)
	#contains = c("EcfReader")
)

setGeneric("setQQPLOT", function(object) standardGeneric("setQQPLOT"))
setMethod("setQQPLOT", signature = (object = "QQPLOT"), function(object) {
	
	aEqcSlotNamesIn = c("acolQQPlot", 
						"astrColour", 
						"anumSymbol", 
						"anumCex", 
						"blnCombined", 
						"blnLegend", 
						"astrLegendText", 
						"strMode",
						"arcdExclude", 
						"numPvalOffset",
						"strAxes",
						"strXlab",
						"strYlab",
						"strTitle",
						"arcdAdd2Plot",
						"strFormat",
						"numCexAxis",
						"numCexLab",
						"numWidth",
						"numHeight",
						"anumParMar",
						"anumParMgp",
						"strParBty",
						"strPlotName",
						"blnYAxisBreak",
						"numYAxisBreak",
						"blnLogPval",
						"blnPlotCI",
						"fileRemove",
						"numRemovePosLim",
						"strRemovedColour",
						"numRemovedSymbol",
						"numRemovedCex",
						"colInChr",
						"colInPos"
						)
	#aEcfSlotNamesIn = c("arcdAddCol", "astrAddColNames")

	objEqcReader <- EqcReader(object@strEqcCommand,aEqcSlotNamesIn)
	
	if(length(objEqcReader@lsEqcSlotsOut) > 0) {
		for(i in 1:length(objEqcReader@lsEqcSlotsOut)) {
			tmpSlot <- names(objEqcReader@lsEqcSlotsOut)[i]
			tmpSlotVal <- objEqcReader@lsEqcSlotsOut[[i]]
			
			if(tmpSlot == "arcdExclude") tmpSlotVal[is.na(tmpSlotVal)]=""
			if(all(!is.na(tmpSlotVal))) slot(object, tmpSlot) <- tmpSlotVal
		}
	}
	isDefaultPngSize = object@numWidth == 640 & object@numHeight == 640
	if(object@strFormat == "pdf" & isDefaultPngSize) {
		object@numWidth <- 6
		object@numHeight <- 6
	}
	
	return(object)
})

#############################################################################################################################
validQQPLOT <- function(objQQPLOT) {
	## Paste validity checks specific for QQPLOT
	## Pval <-> se
	## fileGcSnps exists
	## colMarkerGcSnps exists
	
	if(!(objQQPLOT@strMode == "singleplot" | objQQPLOT@strMode == "subplot"))
		stop(paste("EASY ERROR:QQPLOT\n Wrong strMode defined. PLease use 'singleplot' or ''subplot'. \n", sep=""))
	
	if(objQQPLOT@numPvalOffset<=0 | objQQPLOT@numPvalOffset>1)
		stop(paste("EASY ERROR:QQPLOT\n numPvalOffset must be within ]0;1] \n !!!" ,sep=""))
	
	
	
	if(!(objQQPLOT@strFormat == "png" | objQQPLOT@strFormat == "pdf"))
		stop(paste("EASY ERROR:SPLOT\n Wrong strFormat defined, PLease use 'png' or 'pdf' \n !!!", sep=""))
	
	if(objQQPLOT@numCexAxis <= 0) 
		stop(paste("EASY ERROR:SPLOT\n Wrong numCexAxis defined, PLease use numCexAxis > 0 \n !!!", sep=""))
		
	if(objQQPLOT@numCexLab <= 0) 
		stop(paste("EASY ERROR:SPLOT\n Wrong numCexLab defined, PLease use numCexLab > 0 \n !!!", sep=""))				
	
	if(length(objQQPLOT@anumParMar) != 4) 
		stop(paste("EASY ERROR:SPLOT\n Wrong anumParMar defined, Please use anumParMar of length 4 (Default=5.1;4.1;4.1;2.1) \n !!!", sep=""))				
	
	if(length(objQQPLOT@anumParMgp) != 3) 
		stop(paste("EASY ERROR:SPLOT\n Wrong anumParMgp defined, Please use anumParMgp of length 3 (Default=3;1;0) \n !!!", sep=""))				
	
	if(!(objQQPLOT@strParBty%in%c("o","l","7","c","u","]")))
		stop(paste("EASY ERROR:SPLOT\n Wrong strParBty defined, Please use 'o','l','7','c','u' or ']' \n !!!", sep=""))				
	
	
	if(objQQPLOT@fileRemove != "") {
		if(!file.exists(objQQPLOT@fileRemove))
			stop(paste("EASY ERROR:QQPLOT\n File fileRemove\n ",objQQPLOT@fileRemove,"\n does not exist.", sep=""))
		### Cols exist?
		tblAnnot10<-read.table(objQQPLOT@fileRemove,header=T, sep="",  stringsAsFactors=FALSE, nrows = 10)
	
		isAv <- "Chr" %in% names(tblAnnot10)
		if(!isAv)
			stop(paste(" EASY ERROR:QQPLOT\n Column 'Chr' is not available in fileAnnot. PLease correct file header.", sep=""))
		isAv <- "Pos" %in% names(tblAnnot10)
		if(!isAv)
			stop(paste(" EASY ERROR:QQPLOT\n Column 'Pos' is not available in fileAnnot. PLease correct file header.", sep=""))
	}
	
	
	return(TRUE)
}

QQPLOT.GWADATA.valid <- function(objQQPLOT, objGWA) {
	
	aisMatch = objQQPLOT@acolQQPlot %in% objGWA@aHeader
	
	if(!all(aisMatch))
		stop(paste("EASY ERROR:QQPLOT\n Column \n",paste(objQQPLOT@acolQQPlot[which(!aisMatch)],collapse="\n")," does not exist in file\n",objGWA@fileInShortName,"\n !!!" ,sep=""))
	
	if(objQQPLOT@fileRemove != "") {
		if(objQQPLOT@colInChr == "")
			stop(paste("EASY ERROR:QQPLOT\n Required parameter 'colInChr' is not set. This is needed to remove loci based on position !!!" ,sep=""))
		if(objQQPLOT@colInPos == "")
			stop(paste("EASY ERROR:QQPLOT\n Required parameter 'colInPos' is not set. This is needed to remove loci based on position !!!" ,sep=""))
		isAv <- objQQPLOT@colInChr %in% objGWA@aHeader
		if(!isAv)
			stop(paste("EASY ERROR:QQPLOT\n Column \n",objQQPLOT@colInChr," does not exist in file\n",objGWA@fileInShortName,"\n !!!" ,sep=""))
		isAv <- objQQPLOT@colInPos %in% objGWA@aHeader
		if(!isAv)
			stop(paste("EASY ERROR:QQPLOT\n Column \n",objQQPLOT@colInPos," does not exist in file\n",objGWA@fileInShortName,"\n !!!" ,sep=""))
	}
	
}


QQPLOT.run <- function(objQQPLOT, objGWA) {
	#function(tblPlot,acolQQPlot,astrColour,blnCombined,rcdNumCIBounds,fileOut)
	
	tblPlot 			<- objGWA@tblGWA
	acolQQPlot 			<- objQQPLOT@acolQQPlot
	astrColour 			<- objQQPLOT@astrColour
	anumSymbol 			<- objQQPLOT@anumSymbol
	anumCex 			<- objQQPLOT@anumCex
	arcdExclude 		<- objQQPLOT@arcdExclude
	numPvalOffset 		<- objQQPLOT@numPvalOffset
	blnCombined 		<- objQQPLOT@blnCombined
	astrLegendText		<- objQQPLOT@astrLegendText
	strAxes				<- objQQPLOT@strAxes
	strXlab				<- objQQPLOT@strXlab
	strYlab				<- objQQPLOT@strYlab
	strTitle			<- objQQPLOT@strTitle
	arcdAdd2Plot		<- objQQPLOT@arcdAdd2Plot
	blnLegend			<- objQQPLOT@blnLegend
	numCexAxis			<- objQQPLOT@numCexAxis
	numCexLab			<- objQQPLOT@numCexLab
	blnYAxisBreak		<- objQQPLOT@blnYAxisBreak
	numYAxisBreak		<- objQQPLOT@numYAxisBreak
	blnLogPval			<- objQQPLOT@blnLogPval
	blnPlotCI			<- objQQPLOT@blnPlotCI
	fileRemove 			<- objQQPLOT@fileRemove
	numRemovePosLim		<- objQQPLOT@numRemovePosLim
	strRemovedColour	<- objQQPLOT@strRemovedColour
	numRemovedSymbol	<- objQQPLOT@numRemovedSymbol
	numRemovedCex		<- objQQPLOT@numRemovedCex
	colInChr			<- objQQPLOT@colInChr
	colInPos			<- objQQPLOT@colInPos
	############################
	
	#if(objQQPLOT@strMode == "singleplot") plot.title = objGWA@fileInShortName
	#else plot.title = ""
	if(strTitle == "") plot.title = objGWA@fileInShortName
	else plot.title = strTitle
	
	if(length(astrColour) == 1 & astrColour[1] == "black")
		astrColour 	= rep(astrColour , length(acolQQPlot))
		
	if(length(anumSymbol) == 1 & anumSymbol[1] == 20)
		anumSymbol 	= rep(anumSymbol , length(acolQQPlot))
		
	if(length(anumCex) == 1 & anumCex[1] == 1)
		anumCex = rep(anumCex , length(acolQQPlot))
	
	### Don't scream if length is 1
	
	if(length(acolQQPlot) != length(astrColour)) {
		astrColour = rep(astrColour[1] , length(acolQQPlot)) 
		warning(paste("EASY WARNING:QQPLOT\n Length of astrColour does not match length of acolQQPlot. \n Firts Colour ",astrColour[1]," will be used for all P-Values." ,sep="" ))
	}
	if(length(acolQQPlot) != length(anumSymbol)) {
		anumSymbol = rep(anumSymbol[1] , length(acolQQPlot)) 
		warning(paste("EASY WARNING:QQPLOT\n Length of anumSymbol does not match length of acolQQPlot. \n First Symbol ",anumSymbol[1]," will be used for all P-Values." ,sep="" ))
	}
	if(length(acolQQPlot) != length(anumCex)) {
		anumCex = rep(anumCex[1] , length(acolQQPlot)) 
		warning(paste("EASY WARNING:QQPLOT\n Length of anumCex does not match length of acolQQPlot. \n First symbol size ",anumCex[1]," will be used for all P-Values." ,sep="" ))
	}
	
	
	### Get legend params
	

	if(astrLegendText[1] == "") astrLegendText = acolQQPlot
	astrLegendColour = astrColour
	anumLegendPch = anumSymbol
	###
	
	# xlab = expression("Expected "* -log[10] * " pvalue")
	#print("QQ: Set labels")
	## Set axes labels:
	if(strXlab=="") {
		if(length(acolQQPlot) == 1)
			strXlab = paste("expected -log10 (",acolQQPlot[1],")",sep="")
		else 
			strXlab = expression("expected "* -log[10] * " P-value")
	}
	if(strYlab=="") {
		if(length(acolQQPlot) == 1) 
			strYlab = paste("observed -log10 (",acolQQPlot[1],")",sep="")
		else 
			strYlab = expression("observed "* -log[10] * " P-value")
	}
	
	#print("QQ: Get axis lims")
	## Get axis limits:
	aicolQQ = match(acolQQPlot, names(tblPlot))
	#print(aicolQQ)
	for(icolQQ in aicolQQ) {
		tblPlot[,icolQQ] <- as.numeric(tblPlot[,icolQQ])
		#print(class(tblPlot[,icolQQ]))
	}
	
		
	### Get axis lims strAxes
	plot.xlims <- plot.ylims <- rep(NA,2)
	# old:
	numYlimMax = 0
	for(icolQQ in aicolQQ) {
		pTmp = tblPlot[,icolQQ]
		if(blnLogPval) {
			isinf = pTmp == Inf
			pTmp[isinf] = NA
			numYlimMax = max( numYlimMax, max(pTmp, na.rm=TRUE) )
		} else {
			is0 = pTmp == 0
			pTmp[is0] = NA
			numYlimMax = max( numYlimMax, max(-log10(pTmp), na.rm=TRUE) )
		}
	}
	plot.ylims = c(0, numYlimMax)
	numXlimMax = 0
	if(blnCombined) {
		numClean = 0
		for(icolQQ in aicolQQ) numClean = numClean + length(which(!is.na(tblPlot[,icolQQ])))
		numXlimMax = -log10(1/numClean)
	}
	else {
		for(icolQQ in aicolQQ) numXlimMax = max( numXlimMax, -log10(1/length(which(!is.na(tblPlot[,icolQQ])))))
	}
	plot.xlims = c(0, numXlimMax)
	# old!
		
	if(length(grep("lim",strAxes)) > 0)  {
		## strSPAxes = "lim(x0,x1,y0,y1)"
		str_lim = substr(strAxes,5,nchar(strAxes)-1)
		## str_lim = "x0,x1,y0,y1"
		arr_str_lim = strsplit(str_lim,",")[[1]]
		plot.xlims[1] = ifelse(arr_str_lim[1]!="NULL",as.numeric(arr_str_lim[1]),plot.xlims[1])
		plot.xlims[2] = ifelse(arr_str_lim[2]!="NULL",as.numeric(arr_str_lim[2]),plot.xlims[2])
		plot.ylims[1] = ifelse(arr_str_lim[3]!="NULL",as.numeric(arr_str_lim[3]),plot.ylims[1])
		plot.ylims[2] = ifelse(arr_str_lim[4]!="NULL",as.numeric(arr_str_lim[4]),plot.ylims[2])
	} 
	if(strAxes == "equal")  {
		plot.xlims[1] <- plot.ylims[1] <- 0
		plot.xlims[2] <- plot.ylims[2] <- max(plot.xlims[2], plot.ylims[2])
	} 
		
	#print("QQ: Start plotting")
	
	##########################################################################################################################################

	plot.cex.main = 1.2
	titleWidth = nchar(plot.title)*par()$cin[1]
	if(titleWidth*1.2 > par()$fin[1]) {
		## Zeilenumbruch in Haelfte einfuegen
		iBreakUp = floor(nchar(plot.title)/2)+1
		plot.title = paste(substr(plot.title, 1, iBreakUp), substr(plot.title, iBreakUp+1, nchar(plot.title)), sep="\n")
		if((iBreakUp+1)*1.2*par()$cin[1] > par()$pin[1]) {
		# cex.main = 1.2 (default)
		# par()$cin[1] default character width [inches]
		# par()$pin[1] default plot width [inches]
		# par()$fin[1] default figure width [inches]
		# Rescale
		#plot.cex.main = par()$pin[1]/(nchar(plot.title)*par()$cin[1]*1.2)
		#titleWidth = nchar(plot.title)*par()$cin[1]
		#plot.cex.main = 1.2*par()$fin[1]/titleWidth
		titleWidth = (iBreakUp+1)*par()$cin[1]
		plot.cex.main = 1.2*par()$fin[1]/titleWidth
		}
	}

	
	iPobsUsed = c()
	yAxisTickPos = NULL
	yAxisTickLab = TRUE
	numSlopeID = 1
	#blnYAxisBreak = TRUE
	
	for(i in 1:length(acolQQPlot)) {
	
		iPobs = match(acolQQPlot[i],names(tblPlot))
		aPobs = as.numeric(tblPlot[,iPobs])
		
		###dim(tblplot)==dim(objGWA@tblGWA)!
		
		rcdExclude <- arcdExclude[i]
		if(is.na(rcdExclude)) rcdExclude <- ""
		if(rcdExclude!="") {
			objRCD 	<- RCD(rcdExclude)
			out 	<- RCD.eval(objRCD, objGWA)
			if(any(out)) aPobs <- aPobs[!out]
			else warning(paste("EASY WARNING:\n No SNPs fullfilled the criterion \n",rcdExclude, "\n. in  arcdExclude for file \n ", objGWA@fileIn ,sep="" ))
		}
		###
		
		isMiss = is.na(aPobs)
		if(blnLogPval) {
			isInf = aPobs == Inf
			if(any(isInf & !isMiss)) warning(paste("EASY WARNING:\n There are ",length(which(isInf)), " SNPs with -log10_P=Inf being excluded from the QQ plot for file \n ", objGWA@fileIn ,sep="" ))
			isOutRange = aPobs < 0
			if(any(isOutRange & !isMiss)) warning(paste("EASY WARNING:\n There are ",length(which(isOutRange)), " SNPs with -log10_P<0 being excluded from the QQ plot for file \n ", objGWA@fileIn ,sep="" ))
			aPobsClean = aPobs[!(isMiss | isInf | isOutRange)]
			aPobsCleanSort = sort(aPobsClean,decreasing=TRUE)
			aPobsCleanSortLog = aPobsCleanSort
		} else {
			isZero = aPobs == 0
			if(any(isZero & !isMiss)) warning(paste("EASY WARNING:\n There are ",length(which(isZero)), " SNPs with P=0 being excluded from the QQ plot for file \n ", objGWA@fileIn ,sep="" ))
			isOutRange = aPobs < 0 | aPobs > 1
			if(any(isOutRange & !isMiss)) warning(paste("EASY WARNING:\n There are ",length(which(isOutRange)), " SNPs with P<0 or P>1 being excluded from the QQ plot for file \n ", objGWA@fileIn ,sep="" ))
			aPobsClean = aPobs[!(isMiss | isZero | isOutRange)]
			aPobsCleanSort = sort(aPobsClean)
			aPobsCleanSortLog = -log10(aPobsCleanSort)
		}
		
		rm(aPobs)
		rm(aPobsClean)
		rm(aPobsCleanSort)
		
		aPexpLog = -log10(seq(0,1,length=length(aPobsCleanSortLog)+1)[-1])
		
		N_SNPS <- length(aPexpLog)
		
		if(i == 1 & blnPlotCI) {
			yCiUppLog <- yCiDownLog <- rep(0,N_SNPS)
			for(iCi in 1:N_SNPS){
				yCiUppLog[iCi] <- -log10(qbeta(0.95,iCi,N_SNPS-iCi+1))
				yCiDownLog[iCi] <- -log10(qbeta(0.05,iCi,N_SNPS-iCi+1))
			}
		}
		
		isPlotted = aPobsCleanSortLog >= -log10(numPvalOffset)
		
		if(any(isPlotted)) {
			aPobsCleanSortLog <- aPobsCleanSortLog[isPlotted]
			aPexpLog <- aPexpLog[isPlotted]
			if(i == 1 & blnPlotCI) {
				yCiUppLog <- yCiUppLog[isPlotted]
				yCiDownLog <- yCiDownLog[isPlotted]
			}
		} else {
			warning(paste("EASY WARNING:QQPLOT:\n There are no SNPs left after applying the plotting threshold numPvalOffset for file \n ", objGWA@fileIn ,"\n Instead all SNPs are plotted!!!",sep="" ))
		}
		
		aColour <- rep(astrColour[i],length(aPexpLog))
		aSymbol <- rep(anumSymbol[i],length(aPexpLog))
		aCex 	<- rep(anumCex[i],length(aPexpLog))
		
		x = rev(aPexpLog)
		y = rev(aPobsCleanSortLog)
		if(i == 1 & blnPlotCI) {
			yCiUpp <- rev(yCiUppLog)
			yCiDown <- rev(yCiDownLog)
		}
		
		if(i == 1) yMax = plot.ylims[2]
		if(blnYAxisBreak & numYAxisBreak<yMax) {
			
			#numYAxisBreak = 11
			#yMax = max(y)
			# yu = y[y<=numYAxisBreak] # [0;22] -> [0;80]
			# yo = y[y>numYAxisBreak] # [22;yMax] -> [80;100]
			# yus = yu*80/numYAxisBreak
			# yos = 80 - 20*(numYAxisBreak-yo)/(yMax-numYAxisBreak)
			# y = c(yus,yos)
			
			y <- ifelse(y<=numYAxisBreak, y*80/numYAxisBreak, 80 - 20*(numYAxisBreak-y)/(yMax-numYAxisBreak))
			
			numSlopeID = 80/numYAxisBreak
			
			plot.ylims = c(0,100)
			
			if(i == 1) {
				yAxisTickLabU = pretty(c(0,numYAxisBreak),n=4,min.n=0)
				yAxisTickLabU = yAxisTickLabU[-which(yAxisTickLabU>=numYAxisBreak)]
				yAxisTickPosU = yAxisTickLabU*80/numYAxisBreak
				
				yAxisTickLabO = pretty(c(numYAxisBreak,yMax),n=4,min.n=0)
				yAxisTickLabO = yAxisTickLabO[-which(yAxisTickLabO<=numYAxisBreak)]
				yAxisTickPosO = 80 - 20*(numYAxisBreak-yAxisTickLabO)/(yMax-numYAxisBreak)
				
				yAxisTickLab = c(yAxisTickLabU,yAxisTickLabO)
				yAxisTickPos = c(yAxisTickPosU,yAxisTickPosO)
				
				if(blnPlotCI) {
					### rescale CIs
					yCiUpp <- ifelse(yCiUpp<=numYAxisBreak, yCiUpp*80/numYAxisBreak, 80 - 20*(numYAxisBreak-yCiUpp)/(yMax-numYAxisBreak))
					yCiDown <- ifelse(yCiDown<=numYAxisBreak, yCiDown*80/numYAxisBreak, 80 - 20*(numYAxisBreak-yCiDown)/(yMax-numYAxisBreak))
				}
			}
		}
		
		if(i ==1) 	{
			
			plot(	NULL,NULL, xlab = strXlab,ylab = strYlab, xlim = plot.xlims, ylim = plot.ylims, 
					main = plot.title,cex.main=plot.cex.main,cex.axis=numCexAxis,
					cex.lab=numCexLab,yaxt="n",	las=objQQPLOT@numParLas
				)
			if(blnPlotCI) {
				points(x,yCiUpp,col="grey",type="l")
				points(x,yCiDown,col="grey",type="l")
			}
			points(x,y,pch=aSymbol,col=aColour,cex = aCex)
			# plot(x,y,pch=aSymbol,col=aColour,cex = aCex, 
					# xlab = strXlab,ylab = strYlab, 
					# xlim = plot.xlims, ylim = plot.ylims, 
					# main = plot.title,
					# cex.main=plot.cex.main,cex.axis=numCexAxis,cex.lab=numCexLab,
					# yaxt="n"
			# )
			axis(2,at=yAxisTickPos,labels = yAxisTickLab, cex.axis = numCexAxis, cex.lab=numCexLab)
			abline(a=0,b=numSlopeID,col="black",lty=1,xpd=FALSE)
			if(blnYAxisBreak & numYAxisBreak<yMax) {
				axis.break(2, 80, style = "zigzag")
				abline(80,0,col = "grey",lty = 1)	
			}			
		} else	{
			points(x,y,col=aColour,pch=aSymbol,cex=aCex)
		}
		
		iPobsUsed = c(iPobsUsed, iPobs)
		rm(aPexpLog)
		rm(aPobsCleanSortLog)
	}
	
	
	if(fileRemove != "") {
	
		tblRemove<-read.table(fileRemove,header=T, sep="",  stringsAsFactors=FALSE)
		isRemove <- rep(FALSE, nrow(tblPlot))
		for(iLocus in 1:nrow(tblRemove)) {
			isRemove <- isRemove | tblPlot[,colInChr] == tblRemove$Chr[iLocus] & abs(tblPlot[,colInPos]-tblRemove$Pos[iLocus])<numRemovePosLim
		}
		## Remove SNPs from very first P-Val mentioned
		iPobs = match(acolQQPlot[1],names(tblPlot))
		aPobs = as.numeric(tblPlot[,iPobs])
		
		aPobs <- aPobs[!isRemove]
		aPobs <- aPobs[!is.na(aPobs)]

		if(blnLogPval) aPobsLog = sort(aPobs,decreasing=TRUE)
		else aPobsLog = -log10(sort(aPobs))
		rm(aPobs)
		aPexpLog = -log10(seq(0,1,length=length(aPobsLog)+1)[-1])
		
		xr = rev(aPexpLog)
		yr = rev(aPobsLog)

		if(blnYAxisBreak & numYAxisBreak<yMax) {
			yr <- ifelse(yr<=numYAxisBreak, yr*80/numYAxisBreak, 80 - 20*(numYAxisBreak-yr)/(yMax-numYAxisBreak))
		}
		points(xr, yr, col = strRemovedColour , pch = numRemovedSymbol, cex = numRemovedCex)
		rm(xr)
		rm(yr)
	}
	
	if(blnCombined) {
		
		aPobsCombined = c()
		for(iPobs in iPobsUsed) aPobsCombined = c(aPobsCombined,as.numeric(tblPlot[,iPobs]))
		
		aPobsCombinedClean = aPobsCombined[!is.na(aPobsCombined)]
		rm(aPobsCombined)
		
		if(blnLogPval) aPobsCombinedCleanSortLog = sort(aPobsCombinedClean,decreasing=TRUE)
		else aPobsCombinedCleanSortLog = -log10(sort(aPobsCombinedClean))
		rm(aPobsCombinedClean)
		aPexpLog = -log10(seq(0,1,length=length(aPobsCombinedCleanSortLog)+1)[-1])
		
		xc = rev(aPexpLog)
		yc = rev(aPobsCombinedCleanSortLog)

		if(blnYAxisBreak & numYAxisBreak<yMax) {
			yc <- ifelse(yc<=numYAxisBreak, yc*80/numYAxisBreak, 80 - 20*(numYAxisBreak-yc)/(yMax-numYAxisBreak))		
			# ycu = yc[yc<=numYAxisBreak] # [0;22] -> [0;80]
			# yco = yc[yc>numYAxisBreak] # [22;yMax] -> [80;100]
			# ycu = ycu*80/numYAxisBreak
			# yco = 80 - 20*(numYAxisBreak-yco)/(yMax-numYAxisBreak)
			# yc = c(ycu,yco)
		}
		points(xc,yc, col = "black" , pch = 2)
		astrLegendText = c(astrLegendText, "COMBINED")
		astrLegendColour = c(astrLegendColour, "black")
		anumLegendPch = c(anumLegendPch, 2)
	}
	
	if(arcdAdd2Plot[1] != "") {
		for(rcdTmp in arcdAdd2Plot) {
			for(colTmp in names(tblPlot)) {
				if(length(grep(colTmp,rcdTmp)) > 0) {
					icolTmp=which(names(tblPlot)==colTmp)
					assign(colTmp,as.numeric(as.character(tblPlot[,icolTmp])))
					# create arrays, that are used within criterion and name them accordingly
					# tblIn$pmen -> create array pmen if criterion e.g. pmen>0.3
				}
			}
			
			expTmp=parse(text=rcdTmp)	
			
			eval(expTmp)
		}
	}
	
	if(blnLegend) {
		
		ausr = par()$usr
		xleg = ausr[2]
		yleg = ausr[4]

		legend(	x=xleg, 
				y=yleg,
				inset=0.005,
				legend=astrLegendText,
				col =astrLegendColour,
				pch = anumLegendPch, 
				#cex = legendcex,
				bg="skyblue",
				xpd=TRUE
				)
	}
}


QQPLOT <- function(strEqcCommand){ 
	## Wrapper for class definition
	
	
	
	QQPLOTout <- setQQPLOT(new("QQPLOT", strEqcCommand = strEqcCommand))
	validQQPLOT(QQPLOTout)
	#ADDCOLout.valid <- validADDCOL(ADDCOLout)
	return(QQPLOTout)
	#validECF(ECFout)
	#return(ECFout)
	
	## Identical:
	# ECFin <- new("ECF5", fileECF = fileECFIn) 
	# ECFout <- setECF5(ECFin)
	# return(ECFout)
}

