setClass("SPLOT",
	representation = representation(
						strEqcCommand		=	"character",
						rcdSPlotX			=	"character",
						rcdSPlotY			=	"character",
						strDefaultColour	=	"character",
						numDefaultSymbol	=	"numeric",
						numDefaultCex		=	"numeric",
						arcdColourCrit		=	"character",
						astrColour			=	"character",
						astrColourLegend	=	"character",
						arcdSymbolCrit		=	"character",
						anumSymbol			=	"numeric",
						astrSymbolLegend	=	"character",
						arcdCexCrit			=	"character",
						anumCex				=	"numeric",
						#astrSPlotCexLegend	=	"character",
						strAxes				=	"character",
						strXlab				=	"character",
						strYlab				=	"character",
						strTitle			=	"character",
						blnLegend			=	"logical",
						arcdAdd2Plot		=	"character",
						strMode				=	"character",
						strFormat			=	"character",
						rcdExclude			=	"character",
						numCexAxis			=	"numeric",
						numCexLab			=	"numeric",
						numWidth			=	"numeric",
						numHeight			=	"numeric",
						anumParMar			=	"numeric",
						anumParMgp			=	"numeric",
						strParBty			=	"character",
						blnGrid				=	"logical",
						strPlotName			= 	"character",
						blnPlotCI			= 	"logical",
						rcdCIlengthX		=	"character",
						rcdCIlengthY		=	"character",
						strCIColour			=	"character",
						blnAxesLines		=	"logical",
						numParLas			=	"numeric"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						rcdSPlotX			=	"",
						rcdSPlotY			=	"",
						strDefaultColour	=	"black",
						numDefaultSymbol	=	20,
						numDefaultCex		=	1,
						arcdColourCrit		=	"",
						astrColour			=	"",
						astrColourLegend	=	"",
						arcdSymbolCrit		=	"",
						anumSymbol			=	-1,
						astrSymbolLegend	=	"",
						arcdCexCrit			=	"",
						anumCex				=	-1,
						#astrSPlotCexLegend	=	"",
						strAxes				=	"",
						strXlab				=	"",
						strYlab				=	"",
						strTitle			=	"",
						blnLegend			=	FALSE,
						arcdAdd2Plot		=	"",
						strMode				=	"singleplot",
						strFormat			=	"png",
						rcdExclude			=	"",
						numCexAxis			=	1,
						numCexLab			=	1,
						numWidth			=	640,	# pixel px # needs change for pdf!
						numHeight			=	640,	# pixel px # needs change for pdf!
						anumParMar			=	c(5.1, 4.1, 4.1, 2.1),
						anumParMgp			=	c(3, 1, 0),
						strParBty			=	"o",
						blnGrid				=	TRUE,
						strPlotName			= 	"sp",
						blnPlotCI			= 	FALSE,
						rcdCIlengthX		=	"",
						rcdCIlengthY		=	"",
						strCIColour			=	"grey",
						blnAxesLines		=	TRUE,
						numParLas			=	0
						)
	#contains = c("EcfReader")
)

setGeneric("setSPLOT", function(object) standardGeneric("setSPLOT"))
setMethod("setSPLOT", signature = (object = "SPLOT"), function(object) {
	
	aEqcSlotNamesIn = c("rcdSPlotX", 
						"rcdSPlotY", 
						"strDefaultColour",
						"numDefaultSymbol",
						"numDefaultCex",
						"arcdColourCrit", 
						"astrColour",
						"astrColourLegend",
						"arcdSymbolCrit",
						"anumSymbol",
						"astrSymbolLegend",
						"arcdCexCrit",
						"anumCex",
						#"astrSPlotCexLegend",
						"strAxes",
						"strXlab",
						"strYlab",
						"strTitle",
						"blnLegend",
						"arcdAdd2Plot",
						"strMode",
						"strFormat",
						"rcdExclude",
						"numCexAxis",
						"numCexLab",
						"numWidth",
						"numHeight",
						"anumParMar",
						"anumParMgp",
						"strParBty",
						"blnGrid",
						"strPlotName",
						"blnPlotCI",
						"rcdCIlengthX",
						"rcdCIlengthY",
						"strCIColour",
						"blnAxesLines",
						"numParLas"
						)
	#aEcfSlotNamesIn = c("arcdAddCol", "astrAddColNames")

	objEqcReader <- EqcReader(object@strEqcCommand,aEqcSlotNamesIn)
	
	if(length(objEqcReader@lsEqcSlotsOut) > 0) {
		for(i in 1:length(objEqcReader@lsEqcSlotsOut)) {
			tmpSlot <- names(objEqcReader@lsEqcSlotsOut)[i]
			tmpSlotVal <- objEqcReader@lsEqcSlotsOut[[i]]
			
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
validSPLOT <- function(objSPLOT) {

	if(!(objSPLOT@strMode == "singleplot" | objSPLOT@strMode == "subplot"))
		stop(paste("EASY ERROR:SPLOT\n Wrong strMode defined, PLease use 'singleplot' or ''subplot' \n !!!", sep=""))
	
	if(!(objSPLOT@numDefaultSymbol%in%c(0:25)))
		stop(paste("EASY ERROR:SPLOT\n Wrong numDefaultSymbol defined. PLease use integer value between 0 and 25.", sep=""))

	if(objSPLOT@numDefaultCex <= 0)
		stop(paste("EASY ERROR:SPLOT\n Wrong numDefaultCex defined. PLease use positive numeric value.", sep=""))
	
	strColour = objSPLOT@strDefaultColour
	isOkColors = strColour %in% colors()
	
	isOkColorHex =  nchar(strColour)==7 & substring(strColour,1,1)=="#"
	sTemp = substring(strColour,2,7) 
	isOkColorHex = isOkColorHex & all(strsplit(sTemp,"",fixed=TRUE)[[1]]%in%c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"))
	
	if(!isOkColors & !isOkColorHex)
		stop(paste("EASY ERROR:SPLOT\n Wrong strDefaultColour defined. PLease use color from R colors() or hexadecimal nomenclature, eg #FFFFFF.", sep=""))
	
	if(objSPLOT@astrColour[1] != "") {
		for(i in 1:length(objSPLOT@astrColour)) {
			
			strColour = objSPLOT@astrColour[i]
			isOkColors = strColour %in% colors()
			
			isOkColorHex =  nchar(strColour)==7 & substring(strColour,1,1)=="#"
			sTemp = substring(strColour,2,7) 
			isOkColorHex = isOkColorHex & all(strsplit(sTemp,"",fixed=TRUE)[[1]]%in%c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"))
			
			if(!isOkColors & !isOkColorHex)
				stop(paste("EASY ERROR:SPLOT\n Wrong astrColour defined. PLease use color from R colors() or hexadecimal nomenclature, eg #FFFFFF.", sep=""))
		}
	}
	
	if(objSPLOT@anumSymbol[1] != -1) {
		for(i in 1:length(objSPLOT@anumSymbol)) {
			
			numSymbol = objSPLOT@anumSymbol[i]
			if(!(numSymbol%in%c(0:25)))
				stop(paste("EASY ERROR:SPLOT\n Wrong anumSymbol defined. PLease use integer value between 0 and 25.", sep=""))
				
		}
	}
	
	if(objSPLOT@anumCex[1] != -1) {
		for(i in 1:length(objSPLOT@anumCex)) {
			
			numCex = objSPLOT@anumCex[i]
			if(numCex <= 0)
				stop(paste("EASY ERROR:SPLOT\n Wrong anumCex defined. PLease use positive numeric value.", sep=""))			
		}
	}
	
	if(!(objSPLOT@strFormat == "png" | objSPLOT@strFormat == "pdf"))
		stop(paste("EASY ERROR:SPLOT\n Wrong strFormat defined, PLease use 'png' or 'pdf' \n !!!", sep=""))
	
	if(objSPLOT@numCexAxis <= 0) 
		stop(paste("EASY ERROR:SPLOT\n Wrong numCexAxis defined, PLease use numCexAxis > 0 \n !!!", sep=""))
		
	if(objSPLOT@numCexLab <= 0) 
		stop(paste("EASY ERROR:SPLOT\n Wrong numCexLab defined, PLease use numCexLab > 0 \n !!!", sep=""))				
	
	if(length(objSPLOT@anumParMar) != 4) 
		stop(paste("EASY ERROR:SPLOT\n Wrong anumParMar defined, Please use anumParMar of length 4 (Default=5.1;4.1;4.1;2.1) \n !!!", sep=""))				
	
	if(length(objSPLOT@anumParMgp) != 3) 
		stop(paste("EASY ERROR:SPLOT\n Wrong anumParMgp defined, Please use anumParMgp of length 3 (Default=3;1;0) \n !!!", sep=""))				
	
	if(!(objSPLOT@strParBty%in%c("o","l","7","c","u","]")))
		stop(paste("EASY ERROR:SPLOT\n Wrong strParBty defined, Please use 'o','l','7','c','u' or ']' \n !!!", sep=""))				
	
	return(TRUE)
}

SPLOT.GWADATA.valid <- function(objSPLOT, objGWA) {
	
	objRCD 	<- RCD(objSPLOT@rcdSPlotX)
	x 		<- RCD.eval(objRCD, objGWA)
	
	#if(all(is.na(x))) 
	#	stop(paste("EASY ERROR:SPLOT\n rcdSPlotX \n",objSPLOT@rcdSPlotX," is not valid with file\n",objGWA@fileInShortName,"\n All returned values are NA. PLease check rcd or file.", sep=""))
	
	objRCD 	<- RCD(objSPLOT@rcdSPlotY)
	y 		<- RCD.eval(objRCD, objGWA)
	
	#if(all(is.na(x))) 
	#	stop(paste("EASY ERROR:SPLOT\n rcdSPlotY \n",objSPLOT@rcdSPlotY," is not valid with file\n",objGWA@fileInShortName,"\n All returned values are NA. PLease check rcd or file.", sep=""))
		
	if(length(x) != length(y))
		stop(paste("EASY ERROR:SPLOT\n rcdSPlotX and rcdSPlotY are not valid with file\n",objGWA@fileInShortName,"\n The two arrays returned are of different size. Please check rcds or file.", sep=""))
		
	if(class(x) != "numeric" & class(x) != "integer" & class(x) != "double")
		stop(paste("EASY ERROR:SPLOT\n rcdSPlotX gives non-numeric values for \n",objGWA@fileInShortName,"\n Please check rcdSPlotX or file.", sep=""))
		
	if(class(y) != "numeric" & class(y) != "integer" & class(y) != "double")
		stop(paste("EASY ERROR:SPLOT\n rcdSPlotY gives non-numeric values for \n",objGWA@fileInShortName,"\n Please check rcdSPlotY or file.", sep=""))
	
	
	if(objSPLOT@blnPlotCI) {
		objRCD 	<- RCD(objSPLOT@rcdCIlengthX)
		cix		<- RCD.eval(objRCD, objGWA)
		objRCD 	<- RCD(objSPLOT@rcdCIlengthY)
		ciy 	<- RCD.eval(objRCD, objGWA)
	}
	
}

SPLOT.run <- function(objSPLOT, objGWA) {
	
	tblPlot					=	objGWA@tblGWA
	
	rcdSPlotX				=	objSPLOT@rcdSPlotX
	rcdSPlotY				=	objSPLOT@rcdSPlotY
	strDefaultColour		=	objSPLOT@strDefaultColour
	numDefaultSymbol		=	objSPLOT@numDefaultSymbol
	numDefaultCex			=	objSPLOT@numDefaultCex
	arcdColourCrit			=	objSPLOT@arcdColourCrit
	astrColour				=	objSPLOT@astrColour
	astrColourLegend		=	objSPLOT@astrColourLegend
	arcdSymbolCrit			=	objSPLOT@arcdSymbolCrit
	anumSymbol				=	objSPLOT@anumSymbol
	astrSymbolLegend		=	objSPLOT@astrSymbolLegend
	arcdCexCrit				=	objSPLOT@arcdCexCrit
	anumCex					=	objSPLOT@anumCex
	#astrSPlotCexLegend		=	objSPLOT@astrSPlotCexLegend
	strAxes					=	objSPLOT@strAxes
	strXlab					=	objSPLOT@strXlab
	strYlab					=	objSPLOT@strYlab
	strTitle				=	objSPLOT@strTitle
	arcdAdd2Plot			=	objSPLOT@arcdAdd2Plot
	blnLegend				=	objSPLOT@blnLegend
	rcdExclude				=	objSPLOT@rcdExclude
	numCexAxis				=	objSPLOT@numCexAxis
	numCexLab				=	objSPLOT@numCexLab
	blnGrid					=	objSPLOT@blnGrid
				
	blnPlotCI 				=	objSPLOT@blnPlotCI
	rcdCIlengthX 			=	objSPLOT@rcdCIlengthX
	rcdCIlengthY 			=	objSPLOT@rcdCIlengthY
	strCIColour				=	objSPLOT@strCIColour
	
	blnAxesLines 			=	objSPLOT@blnAxesLines
	
	#if(objSPLOT@strMode == "singleplot") strTitle = objGWA@fileInShortName
	if(strTitle == "") strTitle = objGWA@fileInShortName
	if(strTitle == "omit") strTitle = ""
	
	### Set Plot colours, symbols and symbol sizes
	
	plot.colour = rep(strDefaultColour, dim(tblPlot)[1])
	plot.colour.order = rep(1, dim(tblPlot)[1])
	plot.symbol = rep(numDefaultSymbol, dim(tblPlot)[1])
	plot.cex 	= rep(numDefaultCex, dim(tblPlot)[1])
	
	if(arcdColourCrit[1] != "") {
		for(iColourCrit in 1:length(arcdColourCrit)) {
			
			rcdColourCrit 		<- arcdColourCrit[iColourCrit]
			strColour 			<- astrColour[iColourCrit]
			
			objRCD 	<- RCD(rcdColourCrit)
			out 	<- RCD.eval(objRCD, objGWA)
			
			plot.colour[out] <- strColour
			plot.colour.order[out] <- iColourCrit + 1
		}
	} 
	
	if(arcdSymbolCrit[1] != "") {
		for(iSymbolCrit in 1:length(arcdSymbolCrit)) {
			
			rcdSymbolCrit 		<- arcdSymbolCrit[iSymbolCrit]
			numSymbol 			<- anumSymbol[iSymbolCrit]
				
			objRCD 	<- RCD(rcdSymbolCrit)
			out 	<- RCD.eval(objRCD, objGWA)
			
			plot.symbol[out] <- numSymbol
			
		}
	}
	
	if(arcdCexCrit[1] != "") {
		for(iCexCrit in 1:length(arcdCexCrit)) {
			
			rcdCexCrit 		<- arcdCexCrit[iCexCrit]
			numCex 			<- anumCex[iCexCrit]
		
			objRCD 	<- RCD(rcdCexCrit)
			out 	<- RCD.eval(objRCD, objGWA)
			
			plot.cex[out] <- numCex

		}
	}
	
	### Get x,y 
	
	objRCD 	<- RCD(rcdSPlotX)
	x 		<- RCD.eval(objRCD, objGWA)
	objRCD 	<- RCD(rcdSPlotY)
	y 		<- RCD.eval(objRCD, objGWA)
	
	if(blnPlotCI) {
		objRCD 		<- RCD(rcdCIlengthX)
		ci_x 		<- RCD.eval(objRCD, objGWA)
		objRCD 		<- RCD(rcdCIlengthY)
		ci_y 		<- RCD.eval(objRCD, objGWA)
	}
	
	### Exclude
	if(rcdExclude!="") {
		objRCD 	<- RCD(rcdExclude)
		out 	<- RCD.eval(objRCD, objGWA)
		if(any(out)) {
			x 	= x[which(!out)]
			y 	= y[which(!out)]
			plot.colour 		= plot.colour[which(!out)]
			plot.colour.order 	= plot.colour.order[which(!out)]
			plot.symbol 		= plot.symbol[which(!out)]
			tblPlot 			= tblPlot[which(!out),]
			if(blnPlotCI) {
				ci_x = ci_x[which(!out)]
				ci_y = ci_y[which(!out)]
			}
		}
	}

	
	### Clean x, y
	ablnIsNA = is.na(x) | is.na(y)
	
	x 					= x[which(!ablnIsNA)]
	y 					= y[which(!ablnIsNA)]
	plot.colour 		= plot.colour[which(!ablnIsNA)]
	plot.colour.order 	= plot.colour.order[which(!ablnIsNA)]
	plot.symbol 		= plot.symbol[which(!ablnIsNA)]
	tblPlot 			= tblPlot[which(!ablnIsNA),]
	if(blnPlotCI) {
		ci_x = ci_x[which(!ablnIsNA)]
		ci_y = ci_y[which(!ablnIsNA)]
	}
			
	## Remove infinity values
	isInfX = abs(x) == Inf
	if(any(isInfX)) warning(paste("EASY WARNING:\n There are ",length(which(isInfX)), " SNPs with Inf values on the x-axis being excluded for from the SPLOT for the rcd \n", rcdSPlotX, " \n and file \n ", objGWA@fileIn ,sep="" ))
	isInfY = abs(y) == Inf
	if(any(isInfY)) warning(paste("EASY WARNING:\n There are ",length(which(isInfY)), " SNPs with Inf values on the y-axis being excluded for from the SPLOT for the rcd \n", rcdSPlotY, " \n and file \n ", objGWA@fileIn ,sep="" ))
	isExclude = isInfX | isInfY
	
	x 	= x[which(!isExclude)]
	y 	= y[which(!isExclude)]
	plot.colour 	= plot.colour[which(!isExclude)]
	plot.colour.order = plot.colour.order[which(!isExclude)]
	plot.symbol 	= plot.symbol[which(!isExclude)]
	tblPlot 	= tblPlot[which(!isExclude),]
	if(blnPlotCI) {
		ci_x = ci_x[which(!isExclude)]
		ci_y = ci_y[which(!isExclude)]
	}
	
	if(dim(tblPlot)[1]==0) {
		warning(paste("EASY WARNING:\n There are no SNPs with non-missing values for the SPLOT of rcdSPlotX '", rcdSPlotX, "' rcdSPlotY '", rcdSPlotY, "' and the file \n ", objGWA@fileIn ,sep="" ))
		plot(0,0,col="white")
		text(0,0,"Empty data!")
		return()
	}
	
	## Change plot order according to defined colouring criterion
	iSort = order(plot.colour.order)
	x = x[iSort]
	y = y[iSort]
	plot.colour 	= plot.colour[iSort]
	plot.colour.order = plot.colour.order[iSort]
	plot.symbol 	= plot.symbol[iSort]
	tblPlot 	= tblPlot[iSort,]
	if(blnPlotCI) {
		ci_x = ci_x[iSort]
		ci_y = ci_y[iSort]
	}
	
	### Get axis lims strAxes
	
	plot.xlims <- plot.ylims <- rep(NA,2)

	if(strAxes == "equal")  {
		plot.xlims[1] <- plot.ylims[1] <- min(min(x), min(y))
		plot.xlims[2] <- plot.ylims[2] <- max(max(x), max(y))
	} else if(strAxes == "4quad") {
		plot.xlims[1] <- -max(abs(x))
		plot.xlims[2] <- max(abs(x))
		plot.ylims[1] <- -max(abs(y))
		plot.ylims[2] <- max(abs(y))
	} else if(strAxes == "4quadequal") {
		plot.xlims[1] <- plot.ylims[1] <- -max(max(abs(x)), max(abs(y)))
		plot.xlims[2] <- plot.ylims[2] <- max(max(abs(x)), max(abs(y)))
	} else if(strAxes == "zeroequal") {
		plot.xlims[1] <- plot.ylims[1] <- 0
		plot.xlims[2] <- plot.ylims[2] <- max(max(abs(x)), max(abs(y)))
	} else if(length(grep("lim",strAxes)) > 0)  {
		## strSPAxes = "lim(x0,x1,y0,y1)"
		str_lim = substr(strAxes,5,nchar(strAxes)-1)
		## str_lim = "x0,x1,y0,y1"
		arr_str_lim = strsplit(str_lim,",")[[1]]

		plot.xlims[1] = ifelse(arr_str_lim[1]!="NULL",as.numeric(arr_str_lim[1]),min(x))
		plot.xlims[2] = ifelse(arr_str_lim[2]!="NULL",as.numeric(arr_str_lim[2]),max(x))
		plot.ylims[1] = ifelse(arr_str_lim[3]!="NULL",as.numeric(arr_str_lim[3]),min(y))
		plot.ylims[2] = ifelse(arr_str_lim[4]!="NULL",as.numeric(arr_str_lim[4]),max(y))
	} else {
		#plot.xlims = plot.ylims = NULL
		plot.xlims = c(min(x),max(x))
		plot.ylims = c(min(y),max(y))
	}
	
	###
	plot.xlab = ifelse(strXlab == "", rcdSPlotX, strXlab)
	plot.ylab = ifelse(strYlab == "", rcdSPlotY, strYlab)
	plot.title = strTitle
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
	
	# print(plot.xlims)
	# print(plot.ylims)
	
	plot( NULL, NULL,
			xlim = plot.xlims, 
			ylim = plot.ylims,
			cex.main=plot.cex.main,
			cex.axis=numCexAxis,
			cex.lab=numCexLab,
			main=plot.title,
			xlab=plot.xlab,
			ylab=plot.ylab,
			las=objSPLOT@numParLas
			#lwd=3
	)
	
	### Start plotting
	# plot( x, y,
			# xlim = plot.xlims, 
			# ylim = plot.ylims,
			# col=plot.colour,
			# pch=plot.symbol,
			# cex=plot.cex,
			# cex.main=plot.cex.main,
			# cex.axis=numCexAxis,
			# cex.lab=numCexLab,
			# main=plot.title,
			# xlab=plot.xlab,
			# ylab=plot.ylab
			# #lwd=3
	# )
	if(blnPlotCI) {
		arrows(x0=x,x1=x,y0=y,y1=y+ci_y,angle=90,length=0.05,col=strCIColour)
		arrows(x0=x,x1=x,y0=y,y1=y-ci_y,angle=90,length=0.05,col=strCIColour)
		arrows(x0=x,x1=x+ci_x,y0=y,y1=y,angle=90,length=0.05,col=strCIColour)
		arrows(x0=x,x1=x-ci_x,y0=y,y1=y,angle=90,length=0.05,col=strCIColour)
		
		# overwrite arows:
		points( x, y, col=plot.colour,pch=plot.symbol,cex=plot.cex)
		
		# for (i in 1:length(x)) {
			# polygon(x=c(x[i]-ci_x[i],x[i],x[i]+ci_x[i],x[i]),
					# y=c(y[i],y[i]+ci_y[i],y[i],y[i]-ci_y[i]),
					# col=plot.colour[i])
		# }		
	} else {
		points(x,y,col=plot.colour,pch=plot.symbol,cex=plot.cex)
	}
	
	# Achsenlinien:
	#abline(a=0,b=0,col="gray3",lty=6)
	#abline(v=0,col="gray3",lty=6)
	if(blnGrid) grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")
	#box(lty = 1)
	
	if(blnAxesLines) {
		abline(v=0)
		abline(h=0)
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
	
	#if(blnLegend) legend(x="bottomright",inset=0.005,legend=astrColourLegend,col=astrColour,pch = anumLegendPch,bg="skyblue")
	#par(xpd=TRUE)
	if(blnLegend) {
		
		ausr = par()$usr
		xleg = ausr[2]
		yleg = ausr[4]
	
		#legend(x=xmax, y=ymax,inset=0.005,legend=astrColourLegend,col=astrColour,pch = 20, bg="skyblue",xpd=TRUE)
		#legend(x=xleg, y=yleg,inset=0.005,legend=astrColourLegend,col=astrColour,pch = 20, bg="skyblue",xpd=TRUE)
		
		legendcol <- legendsym <- c()
		legendtext <- c()
		if(length(astrColourLegend)>0 & all(!is.na(astrColourLegend))) {
			legendcol = c(legendcol,astrColour)
			legendsym = c(legendsym, rep(20, length(astrColour)))
			#legendcex = c(legendcex, rep(1, length(astrColour)))
			
			legendtext = c(legendtext, astrColourLegend)
		}
		if(length(astrSymbolLegend)>0 & all(!is.na(astrSymbolLegend))) {
			legendcol = c(legendcol,rep("black", length(anumSymbol)))
			legendsym = c(legendsym, anumSymbol)
			#legendcex = c(legendcex, rep(1, length(anumSymbol)))
			
			legendtext = c(legendtext, astrSymbolLegend)
		}
		# if(length(astrSPlotCexLegend)>0 & all(!is.na(astrSPlotCexLegend))) {
			# legendcol = c(legendcol,rep("black", length(anumSymbol)))
			# legendsym = c(legendsym, rep(20, length(astrColour)))
			# legendcex = c(legendcex, anumCex)
			
			# legendtext = c(legendtext, astrSPlotCexLegend)
		# }
		if(length(legendtext) > 0) {
			legend(	x=xleg, 
					y=yleg,
					inset=0.005,
					legend=legendtext,
					col=legendcol,
					pch = legendsym, 
					#cex = legendcex,
					bg="skyblue",
					xpd=TRUE
					)
		}
	}
}


SPLOT <- function(strEqcCommand){ 
	## Wrapper for class definition
	
	SPLOTout <- setSPLOT(new("SPLOT", strEqcCommand = strEqcCommand))
	validSPLOT(SPLOTout)
	#ADDCOLout.valid <- validADDCOL(ADDCOLout)
	return(SPLOTout)
	#validECF(ECFout)
	#return(ECFout)
	
	## Identical:
	# ECFin <- new("ECF5", fileECF = fileECFIn) 
	# ECFout <- setECF5(ECFin)
	# return(ECFout)
}


