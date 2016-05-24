setClass("RPLOT",
	representation = representation(
						strEqcCommand		=	"character",
						rcdRPlotX			=	"character",
						rcdRPlotY			=	"character",
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
						#astrRPLOTCexLegend	=	"character",
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
						strPlotName			= 	"character"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						rcdRPlotX			=	"",
						rcdRPlotY			=	"",
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
						#astrRPLOTCexLegend	=	"",
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
						strPlotName			= 	"rp"
						)
	#contains = c("EcfReader")
)

setGeneric("setRPLOT", function(object) standardGeneric("setRPLOT"))
setMethod("setRPLOT", signature = (object = "RPLOT"), function(object) {
	
	aEqcSlotNamesIn = c("rcdRPlotX", 
						"rcdRPlotY", 
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
						#"astrRPLOTCexLegend",
						"strAxes",
						"strXlab",
						"strYlab",
						"strTitle",
						"blnLegend",
						"arcdAdd2Plot",
						#"strMode", ## May not be changed for report
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
						"strPlotName"
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
validRPLOT <- function(objRPLOT) {

	if(!(objRPLOT@strMode == "singleplot" | objRPLOT@strMode == "subplot"))
		stop(paste("EASY ERROR:RPLOT\n Wrong strMode defined, PLease use 'singleplot' or ''subplot' \n !!!", sep=""))
	
	if(!(objRPLOT@numDefaultSymbol%in%c(0:25)))
		stop(paste("EASY ERROR:RPLOT\n Wrong numDefaultSymbol defined. PLease use integer value between 0 and 25.", sep=""))

	if(objRPLOT@numDefaultCex <= 0)
		stop(paste("EASY ERROR:RPLOT\n Wrong numDefaultCex defined. PLease use positive numeric value.", sep=""))
	
	strColour = objRPLOT@strDefaultColour
	isOkColors = strColour %in% colors()
	
	isOkColorHex =  nchar(strColour)==7 & substring(strColour,1,1)=="#"
	sTemp = substring(strColour,2,7) 
	isOkColorHex = isOkColorHex & all(strsplit(sTemp,"",fixed=TRUE)[[1]]%in%c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"))
	
	if(!isOkColors & !isOkColorHex)
		stop(paste("EASY ERROR:RPLOT\n Wrong strDefaultColour defined. PLease use color from R colors() or hexadecimal nomenclature, eg #FFFFFF.", sep=""))
	
	if(objRPLOT@astrColour[1] != "") {
		for(i in 1:length(objRPLOT@astrColour)) {
			
			strColour = objRPLOT@astrColour[i]
			isOkColors = strColour %in% colors()
			
			isOkColorHex =  nchar(strColour)==7 & substring(strColour,1,1)=="#"
			sTemp = substring(strColour,2,7) 
			isOkColorHex = isOkColorHex & all(strsplit(sTemp,"",fixed=TRUE)[[1]]%in%c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"))
			
			if(!isOkColors & !isOkColorHex)
				stop(paste("EASY ERROR:RPLOT\n Wrong astrColour defined. PLease use color from R colors() or hexadecimal nomenclature, eg #FFFFFF.", sep=""))
		}
	}
	
	if(objRPLOT@anumSymbol[1] != -1) {
		for(i in 1:length(objRPLOT@anumSymbol)) {
			
			numSymbol = objRPLOT@anumSymbol[i]
			if(!(numSymbol%in%c(0:25)))
				stop(paste("EASY ERROR:RPLOT\n Wrong anumSymbol defined. PLease use integer value between 0 and 25.", sep=""))
				
		}
	}
	
	if(objRPLOT@anumCex[1] != -1) {
		for(i in 1:length(objRPLOT@anumCex)) {
			
			numCex = objRPLOT@anumCex[i]
			if(numCex <= 0)
				stop(paste("EASY ERROR:RPLOT\n Wrong anumCex defined. PLease use positive numeric value.", sep=""))			
		}
	}
	
	if(!(objRPLOT@strFormat == "png" | objRPLOT@strFormat == "pdf"))
		stop(paste("EASY ERROR:RPLOT\n Wrong strFormat defined, PLease use 'png' or 'pdf' \n !!!", sep=""))
	
	if(objRPLOT@numCexAxis <= 0) 
		stop(paste("EASY ERROR:RPLOT\n Wrong numCexAxis defined, PLease use numCexAxis > 0 \n !!!", sep=""))
		
	if(objRPLOT@numCexLab <= 0) 
		stop(paste("EASY ERROR:RPLOT\n Wrong numCexLab defined, PLease use numCexLab > 0 \n !!!", sep=""))				
	
	if(length(objRPLOT@anumParMar) != 4) 
		stop(paste("EASY ERROR:RPLOT\n Wrong anumParMar defined, Please use anumParMar of length 4 (Default=5.1;4.1;4.1;2.1) \n !!!", sep=""))				
	
	if(length(objRPLOT@anumParMgp) != 3) 
		stop(paste("EASY ERROR:RPLOT\n Wrong anumParMgp defined, Please use anumParMgp of length 3 (Default=3;1;0) \n !!!", sep=""))				
	
	if(!(objRPLOT@strParBty%in%c("o","l","7","c","u","]")))
		stop(paste("EASY ERROR:RPLOT\n Wrong strParBty defined, Please use 'o','l','7','c','u' or ']' \n !!!", sep=""))				
	
	return(TRUE)
}

RPLOT.REPORT.valid <- function(objRPLOT, objREPORT) {
	
	objRCD 	<- RCD(objRPLOT@rcdRPlotX)
	x 		<- RCD.eval.report(objRCD, objREPORT)
	if(all(is.na(x)&class(x)=="character")) x <- as.numeric(x)
	
	objRCD 	<- RCD(objRPLOT@rcdRPlotY)
	y 		<- RCD.eval.report(objRCD, objREPORT)
	if(all(is.na(y)&class(y)=="character")) y <- as.numeric(y)

	#if(all(is.na(x))) 
	#	stop(paste("EASY ERROR:RPLOT\n rcdRPlotY \n",objRPLOT@rcdRPlotY," is not valid with file\n",objREPORT@fileInShortName,"\n All returned values are NA. PLease check rcd or file.", sep=""))
		
	if(length(x) != length(y))
		stop(paste("EASY ERROR:RPLOT\n rcdRPlotX and rcdRPlotY are not valid with file\n",objREPORT@fileReport,"\n The two arrays returned are of different size. Please check rcds or file.", sep=""))
		
	if(class(x) != "numeric" & class(x) != "integer" & class(x) != "double")
		stop(paste("EASY ERROR:RPLOT\n rcdRPlotX gives non-numeric values for \n",objREPORT@fileReport,"\n Please check rcdRPlotX or file.", sep=""))
		
	if(class(y) != "numeric" & class(y) != "integer" & class(y) != "double")
		stop(paste("EASY ERROR:RPLOT\n rcdRPlotY gives non-numeric values for \n",objREPORT@fileReport,"\n Please check rcdRPlotY or file.", sep=""))
	

}

RPLOT.run <- function(objRPLOT, objREPORT) {
	
	tblPlot					=	objREPORT@tblReport
	
	rcdRPlotX				=	objRPLOT@rcdRPlotX
	rcdRPlotY				=	objRPLOT@rcdRPlotY
	strDefaultColour		=	objRPLOT@strDefaultColour
	numDefaultSymbol		=	objRPLOT@numDefaultSymbol
	numDefaultCex			=	objRPLOT@numDefaultCex
	arcdColourCrit			=	objRPLOT@arcdColourCrit
	astrColour				=	objRPLOT@astrColour
	astrColourLegend		=	objRPLOT@astrColourLegend
	arcdSymbolCrit			=	objRPLOT@arcdSymbolCrit
	anumSymbol				=	objRPLOT@anumSymbol
	astrSymbolLegend		=	objRPLOT@astrSymbolLegend
	arcdCexCrit				=	objRPLOT@arcdCexCrit
	anumCex					=	objRPLOT@anumCex
	#astrRPLOTCexLegend		=	objRPLOT@astrRPLOTCexLegend
	strAxes					=	objRPLOT@strAxes
	strXlab					=	objRPLOT@strXlab
	strYlab					=	objRPLOT@strYlab
	strTitle				=	objRPLOT@strTitle
	arcdAdd2Plot			=	objRPLOT@arcdAdd2Plot
	blnLegend				=	objRPLOT@blnLegend
	rcdExclude				=	objRPLOT@rcdExclude
	numCexAxis				=	objRPLOT@numCexAxis
	numCexLab				=	objRPLOT@numCexLab
	blnGrid					=	objRPLOT@blnGrid
						
	#if(objRPLOT@strMode == "singleplot") strTitle = objREPORT@fileInShortName
	if(strTitle == "") strTitle = strsplit(objREPORT@fileReportBody,"/",fixed=TRUE)[[1]][length(strsplit(objREPORT@fileReportBody,"/",fixed=TRUE)[[1]])]
	
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
			out 	<- RCD.eval.report(objRCD, objREPORT)
			
			plot.colour[out] <- strColour
			plot.colour.order[out] <- iColourCrit + 1
		}
	} 
	
	if(arcdSymbolCrit[1] != "") {
		for(iSymbolCrit in 1:length(arcdSymbolCrit)) {
			
			rcdSymbolCrit 		<- arcdSymbolCrit[iSymbolCrit]
			numSymbol 			<- anumSymbol[iSymbolCrit]
				
			objRCD 	<- RCD(rcdSymbolCrit)
			out 	<- RCD.eval.report(objRCD, objREPORT)
			
			plot.symbol[out] <- numSymbol
			
		}
	}
	
	if(arcdCexCrit[1] != "") {
		for(iCexCrit in 1:length(arcdCexCrit)) {
			
			rcdCexCrit 		<- arcdCexCrit[iCexCrit]
			numCex 			<- anumCex[iCexCrit]
		
			objRCD 	<- RCD(rcdCexCrit)
			out 	<- RCD.eval.report(objRCD, objREPORT)
			
			plot.cex[out] <- numCex

		}
	}
	
	### Get x,y 
	
	objRCD 	<- RCD(rcdRPlotX)
	x 		<- RCD.eval.report(objRCD, objREPORT)
	objRCD 	<- RCD(rcdRPlotY)
	y 		<- RCD.eval.report(objRCD, objREPORT)
	
	### Exclude
	if(rcdExclude!="") {
		objRCD 	<- RCD(rcdExclude)
		out 	<- RCD.eval.report(objRCD, objREPORT)
		if(any(out)) {
			x 	= x[which(!out)]
			y 	= y[which(!out)]
			plot.colour 		= plot.colour[which(!out)]
			plot.colour.order 	= plot.colour.order[which(!out)]
			plot.symbol 		= plot.symbol[which(!out)]
			tblPlot 			= tblPlot[which(!out),]
		}
		else warning(paste("EASY WARNING:\n No lines fullfilled the criterion \n",rcdExclude, "\n. in  arcdQQPlotExclude for file \n ", objREPORT@fileIn ,sep="" ))
	}

	
	### Clean x, y
	ablnIsNA = is.na(x) | is.na(y)
	
	x 			= x[which(!ablnIsNA)]
	y 			= y[which(!ablnIsNA)]
	
	plot.colour = plot.colour[which(!ablnIsNA)]
	plot.colour.order = plot.colour.order[which(!ablnIsNA)]
	plot.symbol = plot.symbol[which(!ablnIsNA)]
	tblPlot 	= tblPlot[which(!ablnIsNA),]
	
	
	## Remove infinity values
	isInfX = abs(x) == Inf
	if(any(isInfX)) warning(paste("EASY WARNING:\n There are ",length(which(isInfX)), " lines with Inf values on the x-axis being excluded for from the RPLOT for the rcd \n", rcdRPlotX, sep="" ))
	isInfY = abs(y) == Inf
	if(any(isInfY)) warning(paste("EASY WARNING:\n There are ",length(which(isInfY)), " lines with Inf values on the y-axis being excluded for from the RPLOT for the rcd \n", rcdRPlotY,sep="" ))
	isExclude = isInfX | isInfY
	
	x 	= x[which(!isExclude)]
	y 	= y[which(!isExclude)]
	plot.colour 	= plot.colour[which(!isExclude)]
	plot.colour.order = plot.colour.order[which(!isExclude)]
	plot.symbol 	= plot.symbol[which(!isExclude)]
	tblPlot 	= tblPlot[which(!isExclude),]
	
	if(dim(tblPlot)[1]==0) {
		warning(paste("EASY WARNING:\n There are no SNPs with non-missing values for the RPLOT of rcdRPlotX '", rcdRPlotX, "' and rcdRPlotY '", rcdRPlotY, "' !",sep="" ))
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
		plot.xlims = plot.ylims = NULL
	}
	
	###
	plot.xlab = ifelse(strXlab == "", rcdRPlotX, strXlab)
	plot.ylab = ifelse(strYlab == "", rcdRPlotY, strYlab)
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
	
	### Start plotting
	plot( x, y,
			xlim = plot.xlims, 
			ylim = plot.ylims,
			col=plot.colour,
			pch=plot.symbol,
			cex=plot.cex,
			cex.main=plot.cex.main,
			cex.axis=numCexAxis,
			cex.lab=numCexLab,
			main=plot.title,
			xlab=plot.xlab,
			ylab=plot.ylab
			#lwd=3
	)
	
	# Achsenlinien:
	#abline(a=0,b=0,col="gray3",lty=6)
	#abline(v=0,col="gray3",lty=6)
	if(blnGrid) grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")
	#box(lty = 1)
	
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
		# if(length(astrRPLOTCexLegend)>0 & all(!is.na(astrRPLOTCexLegend))) {
			# legendcol = c(legendcol,rep("black", length(anumSymbol)))
			# legendsym = c(legendsym, rep(20, length(astrColour)))
			# legendcex = c(legendcex, anumCex)
			
			# legendtext = c(legendtext, astrRPLOTCexLegend)
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


RPLOT <- function(strEqcCommand){ 
	## Wrapper for class definition
	
	RPLOTout <- setRPLOT(new("RPLOT", strEqcCommand = strEqcCommand))
	validRPLOT(RPLOTout)
	#ADDCOLout.valid <- validADDCOL(ADDCOLout)
	return(RPLOTout)
	#validECF(ECFout)
	#return(ECFout)
	
	## Identical:
	# ECFin <- new("ECF5", fileECF = fileECFIn) 
	# ECFout <- setECF5(ECFin)
	# return(ECFout)
}


