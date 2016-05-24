fnOpenPlot <- function(objIn, objPlot, strSuffix="") {
	
	strPlotName = objPlot@strPlotName
	strFormat = objPlot@strFormat
	
	if(class(objIn) == "GWADATA") {
		pathOut <- objIn@pathOut
		if(strsplit(pathOut,"")[[1]][nchar(pathOut)] != "/") pathOut <- paste(pathOut,"/",sep="")	
		fileOutBody	<- paste(pathOut,objIn@fileInShortName,sep="")
		fileOut	<- paste(fileOutBody,".",strPlotName,".",strFormat,sep="")
		i = 1
		while(file.exists(fileOut)) {
			fileOut <- paste(fileOutBody,".",strPlotName,i,".",strFormat,sep="")
			i = i + 1
		}
	}
	if(class(objIn) == "REPORT") {
		fileOutBody <- objIn@fileReportBody
		#fileOut	<- paste(fileOutBody,".",strFormat,sep="")
		fileOut	<- paste(fileOutBody,".",strPlotName,".",strFormat,sep="")
		i = 1
		while(file.exists(fileOut)) {
			#fileOut <- paste(fileOutBody,i,".",strFormat,sep="")
			fileOut <- paste(fileOutBody,".",strPlotName,i,".",strFormat,sep="")
			i = i + 1
		}
	}
	
	if(strFormat == "png")
		CairoPNG(filename = fileOut,width=objPlot@numWidth,height=objPlot@numHeight)
	else 
		CairoPDF(file = fileOut, width=objPlot@numWidth,height=objPlot@numHeight)
	
	par(xpd=F, 
		bty=objPlot@strParBty, 
		mar=objPlot@anumParMar, 
		mgp=objPlot@anumParMgp
		)
	
	## GRAPHICAL PARAMETERS USING par():
	## mar = Margin sizes c(bottom, left, top, right) = c(5, 4, 4, 2) + 0.1 default
	##		in LINES
	## xpd = T: PLottet über Rand hinaus
	## mgp = Distance in LINES 
	##		of c(axis labels or titles from the axes; tick mark labels from the axes; tick mark symbols from the axes) 
	##		= default = c(3, 1, 0)
	## bty = Box around plot:  If âbtyâ is one of "o" (the default), "l", "7", "c", "u", or "]" the resulting
    ## A value of "n" suppresses the box.
	
}
#fnOpenMultiPlot <- function(fileOutBody, strSuffix, numFiles) {
fnOpenMultiPlot <- function(fileOutBody, objPlot, numFiles, afileOutMultiPlot) {
	
	strPlotName = objPlot@strPlotName
	strFormat = objPlot@strFormat
	
	fileOut	<- paste(fileOutBody,".multi.",strPlotName,".",strFormat,sep="")
	i = 1
	while(file.exists(fileOut) | fileOut%in%afileOutMultiPlot) {
		fileOut <- paste(fileOutBody,".multi.",strPlotName,i,".",strFormat,sep="")
		i = i + 1
	}
	
	numX = ceiling(sqrt(numFiles))
	numY = ceiling(numFiles/ceiling(sqrt(numFiles)))

	#png(file=fileOut,width=(250*numX),height=(250*numY))
	if(strFormat == "png")
		CairoPNG(filename=fileOut,width=(250*numX),height=(250*numY))
	else 
		CairoPDF(file = fileOut, width=(3*numX),height=(3*numY))
	
	par(mfrow=c(numY,numX)) 
	par(mar=c(3,3,3,1)) 
	par(mgp=c(2,0.7,0))
	
	return(fileOut)
}
fnOpenPng.MH <- function(objGWA, strSuffix) {
	
	pathOut <- objGWA@pathOut
	if(strsplit(pathOut,"")[[1]][nchar(pathOut)] != "/") pathOut <- paste(pathOut,"/",sep="")	
	fileOutBody	<- paste(pathOut,objGWA@fileInShortName,sep="")
	
	fileOut	<- paste(fileOutBody,".",strSuffix,".png",sep="")
	i = 1
	while(file.exists(fileOut)) {
		fileOut <- paste(fileOutBody,".",strSuffix,i,".png",sep="")
		i = i + 1
	}
	
	#png(file = fileOut,width=900,height=800)
	CairoPNG(filename=fileOut, width=1600, height=600, res=600)
	#CairoPNG(filename = fileOut,width=900,height=800)
	par(mar=c(7,6,4,10)) 
	#par(xpd=F, mar=par()$mar+c(0,0,0,8))	
	## xpd = T: PLottet über Rand hinaus
}
fnOpenPng.Miami <- function(objGWA, strSuffix) {
	
	pathOut <- objGWA@pathOut
	if(strsplit(pathOut,"")[[1]][nchar(pathOut)] != "/") pathOut <- paste(pathOut,"/",sep="")	
	fileOutBody	<- paste(pathOut,objGWA@fileInShortName,sep="")
	
	fileOut	<- paste(fileOutBody,".",strSuffix,".png",sep="")
	i = 1
	while(file.exists(fileOut)) {
		fileOut <- paste(fileOutBody,".",strSuffix,i,".png",sep="")
		i = i + 1
	}
	
	#png(file = fileOut,width=900,height=800)
	CairoPNG(filename=fileOut, width=1600, height=800, res=600)
	#CairoPNG(filename = fileOut,width=900,height=800)
	par(mar=c(7,6,4,10)) 
	#par(xpd=F, mar=par()$mar+c(0,0,0,8))	
	## xpd = T: PLottet über Rand hinaus
}

fnOpenPng.Boxplot <- function(fileOutBody, numFiles) {

	fileOut	<- paste(fileOutBody,".boxplot.png",sep="")
	i = 1
	while(file.exists(fileOut)) {
		fileOut <- paste(fileOutBody,".boxplot",i,".png",sep="")
		i = i + 1
	}
	
	#png(fileOut, height=(100+50*numFiles), width=1000)
	CairoPNG(filename=fileOut, height=(100+50*numFiles), width=1000)
	par(mar=c(8,25,2,1))
	
}
fnClosePlot <- function() {	
	dev.off()
}
fnAddPlot <- function(iMultiPlot) {
	dev.set(dev.list()[iMultiPlot])
}
fnOpenMultiPng <- function(fileOutBody, strSuffix, numFiles) {
	fileOut	<- paste(fileOutBody,".multi_",strSuffix,".png",sep="")
	i = 1
	while(file.exists(fileOut)) {
		fileOut <- paste(fileOutBody,".multi_",strSuffix,i,".png",sep="")
		i = i + 1
	}
	
	numX = ceiling(sqrt(numFiles))
	numY = ceiling(numFiles/ceiling(sqrt(numFiles)))

	#png(file=fileOut,width=(250*numX),height=(250*numY))
	CairoPNG(filename=fileOut,width=(250*numX),height=(250*numY))
	par(mfrow=c(numY,numX)) 
	par(mar=c(3,3,3,1)) 
	par(mgp=c(2,0.7,0))
}
fnAddMultiPng <- function(iMultiPlot) {
	dev.set(dev.list()[iMultiPlot])
}
fnCloseMultiPng <- function() {
	dev.off()
}