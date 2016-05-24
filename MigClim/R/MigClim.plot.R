#
# The R functions for the MigClim package.
# Robin Engler & Wim Hordijk
# MigClim.plot: Last modified: 26 April 2012
#
# MigClim.plot: displays ascii files produced by the MigClim.migrate() function.
#
MigClim.plot <- function(asciiFile, outDir="", fileFormat="jpeg", fullOutput=FALSE){
	
	### Verify user input
	if(substr(asciiFile, nchar(asciiFile)-3, nchar(asciiFile))!=".asc") asciiFile <- paste(asciiFile, ".asc", sep="")
	if(fileFormat!="jpeg" & fileFormat!="png" & fileFormat!="inR") stop("Input error: 'fileFormat' must be one of 'jpeg' or 'png'.\n")
	if(!file.exists(asciiFile)) stop("Input error: ", asciiFile, " could not be found.\n")
	if(outDir!="") if(!file.exists(outDir)) stop("Input error: 'outDir' directory could not be found.\n")
	
	### load raster library (this is no longer needed, R does this automatically).
	#if(require(raster, quietly=T)==F) stop("This function requires the 'raster' package. Please install 'raster' on your computer and try again.")
	
	### get root name and output directory
	baseName <- substr(basename(asciiFile), 1, nchar(basename(asciiFile))-11)  # The "baseName" of the MigClim simulation.
	if(outDir=="" & dirname(asciiFile)!=".") outDir <- dirname(asciiFile)
	outDir <- paste(outDir,"/",sep="")
	inDir <- paste(dirname(asciiFile),"/",sep="")
	if(inDir=="./") inDir <- ""
	
	### get all files
	if(fullOutput){
		stepName <- paste(baseName, "_step_", sep="")
		fileList <- list.files(dirname(asciiFile))
		fileList <- fileList[which(substr(fileList, 1, nchar(stepName))==stepName)]
		fileList <- c(basename(asciiFile), fileList)
		rm(stepName)
	} else fileList <- basename(asciiFile)
	fileList <- substr(fileList,1,nchar(fileList)-4) # remove ".asc" extension.
	
	for(fileName in fileList){
		
		#fileName <- fileList[1]
		cat("plotting data for", paste(fileName,".asc",sep=""), "\n")
		Rst <- raster(paste(inDir,fileName,".asc",sep=""))
		
		### Create color ramp
		rstVals <- sort(raster::unique(Rst))
		negativeNb <- length(which(rstVals<0))
		positiveNb <- length(which(rstVals>1 & rstVals<30000))
		zeroExists <- any(rstVals==0)
		oneExists <- any(rstVals==1)
		unilimtedExists <- any(rstVals==30000)
		Colors <- rep("yellow", negativeNb)
		if(zeroExists) Colors <- c(Colors, "grey94")
		if(oneExists) Colors <- c(Colors, "black")
		Colors <- c(Colors, rainbow(positiveNb, start=0, end=0.4))
		if(unilimtedExists) Colors <- c(Colors, "pink")
		
		### Create a graphic window that has the correct ratio (to match the map)
		#x11(width=7, height=7*((ymax(Rst)-ymin(Rst))/(xmax(Rst)-xmin(Rst))))
		#image(Rst, col=Colors, breaks=c(min(rstVals)-1,rstVals))
		#savePlot(filename= "mcTest_raster", type = "jpeg", device = dev.cur(), restoreConsole = TRUE)
		
		### Plot
		if(fileFormat=="jpeg"){
			jpeg(filename=paste(outDir,fileName,".jpg",sep=""), width = 2000, height = 2000*((ymax(Rst)-ymin(Rst))/(xmax(Rst)-xmin(Rst))), quality=90, res=300)
			image(Rst, col=Colors, breaks=c(min(rstVals)-1,rstVals))
			dev.off()
		}
		if(fileFormat=="png"){
			png(filename=paste(outDir,fileName,".png",sep=""), width = 2000, height = 2000*((ymax(Rst)-ymin(Rst))/(xmax(Rst)-xmin(Rst))), res=300)
			image(Rst, col=Colors, breaks=c(min(rstVals)-1,rstVals))
			dev.off()
		}
		if(fileFormat=="inR"){
			dev.new(width=7, height=7*((ymax(Rst)-ymin(Rst))/(xmax(Rst)-xmin(Rst))))
			#image(Rst, col=Colors, breaks=c(min(rstVals)-1,rstVals)) # this function does not work 100% anymore (axes are missing)
			plot(Rst, col=Colors, breaks=c(min(rstVals)-1,rstVals), legend=FALSE)
		}
		
		### Release memory
		rm(Rst,rstVals,negativeNb,positiveNb,zeroExists,oneExists,unilimtedExists,Colors)
	}
	
}

