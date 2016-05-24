data.main <- function(extraction.done=TRUE, out2File=FALSE, 
			grid.related.to.temperaturefiles, valid.years=1952:2009,
			combine.time.series=TRUE, range=10000, 
			alt.range=50, temperature.scale.factor=1, 
			dataPath=getwd(),
			temperature.filenames,
			temperature.matrix, 
			pathForTmpFiles=getwd(),
			pathToSave=getwd(), plant="beech"){
	if (!extraction.done){
		cat("Reading ",plant,"-budburst TSV-file:\n",sep="")
		data.budburst <- data.read.phase(path=dataPath, 
				filename=paste(plant,"_budburst.tsv",sep=""))
		save(data.budburst, 
			file=paste(pathForTmpFiles,"/",plant,"-budburst.RData",sep=""))
		cat("Done!\n")

		cat("Reading ",plant,"-leafcolouring TSV-file:\n",sep="")
		data.leafcolouring <- data.read.phase(path=dataPath, 
				filename=paste(plant,"_leafcolouring.tsv",sep=""))
		save(data.leafcolouring, 
			file=paste(pathForTmpFiles,"/",plant,"-leafcolouring.RData",sep=""))
		cat("Done!\n")

		cat("Extracting essential data:\n")
		data.extracted <- data.extract(data.budburst=data.budburst, 
						data.leafcolouring=data.leafcolouring,
						out2File=out2File, silent=FALSE)
		save(data.extracted, 
			file=paste(pathForTmpFiles,"/dataset-extracted-",plant,".RData",sep=""))
		cat("Done!\n")
	} else {
		load(paste(pathForTmpFiles,"/dataset-extracted-",plant,".RData",sep=""))
	}

	if (combine.time.series){
		data.combined <- data.combine(dataset=data.extracted, 
				range=range, alt.range=alt.range, 
				shuffle=TRUE, tries=100, silent=FALSE, 
				out2File, clusters.tmp.file=paste(
					pathForTmpFiles,"/clusters.RData",sep=""))
		save(data.combined, 
			file=paste(pathForTmpFiles,"/dataset-combined-",plant,".RData",sep=""))
		data.modified <- data.combined
	} else {
		data.modified <- data.extracted
	}

	cat("Adding temperatures:\n")
	dataset <- data.addTemperatures(dataset=data.modified, out2File=out2File,
				grid.related.to.Temperatures=grid.related.to.temperaturefiles,
				temperature.filenames=temperature.filenames,
				temperature.matrix=NULL,  
				temperature.scale.factor=temperature.scale.factor, silent=FALSE)
	
	save(dataset, file=paste(pathToSave,"/",plant,"-dataset",
		ifelse(combine.time.series,"-cts",""),".RData",sep=""))

	cat("Done!\n")
}
