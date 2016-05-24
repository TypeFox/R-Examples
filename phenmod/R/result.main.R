result.main <- function(mask.grid, result.grid, plant="beech", 
				model="pim", year=1954, picPath=getwd(), 
				picName="beech-budburst", createFiles=TRUE,
				rsquarePath=getwd(), rsquareFile="rsquare.RData", 
				rsquare.type="cod", silent=FALSE, 
				withOutliers=FALSE){

	result.values <- result.extract.main(mask.grid, 
				result.grid, model=model, 
				interpolate=FALSE, silent=silent, 
				withOutliers=withOutliers)

	#add model name to picName
	picName <- paste(picName, "-", model, sep="")

	#R-Square
	rsquare <- result.rsquare(result.values, type=rsquare.type)
	if (createFiles){
		if (length(list.files(rsquarePath, pattern=rsquareFile))==0){
			rsquare.data <- data.frame(plant=plant, 
						year=year, rsquare=rsquare, 
						stringsAsFactors=FALSE)
		} else {
			load(paste(rsquarePath,"/", rsquareFile, sep=""))
			rsquare.data[dim(rsquare.data)[1]+1,] <- c(plant, year, rsquare)
		}
		save(rsquare.data, file=paste(rsquarePath,"/", rsquareFile,sep=""))
	}
	if (!silent){ cat("R^2=",rsquare,"\n",sep="") }

	# histogramm only for non-interpolated values!
	result.pic.histogramm(values=result.values, 
		picPath=picPath, picName=picName, 
		silent=silent, createFile=createFiles)	
	
	# scatterplots only for non-interpolated values!
	result.pic.scatterplot(values=result.values, 
			picPath=picPath, picName=picName, 
			createFile=createFiles)

	# interpolate
	if (!silent) { cat("Interpolate values:\n",sep="") }
	result.values$doy.model <- result.extract.interpolate(mask.grid=mask.grid, 
					values=result.values$doy.model, 
					alt=mask.grid$alt, x=mask.grid$x, y=mask.grid$y)
	
	result.values$doy.observed <- result.extract.interpolate(mask.grid=mask.grid, 
					values=result.values$doy.observed, 
					alt=mask.grid$alt, x=mask.grid$x, y=mask.grid$y)
		
	result.values$doy.dif <- result.extract.interpolate(mask.grid=mask.grid, 
					values=result.values$doy.dif, 
					alt=mask.grid$alt, x=mask.grid$x, y=mask.grid$y)
	if (!silent){ cat("Interpolation done!\n") }

	result.pic.maps(values=result.values, picPath=picPath, picName=picName, 
			silent=silent, createFile=createFiles)
}
