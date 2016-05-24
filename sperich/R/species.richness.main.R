species.richness.main <-
function(dataset.all.species, dataset.landwater, dataset.height, distances=1:10, weight=0.5, 
		resolution=1, narrow.endemic=FALSE, narrow.endemic.limit=5, 
		upperbound=5, cross.validation=FALSE, fold=5, loocv.limit=10,
		create.image=FALSE, image.title="Interpolated Species Richness", 
		directory=getwd(), filename="species.richness.png",
		evaluation=FALSE, eval.title="Histogramm", adjust=FALSE, 
		clusterlimit=100, predefinedClusterlist=NULL, all.species=-1, 
		export=FALSE, drivername="GTiff", exportname="species.richness.tif", 
		noninterpolatedgrid=NULL, silent=TRUE, do.parallel=FALSE){

	#calculating coordinates and dimension of grid
	if (!silent)
		cat("Preparing grid.. ")
	dimension <- getDimension(dataset.all.species, resolution)
	shift <- getShift(dataset.all.species)
	if (!silent)
		cat("Done!\n")

	#preparing landwatermask
	if (!silent)
		cat("Preparing land-water-mask.. ")
	landwatermask <- createLandwatermask(dataset.landwater, dimension, shift, resolution)
	if (!silent)	
		cat("Done!\n")
	
	if (!silent)
		cat("Preparing height-information.. ")
	height.matrix <- createHeightmask(dataset.height, dimension, shift, resolution)
	if (!silent)	
		cat("Done!\n")

	#add height to landwatermask
	height.matrix[which(height.matrix < 0)] <- 0
	landwatermask[which(landwatermask >= 0)] <- height.matrix[which(landwatermask >= 0)]

	if (narrow.endemic){
		if (!silent)
			cat("Extracting narrow endemic species.. ",sep="")
		# search narrow endemic species out of dataset
		all.species <- getNarrowEndemics(dataset.all.species, all.species, 
							narrow.endemic.limit, dimension, shift, 
							resolution)
		if (!silent)	
			cat("Done!\n")
	}

	if (cross.validation){
		if (!silent)
			cat("Calculating cross-validated species richness.. \n")
		species.richness.weighted <- species.richness.cv(dataset.all.species, landwatermask, fold, 
									loocv.limit, distances, weight, dimension, 
									shift, resolution, upperbound,  all.species, 
									silent, do.parallel)
	} else {
		if (!silent)
			cat("Calculating species richness.. \n")
		species.richness.weighted <- species.richness(dataset.all.species, landwatermask, distances, 
									weight, dimension, shift, resolution,
									upperbound, all.species, silent, do.parallel)
	}

	if (adjust){
		if (!silent)
			cat("Adjusting grid .. ")
		if (is.null(noninterpolatedgrid)){
			noninterpolatedgrid <- createNonInterpolatedGrid(dataset.all.species, 
							dimension, shift, resolution, all.species)
		}
		
		if (is.null(predefinedClusterlist)){
			clusterlist <- searchClusters(species.richness.weighted, dimension, shift, 
							resolution=1, clusterlimit)
		} else {
			clusterlist <- predefinedClusterlist
		}

		species.richness.weighted <- adjustment(species.richness.weighted, noninterpolatedgrid,
								clusterlist)
		if (!silent)
			cat("Done!\n")
	}

	if (create.image){
		if (!silent)
			cat("Creating Map (PNG-File).. ")
		if(createImage(species.richness.weighted, landwatermask, image.title, directory, 
				filename, shift, parts=10, resolution)){
			if (!silent)
				cat("PNG-File successfully created!\n")
		} else {
			cat("PNG-File-Creation failed!\n")
		}
	}

	if (evaluation) {
		if (!silent)
			cat("Creating Evaluation (PNG-File).. ")
		if(evaluate(result.grid.one=species.richness.weighted, result.grid.two=NULL, 
				title.one=eval.title, title.two=NULL, xmax=400, directory=directory, 
				filename=paste(substring(filename, 1, nchar(filename)-4),".histogramm.png",sep=""))){
			if (!silent)
				cat("PNG-File successfully created!\n")
		} else {
			if (!silent)
				cat("PNG-File-Creation failed!\n")
		}
	}
	
	if (export) {
		if (!silent)
			cat("Export .. ")
		exportAsGDAL(species.richness.weighted, shift, resolution, 
				directory=directory, filename=exportname, 
				drivername=drivername)
		if (!silent)
			cat("Done!\n")
	}

	return(species.richness.weighted)
}
