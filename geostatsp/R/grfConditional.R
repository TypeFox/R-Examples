grfConditional = function(data, y=1, 
			param, locations, Nsim,
		 	fun=NULL, nuggetInPrediction=TRUE,
			mc.cores=2){
		
		
	if(is.numeric(locations)){
		# locations is number of cells in the x direction
		Nx = locations[1]
		myExtent = 	extent(data)
		Ny = round(locations*diff(data@bbox[2,])/diff(data@bbox[1,]))
		myExtent@ymax = myExtent@ymin + Ny * diff(data@bbox[1,])/Nx
		locations = raster(myExtent, Ny, Nx,
				crs=projection(data))	
	}
	if(nrow(locations) * ncol(locations) > 10^7) warning("there are lots of cells in the prediction raster,\n this might take a very long time")
	
	if(length(y)==1 & length(names(data)) > 1) {
		y = data.frame(data)[,y]
	}

	
	if(length(dim(y))==3) {
		y = matrix(y, nrow=prod(dim(y)[1:2]), ncol=dim(y)[3])
	}

	# convert param to a matrix if it's a list
	# fill in missing parameters with defaults
	param = fillParam(param)
	if(class(param)%in%c("integer","numeric"))
		param = t(as.matrix(param))
	if(class(y)%in%c("integer","numeric")) {
		y = t(as.matrix(y))
	
	Nsamples = unique(c(dim(param)[1], dim(y)[1]))
	Nsamples = Nsamples[Nsamples!=1]
	if(length(Nsamples)>1.5)
		warning("number of samples in y and param is different")
	if(!length(Nsamples))
		Nsamples = 1
	
	param = param[round(seq(1,dim(param)[1], len=Nsim)),,drop=FALSE]
	y = y[round(seq(1,dim(y)[1], len=Nsim)),,drop=FALSE]
}	



rasterCopy = raster(locations)

simFun = function(D) {

	if(param[D,"nugget"] > 0) {
		err.model = "nugget"
		err.param=c(1,param[D,"nugget"],0)
	} else {
		err.model = err.param = NULL
	}
	
	modelv = modelRandomFields(param[D,])
	
	given = SpatialPointsDataFrame(
			data=data.frame(y=y[D,]),
			coords=coordinates(data)
	)
	given = as(given, "RFspatialPointsDataFrame")
	given@.RFparams = list(n=1, vdim=1)
	
	
	res = raster(RandomFields::RFsimulate(
			model=modelv,
			x=as(locations, "GridTopology"),
			data=given			
			))
	
#	res = RandomFields::CondSimu(krige.method="O",
#			x=xseq, y = yseq, grid=TRUE, gridtriple=TRUE,
#			param=NULL, model=modelv,
#			given=data@coords,
#			data=y[D,],
#			err.model=err.model,
#			err.param=err.param, method="direct decomp."
#	)		
	
	
	if(nuggetInPrediction){
		values(rasterCopy) =  
				rnorm(ncell(res), sd=sqrt(param[,"nugget"]))
		res = res + rasterCopy
	}		
	if(!is.null(fun)) {
		res = fun(res)
	}

	res
	}		
	if(requireNamespace("parallel", quietly=TRUE)){
	result = parallel::mcmapply(simFun, 1:Nsim, SIMPLIFY=TRUE, 
			mc.cores=mc.cores)
} else {
	result = mapply(simFun, 1:Nsim, SIMPLIFY=TRUE)
	
}
	if(all(sapply(result, class)=="RasterLayer"))
		result = do.call(brick, result)

result	
}


