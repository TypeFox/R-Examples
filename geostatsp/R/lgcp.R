lgcp = function(formula=NULL, data,  grid, covariates=NULL, 
		border,
		...) {
	
	if(!missing(border)){
		if(!.compareCRS(data, border))
			border = spTransform(border, CRS(proj4string(data)))
	}
	
	if(is.numeric(grid)) {
		if(!missing(border)){
			cells = squareRaster(border,grid)
		} else {
			cells = squareRaster(data,grid)
		}
	} else {
		cells = squareRaster(grid)
	}
	
# create data
	
	if(!missing(border)) {
		inBorder = over(
				data, 
				as(border, 'SpatialPolygons')
		)
		data = data[!is.na(inBorder),]
	}
	
	counts = rasterize(
			SpatialPoints(data), 
			cells, fun="count")
	names(counts) = "count"
	counts[is.na(counts)] = 0
	
	if(!missing(border)) {
		counts = raster::mask(counts, border)
	}
	
# the formula	
	if(is.null(formula)) {
		formula = as.formula(
				paste(c("count ~ 1", names(covariates)), collapse="+")
		)
	}

	formula	= update.formula(formula,
			.~.+offset(logCellSize) 
	)
	formula = update.formula(formula, count ~ .)

	
	# cell size offset
	logCellSize = cells
	names(logCellSize) = "logCellSize"
	values(logCellSize) =  sum(log(res(cells)) )

  	if(class(covariates)=="RasterLayer") {
		covariates = list( logCellSize, covariates)
		names(covariates) = unlist(lapply(covariates, names))
	} else {
    if(length(covariates)){
		  covariates = c(covariates, logCellSize=logCellSize)
    } else {
      covariates = logCellSize
    }
	}
	
	dots = list(...)
	if(!any(names(list(...))=='family')) {
		dots$family='poisson'
	}
	dots$formula = formula
	dots$data = counts
	dots$grid = cells
	dots$covariates=covariates
	

	result = do.call(glgm, dots)

result

}


