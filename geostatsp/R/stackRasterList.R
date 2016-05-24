stackRasterList = function(x, template=x[[1]],method='ngb',mc.cores=NULL) {

	if(class(x)=="SpatialPolygonsDataFrame")
		x = list(x)
	
	if(class(x)=="RasterLayer") {
		x = list(x)
		names(x) =names(x[[1]])
	}
	
	
	if(is.list(x)) {
		if(is.null(names(x)))
			names(x) = paste("c", seq(1, length(x)),sep="")	
	}
	
	Nlayers = length(names(x))
	
	if(length(method)==Nlayers) {
		if(length(names(x)) & all( names(x)%in%names(method)))
			method = method[names(x)]
	} else {
		method = rep(method, Nlayers)
	}
 	
	modefun = function(qq, na.rm=NULL) c(as.numeric(names(which.max(table(qq)))), NA)[1]
	funList = list(ngb=modefun,
			bilinear=mean)
	
	
	template = as(template, "BasicRaster")
	template2 = raster(template)
	
	# function to reproject rasters
	projfun = function(D) {
		if(class(x[[D]])=="SpatialPolygonsDataFrame"){
			if(length(names(x[[D]]))!=1)
				warning("polygon ", D, "has more than one data column, using the first" )
			
			requireNamespace('rgdal', quietly=TRUE ) 
			
			toAdd =  
					rasterize(
							spTransform(x[[D]], CRS(projection(template))), 
							raster(template), field=names(x[[D]][1])
					)
 
		} else { # not a spdf
			if(as(x[[D]], 'BasicRaster')==template) {
				# same projection, same resolution
				toAdd =  x[[D]]			
			}	 else { # different projection or resolution
				# check to see if it's a categorical variable
				if(x[[D]]@data@isfactor) {
					method[D] = "ngb"
				} 
				
				thelevels = levels(x[[D]])
				
				# same projection, different resolution
				testcrs =compareRaster(template, x[[D]],
					extent=FALSE,rowcol=FALSE,crs=TRUE,stopiffalse=FALSE)				
				if(is.na(testcrs)) testcrs = TRUE
				if(testcrs) { # same resolution
					# should we aggregate?
					toAgg = floor(min((dim(x[[D]])/dim(template2))[1:2]))
					if(toAgg > 1) {
						aggFun = funList[[method[D]]]
						xagg = raster::aggregate(x[[D]], fact=toAgg,
								fun=aggFun)
					} else {
						xagg = x[[D]]
					}
					
					toAdd = raster::resample(xagg, template2, method=method[D])
				} else { # differenet resolution
					# different resolution
					toAdd = projectRaster(x[[D]], template,method=method[D])
				} # end different resolution
				
				if(!is.null(thelevels))
					levels(toAdd) = thelevels
			} # end different projection or resolution
		} # end not SPDF
		toAdd
	} # end projfun
	
	# reproject all the rasters

	
	if(!is.null(mc.cores)) {
		resultList = parallel::mcmapply(
				projfun, D=1:Nlayers,
				mc.cores=mc.cores)
	} else {
		resultList = mapply(projfun, D=1:Nlayers)
	}
	
	

	result = resultList[[1]]
	if(Nlayers >1) {
		for(D in 2:Nlayers )
			result = stack(result, resultList[[D]])
	}
	
	
	if(Nlayers == (dim(result)[3]-1) )
		result = result[[-1]]
	names(result) = names(x)
	result
}