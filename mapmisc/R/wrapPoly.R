wrapPoly = function(x, crs){
	
	if(is.null(attributes(crs)$crop)) {
		attributes(crs)$crop = llCropBox(crs)$crop
	}
	
	if(requireNamespace('rgeos', quietly=TRUE) & 
			requireNamespace('rgdal', quietly=TRUE)) {	
		toCropX = spTransform(attributes(crs)$crop, crs(x))
		xCrop = rgeos::gDifference(x, toCropX, byid=TRUE)

		row.names(xCrop) = gsub(" (buffer|[[:digit:]]+)$","", row.names(xCrop))
		
		
		# remove short line segments
		if(length(grep("^SpatialLines", class(xCrop)))){
		
		NsegPerLine = unlist(lapply(xCrop@lines, 
				function(qq) length(qq@Lines)))

		for(D in which(NsegPerLine > 3)) {
			# retain only three biggest segments
			lineD = xCrop@lines[[D]]@Lines
			NpointsPerSeg = unlist(lapply(lineD, 
					function(qq) nrow(qq@coords)))
			lineD = lineD[which(NpointsPerSeg > 2)]
			xCrop@lines[[D]]@Lines = lineD
		}
	}
		
		if(any(slotNames(x)=='data')) {
		
			xCropData = x@data[match(
						row.names(xCrop),
						rownames(x@data)
				),]

		rownames(xCropData) = names(xCrop)
		
		xCrop = SpatialPolygonsDataFrame(
				xCrop,
				data=xCropData
		)
		
	}
	
		xTcrop = spTransform(xCrop, crs)
		
	} else {
		xTcrop = NULL
	}
	
	xTcrop
}

llCropBox = function(crs, 
		res=0.5, keepInner=FALSE) {
	
	polarLimit = 90-res
	extentLL = extent(-180,180,-polarLimit,polarLimit)
	
	if(length(grep("proj=moll", as.character(crs)))){
		
		projMoll = CRS("+proj=moll +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0,0,0,0,0 ")
		
		N = 2*round(180/res)
		buffer = res
		
		lonSeq = exp(seq(0, log(90), len=N))
		lonSeq = sort(unique(c(lonSeq, -lonSeq)))


		lonMat = cbind(
				rep(180, N),
				90-c(0, exp(seq(0,log(90),len=N-1)))
		)
		edgePoints = SpatialPoints(lonMat, proj4string=crsLL)
		edgePointsT = spTransform(edgePoints, projMoll)
		
		edgeCoords = abs(edgePointsT@coords)
		edgeCoords = approx(edgeCoords[,1], edgeCoords[,2], 
				xout = max(edgeCoords[,1]) * sin(seq(0, pi/2, len=N)[-N]),
				rule=2)
		
		edgeCoords = cbind(
				c(edgeCoords$x[-1], rev(edgeCoords$x)[-1], -edgeCoords$x[-1], -rev(edgeCoords$x)[-1]),
				c(edgeCoords$y[-1], -rev(edgeCoords$y)[-1], -edgeCoords$y[-1], rev(edgeCoords$y)[-1])
				)
				
				
		toCropPoly = SpatialPolygons(list(
				Polygons(list(
								Polygon(edgeCoords, hole=FALSE)
						), 1)
		), proj4string = crs)

	edgePointsLL = spTransform(
			SpatialPoints(edgeCoords, proj4string=crs),
			crsLL)
	edgePointsLL = raster::crop(edgePointsLL, 
			extent(-181, 181, -polarLimit,polarLimit))
	crs(edgePointsLL) = NA
	toCropLL = rgeos::gBuffer(edgePointsLL, width=res)
	crs(toCropLL) = crsLL

} else {

		toCropPoly = NULL
		
		rasterLLorig = raster(
  			extentLL,
				res=res, crs=mapmisc::crsLL
		)
		
		rasterTorig = projectExtent(rasterLLorig, crs)
		rasterTorig = disaggregate(rasterTorig, 2)
		
		# put 0's around the border
		rasterTsmall = crop(rasterTorig, extend(extent(rasterTorig), -6*res(rasterTorig)))
		values(rasterTsmall) = 1
		rasterT = extend(rasterTsmall, extend(extent(rasterTsmall), 5*res(rasterTsmall)), value=0)

		
		rasterLL = projectRaster(
				from=rasterT, 
				crs=mapmisc::crsLL, 
				res = res, method='ngb')
		rasterLL = crop(rasterLL, extentLL)
		if(keepInner){
			values(rasterLL)[is.na(values(rasterLL))] = 1
		} else {
			values(rasterLL)[is.na(values(rasterLL))] = 0
		}
  	
		borderLL = rasterToPoints(rasterLL, fun=function(x) {x<1}, spatial=TRUE)
		if(requireNamespace('rgeos', quietly=TRUE)) {
			crs(borderLL) = NA
			toCropLL = rgeos::gBuffer(borderLL, width=mean(res*1.5))
			crs(toCropLL) = mapmisc::crsLL
		} else {
			toCropLL = NULL
		}
		
		
	}
	
	if(!is.null(toCropLL)){
		toCropLL = raster::crop(toCropLL, extentLL)
	}
	
	list(crop=toCropLL, poly=toCropPoly)
}
