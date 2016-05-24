
#+ brickFunction
spdfToBrick = function(x, 
    template,
    logSumExpected=FALSE,
    pattern = '^expected_[[:digit:]]+$'
) {
  
  
  
  if(class(x)=='SpatialPolygonsDataFrame'){
    x = list(x) 
  }
  if(is.null(names(x)))
    names(x) = as.character(1:length(x))
  
  if(is.numeric(template)) {
    template = squareRaster(
        x[[1]], template[1])
  }
  
  forRaster = NULL
  
  haveRgdal = requireNamespace('rgdal', quietly=TRUE)

  for(Dcensus in rev(names(x))){
    
    if(haveRgdal &
        !identical(
            projection(x[[Dcensus]]),
            projection(template)
            )){
      x[[Dcensus]] = spTransform(
          x[[Dcensus]],
          CRS(projection(template))
      )
    }
    
		Sx = 1:length(x[[Dcensus]])
		
    Sid = rasterize(
        x[[Dcensus]],
        template,
        field=Sx
    )
		
		# number of raster cells assigned to each polygon
		nCellPerPoly = rep(1, length(Sx))
		nCellTable = table(values(Sid))
		nCellPerPoly[as.numeric(names(nCellTable))] = nCellTable
		
		dataHere = as.matrix(x[[Dcensus]]@data[,
        grep(pattern, names(x[[Dcensus]])), drop=FALSE
    ]) / nCellPerPoly

		dataHere[is.na(dataHere)] = 0

		# assign expected counts to raster cells
		forRasterHere = dataHere[values(Sid), ]
		colnames(forRasterHere) = colnames(dataHere)
		forRasterHere[is.na(forRasterHere)] = 0
		
		# polygons not assigned to cells
		notInRaster = which(! Sx %in% values(Sid))
		if(length(notInRaster)){
			polyCentres = SpatialPoints(x[[Dcensus]][notInRaster,])
			polyCell = cellFromXY(template, polyCentres@coords)
			
			dataNotInRaster = aggregate(
					dataHere[notInRaster,], 
					list(cell=polyCell), 
					FUN=sum, na.rm=TRUE)
			forRasterHere[dataNotInRaster$cell, ] = 
					forRasterHere[dataNotInRaster$cell, ] + 
					as.matrix(dataNotInRaster[,-1])
		}
		
    forRaster = cbind(
				forRasterHere,
        forRaster)
  }
  

	
  if(logSumExpected){
    forRaster = apply(forRaster, 1, sum, na.rm=TRUE)
    forRaster[forRaster<=0] = NA	
		# divide by cell size to get intensity
    forRaster = matrix(log(forRaster) - sum(log(res(template))))
		colnames(forRaster) = 'logExpected'
  } else {
		# divide by cell size to get intensity
    forRaster = forRaster / prod(res(template))
	}
  
  if(is.matrix(forRaster)){
    result = brick(template, nl=ncol(forRaster))
    values(result) = forRaster
  } else if(is.vector(forRaster)){
    result = template
    values(result) = forRaster
  } else {
    warning("no data extracted")
    result = template
  }  
  
  result
}

#'

