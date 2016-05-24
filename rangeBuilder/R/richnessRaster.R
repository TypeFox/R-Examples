# Function to receive either a list of SpatialPolygons or a rasterStack and generate a richness raster
# add option to allow species that have ranges smaller than grid size to get enlarged so as to count


# ranges is list of SpatialPolygons or SpatialPolygonsDataFrames, or a rasterstack
# resolution: size of cells
# resUnits: if 'degrees', then raster will be unprojected, if 'meters' then raster will be projected to equal area Behrmann projection
# extent: if 'auto', then the maximal extent of the polygons will be used, if input is raster, this will be ignored as the extent of the rasterstack will be used. If not auto, must be a numeric vector of length 4 with minLong, maxLong, minLat, maxLat

# if a rasterstack is provided, extent, resolution and projection is taken from that

richnessRaster <- function(ranges, resolution = 1, resUnits = 'degrees', extent = 'auto', speciesByCell = FALSE, coverCutoff = 0.5, nthreads = 1) {
	
	if (!resUnits %in% c('degrees', 'meters')) {
		stop('resUnits must either be degrees or meters')
	}
	
	if (nthreads > 1) {
		if (!"package:parallel" %in% search()) {
			stop("Please load package 'parallel' for using the multi-thread option\n");
		}
	}

	if (class(ranges) == 'list') {
		if (class(ranges[[1]]) != 'SpatialPolygons' & class(ranges[[1]]) != 'SpatialPolygonsDataFrame') {
			stop('Input must be a list of SpatialPolygons or a RasterStack.')
		}
	}
	
	# if spatialpolygons
	if (class(ranges) == 'list') {
		if (class(ranges[[1]]) == 'RasterLayer') {
			stop('Rasters must be provided as a RasterStack.')
		}
		if (class(ranges[[1]]) == 'SpatialPolygons' | class(ranges[[1]]) == 'SpatialPolygonsDataFrame') {
			
			# test that all have same CRS
			if (length(unique(sapply(ranges, proj4string))) != 1) {
				stop('proj4string of all polygons must match.')
			}

			if ('auto' %in% extent) {
				#get overall extent
				masterExtent <- getExtentOfList(ranges)
			} else if (is.numeric(extent) & length(extent) == 4) {
				masterExtent <- list(minLong = extent[1], maxLong = extent[2], minLat = extent[3], maxLat = extent[4])
			} else {
				stop("extent must be 'auto' or a vector with minLong, maxLong, minLat, maxLat.")
			}
			
			if (resUnits == 'degrees') {
				proj <- '+proj=longlat +datum=WGS84'
			} else {
				#Behrmann equal area projection
				proj <- '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
			}
			
			# if proj does not match CRS of polygons, transform
			if (proj != proj4string(ranges[[1]])) {
				ranges <- lapply(ranges, function(x) sp::spTransform(x, CRS(proj)))
			}
			
			#create template raster
			ras <- raster(xmn = masterExtent$minLong, xmx = masterExtent$maxLong, ymn = masterExtent$minLat, ymx = masterExtent$maxLat, resolution = rep(resolution, 2), crs = proj)
			
			# get the percent cover for each polygon, for each cell
			if (nthreads > 1) {
				cl <- parallel::makePSOCKcluster(nthreads)
				parallel::clusterExport(cl = cl, varlist = c('ranges', 'ras', 'rasterize'), envir = environment())
				polycover <- parallel::parLapply(cl, ranges, function(x) rasterize(x, ras, getCover = TRUE))
				parallel::stopCluster(cl)
			} else {
				polycover <- lapply(ranges, function(x) rasterize(x, ras, getCover = TRUE))
			}
			
			polycover <- lapply(polycover, function(x) x / 100)
			
			#identify cells for each polygon that satisfy cutoff
			cellInd <- lapply(polycover, function(x) which(values(x) > coverCutoff))
			
			if (!speciesByCell) {
				# create richness raster
				values(ras) <- 0
				pb <- txtProgressBar(min = 0, max = length(cellInd), style = 3)
				for (i in 1:length(cellInd)) {
					setTxtProgressBar(pb, i)
					values(ras)[cellInd[[i]]] <- values(ras)[cellInd[[i]]] + 1
				}
				ras[ras == 0] <- NA
				names(ras) <- 'spRichness'
				return(ras)
			} else {
				#create list of cells with species names
				spByCell <- vector('list', length = ncell(ras))
				for (i in 1:length(cellInd)) {
					if (length(cellInd[[i]]) > 0) {
						for (j in 1:length(cellInd[[i]])) {
							spByCell[[cellInd[[i]][j]]] <- c(spByCell[[cellInd[[i]][j]]], names(cellInd)[i])
						}
					}
				}
				ind <- which(sapply(spByCell, is.null) == TRUE)
				spByCell[ind] <- lapply(spByCell[ind], function(x) vector('character', length = 0))
				return(spByCell)
			}
		}
	}
	
	
	# if rasterstack as input
	if (class(ranges) == 'RasterStack') {
		
		#check that all rasters have values
		valCheck <- minValue(ranges)
		badEntries <- which(is.na(valCheck))
		if (length(badEntries) > 0) {
			badEntries <- paste(which(is.na(valCheck)), collapse = ', ')
			stop(paste0('The following rasters have no non-NA cells: ', badEntries, '.'))
		}
		
		# rasterstack calculations only
		# convert all rasters to presence/absence: create matrix of cells (rows) x raster (cols)
		if (nthreads > 1) {
			cl <- parallel::makePSOCKcluster(nthreads)
			parallel::clusterExport(cl = cl, varlist = c('ranges', 'values'), envir = environment())
			vals <- parallel::parLapply(cl, 1:nlayers(ranges), function(x) values(ranges[[x]]))
			vals <- parallel::parLapply(cl, vals, function(x) ifelse(is.na(x), 0, 1))
			vals <- do.call(cbind, vals)
			parallel::stopCluster(cl)
		} else {
			vals <- lapply(1:nlayers(ranges), function(x) values(ranges[[x]]))
			vals <- lapply(vals, function(x) ifelse(is.na(x), 0, 1))
			vals <- do.call(cbind, vals)
		}
		
		if (!speciesByCell) {
			#get sum values
			cellSums <- rowSums(vals)
			
			ras <- ranges[[1]]
			values(ras) <- cellSums
			
			#remove zero cells
			ras[ras == 0] <- NA
			names(ras) <- 'spRichness'
			return(ras)

		} else {
			# create list where each item is a cell that contains a vector of species names
			# in vals, rows are cells, columns are species
			if (nthreads > 1) {
				cl <- parallel::makePSOCKcluster(nthreads)
				parallel::clusterExport(cl = cl, varlist = 'vals', envir = environment())
				spByCell <- parallel::parApply(cl, vals, 1, function(x) which(x == 1))
				spByCell <- parallel::parLapply(cl, spByCell, function(x) names(ranges)[x])
				parallel::stopCluster(cl)
			} else {
				spByCell <- apply(vals, 1, function(x) which(x == 1))
				spByCell <- lapply(spByCell, function(x) names(ranges)[x])
			}
			return(spByCell)
		}
	}	
	
}


