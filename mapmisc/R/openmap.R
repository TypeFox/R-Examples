osmTiles = function(name) {
	result = c(
			osm = "http://tile.openstreetmap.org",
			"osm-no-labels"="http://a.tiles.wmflabs.org/osm-no-labels/",
      "osm-de"="http://c.tile.openstreetmap.de/tiles/osmde",
			"osm-transport"="http://tile2.opencyclemap.org/transport/",
			"bw-mapnik"="http://b.tiles.wmflabs.org/bw-mapnik",
			mapquest="http://otile1.mqcdn.com/tiles/1.0.0/osm/",
			"mapquest-sat"="http://otile1.mqcdn.com/tiles/1.0.0/sat",
      "mapquest-labels"='http://otile3.mqcdn.com/tiles/1.0.0/hyb/',
      'osm-cyclemap' = 'http://a.tile.opencyclemap.org/cycle/',
      'osm-seamap' = 'http://tiles.openseamap.org/seamark/',
      'osm-fr' = 'http://a.tile.openstreetmap.fr/osmfr/',
#      'osm-rail' = 'http://a.tiles.openrailwaymap.org/standard/',
# rail is 512 insstead of 256 tiles
      #			hill="http://www.toolserver.org/~cmarqu/hill/",
			'landscape'="http://tile.opencyclemap.org/landscape/",
#		"osm-retina"="http://tile.geofabrik.de/osm_retina/",
		"opentopomap" = "http://opentopomap.org/tiles/",
#		"osm2world"="http://tiles.osm2world.org/osm/pngtiles/n/",
#		bvg="http://mobil.bvg.de/tiles/",
#	landshaded="http://tiles.openpistemap.org/landshaded/",
	"maptoolkit"="http://tile2.maptoolkit.net/terrain/",
#	skobbler="http://tiles3.skobbler.net/osm_tiles2/",	
#  'sputnik' = 'http://a.tiles.maps.sputnik.ru/tiles/kmt2/',
	waze="http://tilesworld.waze.com/tiles/",
  'waze-us'='http://livemap-tiles1.waze.com/tiles/',
#	eu="http://alpha.map1.eu/tiles/"#,
	humanitarian="http://a.tile.openstreetmap.fr/hot/",
 cartodb='http://c.basemaps.cartocdn.com/light_all/',
'cartodb-dark'='http://c.basemaps.cartocdn.com/dark_all/',
historical='http://www.openhistoricalmap.org/ohm_tiles/'#,
#'esri' = 'http://services.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/',
#' 'esri-grey' = 'http://services.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/',
#' 'esri-transport'='http://server.arcgisonline.com/ArcGIS/rest/services/Reference/World_Transportation/MapServer/tile/',
#' 'esri-topo' = 'http://services.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/'
	)
	
	
  

	# language labels don't appear to be working
	languages = c("en","fr","de", "it","es","ru")
	toadd =	paste("http://a.www.toolserver.org/tiles/osm-labels-", languages,"/", sep="")
	names(toadd) = paste("osm-labels-", languages, sep="")
#	result = c(result, toadd)
	
	stamen = c("toner","watercolor")#,"terrain","terrain-background")
	toadd = paste("http://tile.stamen.com/", stamen, "/",sep="")
	names(toadd) = paste("stamen-", stamen, sep="")
	
	result = c(result, toadd)
	
	
	if(!missing(name)) {
		if(name %in% names(result)) {
			result = result[name]
		} else {
			warning("name ", name, " is not a tile path, returning vector of all tile paths")
		}
	}
	
	result
	
}

openmap = function(x, zoom, 
	path="http://tile.openstreetmap.org/",
	maxTiles = 9,
	crs=projection(x),
  buffer=0, fact=1,
	verbose=FALSE) {


	alltiles = osmTiles()
	pathOrig = path
	pathisname = path %in% names(alltiles)
	path[pathisname] = alltiles[path[pathisname]]
	
	if(length(grep("/$", path, invert=TRUE)))
		path[ grep("/$", path, invert=TRUE)] =
				paste(path[ grep("/$", path, invert=TRUE)], "/", sep="")

	if(length(grep("^http[s]*://", path, invert=TRUE)))
		path[ grep("^http[s]*://", path, invert=TRUE)] = 
				paste("http://", 
						path[ grep("^http[s]*://", path, invert=TRUE)], sep="")
	names(pathOrig) = path

	crsOut=crs
	crsIn = crs(x)
	if(all(is.na(crsIn))) {
		if(is.vector(x)){
			crsIn=crsLL
		} else{
			crsIn = crs	
		}
	}
	
	
	extMerc = .getExtent(x,crsIn, buffer, crsMercSphere)
  extMerc = cropExtent(extMerc, openmapExtentMercSphere)
  
	if(missing(zoom)) {
	zoom = 1
	while(nTilesMerc(extMerc, zoom) <= maxTiles & zoom <= 18) {
		zoom = zoom + 1
	}
	zoom = min(c(18,max(c(1, zoom-1))))
	}
	if(verbose) cat("zoom is ", zoom, ", ", nTilesMerc(extMerc, zoom), "tiles\n")

	result = NULL
  
	for(Dpath in rev(path)) {
		thistile = try(
				getTilesMerc(extMerc, zoom=zoom,
				path=Dpath,
				verbose=verbose),
		silent=TRUE	)

		if(class(thistile)=="try-error"){
			message(paste(Dpath, "not accessible"))
      thistile=NULL
		}	else {
			if(length(names(thistile))) {
				theprefix=strsplit(names(thistile), "([rR]ed|[gG]reen|[bB]lue)$",fixed=FALSE)[[1]]
				names(thistile) = gsub(theprefix, paste(pathOrig[Dpath], "",sep=""), 
					names(thistile),fixed=TRUE)		
			}

		ctable = NULL
		if(!is.null(thistile)) {
			if(nlayers(thistile)==1)
				ctable = thistile@legend@colortable
      result =  stack(thistile, result)	
		}
		
		if(length(ctable))
				result[[1]]@legend@colortable = ctable
		
		} # end not try-error
	} # end loop through path	

	
	if(is.null(result)) {
		result = raster(openmapExtentMercSphere,1,1,crs=crsMercSphere)
		values(result) = NA
    attributes(result)$openmap = list(
        tiles=NA,
        message=thistile,
        path=path,
        zoom=zoom
        )
	} 

	
	if(!is.na(crsOut)  ){
		oldColorTable = list()
		for(D in names(result))
			oldColorTable[[D]] = result[[D]]@legend@colortable
		
		if(verbose) cat("reprojecting ", ncell(result), " cells...")

    toRaster = projectExtent(result, crsOut)
		
		
		if(any(fact > 1)){
			res(toRaster) = res(toRaster) / rep_len(fact,2)
		}
		
		resultProj = stack(projectRaster(result, toRaster, method="ngb"))

		for(D in names(resultProj))
			resultProj[[D]]@legend@colortable = oldColorTable[[D]]

		
    if(verbose) cat("done\n")
    
	} else {
		resultProj = stack(result)
	}


		#	resultProj@legend@colortable = result@legend@colortable


	for(D in names(resultProj)) {
      if(length(result[[D]]@legend@colortable)) {
			resultProj[[D]]@legend@colortable =
					result[[D]]@legend@colortable
      if(any(values(resultProj[[D]])==0,na.rm=TRUE)) {
      # set NA's to transparent
      resultProj[[D]]@legend@colortable =
          c('#FFFFFF00',
              resultProj[[D]]@legend@colortable 
          )
      values(resultProj[[D]]) = 1+values(resultProj[[D]])
	}
}
} # end D in names(resultProj)

	if(nlayers(resultProj)==1) 
		resultProj = resultProj[[1]]
		
  attributes(resultProj)$tiles = attributes(thistile)$tiles
  attributes(resultProj)$tiles$path = path
	
	resultProj
}

