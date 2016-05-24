# Purpose        : NEW plotKML CLASSES;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : Pierre Roudier (pierre.roudier@landcare.nz); Dylan Beaudette (debeaudette@ucdavis.edu); 
# Dev Status     : Pre-Alpha
# Note           : for more info see [http://cran.r-project.org/doc/manuals/R-exts.html];


## Color palette:
setClass("sp.palette", representation(type = 'character', bounds = 'vector', color = 'character', names = 'character', icons = 'character'), validity = function(object) {
   if(!class(object@bounds)=="numeric")
      return('Vector with upper and lower limits required')
   if((length(object@bounds)-1)!=length(object@color)|(length(object@bounds)-1)!=length(object@names))
      return('Size of bounds (-1), colors and element names must be equal')
   if(any(nchar(object@color)<7|nchar(object@color)>9))
      return('Colors in the hex system required') 
})

## A new class for SpatialMetadata:
setClass("SpatialMetadata", representation(xml = "XMLInternalDocument", field.names = "character", palette = "sp.palette", sp = "Spatial"), validity = function(object) {
    if(!(xmlName(xmlRoot(object@xml))=="metadata"|xmlName(xmlRoot(object@xml))=="MD_Metadata"))
      return("XML file tagged 'metadata' or 'MD_Metadata' not found")
    ## check the metadata names:
    ny <- unlist(xmlToList(object@xml))
    met <- data.frame(metadata=names(ny), value=paste(ny))
    mdnames <- read.csv(system.file("mdnames.csv", package="plotKML"))
    field_names <- merge(met, mdnames[,c("metadata","field.names")], by="metadata", all.x=TRUE, all.y=FALSE)[,"field.names"]    
    if(!any(field_names %in% object@field.names))
      return("Field names do not match the column names in the xml slot")
    if(!class(object@field.names)=="character")
      return("Field names as character vector required")      
})

## A new class for SpatialPhoto:
setClass("SpatialPhotoOverlay", representation(filename = "character", pixmap = "pixmapRGB", exif.info = "list", PhotoOverlay = "list", sp = "SpatialPoints"), validity = function(object) {
    if(length(object@filename)==0&is.null(object@pixmap)){
      return("Either 'pixmap' slot or 'filename' need to be specified.")
    }
    if(length(object@sp)>1){
      return("'SpatialPoints' object of length 1 required.")
    }
    
    # minimum info [http://www.sno.phy.queensu.ca/~phil/exiftool/TagNames/EXIF.html]:
    exif.sel <- c("DateTime", "ExposureTime", "FocalLength", "Flash")
    for(i in 1:length(exif.sel)){
    if(length(which(names(object@exif.info) %in% exif.sel[i]))==0)
      return(paste("'", exif.sel[i], "' field required.", sep=""))
    }
    
    # minimum info [http://code.google.com/apis/kml/documentation/kmlreference.html#photooverlay]:
    geo.sel <- c("rotation", "leftFov", "rightFov", "bottomFov", "topFov", "near", "shape", "range", "tilt", "heading")
    for(i in 1:length(geo.sel)){
    if(length(which(names(object@PhotoOverlay) %in% geo.sel[i]))==0)
      return(paste("'", geo.sel[i], "' field required.", sep=""))
    }
    
    if(object@PhotoOverlay$rotation < 0 | object@PhotoOverlay$rotation > 180){
      return("Check KML validity: kml:angle180")
    }
    if(object@PhotoOverlay$leftFov < -180 | object@PhotoOverlay$leftFov > 0){
      return("Check KML validity: kml:angle180")
    }
    if(object@PhotoOverlay$rightFov < 0 | object@PhotoOverlay$rightFov > 180){
      return("Check KML validity: kml:angle180")
    }
    if(object@PhotoOverlay$bottomFov < -90 | object@PhotoOverlay$bottomFov > 0){
      return("Check KML validity: kml:angle90")
    }
    if(object@PhotoOverlay$topFov < 0 | object@PhotoOverlay$topFov > 90){
      return("Check KML validity: kml:angle90")
    }
    if(!(object@PhotoOverlay$shape %in% c("rectangle", "cylinder", "sphere"))){
      return("Shape can be only one of the following: 'rectangle', 'cylinder', 'sphere'.")
    }
    if(object@PhotoOverlay$range < 0){
      return("Check KML validity: positive value required")
    }
    if(object@PhotoOverlay$tilt < 0 | object@PhotoOverlay$tilt > 90){
      return("Check KML validity: kml:angle90")
    }     
    if(object@PhotoOverlay$heading < 0 | object@PhotoOverlay$heading > 360){
      return("Check KML validity: kml:angle360")
    }      
})

## A new class for SpatialPredictions:
setClass("SpatialPredictions", representation(variable = "character", observed = "SpatialPointsDataFrame", regModel.summary = "ANY", vgmModel = "data.frame", predicted = "SpatialPixelsDataFrame", validation = "SpatialPointsDataFrame"), validity = function(object) {
    if(any(!(object@variable %in% names(object@observed@data)))){
      return("Variable name not available in the 'data' slot")
    }
    if(any(!(object@variable %in% names(object@predicted@data)))){
      return("Variable name not available in the 'predicted' slot")
    }
    if(length(object@validation) <50){
      return("Validation data critically small (<50) for reliable validation") 
    }
})

## New classes for SpatialSimulations:
setClass("SpatialVectorsSimulations", representation(realizations = "list", summaries = "SpatialGridDataFrame"), validity = function(object) {
   require(sp)
   object.ov <- over(y=object@summaries, x=as(object@realizations[[1]], "SpatialPoints"))
    if(length(object.ov)==0)
      return("'Realizations' and 'summaries' objects do not overlap spatially")
    if(length(names(object@summaries))<2)
      return("The 'summaries' slot should contain at least two layers (aggregate values and information entropy)")
})

setClass("RasterBrickSimulations", representation(variable = "character", sampled = "SpatialLines", realizations = "RasterBrick"), validity = function(object) {
   require(raster)
   object.ov <- extract(object@realizations, as(object@sampled, "SpatialPoints"))
    if(length(object.ov)==0)
      return("'Realizations' and 'sampled' objects do not overlap spatially")
})


## A new class for SamplingPatterns:
setClass("SpatialSamplingPattern", representation(method = "character", pattern = "SpatialPoints", sp.domain = "SpatialPolygonsDataFrame"), validity = function(object) {
    ov <- over(y=object@sp.domain, x=object@pattern)
    if(length(ov)==0)
      return("'Pattern' and 'sp.domain' do not overlap spatially")
})

## A new class for RasterBrickTimeSeries:
setClass("RasterBrickTimeSeries", representation(variable = "character", sampled = "SpatialPointsDataFrame", rasters = "RasterBrick", TimeSpan.begin = "POSIXct", TimeSpan.end = "POSIXct"), validity = function(object) {
    sel <- !is.na(object@TimeSpan.begin)&!is.na(object@TimeSpan.end)
    if(any(object@TimeSpan.begin[sel] > object@TimeSpan.end[sel]))
      return("'TimeSpan.begin' must indicate time before or equal to 'TimeSpan.end'")
    if(!(length(object@TimeSpan.begin)==length(object@TimeSpan.end)&length(object@TimeSpan.begin)==nlayers(object@rasters)))
      return("Length of the 'TimeSpan.begin' and 'TimeSpan.end' slots and the total number of rasters do not match")
    ov <- extract(object@rasters, object@sampled)
    if(!nrow(ov)==length(object@sampled)){
      return("Not all points can be overlaid using the data in the @rasters slot")
    }
})

### A new class for SpeciesDistributionMap:
setClass("SpatialMaxEntOutput", representation(sciname = "character", occurrences = "SpatialPoints", TimeSpan.begin = "POSIXct", TimeSpan.end = "POSIXct", maxent = "MaxEnt", sp.domain = "Spatial", predicted = "RasterLayer"), validity = function(object) {
    if(object@TimeSpan.begin > object@TimeSpan.end)
      return("'TimeSpan.begin' must indicate time before or equal to 'TimeSpan.end'")    
    if(length(object@occurrences) <5)
      return("Occurences critically small (<5) for reliable validation")  
    object.ov <- extract(x=object@predicted, y=object@occurrences)
    if(length(object.ov)==0)
      return("'Occurences' and 'rasters' do not overlap spatially")
    if(is.na(object@predicted@crs@projargs))
      return("Proj4 string must be specified")  
})



################## generic functions ##############


if (!isGeneric("spMetadata")){
  setGeneric("spMetadata", function(obj, ...){standardGeneric("spMetadata")})
}

if (!isGeneric("GetNames")){
  setGeneric("GetNames", function(obj){standardGeneric("GetNames")})
}

if (!isGeneric("GetPalette")){
  setGeneric("GetPalette", function(obj){standardGeneric("GetPalette")})
}

if (!isGeneric("metadata2SLD")){
  setGeneric("metadata2SLD", function(obj, ...){standardGeneric("metadata2SLD")})
}

if (!isGeneric("kml_layer")){
  setGeneric("kml_layer", function(obj, ...){standardGeneric("kml_layer")})
}

if (!isGeneric("kml_metadata")){
  setGeneric("kml_metadata", function(obj, ...){standardGeneric("kml_metadata")})
}

if (!isGeneric("kml")){
  setGeneric("kml", function(obj, ...){standardGeneric("kml")})
}

if (!isGeneric("getCRS")){
  setGeneric("getCRS", function(obj, ...){standardGeneric("getCRS")})
}

if (!isGeneric("reproject")){
  setGeneric("reproject", function(obj, ...){standardGeneric("reproject")})
}

if (!isGeneric("vect2rast")){
  setGeneric("vect2rast", function(obj, ...){standardGeneric("vect2rast")})
}

if (!isGeneric("plotKML")){
  setGeneric("plotKML", function(obj, ...){standardGeneric("plotKML")})
 }
 
# end of script;