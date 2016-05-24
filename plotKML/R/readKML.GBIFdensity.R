# Purpose        : Reads GBIF celldensity records from KML to RasterBrick;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : Dylan Beaudette (debeaudette@ucdavis.edu); 
# Status         : pre-alpha
# Note           : Description of the GBIF celldensity records is available at: [http://data.gbif.org/ws/rest/density];


readKML.GBIFdensity <- function(kml.file, gbif.url = FALSE, silent = FALSE){

  if(is.na(file.info(kml.file)$size)){
    stop(paste(kml.file, "does not exist in the working directory"))
  }
  else{
  if(file.info(kml.file)$size < 1024){
    warning("File size < 1KB... possibly an empty object")
  }
  }
  
  if(silent==FALSE){
  warning("All usage of these data must be in accordance with the GBIF Data Use Agreement - see http://www.gbif.org/DataProviders/Agreements/DUA", immediate. = TRUE)
  }
  
  ret <- xmlTreeParse(kml.file, useInternalNodes = TRUE)
  # top structure: 
  top <- xmlRoot(ret)
  doc.name <- xmlValue(top[["Document"]][["name"]])
  # strip the latin name:
  lname <- strsplit(doc.name, "GBIF Data Portal Occurrence Density Layer for species: ")[[1]][2]
  # get coordinates and counts:
  count.lst <- NULL
  crds.lst <- NULL
  desc.lst <- NULL
  # TH: xpath commands e.g. xpathApply(top, "//Placemark/name", xmlValue) would not produce anything;
  nf <- which(names(top[["Document"]]) %in% "Folder")
  if(length(nf)>0){
  if(silent == FALSE){
  pb <- txtProgressBar(min=0, max=length(nf), style=3)
  }
  for(i in seq_along(nf)){
  count.lst[[i]] <- xmlValue(top[["Document"]][[nf[i]]][["Placemark"]][["name"]])
  crds.lst[[i]] <- xmlValue(top[["Document"]][[nf[i]]][["Placemark"]][["Point"]][["coordinates"]])  
  if(gbif.url == TRUE){
  desc.lst[[i]] <- xmlValue(top[["Document"]][[nf[i]]][["Placemark"]][["description"]])
  }
  if(silent == FALSE){ setTxtProgressBar(pb, i) }
  }
  if(silent == FALSE){ close(pb) }
  # number of records:
  nrecs <- as.numeric(sapply(count.lst, function(x){strsplit(strsplit(x, " - ")[[1]][2], "record")[[1]][1]}))
  # coordinates of the cell:
  lon <- as.numeric(sapply(crds.lst, function(x){strsplit(x, ",")[[1]][1]}))
  lat <- as.numeric(sapply(crds.lst, function(x){strsplit(x, ",")[[1]][2]})) 
  # cell id and taxon concept keys:
  if(gbif.url == TRUE){
  cellid <- paste(sapply(sapply(desc.lst, function(x){strsplit(x, "cellid=")[[1]][2]}), function(x){strsplit(x, "&")[[1]][1]}))
  taxonconceptkey <- paste(sapply(sapply(desc.lst, function(x){strsplit(x, "taxonconceptkey=")[[1]][2]}), function(x){strsplit(x, '">')[[1]][1]}))
  out <- data.frame(lon, lat, nrecs, cellid, taxonconceptkey)
  }
  else{ 
  out <- data.frame(lon, lat, nrecs)
  }
  ## Rasterize:
  coordinates(out) <- ~lon+lat
  # create an empty grid:
  crs = CRS("+proj=longlat +datum=WGS84")
  r.sp <- SpatialGrid(GridTopology(cellcentre.offset=c(-179.5,-89.5), cellsize=c(1,1), cells.dim=c(360,180)), proj4string = crs)
  r <- raster(r.sp)
  # rasterize - convert points to a raster map:
  out.r <- rasterize(out, r, "nrecs", fun=mean)
  names(out.r) <- lname

  if(gbif.url == TRUE){
  out2.r <- rasterize(out, r, "cellid")
  out2.r@data@isfactor <- TRUE
  out2.r@legend@values <- levels(as.factor(cellid))
  names(out.r) <- "cellid"
  out3.r <- rasterize(out, r, "taxonconceptkey")
  out3.r@data@isfactor <- TRUE
  out3.r@legend@values <- levels(as.factor(taxonconceptkey))
  names(out3.r) <- "taxonconceptkey"  
  out.r <- brick(out.r, out2.r, out3.r)
  }

  out.r@title <- paste("GBIF cell density count for:", lname)
  out.r@history <- paste(kml.file)
  out.r@zname <- "Last update"
  out.r <- setZ(out.r, as.character(file.info(kml.file)$mtime))
  
  return(out.r)
  
  }
  else 
  { 
  warning("No nodes of type 'Folder' were found. Not a valid GBIF celldensity KML.") 
  return(NULL)
  }
}
  
# end of script;