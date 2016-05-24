# Purpose        : Generic methods to plot large vectors;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : ;
# Dev Status     : Alpha
# Note           : Implemented for parallel processing;

kml.tiles <- function(obj, 
  folder.name,
  file.name, 
  block.x,
  kml.logo, 
  cpus, 
  home.url=".", 
  desc=NULL, 
  open.kml=TRUE, 
  return.list=FALSE, 
  ...){
  
  if(missing(folder.name)){ folder.name <- normalizeFilename(deparse(substitute(obj))) }
  if(missing(file.name)){ file.name <- paste(normalizeFilename(deparse(substitute(obj))), ".kml", sep="") }
  
  ## check class of object:
  if(any(!(class(obj)=="SpatialPointsDataFrame"|class(obj)=="SpatialLinesDataFrame"|class(obj)=="SpatialPolygonsDataFrame"))){
    stop("Object of class SpatialPoints*, SpatialLines*, SpatialPolygons* expected")
  }
  ## reproject if necessary:
  prj.check <- check_projection(obj, control = TRUE)
  if(!prj.check) { suppressMessages( obj <- reproject(obj) ) }
  ## tile object:
  if(requireNamespace("GSIF", quietly = TRUE)){
    if(any(class(obj)=="SpatialPointsDataFrame")){
      obj.lst <- GSIF::tile(obj, block.x=block.x)
    } else {
      obj.lst <- GSIF::tile(obj, block.x=block.x, tmp.file = TRUE)
    }
    ## some tiles might be empty and need to be removed...
    obj.lst <- obj.lst[sapply(obj.lst, length)>0]
  } else {
    stop("Install and load package 'GSIF'")
  }
  ## list of bounding boxes:
  file.lst <- sapply(1:length(obj.lst), function(j){paste0(folder.name, "_T", j, ".kml")})
  folder.lst <- sapply(1:length(obj.lst), function(j){paste0(folder.name, "_T", j)})
  if(any(class(obj)=="SpatialPointsDataFrame")){
    north <- sapply(obj.lst, function(j){max(slot(j, "coords")[,2])})
    south <- sapply(obj.lst, function(j){min(slot(j, "coords")[,2])})
    east <- sapply(obj.lst, function(j){max(slot(j, "coords")[,1])})
    west <- sapply(obj.lst, function(j){min(slot(j, "coords")[,1])})
  } else {
    north <- sapply(obj.lst, function(j){j@bbox[2,2]})
    south <- sapply(obj.lst, function(j){j@bbox[2,1]})
    east <- sapply(obj.lst, function(j){j@bbox[1,2]})
    west <- sapply(obj.lst, function(j){j@bbox[1,1]})  
  }
  ## write all tiles to KML:
  if(requireNamespace("snowfall", quietly = TRUE)&requireNamespace("parallel", quietly = TRUE)){  
    if(missing(cpus)){ cpus <- parallel::detectCores(all.tests = FALSE, logical = FALSE) }
    snowfall::sfInit(parallel=TRUE, cpus=cpus)
    ## this might take a lot of RAM...
    snowfall::sfExportAll()
    snowfall::sfLibrary(package="rgdal", character.only=TRUE)
    snowfall::sfLibrary(package="sp", character.only=TRUE)
    snowfall::sfLibrary(package="plotKML", character.only=TRUE)
    snowfall::sfLibrary(package="XML", character.only=TRUE)
    x <- snowfall::sfLapply(1:length(obj.lst), function(j){kml(obj.lst[[j]], file.name=file.lst[j], ...)})
    snowfall::sfStop()
  } else {
    x <- lapply(1:length(obj.lst), function(j){kml(obj.lst[j], folder.name=folder.lst[j], file.name=file.lst[j], ...)})
  }
  lst <- data.frame(kml.tile=file.lst, north=north, south=south, east=east, west=west)
  
  kml_open(file.name)
  kml.out <- get("kml.out", envir=plotKML.fileIO)
  
  ## add description:
  if(!is.null(desc)){ 
    description_txt <- sprintf('<description><![CDATA[%s]]></description>', desc)
    parseXMLAndAdd(description_txt, parent=kml.out[["Document"]]) 
  }
  ## Region and network link section:
  network_txt <- sprintf('
    <NetworkLink>
        <name>%s</name>
        <Region>
          <Lod>
            <minLodPixels>128</minLodPixels><maxLodPixels>-1</maxLodPixels>
          </Lod>
          <LatLonAltBox>
            <north>%.5f</north><south>%.5f</south>
            <east>%.5f</east><west>%.5f</west>
          </LatLonAltBox>
        </Region>
        <Link>
          <href>%s</href>
          <viewRefreshMode>onRegion</viewRefreshMode>
        </Link>
      </NetworkLink>', unlist(lst[["kml.tile"]]), unlist(lst[["north"]]), unlist(lst[["south"]]), unlist(lst[["east"]]), unlist(lst[["west"]]), paste(home.url, unlist(lst[["kml.tile"]]), sep="/"))   
  parseXMLAndAdd(network_txt, parent=kml.out[["Document"]])
  assign('kml.out', kml.out, envir=plotKML.fileIO)
  if(!missing(kml.logo)){ kml_screen(image.file = kml.logo, position = "UR", sname = "logo") }
  kml_close(file.name)
  if(open.kml==TRUE){
    kml_View(file.name)
  } else {
    message(paste("Object written to:", file.name))
  }
  if(return.list==TRUE){
    return(obj.lst)
  }
}