# Purpose        : Generic methods to plot large rasters;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : ;
# Dev Status     : Alpha
# Note           : Implemented for parallel processing. SuperOverlay file (https://developers.google.com/kml/documentation/kml_21tutorial?csw=1#superoverlays);


plotKML.GDALobj <- function(obj, file.name, block.x, tiles=NULL, tiles.sel=NULL, altitude=0, altitudeMode="relativeToGround", colour_scale, z.lim=NULL, breaks.lst=NULL, kml.logo, overwrite=TRUE, cpus, home.url=".", desc=NULL, open.kml=TRUE, CRS=attr(obj, "projection"), plot.legend=TRUE){

  if(!class(obj)=="GDALobj"){
    stop("Object of class \"GDALobj\" required.")
  }
  if(missing(colour_scale)){ 
    colour_scale <- SAGA_pal[[1]] 
  }
  if(missing(z.lim)&is.null(breaks.lst)){
    try( z.lim <- c(attr(obj, "df")$Bmin, attr(obj, "df")$Bmax) )
    if(!is.numeric(z.lim)&length(z.lim)==2){ stop("Please specify 'z.lim' argument") }
  }
  if(!is.null(breaks.lst)&length(breaks.lst)<15){
    stop("'breaks.lst' must contain at least 15 elements")
  }
  if(is.na(CRS)){
    stop("'projection' missing or 'NA'")
  }

  if(!length(colour_scale)==(length(breaks.lst)-1)&!is.null(breaks.lst)){ stop("'length(colour_scale)' and 'length(breaks.lst)-1' of equal length required") }
  GDALobj.file <- attr(obj, "file")
  if(is.null(tiles)){
    if(requireNamespace("GSIF", quietly = TRUE)){
      ## if missing get tiling system using block.x:
      tiles <- GSIF::getSpatialTiles(obj, block.x = block.x, return.SpatialPolygons = FALSE)
    }
  }
  if(is.null(tiles.sel)){
    ## if missing do all tiles:
    tiles.sel <- 1:nrow(tiles)
  }
  if(min(tiles.sel)<1|max(tiles.sel)>nrow(tiles)|!is.integer(tiles.sel)){ stop("'tiles.sel' must be: integer within the range of tile numbers (see ?getSpatialTiles)") }
  ## write all tiles:
  if(requireNamespace("parallel", quietly = TRUE)&requireNamespace("snowfall", quietly = TRUE)){  
    if(missing(cpus)){ cpus <- parallel::detectCores(all.tests = FALSE, logical = FALSE) }
    snowfall::sfInit(parallel=TRUE, cpus=cpus)
    snowfall::sfExport("GDALobj.file", "tiles", "tiles.sel", "breaks.lst", "altitude", "altitudeMode", "colour_scale", "z.lim", "overwrite", "CRS")
    snowfall::sfLibrary(package="rgdal", character.only=TRUE)
    snowfall::sfLibrary(package="sp", character.only=TRUE)
    snowfall::sfLibrary(package="plotKML", character.only=TRUE)
    snowfall::sfLibrary(package="XML", character.only=TRUE)
    snowfall::sfLibrary(package="RSAGA", character.only=TRUE)
    snowfall::sfLibrary(package="raster", character.only=TRUE)
    lst <- snowfall::sfLapply(tiles.sel, .kml_SpatialGrid_tile, GDALobj.file=GDALobj.file, tiles=tiles, altitude=altitude, altitudeMode=altitudeMode, colour_scale=colour_scale, breaks.lst=breaks.lst, z.lim=z.lim, overwrite=overwrite, CRS=CRS)
    snowfall::sfStop()
    lst <- do.call(rbind, lst)
  } else {
    warning("To speed up writing of tiles install and load package 'snowfall'")
    lst <- list(NULL)
    for(i in 1:length(tiles.sel)){
      lst[[i]] <- .kml_SpatialGrid_tile(i, GDALobj.file=GDALobj.file, tiles=tiles, altitude=altitude, altitudeMode=altitudeMode, colour_scale=colour_scale, breaks.lst=breaks.lst, z.lim=z.lim, overwrite=overwrite)
    }
  }
  
  if(missing(file.name)){
    file.name <- set.file.extension(basename(GDALobj.file), ".kml")
  }
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
  ## add legend and/or logo:
  if(plot.legend==TRUE){
    kml.legend <- paste0(strsplit(file.name, ".kml")[[1]][1], "_legend.png")
    if(is.null(breaks.lst)){
      kml_legend.bar(x=signif(seq(z.lim[1], z.lim[2], length.out=25), 3), legend.file=kml.legend, legend.pal=colour_scale)
    } else {
      breaks.s <- seq(1, length(breaks.lst), length.out=15)
      kml_legend.bar(x=as.factor(signif(breaks.lst[breaks.s], 2)), legend.file=kml.legend, legend.pal=colour_scale[breaks.s])
    }
    kml_screen(image.file = kml.legend, position = "UL", sname = "legend")
  }
  if(!missing(kml.logo)){ kml_screen(image.file = kml.logo, position = "UR", sname = "logo") }
  kml_close(file.name)
  if(open.kml==TRUE){
    kml_View(file.name)
  } else {
    message(paste("Object written to:", file.name))
  }
}

## auxiliary function:
.kml_SpatialGrid_tile <- function(i, GDALobj.file, colour, tiles, breaks.lst, colour_scale, altitude, altitudeMode, z.lim, N.min=4, overwrite, CRS){
  if(any(!c("offset.x", "offset.y", "region.dim.x", "region.dim.y") %in% names(tiles))){ stop("Missing columns in the 'tiles' object. See ?GSIF::tile") }
  kml.tile <- set.file.extension(paste0(strsplit(basename(GDALobj.file), "\\.")[[1]][1], "_T", i), ".kml")
  raster_name = set.file.extension(paste0(strsplit(basename(GDALobj.file), "\\.")[[1]][1], "_T", i), ".png")
  r <- readGDAL(GDALobj.file, offset=c(tiles$offset.y[i], tiles$offset.x[i]), region.dim=c(tiles$region.dim.y[i], tiles$region.dim.x[i]), silent=TRUE)
  proj4string(r) = CRS
  if(sum(!is.na(r@data[,1]))>N.min){
    prj.check <- check_projection(r, control = TRUE)
    if(!prj.check) { suppressMessages( r <- reproject(r) ) }
    if(is.numeric(r@data[,1])){
      if(!is.null(breaks.lst)){
        if(missing(colour)){
          r@data[,"colour"] <- cut(r@data[,1], breaks=breaks.lst, include.lowest=TRUE, right=FALSE)
        }
      } else {
        names(r)[1] <- "colour" 
      }
      pw = r@grid@cells.dim[1]*3
      ph = r@grid@cells.dim[2]*3
      ## get the bounding box:
      north <- r@bbox[2,2]
      south <- r@bbox[2,1]
      east <- r@bbox[1,2]
      west <- r@bbox[1,1]
      if(!(overwrite==FALSE&file.exists(raster_name))){
        kml_open(kml.tile)
        if(!is.null(z.lim)){
          suppressMessages( kml_layer(r, colour=colour, colour_scale=colour_scale, altitude=altitude, altitudeMode=altitudeMode, raster_name=raster_name, plot.legend=FALSE, png.width=pw, png.height=ph, z.lim=z.lim) )
        } else {
          suppressMessages( kml_layer(r, colour=colour, colour_scale=colour_scale, altitude=altitude, altitudeMode=altitudeMode, raster_name=raster_name, plot.legend=FALSE, png.width=pw, png.height=ph) )
        }
        kml_close(kml.tile)
      }
    } else {
      kml.tile <- NULL
      north <- NULL
      south <- NULL
      east <- NULL
      west <- NULL
    }
    out <- data.frame(kml.tile=kml.tile, north=north, south=south, east=east, west=west)
    return(out)
  }
}
