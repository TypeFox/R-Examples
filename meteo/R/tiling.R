tiling <- function(filename="", # path to grid file in SAGA format / raster formats 
                   tilesize=500, # in cells nx= ny
                   overlapping=50, # in cells
                   aspoints= FALSE,
                   asfiles=FALSE,                   
                   tilename="tile", # prexif to be given to tile names
                   format="GTiff",
                   tiles_folder=paste(getwd(),'tiles',sep='/'), # resulting folder
                   parallel.processing=FALSE,
                   cpus=6) {
  
if(class(filename)=="RasterLayer") {r=filename} else{r= raster(filename) }
  
  bb <- extent(r)
  bb <- cbind(c(bb@xmin,bb@xmax),c(bb@ymin,bb@ymax))
  
  cel <- res(r) [1]
  
  bb[1,]=bb[1,]-overlapping*cel
  bb[2,]=bb[2,]+overlapping*cel
  
 
  
  step=tilesize*cel+overlapping*cel
  nlon=ceiling(diff(bb[,1]) / step )
  nlat=ceiling(diff(bb[,2]) / step )
  
  
  p.l <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmin=seq(bb[1,1],bb[1,1]+(nlon-1)*step,   by=step), latmin=seq(bb[1,2],bb[1,2]+(nlat-1)*step  ,by=step)) - overlapping*cel    
  p.l$lonmin[p.l$lonmin<extent(r)@xmin] = extent(r)@xmin
  p.l$latmin[p.l$latmin<extent(r)@ymin] = extent(r)@ymin
  p.u <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmax=seq(bb[1,1]+step, bb[1,1]+nlon*step   ,by=step), latmax=seq(bb[1,2]+step   ,bb[1,2]+nlat*step,  by=step))  +overlapping*cel  
  p.u$lonmax[p.u$lonmax>extent(r)@xmax] = extent(r)@xmax
  p.u$latmax[p.u$latmax> extent(r)@ymax] = extent(r)@ymax
  
  ptiles <- cbind(p.l, p.u)
  
  poligoni=as.list(rep(NA,length(ptiles$lonmin)))
  
  for(i in 1:length(ptiles$lonmin)) {
    
    x <- rbind(ptiles$lonmin[i], ptiles$lonmax[i],ptiles$lonmax[i],ptiles$lonmin[i],ptiles$lonmin[i])
    y <- rbind( ptiles$latmin[i], ptiles$latmin[i], ptiles$latmax[i], ptiles$latmax[i],ptiles$latmin[i])
    
    poligoni[[i]]<-  Polygons( list(Polygon(cbind(x,y))) ,i)  
  }
  
  poll<-SpatialPolygons(poligoni)
  
  tiles=as.list(rep(NA,length(poll@polygons)))
  
  
  
  if(parallel.processing){
    if(!sfParallel()){
      sfInit ( parallel = TRUE, cpus =cpus) 
      sfst=FALSE}else{sfst=TRUE}
    
    sfLibrary(package="raster", character.only=TRUE)
    sfLibrary(package="rgdal", character.only=TRUE)
    
    sfExport( "poll" )
    sfExport( "r" )
    sfExport( "tiles_folder" )
    sfExport( "tilename" )
    sfExport( "format" )
    if(asfiles){
        if(!dir.exists(tiles_folder)){
            dir.create(tiles_folder)
        }
        tiles <- sfLapply(1:length(poll@polygons), function(i)  {
            crop(r, poll[i,]);
            writeRaster(tiles[[1]], filename= paste(tiles_folder,"/",tilename,i,".",format,sep=""), format=format, overwrite=T) } )
    }else{
      tiles <- sfLapply(1:length(poll@polygons), function(i)  crop(r, poll[i,] )  )
    }

  }else{
    if(asfiles){
        if(!dir.exists(tiles_folder)){
            dir.create(tiles_folder)
        }
        for(i in 1:length(poll@polygons)){ tiles[i] <- crop(r, poll[i,]);
                                           writeRaster(tiles[i][[1]], filename= paste(tiles_folder,"/",tilename,i,".",format,sep=""), format=format, overwrite=T)
        }
    }else{
        for(i in 1:length(poll@polygons)){ tiles[i] <- crop(r, poll[i,])}
    }
  }

  
  if(aspoints){
    rem <- sapply(tiles, function(i) all(is.na(i@data@values) ) )
    tiles = tiles[!rem]
    tiles= lapply(tiles, function(i) rasterToPoints(i,spatial=TRUE) )
   
  }

  if(parallel.processing) { sfStop() }
  
return(tiles)
}

  

  

