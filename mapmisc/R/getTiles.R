
getRasterSphere = function(zoom) {
  N = 2^zoom 
  raster(openmapExtentMercSphere, nrows = N, ncols=N, crs=crsMercSphere)
}

getRowCol <- function(extMerc,
    zoom, 
    rasterSphere = getRasterSphere(zoom)){

  Sy=rowFromY(rasterSphere, c(
          ymax(extMerc), ymin(extMerc)
          ))
  Sx=colFromX(rasterSphere, c(
              xmin(extMerc), xmax(extMerc)
      ))
  
  
  list(
  col = seq(Sx[1],Sx[2]),
  row = seq(Sy[1],Sy[2])
) 
}

nTilesMerc <- function(extMerc,zoom){
  
  SrowCol = getRowCol(extMerc, zoom=zoom)
  
  length(SrowCol[[1]])*length(SrowCol[[2]])
  
}


getTilesMerc = function(
    extMerc=openmapExtentMercSphere, 
    zoom=1, 
    path="http://tile.openstreetmap.org/",
    cacheDir=file.path(tempdir(), paste("X",make.names(path),sep="")),
    verbose=FALSE){
  
  rasterSphere = getRasterSphere(zoom)  
  
  SrowCol = getRowCol(extMerc, rasterSphere=rasterSphere)
  
  rasters = list()
  colourtable = NULL
  
  
  for(Dx in SrowCol$col){

    rasterCol= list()
    Dcache = file.path(cacheDir, zoom, Dx-1)
    Dpath = paste(path,zoom,'/',Dx-1,'/',sep='')
    dir.create(Dcache,recursive=TRUE,showWarnings=FALSE)
    for(Dy in SrowCol$row) {
      
      Dcell = cellFromRowCol(rasterSphere, Dy, Dx)
      Dextent = extent(rasterFromCells(rasterSphere, Dcell))
      
      Dtile = paste(Dy-1, '.png',sep='')
      Dfile = file.path(Dcache, Dtile)
      Durl = paste(Dpath, Dtile, sep='')
      
      Dsize = file.info(Dfile)['size']
      if(!any(Dsize > 0,na.rm=TRUE)) {
        try(utils::download.file(
						Durl, Dfile, quiet=!verbose, 
						method='internal', mode = 'wb'
		), silent=TRUE)
      } else {
        if(verbose) cat("tile ", Dfile, " cached\n")
      }

      thisimage = try(raster::brick(Dfile), silent=TRUE)
      
      if(class(thisimage)=='try-error') {
        if(verbose) warning("tile ", Dfile, " cannot be loaded")
        
        thisimage = raster(
            Dextent, nrows=256, ncols=256, crs=crsMercSphere
            )
        values(thisimage) = NA
      } else {
        crs(thisimage) = crsMercSphere
        extent(thisimage) = Dextent
      }
      
      if(nlayers(thisimage)==1) {
        # single layer, there's a colortable
        newcolours = thisimage@legend@colortable
        thisimage=thisimage[[1]]
        if(length(newcolours)) {
          if(any(is.na(values(thisimage)))) {
            newcolours[1] = NA
          }
          thisimage = thisimage + length(colourtable)
          colourtable = c(colourtable, newcolours)
          thisimage@legend@colortable = colourtable
        }
        names(thisimage) = gsub("^http://|/$", "", path)
      } else if (nlayers(thisimage)>1){
        
        cnames = c('Red','Green','Blue','Trans')[1:min(c(4,nlayers(thisimage)))]
        
        names(thisimage)[1:length(cnames)] = paste(
            gsub("^http://|/$", "", path),
            cnames,
            sep="")
        
        transLayer = grep("Trans$", names(thisimage), value=TRUE)
        colLayer = grep("Trans$", names(thisimage), value=TRUE,invert=TRUE)
        if(length(transLayer)) {
          # convert transparent to NA
          transLayer = reclassify(thisimage[[transLayer]], data.frame(-Inf, 10, NA))
          thisimage = raster::mask(thisimage, transLayer)
        }
      } # end more than one layer 
      
      rasterCol[[paste('y', Dy,sep='')]] = thisimage
    } # end Dy
    
    if(length(rasterCol)> 1) {
      names(rasterCol)[1:2] = c('x','y')
      rasters[[paste('x',Dx,sep='')]] = do.call(
          raster::merge, rasterCol
      )
      names(rasters[[paste('x',Dx,sep='')]]) = 
          names( rasterCol[[1]])
    } else {
      rasters[[paste('x',Dx,sep='')]] = rasterCol[[1]]
    }
    } # end Dx
    
    # merge the rows
    if(length(rasters) > 1) {
      thenames = names(rasters[[1]])
      
      names(rasters)[1:2] = c('x','y')
      rastersMerged = do.call(
              raster::merge, rasters
          )
          names(rastersMerged) = names(rasters[[1]])
    } else {
      rastersMerged = rasters[[1]]
    }
    
    # add the colortable
    if( length(colourtable)) {
      # re-order colours so most common colours are first
      newtable = unique(colourtable)
      tomatch = match(colourtable,newtable)
      rastersMerged@legend@colortable = newtable
      values(rastersMerged) = tomatch[values(rastersMerged)+1]-1
    }
    attributes(rastersMerged)$tiles = 
        list(tiles = length(SrowCol[[1]])*length(SrowCol[[2]]), 
            zoom=zoom,
            path=path)
    
    return(rastersMerged)	
    
}


