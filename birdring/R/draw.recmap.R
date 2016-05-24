draw.recmap <-
function(data, file='screen', subset=NULL, group=NULL, bbox=NULL, projection='longlat', zoom=2, margin=1, 
         points=0, pcol='red', lines=0, lcol='red', gcircle=FALSE, density=FALSE, grid.size=300, mask='sea', all.data=FALSE,          
         height=1000, width=0, units='px', col='lemonchiffon', bg='lightblue1', border=NA, lwd=1, 
         legend='none', title=NULL, labels=NULL, lcex=0, alpha=0.2, ...){
  
  if(!is.null(bbox)) bbox <- bbox[c(3,4,1,2)] # change the order to be consistent with the draw.map function  

  # first check lat/long columns exist and subset data
  if( !all(c('lat', 'lon') %in% names(data)) ) 
    stop("data should have  columns labelled lat and long")  
  if ( !is.null(subset) ) # subset data if required
    data <- subset(data, eval(parse(text=subset)))
  if( nrow(data)==0 )
    stop('no data selected')
 
  # set-up the projections
  if( tolower(substr(projection,1,4)) == 'long' | tolower(substr(projection,1,3)) == 'lat')    
    p4s <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  else if( tolower(substr(projection,1,4)) == 'merc' ) # this one used by OpenStreetMap    
    p4s <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"
  else
    stop('unrecognised projection')

  # load the map data
  if( zoom==1 || substr(tolower(zoom), 1, 1)=='c' ) {   
    requireData("rworldmap") #data(countriesCoarse, envir = environment())
    countriesCoarse <- countriesCoarse
    worldmap <- spTransform(countriesCoarse, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  } else if( zoom==2 || substr(tolower(zoom), 1, 1)=='l' ){
    requireData("rworldmap") #data(countriesLow, envir = environment())
    countriesLow <- countriesLow
    worldmap <- spTransform(countriesLow, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  } else if( zoom==3 || substr(tolower(zoom), 1, 1)=='h' ){
    requireData("rworldxtra") #data(countriesHigh, envir = environment())
    countriesHigh <- countriesHigh
    worldmap <- spTransform(countriesHigh, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  } 
  
  # work out the bounding box for the map
  mapdata <- data # may need this for the density contours
  if( is.null(bbox) )
    # bbox format lat_min, lat_max, lon_min, lon_max
    bbox <- c( min(data$lat)-margin, max(data$lat)+margin, min(data$lon)-margin, max(data$lon)+margin )
  else { # only plot data within bbox and warn
    if( length(bbox)!=4 )
      stop('Map limits need to be as c(lat_min, lat_max, lon_min, lon_max)', call.=FALSE)
    outwith <- c(which(data$lat<bbox[1] | data$lat>bbox[2]), 
                 which(data$lon<bbox[3] | data$lon>bbox[4]))
    if( length(outwith)>0 ){
      warning('some data are outside map limits')
      data <- data[-c(outwith), ]
      if( all.data==FALSE )
        mapdata <- data
    }
  }
  corners <- matrix(c(bbox[3], bbox[1], bbox[4], bbox[2]),nrow=2) # x1, y1, x2, y2
  
  # clip the map to the bounding box
  clip.extent <- as(extent(corners[1], corners[3], corners[2], corners[4]), "SpatialPolygons")
  proj4string(clip.extent) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  gI <- gIntersects(worldmap, clip.extent, byid=TRUE )
  out <- lapply(which(gI), function(x){ gIntersection(worldmap[x,], clip.extent) })
  map <- SpatialPolygons(lapply(1:length(out), function(i){ Pol<-slot(out[[i]], "polygons")[[1]]; slot(Pol, "ID")<-as.character(i); Pol}))
  proj4string(map) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  mapt <- spTransform(map, CRS(p4s)) 

  cornerst <- bbox(spTransform(clip.extent, CRS(p4s)))
 
  # now sort out the canvas
  mapasp <- (cornerst[4]-cornerst[2]) / (cornerst[3]-cornerst[1])
    
  if ( height>0 & width>0 ){
    warning('do not give both height and width, height will be ignored')
    height <- 0
  }
  if ( height > 0 )
    width <- floor(height / mapasp)
  else 
    height <- floor(width * mapasp)
  
  if ( !is.null(file) ){ # work file output type using extension
    ftype <- tolower(unlist(strsplit(file,'\\.'))[2])
    if ( is.na(ftype) ) 
      ftype <- 'stdio'
    else if ( ftype=='png' ) 
      png(file, width, height, units)
    else if ( ftype=='jpg' | ftype=='jpeg' )
      jpeg(file, width, height, units)
    else if ( ftype=='tif' | ftype=='tiff' )
      tiff(file, width, height, units)
    else 
      warning(paste('unrecognised file type:', ftype))
  } else 
    ftype <- 'stdio' # is.null(file) so send to screen
  
  op <- par(mfrow=c(1,1), plt=c(0,1,0,1), mai=c(0,0,0,0)) # make sure plot covers frame   
  xlim <- cornerst[c(1,3)] # make sure the map covers the data
  ylim <- cornerst[c(2,4)]

  newborder <- ifelse(density, 'transparent', border) # setting border to transparent now makes the lines slightly neater
  parBB <- ifelse(ftype=='stdio', FALSE, TRUE) # graphics device may not be right aspect ratio, file should be
  plot(mapt, xlim=xlim, ylim=ylim, col=col, bg=bg, border=newborder, pbg=bg, asp=TRUE, setParUsrBB=parBB)
    
  if( density==TRUE ){
       
    ddata <- SpatialPoints(mapdata[ ,c('lon','lat')], CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
    ddata <- spTransform(ddata, CRS(p4s))
    map.density2(ddata, grid.size)

    # now mask the land/sea
    sea_mask <- gDifference(clip.extent, gUnionCascaded(map))
    sea_maskt <- spTransform(sea_mask, CRS(p4s))
    mcol <- ifelse(mask=='sea', bg, col)  # determine mask colour
    
    if( mask=='land' )
      plot(mapt, col=mcol, border='transparent', add=TRUE)
    else
      plot(sea_maskt, col=mcol, border='transparent', add=TRUE)
    # overplot the borders so they look nicer 
    plot(mapt, col='transparent', bg='transparent', border=border, lwd=lwd, add=TRUE)
  }
  
  if( points>0 ){
    pdata <- as.data.frame(data[,c('lat', 'lon', group)])
    
    pcol1 <- pcol
    if( !is.null(group) ){
      pdata[ ,3] <- as.factor(pdata[ ,3]) # ensure grouping var is a factor
      glevels <- levels(pdata[ ,3])

      if( length(pcol) != length(glevels) ){
        warn <- paste('number of groups is', length(glevels),'so pcol recycled to:', paste(pcol,collapse=', '), sep=' ')
        warning(warn, call.=FALSE)
        pcol1 <- rep(pcol, length.out=length(glevels))
      }
      pcol1 <- pcol[as.factor(pdata[ ,3])]
    }
    
    # convert to SpatialPointsDataFrame
    coordinates(pdata) <- ~lon+lat
    proj4string(pdata) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    pdata <- spTransform(pdata, CRS(p4s))

    points(pdata, pch=19, col=pcol1, cex=points)
  }

  if( lines>0 ){
    ldata <- data[ ,c('scheme', 'ring', 'lat' , 'lon', group)]
    ldata$birdid <- factor(paste(ldata$scheme, ldata$ring))
    
    lcol1 <- lcol
    if( !is.null(group) ){
      ldata[ ,5] <- as.factor(ldata[ ,5]) # ensure grouping var is a factor
      glevels <- levels(ldata[ ,5])
      
      if( length(lcol) != length(glevels) ){
        lcol1 <- rep(lcol, length.out=length(glevels))
        warn <- paste('number of groups is', length(glevels),'so lcol recycled to:', paste(lcol1,collapse=', '), sep=' ')
        warning(warn, call.=FALSE)
      }
      lcol1 <- lcol1[as.factor(ldata[unique(ldata$birdid) ,5])] # pulls out the first level with each line
    }

    x1 <- list() # create SpatialLines object
    for( ind in levels(ldata$birdid) ){
      x <- ldata[ldata$birdid==ind, ]
      x1[[ind]] <- Line(cbind(x$lon, x$lat))
      x1[[ind]] <- Lines(list(x1[[ind]]), ID=as.character(ind))   
    }

    if( gcircle ){ # transform to great circle routes
      for( ind in 1:length(x1) ){
        track <- coordinates(x1[[ind]])[[1]]  
        ID <- slot(x1[[ind]],'ID')
        t1 <- matrix(nrow=0, ncol=2)
        for( i in 1:(nrow(track)-1) )
          t1 <- rbind(t1, gcIntermediate(track[i,],track[(i+1),],addStartEnd=TRUE))
        x1[[ind]] <- Lines(list(Line(t1)), ID=ID)   
      }
    }
    
    dataSL <- SpatialLines(x1, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    dataSLm <- spTransform(dataSL, CRS(p4s))
    
    plot(dataSLm, col=lcol1, lwd=lines, add=TRUE)
  }
  
  if( tolower(legend)!='none' ){
    if( !is.null(group) ){
      gp.levels <- levels(as.factor(data[, group]))
      if( length(labels)==length(gp.levels) )
        gp.levels <- labels
    }
    if( is.null(title) ){
      if( is.null(group) )
        title='' # so there is a title
      else
        title=group
    }
     
    if( lcex==0 )
      lcex <- (height*0.025)/par('cra')[2]
    lty <- ifelse(lines>0, 1L, 0L)
    symbol <- ifelse(points>0, 'circle', as.integer(NA))
    if( lines>0 )
      xcol <- lcol
    if( points>0 )
      xcol <- pcol
    map.legend(x=legend, text=gp.levels, col=xcol, symbol=symbol, alpha=alpha, lwd=lwd, lty=lty, cex=lcex, 
               title=title, ...)
  }
  
  par(op)
  if( ftype != 'stdio' )
    dev.off()
  
  invisible(mapt)

}
