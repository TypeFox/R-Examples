
pieSP<-function(SPDF,
                  zcol=1:length(SPDF@data),
                  scalelist=TRUE,  # TRUE proportional, FALSE pie charts same size
                  max.radius=100,  #m
                  do.sqrt = TRUE
){
  
  ####  FUNCTION modified from plotGoogleMaps package
  
  colPalette=rainbow(ncol(SPDF@data[,zcol]))
  fillColor=rainbow(length(zcol))
  strokeColor="red"
  
  createSphereSegment <- function(partOfSP,
                                  max.radius=100,  #m
                                  key.entries = as.numeric(partOfSP@data),
                                  scalelist=1,
                                  do.sqrt = TRUE,
                                  fillColor=rainbow(length(key.entries)),
                                  strokeColor="red",
                                  id=length(key.entries)) {
    
    center=coordinates(partOfSP)
    fillColor<-as.character(substr(fillColor,1,7))
    
    obj = as(partOfSP, "SpatialPointsDataFrame")
    z = as.numeric(partOfSP@data)
    # avoid negative values
    #     if (min(key.entries)<0 ){
    #       ke<-abs(min(key.entries))+ key.entries+mean(key.entries)
    #     }else{ke<-key.entries+min(key.entries)}     # no zeros for radius vecor
    ke<-z
    # creating a vector for subgroupes
    if(do.sqrt){
      scale.level<- sqrt(ke/(sum(ke)) ) }else{scale.level<-ke/(sum(ke))}
    radius.level<-max.radius*scale.level
    # list of radiuses for createSphereCircle
    radius.vector<-   radius.level
    dfi<-cbind(rep(NA,1+length(radius.vector)))
    dfi[1]=0
    
    if(max(scale.level)==0)
    {scale.level=0.1
     scalelist=0.01}
    
    dfi[2:(length(radius.vector)+1)]=360/sum(scale.level)*  scale.level
    
    
    fi= cbind(rep(NA,2*length(radius.vector)) )
    dfi=as.numeric(dfi)
    for (i in (seq(2,length(fi),by=2))){
      fi[i]=sum(dfi[1:(i/2+1)])
    }
    fi[1]=0
    for (i in (seq(3,length(fi),by=2))){
      fi[i]=fi[i-1]
    }
    
    
    radius.vector=scalelist*max.radius/6378000
    
    lat1 <- (pi/180)* center[2];
    lng1 <- (pi/180)* center[1];
    
    paths<-as.list(rep(NA,length(radius.vector)))
    
    for(ik in (seq(2,length(fi),by=2))){
      radius<-  radius.vector
      
      coords<-cbind(rep(NA,11),rep(NA,11))
      coords[1,]  <- c(  center[,1], center[,2])
      j<-2
      for (i in seq(fi[ik-1],fi[ik],length.out=10)) {
        tc <- (pi/180)*i
        y <- asin(sin(lat1)*cos(radius)+cos(lat1)*sin(radius)*cos(tc))
        dlng <- atan2(sin(tc)*sin(radius)*cos(lat1),cos(radius)-sin(lat1)*sin(y))
        x <- ((lng1-dlng+pi) %% (2*pi)) - pi
        
        coords[j,1] <-c(x*(180/pi))
        coords[j,2] <-c(y*(180/pi))
        j<-j+1
      }
      
      
      
      lonlat<-coords
      paths[[ik/2]]<-rbind(lonlat,lonlat[1,])
    }
    pol=paths
    lonlat<-rbind( paths[[1]],paths[[1]][1,])
    pol[[1]]<-Polygon(lonlat,hole=FALSE)
    pol[[1]]<-Polygons(list(pol[[1]]),ID=id-1)
    
    for(i in (2:(length(paths)))){
      lonlat<-rbind( paths[[i]],paths[[i]][1,])
      pol[[i]]<-Polygon(lonlat,hole=FALSE)
      pol[[i]]<-Polygons(list(pol[[i]]),ID=id-i)   }
    
    return(pol)
    
  }
  
  SP <-as(SPDF, "SpatialPointsDataFrame")
  projekcija=SP@proj4string
  SP.ll <- spTransform(SP, CRS("+proj=longlat +datum=WGS84"))
  SP.ll<-SP.ll[,zcol]
  
  if(strokeColor!=""){
    rgbc<-col2rgb(strokeColor)
    strokeColor<-rgb(rgbc[1],rgbc[2],rgbc[3],maxColorValue=255) }
  
  if(!is.null(colPalette)){
    rgbc<-col2rgb(colPalette)
    colPalette<-apply(rgbc,2,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))}
  if(scalelist){
    xdata<-SP@data[,zcol]
    ## Existing code commented out
    #xdata <- apply(xdata, 2L, function(x) (x - min(x, na.rm = TRUE))/diff(range(x, na.rm = TRUE)))
    #xsum <- apply(xdata, 1L,function(x) ( sum(x)))
    
    ## Sum across all columns by row to get totals of all "pie" slice values
    ## Dividing into maximum total value gives scaling factor where largest pie
    ## has a scaling factor of 1
    xsum <- apply(xdata, 1L, function(x) (sum(x,na.rm=TRUE)))
    scalelist<-xsum/max(xsum)
    scalelist<-sqrt(scalelist)} else{ scalelist<-rep(1,length(SP.ll@coords[,1])) }
  
  Pols<-as.list(rep(NA,length(SP.ll[,1])))
  Srl<-Pols
  num=(rep(NA,length(zcol)*length(SP.ll@data[,1])) )
  for(i in 1:length(SP.ll@data[,1])){
    num[i*length(zcol)-(0:(length(zcol)-1))]=i
  }
  
  Pols<- lapply( 1:length(SP.ll@data[,1]), function(i) {
    createSphereSegment(SP.ll[i,zcol],
                        max.radius=max.radius,  #m
                        scalelist= scalelist[i],
                        do.sqrt = do.sqrt,
                        fillColor=colPalette,
                        strokeColor= strokeColor,
                        id=i*length(as.numeric(SP.ll@data[i,zcol]))) } )
  
  id=num  # rep(1:length(zcol),length(SP.ll@data[,1]))
  dat=data.frame(id)
  SP$id=1: length(SP.ll@data[,1])
  dat=merge(dat,SP@data)
  Pols=unlist(Pols)
  SPl<-SpatialPolygons(Pols,proj4string=SP.ll@proj4string) 
  SPldf<-SpatialPolygonsDataFrame(SPl,dat,match.ID = FALSE)
  SPldf <- spTransform(SPldf, projekcija)
  return(SPldf)
}

###########################  END OF FUNCTION pie ###############################
#
