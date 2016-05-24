bubbleSP <-
  function(SPDF,
           zcol=1,
           scale_e=1,
           max.radius=100,
           key.entries = quantile(SPDF@data[,zcol],(1:5)/5, na.rm=TRUE),
           do.sqrt = TRUE,
           radius.vector=NULL){
    SP=SPDF
    
    bools=!is.na(SP@data[,zcol[1]])
    
    SP<-SP[which(!is.na(SP@data[,zcol[1]])) ,]
    
    if(is.null(radius.vector) ){
    
    obj = as(SP, "SpatialPointsDataFrame")
    data = obj@data
    if (NCOL(data) == 1){
      z = as.numeric(data[,1])
    }else {
      z = as.numeric( data[, zcol] )  }
    # avoid negative values
    
    kkk=c(min(z),key.entries)
    
    kkk=sapply(2:length(kkk), function(i) mean( c(kkk[i],kkk[i-1]) ) )
    
    ke<-abs(kkk) + mean(abs(kkk))    # no zeros as max for radius vecor
    # creating a vector for subgroups
    if(do.sqrt){scale.level<- sqrt(ke/(max(ke)) ) }else{scale.level<-ke/(max(ke))}
    radius.level<-max.radius*scale.level
    
    # new
    if (key.entries[1] > min(z) ) {
      breakss <- factor(  c(min(z), key.entries) )
    } else {
      breakss <- factor(  c(key.entries) )
    }
    
    break_unique <- as.numeric(levels(breakss))
    
    # new
    if (break_unique[length(break_unique)] < max(z) ) {
      break_unique[length(break_unique)] <- max(z)
    }
    
    #setting of auxiliary variable				
    breaks_used=1:(length(break_unique))
    
    if (length(unique(z)) == length(key.entries)) {
      zz = factor(z, labels = radius.level)
      radius.vector <- floor(as.numeric(as.vector(zz)))
    }  else {
      break_unique2=break_unique
      break_unique2[length(break_unique2)] <- max(z)+1
      #now find out which categories are used ...
      #findInterval(z, break_unique2)
      #data.frame(z, int=findInterval(z, break_unique2), classes=t(break_unique2))
      breaks_used=sort(unique(findInterval(z, break_unique2)))
      #NEW CODE
      zz = factor(cut(z, break_unique, include.lowest = TRUE), labels = radius.level[breaks_used])        
      radius.vector <- floor(as.numeric(as.vector((zz))))
      
    }
    
    
    }else{
      radius.vector=radius.vector[bools]
    }
    
    proj=SP@proj4string
    SP=spTransform(SP,CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84"))

    Pols=lapply(1:length(SP@data[,1]), function (i){
      ellip(a=radius.vector[i]*scale_e,b=radius.vector[i]*scale_e, alpha = radius.vector[i], 
            loc = c(SP@coords[i,1],SP@coords[i,2]), n =50,Id=paste(i) ) } )
    

    SPl<-SpatialPolygons(Pols,proj4string =CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84"))
    SPl=spTransform(SPl,proj)
    SPldf<-SpatialPolygonsDataFrame(SPl,SP@data,match.ID = FALSE)
    return(SPldf)
  }
    
