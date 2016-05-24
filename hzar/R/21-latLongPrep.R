## Calculate the distance between two points on the earth.

hzar.map.greatCircleDistance <-  function( lat1, long1, lat2, long2,units="Km",degrees=TRUE){
  latLong.d(latLong.p(lat1,long1,degrees),
            latLong.p(lat2,long2,degrees),
            units);
}

latLongRE.NE <- "n(orth)?|e(ast)?"
latLongRE.SW <- "w(est)?|s(outh)?"
latLongRE.NESW <- "n(orth)?|e(ast)?|w(est)?|s(outh)?"

hzar.map.dms2deg <- function(deg,min,sec,dir){
  
  if(!all(grepl(latLongRE.NESW,dir,ignore.case=TRUE)))
    stop(paste("Direction",dir[which(!grepl(latLongRE.NESW,dir,ignore.case=TRUE))],"not recognized."))
  
  res <- deg+min/60+sec/3600
  if(length(dir)==1){
    if(grepl(latLongRE.NE,dir,ignore.case=TRUE))
     return(res)
    return(-res)
  }
  ifelse(grepl(latLongRE.NE,dir,ignore.case=TRUE),res,-res)
}

pdms2deg <- function(x){
  if(length(x)==2)
    return(hzar.map.dms2deg(as.numeric(x[1]),0,0,x[2]))
  if(length(x)==3)
    return(hzar.map.dms2deg(as.numeric(x[1]),as.numeric(x[2]),0,x[3]))
  if(length(x)==4)
    return(hzar.map.dms2deg(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[3]),x[4]))
  stop(paste(paste("[",x,"]",collapse=" "),"not understood."))
}
  

hzar.map.latLong.dms <- function(coordinates){
  ## cat("A\n")
##   print(coordinates)
  res <- lapply(strsplit(coordinates,"(?i)(?<=[NSEW]) +",perl=TRUE), strsplit,"-| +")
  result <- matrix(NA,nrow=length(coordinates),ncol=max(sapply(res,length)))
##   print(result)
  res <- lapply(res,lapply,pdms2deg)
  for(iter in 1:length(coordinates)){
##     print(as.numeric(res[[iter]]))
    result[iter,1:length(res[[iter]])] <- as.numeric(res[[iter]])
  }
  return(result)
  
}
hzar.map.latLongSites <- function(siteIDs,site.lat,site.long,degrees=TRUE){
  cbind(site=siteIDs,latLong.p(site.lat,site.long,degrees))
}
hzar.map.latLongSites.dms <- function(siteIDs,coordinates){
  if(!is.character(coordinates))
    stop("Argument coordinates is not a character vector.")
  if(!all(valid <- grepl("^\\d+(\\.\\d*)?((-| )\\d+(\\.\\d*)?((-| )\\d+(\\.\\d*)?)?)?(-| )[NS]\\w* \\d+(\\.\\d*)?((-| )\\d+(\\.\\d*)?((-| )\\d+(\\.\\d*)?)?)?(-| )[EW]\\w*$",coordinates,ignore.case=TRUE)))
    stop(paste("Line",which(!valid),":",coordinates[which(!valid)],"malformed.",collapse="\n"))
  
##   print(coordinates)
  res <- hzar.map.latLong.dms(coordinates)
##   print(res)
##   if(length(siteIDs)!=nrow(res))
##     stop(
  cbind(site=siteIDs,latLong.p(res[,1],res[,2],TRUE))
}
hzar.map.distanceFromSite <-  function(latLongSites,site0,units="Km"){
  if(is.character(site0)&&site0[[1]] %in% latLongSites$site){
    i=which(latLongSites$site %in% site0[1]);
    return(latLong.d(latLongSites[i[1],],latLongSites,units))
  }else{
    stop(paste(site0,"not listed in {",paste(latLongSites$site,collapse=", "),"}"))
  }
}

## Useful values
earthD.Km   = 6371.009 #kilometers
earthD.mile = 3958.761 #statute miles
earthD.ntMl = 3440.069 #nautical miles
earthSphd.r = 298.257223563 #WGS84
earthSphd.ep= (2*earthSphd.r -1)/(earthSphd.r-1)^2
## latLongCoors <- function(lat,long){
##   return(cbind(lat=lat,long=long));
## }

latLong.p <- function(lat,long,degrees=TRUE){
  if(!degrees){
    
    lat=180*lat/pi
    long=180*long/pi
  }
##   cat(lat,"\n",long,"\n")
  latD=lat%%360;
  longD=long%%360;
  longD=ifelse((latD>=270) |(latD <=90),
    ifelse(longD>180,longD-360,longD),
    longD-180)
  
  latD=ifelse(latD>=270,
    latD-360,
    ifelse(latD>90,180-latD,latD))
  
  lat=pi*lat/180
  long=pi*long/180
  return(data.frame(lat.rad=lat,long.rad=long,lat.deg=latD,long.deg=longD));
  ##}
  
  ##return(data.frame(lat.rad=lat,long.rad=long,lat.deg=latD,long.deg=longD));
}


latLong.d <- function(p1,p2,units="Km"){
  if(identical(units,"Km")){ 
    R=earthD.Km;
  }else{
    if(identical(units,"miles")){
       R=earthD.mile;
     } else {
       if(identical(units,"nautical")){
         R=earthD.ntMl;
       }else{
         stop(units +"not identifiable");
       }
     }
  }
  dLat=p2$lat.rad-p1$lat.rad;
  dLong=p2$long.rad-p1$long.rad;
  dLong=ifelse(
    dLong>pi,
    dLong-2*pi,
    ifelse(-dLong>pi,
           dLong+2*pi,
           dLong));
  
  reLat1=atan((earthSphd.r -1)*tan(p1$lat.rad) /earthSphd.r )
  reLat2=atan((earthSphd.r -1)*tan(p2$lat.rad) /earthSphd.r )
  cenNum=sqrt((cos(reLat2)*sin(dLong))^2+(cos(reLat1)*sin(reLat2)-cos(reLat2)*sin(reLat1)*cos(dLong))^2);
  cenDen=sin(reLat1)*sin(reLat2)+cos(reLat2)*cos(reLat1)*cos(dLong);
  central <- atan2(cenNum,cenDen);
  lFP <- (reLat1+ reLat2)/2 ;
  lFQ <- (-reLat1+ reLat2)/2 ;
  lFX <- (central-sin(central))*sin(lFP)^2*cos(lFQ)^2/cos(central/2)^2;
  lFY <- (central+sin(central))*cos(lFP)^2*sin(lFQ)^2/sin(central/2)^2;
  
  lambertFormulaeD <- ifelse(central==0,0,R*(central-(lFX+lFY)/(2*earthSphd.r)));
  return( lambertFormulaeD);
  
}
