skyarea=function(long=c(129,141),lat=c(-2,3),inunit='deg',outunit='deg2',sep=":"){
  if(inunit %in% c('deg','amin','asec','rad','sex')==FALSE){stop('inunit must be one of deg, amin, asec, rad or sex')}
  if(outunit %in% c('deg2','amin2','asec2','rad2','sr')==FALSE){stop('inunit must be one of deg2, amin2, asec2 or rad2')}
  if(length(long)==1){long=c(0,long)}
  if(length(lat)==1){lat=c(0,lat)}
  fullsky=129600/pi
  if(inunit=='sex'){long=hms2deg(long,sep=sep);lat=dms2deg(lat,sep=sep)}
  if(inunit=='amin'){long=long/60;lat=lat/60}
  if(inunit=='asec'){long=long/3600;lat=lat/3600}
  if(inunit=='rad'){long=long*180/pi;lat=lat*180/pi}
  areafrac=((sin(lat[2]*pi/180)-sin(lat[1]*pi/180))*(long[2]-long[1])/360)/2
  if(areafrac<0){stop('Sky area is non-physical (negative)')}
  if(areafrac>1){stop('Sky area is non-physical (bigger than the full sky)')}
  area=areafrac*fullsky
  if(outunit=='amin2'){area=area*3600}
  if(outunit=='asec2'){area=area*12960000}
  if(outunit=='rad2' | outunit=='sr'){area=area/(180/pi)^2}
  return(c(area=area,areafrac=areafrac))
}