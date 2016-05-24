cosvarcar=function(aside=50, bside=50, cside=50, regions=1){
  aside=abs(aside); bside=abs(bside); cside=abs(cside)
  temp=sort(c(aside,bside),decreasing = TRUE)
  aside=temp[1];bside=temp[2]
  x=((aside/bside)-1.0)**0.5
  cv=as.numeric((1.0-0.03*x)*(219.7-52.4*log10(aside*bside*291.0) + 3.21*(log10(aside*bside*291.0))^2)/sqrt(regions*cside/291.0))
  return(cv)
}

cosvarsph=function(long = c(129, 141), lat = c(-2, 3), zmax=1, zmin=0, regions=1, inunit='deg', sep=":"){
  if(inunit %in% c('deg','amin','asec','rad','sex')==FALSE){stop('inunit must be one of deg, amin, asec, rad or sex')}
  if(length(long)==1){long=c(0,long)}
  if(length(lat)==1){lat=c(0,lat)}
  fullsky=129600/pi
  if(inunit=='sex'){long=hms2deg(long,sep=sep);lat=dms2deg(lat,sep=sep)}
  if(inunit=='amin'){long=long/60;lat=lat/60}
  if(inunit=='asec'){long=long/3600;lat=lat/3600}
  if(inunit=='rad'){long=long*180/pi;lat=lat*180/pi}
  CoDistLow = cosdistCoDist(z=zmin,H0=70,OmegaM=0.3)     
  CoDistHigh = cosdistCoDist(z=zmax,H0=70,OmegaM=0.3)
  cside=CoDistHigh-CoDistLow
  area=skyarea(long = long, lat = lat, inunit = 'deg', outunit='deg2')[1]
  volume=cosvol(area=area, zmax = zmax, zmin=zmin, H0 = 70, OmegaM = 0.3, inunit='deg2')[1]
  aside=cos(mean(lat)*pi/180)*(abs(diff(long))/360)*(CoDistLow+cside/2)
  bside=(abs(diff(long))/180)*(CoDistLow+cside/2)
  scale=sqrt(volume*1e9/(aside*bside*cside))
  aside=aside*scale
  bside=bside*scale
  return(cosvarcar(aside=aside, bside=bside, cside=cside, regions=regions))
}