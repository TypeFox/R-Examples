sinh<-function(data,shw, Ralpha=NULL, Delta=NULL)
{
  data.shw<-filter.shw(data,shw)
  ind<-is.na(data.shw$long)| is.na(data.shw$lat)
  datashw=data.shw[!ind,]

  day.new<-midint(datashw)[,2]
  month.new<-midint(datashw)[,3]
  midtime<-midint(datashw)[,1]
  

  if(is.null(Ralpha) || is.null(Delta)){
    m<-nrow(datashw)
    Ralpha<-rep(0,m)
    Delta<-rep(0,m)
    data(radiant,envir=environment())
    radiant<-get("radiant",envir=environment())
    i<-which(substr(names(radiant),start=1,stop=3)==shw)[1]
    
    for(j in 1:m){
        ind<-day.new[j]==radiant$Day&month.new[j]==radiant$Month
        Ralpha[j]<-radiant[ind,i]
        Delta[j]<-radiant[ind,i+1] }
  }    
  if(any(is.na(Ralpha))) 
    stop("no table values for radiant coordinates")

  t<-hms2rad(ut2ha(yr=datashw$year, mo=month.new, dy=day.new, hr=midtime, ra.sou=paste(as.character(Ralpha*24/360),"h",sep=""), 
                    lon.obs=paste(datashw$EW,paste(datashw$long,"d",sep=""),sep=" ")))
  sine.h<-sin(dms2rad(Delta))*sin(dms2rad(datashw$lat))+cos(dms2rad(Delta))*cos(dms2rad(datashw$lat))*cos(t)
  data.frame(datashw,sine.h=round(sine.h,3))
}
