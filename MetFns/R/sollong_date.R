sollong_date<-function(year,value, month.beg=1,month.end=12,day.beg=1,day.end=31,time.beg=0,time.end=2359)
{
  t1<-86400*(yday(ymd(paste(as.character(year),"-",as.character(month.beg),"-",as.character(day.beg),sep="")))-1)+3600*dec.time(time.beg)
  t2<-86400*(yday(ymd(paste(as.character(year),"-",as.character(month.end),"-",as.character(day.end),sep="")))-1)+3600*dec.time(time.end)

  sol<-function(t){
     x<-as.character(as.POSIXct(t, origin = paste(as.character(year),"-01-01", sep=""),tz="UTC"))
     solar.long(year,as.numeric(substr(x,6,7)),as.numeric(substr(x,9,10)), 
                ifelse(nchar(x)<12,0,dec.time(as.numeric(paste(substr(x,12,13),
                       substr(x,15,16),sep="")))+as.numeric(substr(x,18,19))/3600),prec=11)-value}


  years<-1984:2017
  seconds<-c(6842968,6778871,6801102,6823412,6845580,6781486,6803852,6825734,6847739,
             6783726,6805704,6827620,6849623,6785410,6808050,6830233,6852326,6788426,
             6810351,6832260,6854459,6790090,6812540,6834409,6856218,6792465,6814620,
             6836792,6859193,6794786,6817197,6839165,6860952,6797149)
  if(sol(t2)<0) t2<-seconds[year==years]
  if(sol(t1)>0) t1<-seconds[year==years]+5

  floor_date(as.POSIXct(uniroot(sol,lower=t1, upper=t2)$root, 
             origin = paste(as.character(year),"-01-01", sep=""),tz="UTC"),unit="minute")
}


