#' @export
dateToEpiweek<-function(date,format="%Y-%m-%d",firstday="Sunday"){
  date=strptime(date,format=format)

  if(class(date)[2]!="POSIXt" || is.na(date)){
    print("Wrong format for date!")
    break
  }
  if(!(firstday=="Sunday"|| firstday=="Monday")){
    print("Wrong firstday!")
    break
  }
  year=1900+date$year

  jan4=strptime(paste(year,1,4,sep="-"),format="%Y-%m-%d")
  wday=jan4$wday

  wday[wday==0]=7
  wdaystart=ifelse(firstday=="Sunday",7,1)
  if(wday== wdaystart) weekstart=jan4
  if(wday!= wdaystart) weekstart=jan4-(wday-ifelse(firstday=="Sunday",0,1))*86400

  weeknum=ceiling(as.numeric((difftime(date,weekstart,units="days")+0.1)/7))
  mday=date$mday
  wday=date$wday

  year=ifelse(weeknum==53 & mday-wday>=(ifelse(firstday=="Sunday",29,28)),year+1,year)
  weeknum=ifelse(weeknum==53 & mday-wday>=(ifelse(firstday=="Sunday",29,28)),1,weeknum)
  year.shift=year-1
  jan4.shift=strptime(paste(year.shift,1,4,sep="-"),format="%Y-%m-%d")
  wday=jan4.shift$wday
  wday[wday==0]=7

  wdaystart=ifelse(firstday=="Sunday",7,1)
  if(wday== wdaystart) weekstart=jan4.shift
  if(wday!= wdaystart) weekstart=jan4.shift-(wday-ifelse(firstday=="Sunday",0,1))*86400

  weeknum.shift=ceiling(as.numeric((difftime(date,weekstart)+0.1)/7))
  year=ifelse(weeknum==0,year.shift,year)
  weeknum=ifelse(weeknum==0,weeknum.shift,weeknum)

  return(list("year"=year,"weekno"=weeknum))
}

#' @export
epiweekToDate<-function(year,weekno,firstday="Sunday"){
  if(!(firstday=="Sunday"|| firstday=="Monday")){
    print("Wrong firstday!")
    break
  }
  if(year<0 || weekno<0){
    print("Wrong Input!")
    break
  }

  jan4=strptime(paste(year,1,4,sep="-"),format="%Y-%m-%d")
  wday=jan4$wday
  wday[wday==0]=7
  wdaystart=ifelse(firstday=="Sunday",7,1)
  if(wday== wdaystart) weekstart=jan4
  if(wday!= wdaystart) weekstart=jan4-(wday-ifelse(firstday=="Sunday",0,1))*86400

  jan4_2=strptime(paste(year+1,1,4,sep="-"),format="%Y-%m-%d")

  wday_2=jan4_2$wday
  wday_2[wday_2==0]=7
  wdaystart_2=ifelse(firstday=="Sunday",7,1)
  if(wday_2== wdaystart_2) weekstart_2=jan4_2
  if(wday_2!= wdaystart_2) weekstart_2=jan4_2-(wday_2-ifelse(firstday=="Sunday",0,1))*86400

  if(weekno>((weekstart_2-weekstart)/7)){
    print(paste("There are only ",(weekstart_2-weekstart)/7," weeks in ",year,"!",sep=""))
    break
  }

  d0=weekstart+(weekno-1)*7*86400
  d1=weekstart+(weekno-1)*7*86400+6*86400

  return(list("d0"=strptime(d0,format="%Y-%m-%d"),"d1"=strptime(d1,format="%Y-%m-%d")))
}
