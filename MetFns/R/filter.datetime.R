filter.datetime<-function(data,year,month.beg, month.end=month.beg, day.beg,day.end=day.beg,  time.beg=0, time.end=2359)
{ 
  if(!is.data.frame(data) || 
     !is.numeric(c(year,month.beg,month.end,day.beg,day.end,time.beg,time.end))|| 
     !(year%in%1984:2020) || any(c(month.beg,month.end,day.beg,day.end)<1) || 
     any(c(month.beg,month.end)>12)|| any(c(time.beg,time.end)<0) || any(c(time.beg,time.end)>2359))
     stop("invalid input parameter(s) specification: check 
           data/year/month.beg/month.end/day.beg/day.end/time.beg/time.end ")

 

  max1<-max.day(year,month.beg)
  max2<-max.day(year,month.end)

  if(day.beg>max1 || day.end>max2 || (month.end==month.beg && day.end<day.beg))
     stop("invalid input parameter(s) specification: check month.beg/month.end format")


  filter.sol(data, solar.long(year,month.beg,day.beg,dec.time(time.beg)), 
             solar.long(year,month.end,day.end,dec.time(time.end)))
  
}