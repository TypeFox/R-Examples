midint<-function(data)
{
  if(!is.data.frame(data)) 
     stop("invalid input parameter specification: check data")
    
  time<-cbind(dec.time(data$start),dec.time(data$stop))
  midtime<-apply(time,1,mean)
  ind1<-time[,1]>time[,2]&midtime<12&abs(midtime-12)>0.00001
  ind2<-time[,1]>time[,2]&(midtime>12 | abs(midtime-12)<0.00001)
  
  
  midtime[ind1]<-midtime[ind1]+12
  midtime[ind2]<-midtime[ind2]-12
  midtime=ifelse(midtime<0,0,midtime)
  
  day.new<-data$day
  month.new<-data$month
  maxx<-Vectorize(max.day) 
  maxday<-maxx(data$year,data$month)

  
  day.new[ind2]<-day.new[ind2]+1
  
  ind3<-day.new>maxday
  day.new[ind3]<-day.new[ind3]%%maxday[ind3]
  month.new[ind3]<-month.new[ind3]+1 
  
  cbind(midtime,day.new,month.new)


}
