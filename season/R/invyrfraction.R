# invyrfraction.R
# convert fraction of the year into a date (day and month)
# month on a scale of [1,13)
# type =  monthly/weekly/daily
# Jan 2014

invyrfraction<-function(frac,type='daily',text=TRUE){
  n<-length(frac)
  if(sum(frac<0)+sum(frac>1)>0){stop('Fraction must be between 0 and 1')}
  if (type=='daily'){
     yrlength<-365.25; 
     day<-(frac*yrlength)+1;
     day=day-(365*as.numeric(day>365)); # avoid values > 365
     day=pmax(day,1); # avoid values < 1
     date<-strptime(day,'%j');
     day<-as.numeric(format(date,'%d')); # Day of the month as decimal number (01?31)
     month<-format(date,'%B'); # Month name
     if (text==TRUE){daym<-paste('Month =',month,', day =',day)}
     if (text==FALSE){
        monthnum<-as.numeric(format(date,'%m')); # Month number
        mnthlength<-c(31,28.25,31,30,31,30,31,31,30,31,30,31)
        daym<-monthnum+((day-1)/mnthlength[monthnum])
     } # 
  }
  if (type=='weekly'){
    week<-(frac*52)+1;
    if (text==TRUE){daym<-paste('Week =',round(week,1))}
    if (text==FALSE){daym<-week}
  }
  if (type=='monthly'){
     month<-(frac*12)+1;
     if (text==TRUE){daym<-paste('Month =',round(month,1))}
     if (text==FALSE){daym<-month}
  }
  if (type=='hourly'){
    month<-(frac*24)+1;
    if (text==TRUE){daym<-paste('Hour =',round(month,1))}
    if (text==FALSE){daym<-month}
  }
  return(daym)
}

