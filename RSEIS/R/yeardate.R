`yeardate` <-
function(yr, jd=1, hr=0, mi=0, sec=0)
  {

    if(missing(jd)) jd =1
    if(missing(hr)) hr=0
    if(missing(mi)) mi=0
    if(missing(sec)) sec=0

#####    A year will be a leap year if it is divisible by 4 but not by 100. 
#####    If a year is divisible by 4 and by 100,
#####     it is not a leap year unless it is also divisible by 400.
    
#############   get the leap years

    n = max(c(length(yr), length(jd), length(hr), length(mi), length(sec)))
    YRDAYS = rep(365, n)
    
    YRDAYS[yr%%4 == 0 &  yr%%100 != 0 ] =366

    
    YRDAYS[yr%%4 == 0 &  yr%%100 == 0 & yr%%400 == 0 ] = 366
    
   ######  cbind(yr, YRDAYS)
    
    jd = jd-1
    d  = yr+jd/YRDAYS+hr/(YRDAYS*24)+mi/(YRDAYS*24*60)+sec/(YRDAYS*24*3600)
    return(d)
  }

