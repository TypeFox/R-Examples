`recdate` <-
function(jd=0, hr=0, mi=0, sec=0, yr=0)
{
    if(missing(jd)) jd = 1
    if(missing(hr)) hr = 0
    if(missing(mi)) mi = 0
    if(missing(sec)) sec = 0
    if(missing(yr)) yr = 1972

    N = max(length(jd), length(hr), length(mi), length(sec), length(yr), na.rm=TRUE)
    
    if(length(yr) <= 1 )
      {
        yr = rep(yr, N)
      }

    eday = EPOCHday(yr, jd=jd)

    

    secs = eday$jday*(86400)+hr*(3600)+mi*(60)+sec;

 days = trunc( secs / (86400));
 left1 =  secs - days*(86400);
 hrs = trunc( left1 / (3600));
 left2 =  left1 - hrs*(3600);
 mins = trunc(left2/60.0);
 osec = left2 - mins*60;

    YRJD = EPOCHyear(days, origyr=eday$origyr)

  list( jd=YRJD$jd, hr=hrs, mi=mins, sec=osec, yr=YRJD$yr)
}

