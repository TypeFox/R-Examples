`dateStamp` <-
function(datelist)
  {

if(is.null(datelist$msec)) { datelist$msec = rep(0, length(datelist$sec) ) }
    
    rd = recdate(datelist$jd , datelist$hr , datelist$mi , datelist$sec +datelist$msec /1000, yr=datelist$yr )
    sec = floor(rd$sec)
    msec = 1000*(rd$sec-sec)
    t1 =   (msec-floor(msec))/1000
    msec = floor(msec)
    
    ftime = paste(sep=":", rd$yr,
      formatC(rd$jd, width=3 , flag = "0")  ,
      formatC(rd$hr, width=2 , flag = "0"),
      formatC(rd$mi, width=2 , flag = "0"),
      formatC(sec, width=2 , flag = "0"),
      formatC(msec, width=3 , flag = "0") )
    
    return(ftime)
    
  }

