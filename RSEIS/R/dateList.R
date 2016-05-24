`dateList` <-
function(datevec)
  {
    L = list(yr=datevec[1],
      jd=datevec[2],
      mo=datevec[3],
      dom=datevec[4],
      hr=datevec[5],
      mi=datevec[6],
      sec=datevec[7],
      msec=datevec[8])
    
    return(L)
    
  }


