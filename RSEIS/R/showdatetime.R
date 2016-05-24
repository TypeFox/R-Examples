showdatetime<-function(rd, AMPM=FALSE)
{
  if(missing(AMPM)) { AMPM=FALSE    }

  rd = recdate(jd=rd$jd, hr=rd$hr, mi=rd$mi, sec= rd$sec, yr=rd$yr)

  MODAY = getmoday(rd$jd,  rd$yr)

  gdates = paste(sep="-", formatC(rd$yr, width=3 , flag = "0") ,
    formatC( MODAY$mo , width=2 , flag = "0"),
    formatC( MODAY$dom , width=2 , flag = "0"))

  mysec = floor(rd$sec)
  mymcrsec = floor((rd$sec-mysec)*1000000)

  hours = rd$hr
  A = rep(NULL, length(hours))
  if(AMPM)
    {
      A = rep('AM', length(hours))
      A[hours>=12] = 'PM'
      hours[hours>=13] = hours[hours>=13]-12
    }
  
  gtimes = paste(sep=":",
    formatC(hours, width=2 , flag = "0"),
    formatC(rd$mi, width=2 , flag = "0"),
    formatC(mysec, width=2 , flag = "0") )

  amcrsec =  formatC(as.integer(mymcrsec), width=6, flag = "0")

 

  charvec = paste(gdates, gtimes, amcrsec, A)

  cat(charvec, sep="\n")

  invisible(charvec)

}

