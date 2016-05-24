PCfiledatetime<-function (orgtim, tims) 
{
 ###   orgtime is a vector c(yr, jd, hr, mi, sec)
  if(is.list(orgtim))
    {
      rd = list(yr=orgtim$yr[1],  jd=orgtim$jd[1], hr=orgtim$hr[1], mi=orgtim$mi[1], sec=orgtim$sec[1]+tims)
    }
  else
    {
  rd = list(yr=orgtim[1],  jd=orgtim[2], hr=orgtim[3], mi=orgtim[4], sec=orgtim[5]+tims)
}
  
  rd = RSEIS::recdate(jd = rd$jd, hr = rd$hr, mi = rd$mi, sec = rd$sec, 
    yr = rd$yr)
  MODAY = RSEIS::getmoday(rd$jd, rd$yr)
  gdates = paste(sep = "_", formatC(rd$yr, width = 3, flag = "0"), 
    formatC(MODAY$mo, width = 2, flag = "0"), formatC(MODAY$dom, 
                                 width = 2, flag = "0"))
  mysec = floor(rd$sec)
  mymcrsec = floor((rd$sec - mysec) * 1e+06)
  hours = rd$hr
  gtimes = paste(sep = "_", formatC(hours, width = 2, flag = "0"), 
    formatC(rd$mi, width = 2, flag = "0"), formatC(mysec, width = 2, 
                              flag = "0"))
  amcrsec = formatC(as.integer(mymcrsec), width = 6, flag = "0")
  charvec = paste(sep="_", gdates, gtimes, amcrsec)
  cat(charvec, sep = "\n")
  
  return( charvec[which.min(tims)] )
}




