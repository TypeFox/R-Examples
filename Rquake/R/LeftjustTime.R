LeftjustTime<-function(g1)
  {
#################   rectify all date/times to same minute mark
    ###########  used for earthquake location
    
    d1 = RSEIS::recdate(jd=g1$jd, hr=g1$hr, mi=g1$mi, sec=g1$sec, yr=g1$yr)
    Adate = Qrangedatetime(d1)
    mdate = Adate$min
   YY = RSEIS::YRsecdif(mdate$jd, mdate$hr, mdate$mi, 0, d1$jd, d1$hr, d1$mi, d1$sec, yr1 = mdate$yr, yr2 = d1$yr)
   dm =  RSEIS::getmoday(mdate$jd, mdate$yr)
    ret1 = g1
    n = length(YY)
    ret1$sec = YY
    ret1$yr = rep(mdate$yr, n)
    ret1$jd = rep(mdate$jd, n)
    ret1$hr = rep(mdate$hr, n)
    ret1$mi  = rep(mdate$mi, n)
    ret1$dom  = rep(dm$dom, n)
    ret1$mo  = rep(dm$mo, n)
  
    return( ret1  )
  }
