getCMT<-function(fn, skip=1)
{
  if(missing(skip)) { skip=1 }

  #####  scan the CMT data from the harvard site:
  ####   values are separated by a space
  CMT = scan(file=fn, skip=skip, what=list(lon=0, lat=0, str1=0, dip1=0, rake1=0, str2=0, dip2=0, rake2=0, sc=0, iexp=0, name=""))

######  mechs on this site have two forms for the name
  LENDATE = nchar(CMT$name)
  YR = rep(0, length(CMT$name))
  DY = rep(0, length(CMT$name))
  MO = rep(0, length(CMT$name))
  HR = rep(0, length(CMT$name))
  MI = rep(0, length(CMT$name))
  SE = rep(0, length(CMT$name))
  
  YR[LENDATE==7] = as.numeric(substr(CMT$name[LENDATE==7] , 5,6))
  DY[LENDATE==7] = as.numeric(substr(CMT$name[LENDATE==7] , 3,4))
  MO[LENDATE==7] = as.numeric(substr(CMT$name[LENDATE==7] , 1,2))

  YR[LENDATE>7] = as.numeric(substr(CMT$name[LENDATE>7] , 1,4))
  DY[LENDATE>7] = as.numeric(substr(CMT$name[LENDATE>7] , 7,8))
  MO[LENDATE>7] = as.numeric(substr(CMT$name[LENDATE>7] , 5, 6))

  HR[LENDATE>7] = as.numeric(substr(CMT$name[LENDATE>7] , 9,10))
  MI[LENDATE>7] = as.numeric(substr(CMT$name[LENDATE>7] , 11,12))
 ##  SE[LENDATE>7] = as.numeric(substr(CMT$name[LENDATE>7] , 13, 14))

  YR[YR<70] = YR[YR<70] + 2000
  YR[YR<100] = YR[YR<100] + 1900

  CMT$yr = YR
  CMT$mo = MO
  CMT$dom = DY

  CMT$jd = RSEIS::getjul(CMT$yr,CMT$mo,CMT$dom )

  CMT$hr = HR
  CMT$mi = MI
  CMT$se = SE

  invisible(CMT)

}
