`setWPX` <-
  function(phase=NULL, col=NULL, yr=NULL, jd=NULL, hr=NULL,  mi=NULL, sec=NULL, dur=NULL, name=NULL, comp=NULL, dispcomp=NULL, onoff=NULL)
{

###  s1 = setWPX(phase="Y", yr=2005, jd =5, hr=6, sec=runif(10) )
###  s2 = setWPX(phase="Y", yr=2005, jd =5, hr=6, sec=runif(4) )


###   rpic = recdate(jd=5, hr=6, mi=0, sec=runif(10), yr=2005 )

  N = max(c(length(jd), length(hr), length(mi), length(sec), length(yr))  , na.rm=TRUE)
  
  
  if(missing(phase)) { phase ="P" }
  if(missing(dur)) { dur = 0 }
  if(missing(name)) { name = NA }
  if(missing(comp)) { comp ="V"  }
  if(missing(dispcomp)) { dispcomp = "V"  }
  if(missing(onoff)) { onoff = 1 }
  if(missing(yr)) { yr = 1972 }
  if(missing(jd)) { jd = 1}
  if(missing(hr)) { hr = 1}
  if(missing(mi)) { mi = 1}
  if(missing(sec)) { sec = 1}
  if(missing(col)) { col = 1 }

  
  if(length(phase)==1) { phase =  rep(phase, times=N) }
  if(length(yr)==1) { yr =  rep(yr, times=N) }
  if(length(col)==1) { col  =  rep(col, times=N) }
  if(length(dur)==1) { dur  =  rep(dur, times=N) }
  if(length(jd)==1) { jd  =  rep(jd, times=N) }
  
  if(length(name)==1) { name  =  rep(name, times=N) }
  if(length(comp)==1) { comp  =  rep(comp, times=N) }
  if(length(onoff)==1) { onoff  =  rep(onoff, times=N) }

 ## print(jd)
 ##  print(yr)

  
  rpic = recdate(jd=jd, hr=hr, mi=mi, sec=sec, yr=yr)
  
  gd = getmoday(jd, yr )

  pix= cleanWPX()
  
  pix$yr=rpic$yr
  pix$jd= rpic$jd
  pix$mo= gd$mo
  pix$dom= gd$dom
  pix$hr= rpic$hr
  pix$mi= rpic$mi
  pix$sec= rpic$sec
  pix$dur= dur
  pix$res= dur
  pix$err= dur
  
  pix$phase = phase

  
  pix$name= name

  pix$comp= comp

  pix$flg=rep(0, times=N)
  pix$onoff=rep(1, times=N)
  pix$pol=rep("_", times=N)
  pix$c3=rep(0, times=N)

  pix$tag= name

  pix$col  = col
  
  ##   pix = SORT.pix(pix)

  invisible(pix)

}

