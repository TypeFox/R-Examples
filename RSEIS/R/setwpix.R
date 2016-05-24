`setwpix` <-
function(phase=NULL, col=NULL, yr=NULL, jd=NULL, hr=NULL,  mi=NULL, sec=NULL, dur=NULL, name=NULL, comp=NULL, dispcomp=NULL)
{


  
  if(is.list(sec))
    {
      yr = sec$yr
      jd = sec$jd
      hr =  sec$hr
      mi =  sec$mi
      sec =  sec$sec


    }

 
    N = max(length(jd), length(hr), length(mi), length(sec), length(yr), na.rm=TRUE)
    
  
  if(missing(phase)) { phase = rep("P", times=N) }
  if(missing(dur)) { dur = rep(0, times=N) }
  if(missing(name)) { name = rep(NA, times=N)  }
  if(missing(comp)) { comp = rep("V", times=N)  }
  if(missing(dispcomp)) { dispcomp = rep("V", times=N)  }
    
  if(length(phase)==1) { phase =  rep(phase, times=N) }
  if(length(yr)==1) { yr =  rep(yr, times=N) }
  if(length(col)==1) { col  =  rep(col, times=N) }
  if(length(dur)==1) { dur  =  rep(dur, times=N) }
  if(length(jd)==1) { jd  =  rep(jd, times=N) }
  
  if(length(name)==1) { name  =  rep(name, times=N) }
  if(length(comp)==1) { comp  =  rep(comp, times=N) }
  
  rpic = recdate(jd=jd, hr=hr, mi=mi, sec=sec, yr=yr)
  
  gd = getmoday(jd, yr )

  pix=setypx()
  
  pix$yr=yr
  pix$jd= rpic$jd
  pix$mo= gd$mo
  pix$dom= gd$dom
  pix$hr= rpic$hr
  pix$mi= rpic$mi
  pix$sec= rpic$sec
  pix$dur= dur
  pix$res= dur
    pix$phase = phase

    
  pix$name= name
  pix$comp= comp

  pix$flg=rep(0, times=N)
  pix$onoff=rep(0, times=N)
  pix$pol=rep("_", times=N)
  pix$c3=rep(0, times=N)

  pix$tag=rep(0, times=N)

    pix$col  = col
  
  ##   pix = SORT.pix(pix)
  
  invisible(pix)

}

