cleanWPX<-function()
  {
    ####  return an empty WPX
    WPX = list(
      tag=NA,
      name=NA,
      comp=NA,
      c3=NA,
      phase=NA,
      err=0,
      pol=0,
      flg=0,
      res=0,
      dur=0,
      yr=0,
      mo=0,
      dom=0,
      jd=0,
      hr=0,
      mi=0,
      sec=0,
      col='black',
      onoff =0  )
  ##   WPX = data.frame(WPX, stringsAsFactors = FALSE)
    
    return(WPX)
  }
