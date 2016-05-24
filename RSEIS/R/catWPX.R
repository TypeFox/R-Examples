catWPX<-function(WPX, ppx)
  {
    ####  combine two wpx pickfiles dataframes

 ###    KPX  =vector(mode="list") 

 KPX  = setypx()

    
    KPX$tag= c(WPX$tag,  ppx$tag )
    KPX$name=c(WPX$name,ppx$name )
    KPX$comp=c(WPX$comp,ppx$comp )
    KPX$c3=c(WPX$c3,ppx$c3 )
    KPX$phase=c(WPX$phase,ppx$phase )

    KPX$err=c(WPX$err,ppx$err )
    KPX$pol=c(WPX$pol,ppx$pol )
    KPX$flg=c(WPX$flg,ppx$flg )
    KPX$res= c(WPX$res,ppx$res )
    KPX$dur= c(WPX$dur,ppx$dur )
  
    KPX$yr=c(WPX$yr,ppx$yr )
    KPX$mo= c(WPX$mo,ppx$mo )
    KPX$dom=c(WPX$dom,ppx$dom )
    KPX$jd=c(WPX$jd,ppx$jd )
    KPX$hr=c(WPX$hr,ppx$hr )
    KPX$mi=c(WPX$mi,ppx$mi )
    KPX$sec=c(WPX$sec,ppx$sec )
    KPX$col=c(WPX$col,ppx$col )
    KPX$onoff = c(WPX$onoff,ppx$onoff )


    return(KPX)
    
  }
