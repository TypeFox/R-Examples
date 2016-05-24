addWPX<-function(WPX, ppx)
  {
    ####  add one pick to WPX file


    NPX  = length(WPX$sec)+1
    
    WPX$tag[NPX]= ppx$tag
    WPX$name[NPX]=ppx$name
    WPX$comp[NPX]=ppx$comp
    WPX$c3[NPX]=ppx$c3
    WPX$phase[NPX]=ppx$phase

    WPX$err[NPX]=ppx$err
    WPX$pol[NPX]=ppx$pol
    WPX$flg[NPX]=ppx$flg
    WPX$res[NPX]= ppx$res
    WPX$dur[NPX]= ppx$dur
    
    WPX$yr[NPX]=ppx$yr
    WPX$mo[NPX]= ppx$mo
    WPX$dom[NPX]=ppx$dom
    WPX$jd[NPX]=ppx$jd
    WPX$hr[NPX]=ppx$hr
    WPX$mi[NPX]=ppx$mi
    WPX$sec[NPX]=ppx$sec
    
    WPX$col[NPX]=ppx$col
    WPX$onoff[NPX] = ppx$onoff


    return(WPX)
  }
