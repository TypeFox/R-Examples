INITpickfile <-
function(stas=NULL, src=NULL, WPX=NULL)
  {
##  initialize a pickfile based on picks
    ###  from RSEIS::swig stored in the  WPX

    if(is.null(WPX)) return(NULL)
    
    Apf = RSEIS::EmptyPickfile()


    nona = which( is.na(WPX$tag) )
        
        if(length(nona)>0)
          {
            WPX = RSEIS::deleteWPX(WPX, nona)
          }
        
    A1T = Qrangedatetime(WPX)
    s1 = RSEIS::secdifL(A1T$min,  WPX)

    ords1 = order(s1)

    osta =  WPX$name[ords1]
    lat1 = stas$lat[osta[1]==stas$name]
    lon1 = stas$lon[osta[1]==stas$name]

    spx = data.frame(WPX, stringsAsFactors = FALSE)
    spx = spx[ords1, ]
    
    wfirst = as.list( spx[1, ] )
    
    wfirst$sec = 0
    psecs = RSEIS::secdifL(wfirst, spx)
    
    Apf$LOC$yr =  wfirst$yr
    Apf$LOC$jd =  wfirst$jd
    Apf$LOC$hr =  wfirst$hr
    Apf$LOC$mi =  wfirst$mi
    moday =  RSEIS::getmoday(wfirst$yr, wfirst$jd)
    Apf$LOC$mo =  wfirst$mo
    Apf$LOC$dom =  wfirst$dom
    
    Apf$LOC$lat =  lat1
    Apf$LOC$lon =  lon1
    Apf$LOC$z =  6


    spx =  as.list(spx) 

    Apf$STAS$tag=spx$tag[spx$onoff>=0]
    Apf$STAS$name=spx$name[spx$onoff>=0]
    Apf$STAS$comp=spx$comp[spx$onoff>=0]
    Apf$STAS$c3=spx$c3[spx$onoff>=0]
    Apf$STAS$phase=spx$phase[spx$onoff>=0]
    Apf$STAS$sec=psecs
    Apf$STAS$err=spx$err[spx$onoff>=0]
    Apf$STAS$pol=spx$pol[spx$onoff>=0]
    Apf$STAS$flg=spx$flg[spx$onoff>=0]
    Apf$STAS$res=spx$res[spx$onoff>=0]

    msta = match(spx$name[spx$onoff>=0] , stas$name)

    

    Apf$STAS$lat=stas$lat[msta]
    Apf$STAS$lon= stas$lon[msta]
    
    Apf$STAS$z= stas$z[msta]
    

    return(Apf)
    
  }
