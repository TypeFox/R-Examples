FmakeDB<-function(LF2, kind =1, Iendian=1, BIGLONG=FALSE)
  {
    if(missing(kind)) { kind =1 }
    
    if(missing(Iendian)) { Iendian=1 }
    if(missing(BIGLONG)) { BIGLONG=FALSE }

    ADB = list(fn="", 
      yr=0,
      jd=0,
      hr=0,
      mi=0,
      sec=0,
      dur=0,
      t1=0,
      t2=0,
      sta="", 
      comp="")
    attr(ADB, "origyr")<- 1972
    N = 0

    if(length(kind)==1) kind = rep(kind, times=length(LF2) )
    if(length(Iendian)==1) Iendian = rep(Iendian, times=length(LF2) )
    if(length(BIGLONG)==1) BIGLONG = rep(BIGLONG, times=length(LF2) )
     
    ##   OLD: sinfo  =  getseisinfo(LF2, kind=kind)
    for(i in 1:length(LF2))
      {
        sinfo  = GET.seis(LF2[i], kind=kind[i], Iendian=Iendian[i], BIGLONG=BIGLONG[i] , HEADONLY=TRUE , PLOT=-1)
        
        for(j in 1:length(sinfo))
          {
            REC = sinfo[[j]]
            N = N + 1
            ADB$fn[N] = REC$fn
            ADB$sta[N] = REC$sta
            ADB$comp[N] = REC$comp
            ADB$yr[N] = REC$DATTIM$yr
            ADB$jd[N] = REC$DATTIM$jd
            ADB$hr[N] = REC$DATTIM$hr
            ADB$mi[N] = REC$DATTIM$mi
            ADB$sec[N] = REC$DATTIM$sec+REC$DATTIM$msec/1000
            ADB$dur[N] = REC$DATTIM$dt*REC$N
            
            
          }
        
 
      }
    
       
         origyr = min(ADB$yr)
        eday = EPOCHday(ADB$yr, jd = ADB$jd, origyr = origyr)
        ADB$t1 = eday$jday + ADB$hr/24 + ADB$mi/(24 * 60) + ADB$sec/(24 *
          3600)
        ADB$t2 = ADB$t1 + ADB$dur/(24 * 3600)

    attr(ADB, "origyr")<- origyr
    attr(ADB, "kind")=kind
    attr(ADB, "Iendian")=Iendian
    attr(ADB, "BIGLONG")=BIGLONG

    
    
        invisible(ADB)
        


    
  }
