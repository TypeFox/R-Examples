makeDB<-function(path="." , pattern="R", dirs="", kind =1, Iendian=1, BIGLONG=FALSE)
  {
    if(missing(kind)) { kind =1 }
    if(missing(pattern)) {pattern="R"  }
      if(missing(Iendian)) { Iendian=1 }
       if(missing(BIGLONG)) { BIGLONG=FALSE }
        if(missing(path)) { path="." }
     
if(dirs=="")
  {
  
    if(is.na(pattern) | is.null(pattern)  )
      {
        LF1 = list.files(path=path)
        NDIRS = 1
      }
    else
      {

        LF1 = list.files(path=path, pattern=pattern)
        NDIRS = length(LF1)

      }
    ####
   #### if(dirs=="") { dirs="." }

  }
    else
      {
    if(length(dirs)>=1)
      {
           LF1 = dirs
           NDIRS = length(LF1)
      }
  }
    
    
    
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
    for(i in 1:NDIRS)
      {
        LF2 = list.files(path =paste(sep="/", path, LF1[i]) , full.names=TRUE, recursive = TRUE)



        
        ##   OLD: sinfo  =  getseisinfo(LF2, kind=kind)

        sinfo  = JGET.seis(LF2, kind=kind, Iendian=Iendian, BIGLONG=BIGLONG , HEADONLY=TRUE , PLOT=-1)
        
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
           ADB$dt[N] = REC$DATTIM$dt
           ADB$dur[N] = REC$DATTIM$dt*REC$N
          }

      }



    
    origyr = min(ADB$yr)

 eday = EPOCHday(ADB$yr, jd = ADB$jd, origyr = origyr)
    ADB$t1 = eday$jday + ADB$hr/24 + ADB$mi/(24 * 60) + ADB$sec/(24 *
        3600)
    ADB$t2 = ADB$t1 + ADB$dur/(24 * 3600)
    attr(ADB, "origyr")<- min(origyr,na.rm = TRUE )
    attr(ADB, "kind")=kind
    attr(ADB, "Iendian")=Iendian
    attr(ADB, "BIGLONG")=BIGLONG

    
    invisible(ADB)

    
  }
