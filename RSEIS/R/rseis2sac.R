rseis2sac<-function(GH, sel=1, win=c(0,1), path=".", BIGLONG=FALSE )
  {
##############  convert an RSEIS structure to SAC format
########   if path is provided, the individual files will be
#########    written in that directory
    
    if(missing(sel)) { sel = 1:length(GH$STNS) }

    
    theENDIAN =  .Platform$endian

    aunits = "volts"


    OLDir = getwd()
    RDT = rangedatetime(GH$info)
    newdir  = paste(path,filedatetime(RDT$min), sep="/")

    if(file.exists(newdir))
      {
        setwd(newdir)
      }
    else
      {
        tdir = dir.create(newdir, recursive=TRUE )
        if(tdir==TRUE) { setwd(newdir)  }
        else
          {
            print("ERROR: CANNOT create or write in this directory")
          }
      }
    


    
    if( !identical(path, ".")) 
      {
        OLDir = getwd()
        dir.create(path)
        setwd(path)
      }

    for(j in 1:length(sel))
      {
        i = sel[j]
        
        fn = "GH"
        thesta = GH$STNS[i]
        thecomp = GH$COMPS[i]
        dt = GH$dt[i]
        
        tstart = list(yr=GH$info$yr[i]  ,
          jd=GH$info$jd[i]  ,
          mo=GH$info$mo[i]  ,
          dom=GH$info$dom[i]  ,
          hr=GH$info$hr[i]  ,
          mi=GH$info$mi[i]  ,
          sec=GH$info$sec[i]  ,
          msec=GH$info$msec[i]  ,
          dt=GH$info$dt[i]  ,
          t1=GH$info$t1[i]  ,
          t2=GH$info$t2[i]  ,
          off=GH$info$off[i]  
          )
        N = GH$info$n[i]

        sig = GH$JSTR[[i]]

    sacfn = paste(sep=".",
          formatC(GH$info$yr[i] , width = 4, flag = "0") ,
          formatC(GH$info$jd[i], width = 3, flag = "0"),
          formatC(GH$info$hr[i], width = 2, flag = "0") ,
          formatC(GH$info$mi[i], width = 2, flag = "0") ,
          formatC(GH$info$sec[i], width = 2, flag = "0")  ,
          thesta ,
          thecomp,
          "SAC")

     
        
        a1 = list(fn=fn, sta=thesta,  comp=thecomp, dt=dt, DATTIM=tstart,
          N=N, units=aunits , amp=sig , IO=list(kind=2, Iendian=theENDIAN,
                                          BIGLONG=BIGLONG))

        ######################  write out the SAC file
        write1sac(a1, fn=sacfn, BIGLONG=BIGLONG )
        
      }



        setwd(  OLDir)
    
  }
