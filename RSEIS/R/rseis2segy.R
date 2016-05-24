rseis2segy<-function(GH, sel=1, win=c(0,1), path=".", BIGLONG=FALSE )
    {
        MAX32bit = 2147483647
##############  convert an RSEIS structure to SEGY format
########   if dir is provided, the individual files will be
#########    written in that directory
###  rseis2segy(GH, sel = which(GH$STNS=="FAR" & GH$COMPS %in% c("V", "N", "E")),
    ###################    path = ".", BIGLONG=FALSE)
    SEGYhead.names = c("lineSeq", "reelSeq", "event_number", "channel_number",
      "energySourcePt", "cdpEns", "traceInEnsemble", "traceID",
      "vertSum", "horSum", "dataUse", "sourceToRecDist",
      "recElevation", "sourceSurfaceElevation", "sourceDepth",
      "datumElevRec", "datumElevSource",
      "sourceWaterDepth","recWaterDepth",
      "elevationScale", "coordScale",
      "sourceLongOrX", "sourceLatOrY","recLongOrX", "recLatOrY",
      "coordUnits", "weatheringVelocity", "subWeatheringVelocity",
      "sourceUpholeTime", "recUpholeTime", "sourceStaticCor",
      "recStaticCor", "totalStatic", "lagTimeA", "lagTimeB",
      "delay", "muteStart", "muteEnd",
      "sampleLength", "deltaSample", "gainType", "gainConst",
      "initialGain",
      "correlated", "sweepStart", "sweepEnd", "sweepLength",
      "sweepType", "sweepTaperAtStart", "sweepTaperAtEnd",
      "taperType", "aliasFreq", "aliasSlope",
      "notchFreq","notchSlope", "lowCutFreq", "hiCutFreq",
      "lowCutSlope", "hiCutSlope",
      
      "year", "day", "hour", "minute", "second", "timeBasisCode",
      
      "traceWeightingFactor", "phoneRollPos1", "phoneFirstTrace",
      "phoneLastTrace", "gapSize", "taperOvertravel",
      
      "station_name", "sensor_serial", "channel_name",
      "totalStaticHi", "samp_rate", "data_form", "m_secs",
      "trigyear", "trigday", "trighour", "trigminute", "trigsecond",
      "trigmills", "scale_fac", "inst_no", "not_to_be_used",
      "num_samps", "max", "min")

    initsegy<-function()
      {

        N = length(SEGYhead.names)
        j = rep(0, N)
        HEAD=as.list(j)
        names(HEAD) = SEGYhead.names
        HEAD[72] = "STAXX"
        HEAD[73] ="12345"
        HEAD[74] = "XXX"

        h = list(HEAD=HEAD, amp=0)

        return(h)
      }

        ####   this is a bogus fix, but it is the convention I have used in past
      
    
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


    ############# for back compatability, need something in here if its missing
    if(is.null(  GH$info$scalefac )) GH$info$scalefac = rep(NA, times = length(GH$info$sec) )
        if(is.null(  GH$info$gain ))     GH$info$gain=   rep(1, times = length(GH$info$sec) )

        UComps = unique(GH$COMPS)
    
    for(j in 1:length(sel))
      {
        i = sel[j]
        
        fn = "GH"
        thesta = GH$STNS[i]
        thecomp = GH$COMPS[i]
        dt = GH$dt[i]

         ####  this is artificial
        CHN.NUM = which(thecomp==UComps)


        
        ## N = GH$info$n[i]

        sig = GH$JSTR[[i]]
        N = length(sig)
       
        ##  a1 = list(fn=fn, sta=thesta,  comp=thecomp, dt=dt, DATTIM=tstart,
        ##     N=N, units=aunits , amp=sig , IO=list(kind=2, Iendian=theENDIAN,  BIGLONG=BIGLONG))


####  here we need to take care of the NA's that are either
###  before or after the main signals....

        STARTDATE =  recdate(GH$info$jd[i] , GH$info$hr[i], GH$info$mi[i], GH$info$sec[i]+GH$info$msec[i]/1000,    GH$info$yr[i])

        wawa =  is.na(sig)
        if(any(wawa))
            {
### check if signal is contiguous:
                
                w1 = which(wawa)
                w2 = which(!wawa)
                w1 = which(wawa)
            if(any(diff(w1)>1) ) { print('discontiguous NA sequence') }
                w2 = which(!wawa)
                
                rw1 =  range(w1)
                rw2 =  range(w2)
                
                
####### repair blank spots in signal
                if(  any(w1>rw2[1] & w1<rw2[2]) )
                    { #### need to do some repair work
                        ##  for now replace the NA with a legit value : mean value?
                        sig[w1>rw2[1] & w1<rw2[2]] = mean(sig[w2])
                    }

                sig = sig[w2]
                N = length(sig)
### repair the start of the trace: change the start time
                
                begchunk =  w1<rw2[1]
                if(all(begchunk) )
                    {
                        nsampbeg = length(w1)
                        mw1 = max(w1)
### these should be the same
                        newSTARTsec= GH$info$sec[i]+(mw1)*GH$dt[i]
                        STARTDATE =  recdate(GH$info$jd[i] , GH$info$hr[i], GH$info$mi[i], newSTARTsec,GH$info$yr[i])
                    }
         
                
            }

        
        
       
        
        r = range(sig, na.rm=TRUE)
        scale2 = max(abs(r))
        scale3 = MAX32bit/scale2
        scalefac=scale3
  
         if(any(is.na(sig) ) )
            {
                cat('ERROR!', sep='\n')
                 cat(paste(c(i,scalefac, range(sig, na.rm=TRUE) ), collapse= ' ') , sep='\n')

                break
            }

        a1 =  initsegy()

###  cat(paste("a1$HEAD$",SEGYhead.names,"=",  sep="")  , sep="\n")

        sec = floor(STARTDATE$sec)
        msec = (STARTDATE$sec - sec)*1000

        a1$HEAD$lineSeq=i
        a1$HEAD$year=STARTDATE$yr
        a1$HEAD$day=STARTDATE$jd
        a1$HEAD$hour=STARTDATE$hr
        a1$HEAD$minute=STARTDATE$mi
        a1$HEAD$second=STARTDATE$sec
        a1$HEAD$gainConst= 1
        a1$HEAD$gainType= 1
        
        a1$HEAD$station_name=thesta
        a1$HEAD$sensor_serial="12345"
        a1$HEAD$channel_name=thecomp
        a1$HEAD$channel_number=CHN.NUM
        a1$HEAD$inst_no=999

        a1$HEAD$timeBasisCode = 2

        a1$HEAD$deltaSample = dt*1000000
        
        a1$HEAD$samp_rate=dt*1000000
        a1$HEAD$m_secs=msec


#####  here we need to determine the scale factor
####   turn the doubles (floating point)  to ints 
####        if(is.null(GH$info$scalefac[i]) | is.na(GH$info$scalefac[i]) )
            
                
                ################   GH$info$scalefac[i]=1

        if(is.null(GH$info$gain[i] )|is.na(GH$info$gain[i] )) GH$info$gain[i]=1
        
        a1$HEAD$scale_fac=scalefac
        a1$HEAD$num_samps=N
        
        a1$HEAD$gainConst = GH$info$gain[i]

        

        scalefac =  scalefac / GH$info$gain[i]

        Kounts = round(sig*scalefac)
        
        a1$amp = as.integer(Kounts)

        a1$HEAD$min = min(a1$amp, na.rm=TRUE)
        a1$HEAD$max = max(a1$amp, na.rm=TRUE)
        a1$HEAD$data_form =1

         a1$HEAD$lineSeq= 1 
         a1$HEAD$event_number= 1
         a1$HEAD$reelSeq= 1
         a1$HEAD$traceID = 1
        
        segyfn = paste(sep=".",
          formatC(STARTDATE$yr , width = 4, flag = "0") ,
          formatC(STARTDATE$jd, width = 3, flag = "0"),
          formatC(STARTDATE$hr, width = 2, flag = "0") ,
          formatC(STARTDATE$mi, width = 2, flag = "0") ,
          formatC(STARTDATE$sec, width = 2, flag = "0")  ,
          thesta ,
          thecomp,
          "SEGY")


        cat(segyfn, sep='\n')
        ########cat(paste(c(i, scalefac, range(sig) , range(a1$amp) ), collapse= ' ') , sep='\n')

        ######################  write out the SEGY file
        write1segy(a1, fn=segyfn , BIGLONG=BIGLONG  )
        ###  cat('done write1segy1' , sep='\n')
      }



    
    setwd(  OLDir)
    
  }
