`JSEGY.seis` <-function(fnames, Iendian=1 , HEADONLY=FALSE, BIGLONG=FALSE,  PLOT=-1, RAW=FALSE )
{
  ###  get a bunch of segy files from a directory and store in structure
  ####  
  
  if(missing(PLOT)) { PLOT=-1 }
  
  if(missing(Iendian)) { Iendian=1 }
  if(missing(HEADONLY)) {HEADONLY=FALSE }
  if(missing(BIGLONG)) { BIGLONG=FALSE}
  if(missing(RAW)) { RAW=FALSE }

  GIVE = as.list(1:length(fnames))


 DATAendian =c("little", "big", "swap")
  Kendian =c(1,2,3)

  #### Iendian should now be a vector


 

if(is.character(Iendian))
  {
    endianVEC =Iendian
  }
  else
    {
      endianVEC = DATAendian[match(Iendian , Kendian )]
    }

  if(length(endianVEC)<length(fnames)) { endianVEC = rep(endianVEC, times=length(fnames) ) }


    
  
 #### if(is.character(Iendian))
 ####   {
 ####     Iendian = grep(Iendian, DATAendian)
####      theENDIAN = DATAendian[Iendian]
####    }
 #### else
####    {
####      theENDIAN = DATAendian[Iendian]
      
####    }
  

 ########  on some systems LONG = 8 bytes
   #######   on others LONG = 4 bytes
    #####   there may be other differences that
     ###   I am not aware of but this is enough for now.
  
 if(BIGLONG)
    {
      
    ishort = 2
  iint  = 4
  ilong = 8
  ifloat = 4
  idouble = 8

  }
  else
    {

    ishort = 2
  iint  = 4
  ilong = 4
  ifloat = 4
  idouble = 8



    }

  for(i in 1:length(fnames))
    {

      fn = fnames[i]
      infile = fn
      theENDIAN = endianVEC[i]
####
     ##  print(paste(fn, theENDIAN) );

      
###  if this file does not exist, exit!
      if(file.exists(infile)==FALSE)
        {
          print(paste(sep=' ', "file does not exist", fn) ); 
          next;
        }
      else
        {
###  print(paste(sep=' ', "file exists", fn) );

        }

      SEGYhead.names = c("lineSeq", "reelSeq", "event_number", "channel_number",
        "energySourcePt", "cdpEns", "traceInEnsemble",
        "traceID", "vertSum", "horSum", "dataUse",
        "sourceToRecDist", "recElevation", "sourceSurfaceElevation", "sourceDepth",
        "datumElevRec", "datumElevSource", "sourceWaterDepth","recWaterDepth",
        
        "elevationScale", "coordScale",

        "sourceLongOrX", "sourceLatOrY","recLongOrX", "recLatOrY",
        
        "coordUnits", "weatheringVelocity", "subWeatheringVelocity",
        "sourceUpholeTime", "recUpholeTime",  "sourceStaticCor",
        "recStaticCor", "totalStatic", "lagTimeA", "lagTimeB", "delay", 
        "muteStart", "muteEnd",
        
        "sampleLength",
        "deltaSample",
        "gainType",
        "gainConst", 
        "initialGain",
        
        "correlated", "sweepStart", "sweepEnd", "sweepLength", "sweepType", 
        "sweepTaperAtStart", "sweepTaperAtEnd", "taperType", "aliasFreq", "aliasSlope",
        "notchFreq","notchSlope", "lowCutFreq", "hiCutFreq", "lowCutSlope", "hiCutSlope",
        
        "year", 
        "day",
        "hour",
        "minute",
        "second",
        "timeBasisCode",
        
        "traceWeightingFactor", "phoneRollPos1", "phoneFirstTrace",
        "phoneLastTrace", "gapSize", "taperOvertravel",
        
        "station_name", 
        "sensor_serial",
        "channel_name",
        "totalStaticHi",
        "samp_rate",
        "data_form",
        "m_secs", 
        "trigyear", "trigday", "trighour", "trigminute", "trigsecond", "trigmills", 
        "scale_fac",
        "inst_no",
        "not_to_be_used",
        "num_samps",
        "max", "min")



      
      zz <- file(fn, "rb")

      A1 =  readBin(zz, integer() , n = 7, size = iint, signed = TRUE,
        endian = theENDIAN)
      
      A2 =  readBin(zz, integer() , n = 4, size =ishort , signed = TRUE,
        endian = theENDIAN)
      
      A3 =  readBin(zz, integer() , n = 8, size = iint, signed = TRUE,
        endian = theENDIAN)

     
 ###  short elevationScale;    /*  68 Elevation Scaler: scale = 1 */
 ###   short coordScale;        /*  70 Coordinate Scaler: scale = 1 */
      A4 = readBin(zz, integer() , n = 2 , size =ishort,  endian =theENDIAN, signed = TRUE)
      

      A5 = readBin(zz, integer() , n = 4 , size =iint,  endian = theENDIAN, signed = TRUE)

      coordScale= A4[2]
      
      scfc = 1
      scfe = 1
      
      if (coordScale < 0)
        {
          scfc = -1. / coordScale;
        }
      recLongOrX= scfc* A5[3]
      recLatOrY=scfc*A5[4]
      
      elevationScale=A4[1]
      scfe = elevationScale
      if (elevationScale < 0){ scfe = -1. / elevationScale }
      
      recElevation=scfe*A3[2]

      sourceLongOrX=A5[1]/10000
      sourceLatOrY=A5[2]/10000
      
      
      sourceDepth=scfe*A3[4]
      
      
      coords = list(reclat=recLatOrY, reclon=recLongOrX, recel =recElevation, srclat= sourceLatOrY,
        srclon=sourceLongOrX, srcdep=sourceDepth)
     
      A6 = readBin(zz, integer() , n =13  , size =ishort,  endian = theENDIAN, signed = TRUE)
      coordUnits = A6[1]

      
      sampleLength=readBin(zz, integer() , n = 1 , size =ishort,  endian =theENDIAN, signed = TRUE)
      deltaSample=readBin(zz, integer() , n = 1 , size =ishort,  endian =theENDIAN, signed = TRUE)
      gainType=readBin(zz, integer() , n = 1 , size =ishort,  endian =theENDIAN, signed = TRUE);        
      gainConst=readBin(zz, integer() , n = 1 , size =ishort,  endian =theENDIAN, signed = TRUE);    


      
      initialGain =readBin(zz, integer() , n = 1 , size =ishort,  endian =theENDIAN, signed = TRUE);

      
      

      B6 = readBin(zz, integer() , n =16  , size =ishort,  endian = theENDIAN, signed = TRUE)

   
      year=readBin(zz, integer() , n = 1 , size =ishort,  endian =theENDIAN, signed = TRUE)
      day=readBin(zz, integer() , n = 1 , size =ishort,  endian =theENDIAN, signed = TRUE);
      hour=readBin(zz, integer() , n = 1 , size =ishort,  endian =theENDIAN, signed = TRUE); 
      minute=readBin(zz, integer() , n = 1 , size =ishort,  endian =theENDIAN, signed = TRUE);
      second=readBin(zz, integer() , n = 1 , size =ishort,  endian =theENDIAN, signed = TRUE);

      
      timeBasisCode=readBin(zz, integer() , n = 1 , size =ishort,  endian =theENDIAN, signed = TRUE);

      B7 = readBin(zz, integer() , n =6  , size =ishort,  endian = theENDIAN, signed = TRUE)


 
      stationname=readChar(zz, 6, useBytes = FALSE)
      sensorserial=readChar(zz, 8, useBytes = FALSE)
      channelname=readChar(zz, 4, useBytes = FALSE)
      
      ##  print(paste(sep="", "<",channelname,">"))

##  close(zz)
   ##       cat(paste(sep=" ",paste(sep="", "<",stationname,">"), paste(sep="", "<",channelname,">")), sep="\n") 
      
 ## print(paste(sep=" ", "Name=", stationname,sensorserial,channelname))  

       
      totalStaticHi= readBin(zz, integer() , n = 1 , size =ishort,  endian = theENDIAN, signed = TRUE)


      samprate= readBin(zz, integer() , n = 1 , size =iint,  endian = theENDIAN, signed = TRUE)


      dataform=readBin(zz, integer() , n = 1 , size =ishort ,  endian = theENDIAN, signed = TRUE)

      msecs =readBin(zz, integer() , n = 1 , size =ishort ,  endian = theENDIAN, signed = TRUE)

     ##   print(paste(sep=" ", "time start=", year,day,hour,minute, second,  msecs  ))
      

   ##  print(paste(sep=" ", "gains", gainType, gainConst ))


    
      
      trigyear=readBin(zz, integer() , n = 1 , size =ishort ,  endian = theENDIAN, signed = TRUE)
      trigday=readBin(zz, integer() , n = 1 , size =ishort ,  endian = theENDIAN, signed = TRUE)
      trighour=readBin(zz, integer() , n = 1 , size =ishort ,  endian = theENDIAN, signed = TRUE)
      trigminute=readBin(zz, integer() , n = 1 , size =ishort ,  endian = theENDIAN, signed = TRUE)
      trigsecond=readBin(zz, integer() , n = 1 , size =ishort ,  endian = theENDIAN, signed = TRUE)
      trigmills=readBin(zz, integer() , n = 1 , size =ishort ,  endian = theENDIAN, signed = TRUE)

  ##      print(paste(sep=" ","trigtime=", trigyear,  trigday,trighour , trigminute, trigsecond, trigmills  ))
      
      scalefac=readBin(zz, numeric() , n = 1 , size =ifloat  ,  endian = theENDIAN, signed = TRUE)

      if(gainConst == 0)gainConst = 1
      if(scalefac == 0) scalefac  = 1

      

      scalefac = scalefac / gainConst;

      
      
      instno=readBin(zz, integer() , n = 1 , size =ishort ,  endian = theENDIAN, signed = TRUE)

      nottobeused =readBin(zz, integer() , n = 1 , size =ishort ,  endian = theENDIAN, signed = TRUE)
      numsamps =readBin(zz, integer() , n = 1 , size =iint ,  endian = theENDIAN, signed = TRUE)

      max =readBin(zz, integer() , n = 1 , size =iint ,  endian = theENDIAN, signed = TRUE)

      min =readBin(zz, integer() , n = 1 , size =iint ,  endian = theENDIAN, signed = TRUE)


  ##  print(paste(sep=" ", i, "numsamps=", numsamps, "sampleLength=", sampleLength  ))
  ##  print(paste(sep=" ",  "scalefac=",   scalefac ))
      ####################################      done reading header
       #####    count = 7+4+8+2+4+13+5+16+6+6+4+7+5


SEGYall = c(A1, A2, A3, A4, A5, A6, sampleLength,
  deltaSample,gainType,gainConst,initialGain,B6,
  year, day, hour, minute, second,timeBasisCode,
  B7,stationname,sensorserial,channelname,
  totalStaticHi, samprate, dataform,
  msecs ,trigyear, trigday,trighour,
  trigminute, trigsecond, trigmills,scalefac,
  instno, nottobeused, numsamps, max, min)

      
  
  SEGYH = data.frame(names=SEGYhead.names, values=SEGYall)
  
  


      
      N = numsamps

      if( !is.integer(N))
         {
           ##  this is a problem
           print("ERROR: number of samples is not an integer.")
           print(paste(i, fn))

         }

         
       ####  bugfix from Jake
      ####   dt = as.numeric( deltaSample )/1000000
       
      if(as.numeric( samprate )==0)
        {
          samprate = deltaSample
        }

      
      dt = as.numeric( samprate )/1000000
      ####   you have to have a sample rate
      if(dt==0)
        {
          dt=0.025
        }
      
      DATIM =   c(year, day, hour, minute)

      
      sec = second+  msecs/1000
      
      thesta= stationname
      thecomp= channelname

      
       md = getmoday(DATIM[2], DATIM[1])

      t1 = 0
      t2 = dt*(N-1)
      

      tstart = list(yr=DATIM[1], jd=DATIM[2] , mo=md$mo, dom=md$dom, hr=DATIM[3], mi=DATIM[4], sec=sec, msec=0, dt=dt, t1=t1,
        t2=t2, off=0)

      aunits="unknown"
    ####  if(is.null(thesta))   thesta="XXX"
     #### if(is.null(thecomp))  thecomp="X"
     ####  if(is.na(thesta))   thesta="XXX"
     #### if(is.na(thecomp))  thecomp="X"
           
      

      if(HEADONLY==TRUE)
          {
              
####   print(paste("headonly ", i))
              
              x = NULL
              aunits=NA
              
              
              
          }
      else
          {
####################################   read in integer  samples
####################################   
####################################            
              D1 = readBin(zz, integer() , n = N , size =iint ,  endian = theENDIAN, signed = TRUE)
              
####################################
####################################      
              
###    i1 = unlist( strsplit(split=" ", B4[4]) )


              if(RAW)
                  {
                      
####################################
####################################      
                      
                      x = as.vector(D1)

                      aunits="counts"
###  print("using RAW (counts), range:")
###  print(range(x))
                      
                  }
              else
                  {

####################################
####################################      
                      
#########  convert to floating point numbers with physical units
                      x =   as.vector(D1)/scalefac
####################################
####################################      


                      aunits="volts"
                  }

          }




      
      close(zz)


      
      GIVE[[i]] = list(fn=fn, sta=thesta,  comp=thecomp, dt=dt, DATTIM=tstart, N=N, units=aunits , coords=coords ,  amp=x, HEAD=SEGYH, IO=list(kind=1, Iendian=Iendian,  BIGLONG=BIGLONG) )

      
      if(PLOT>=0)
        {
          plot(x, type='l', main=paste(sep=' ', thesta, thecomp))
          print("left CLICK once in WINDOW for NEXT TRACE:")
           if(PLOT==0) { locator(1) }
          else
            {
              Sys.sleep(PLOT);
            }
          
        }
    }


 ### print("using RAW (counts), range:")
 ### print(range(x))
 ### for(ig in 1:length(GIVE)) { print( range(GIVE[[ig]]$amp ))  }
  
  
  invisible(GIVE)
}

