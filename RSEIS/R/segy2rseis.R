segy2rseis<-function(fnames, Iendian=1 , HEADONLY=FALSE, BIGLONG=FALSE,  PLOT=-1, RAW=FALSE )
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

  

  if(is.character(Iendian))
    {
      endianVEC =Iendian
    }
  else
    {
      endianVEC = DATAendian[match(Iendian , Kendian )]
    }

  if(length(endianVEC)<length(fnames)) { endianVEC = rep(endianVEC, times=length(fnames) ) }

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
        ##   print(paste(sep=' ', "file exists", fn) );

        }


      ###########   NOTE: here the data is counts
      ONEsegy =   read1segy(fn, Iendian=theENDIAN , HEADONLY=FALSE, BIGLONG=BIGLONG ) ;

      gsegy = ONEsegy$HEAD

      gainConst = gsegy$gainConst
      scalefac  = gsegy$scale_fac

      
      if(gainConst == 0)gainConst = 1
      if(scalefac == 0) scalefac  = 1

      

      conversionfac = scalefac / gainConst;



      N = gsegy$num_samps

      
      if( !is.integer(N))
        {
          ##  this is a problem
          print("ERROR: number of samples is not an integer.")
          print(paste(i, fn))

        }

      
####  bugfix from Jake
####   dt = as.numeric( deltaSample )/1000000
      
      samprate = gsegy$samp_rate

      if(as.numeric( samprate )==0)
        {
          samprate = gsegy$deltaSample
        }

      
      dt = as.numeric( samprate )/1000000
####   you have to have a sample rate
      if(dt==0)
        {
          dt=0.025
        }
      
      DATIM =   c(gsegy$year, gsegy$day, gsegy$hour, gsegy$minute)

      
      sec = gsegy$second+  gsegy$m_secs/1000
      
      thesta= gsegy$station_name
      thecomp= gsegy$channel_name

      
      md = getmoday(DATIM[2], DATIM[1])

      t1 = 0
      t2 = dt*(N-1)
      

      tstart = list(yr=DATIM[1], jd=DATIM[2] , mo=md$mo, dom=md$dom,
        hr=DATIM[3], mi=DATIM[4], sec=sec, msec=0, dt=dt, t1=t1,
        t2=t2, off=0)

     
      coords = NULL

      if(RAW==TRUE)
        {
       y = ONEsegy$amp
       aunits="counts"

        }
      else
        {
          
      y = ONEsegy$amp/conversionfac
       aunits="volts"
    }
      
####  if(is.null(thesta))   thesta="XXX"
#### if(is.null(thecomp))  thecomp="X"
####  if(is.na(thesta))   thesta="XXX"
#### if(is.na(thecomp))  thecomp="X"
      
      GIVE[[i]] = list(fn=fn, sta=thesta,  comp=thecomp, dt=dt, DATTIM=tstart, N=N, units=aunits ,
            coords=coords ,  amp=y, HEAD=gsegy, IO=list(kind=1, Iendian=Iendian,  BIGLONG=BIGLONG) )



      if(PLOT>=0)
        {
          tee = seq(from=0, by=dt, length=length(y))
          plot(tee, y, type='l', main=paste(sep=' ', thesta, thecomp), xlab="time, s", ylab=aunits )
         
          if(PLOT==0 | PLOT==TRUE ) {
            print("left CLICK once in WINDOW for NEXT TRACE:")
            locator(1) }
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




