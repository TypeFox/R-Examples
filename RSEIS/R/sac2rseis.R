sac2rseis<-function(fnames, Iendian=1 , HEADONLY=FALSE, BIGLONG=FALSE,  PLOT=-1, RAW=FALSE )
{
###  get a bunch of sac files from a directory and store in structure
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


      ONEsac =   read1sac(fn, Iendian=theENDIAN , HEADONLY=FALSE, BIGLONG=BIGLONG ) ;

      gsac = ONEsac$HEAD

      gainConst = 1
      scalefac  = 1

      N = gsac$npts

      
      if( !is.integer(N))
        {
          ##  this is a problem
          print("ERROR: number of samples is not an integer.")
          print(paste(i, fn))

        }
      
       thecomp1=  gsac$kcmpnm

      sblank  = unlist( strsplit(thecomp1, split="") )

      thecomp=   paste(sblank[which(sblank!=" ")], collapse="")

      dt = gsac$delta
      
      
      
yr = gsac$nzyear

jd = gsac$nzjday
hr =  gsac$nzhour
 mi =  gsac$nzmin
      msec  =   gsac$nzmsec    
 sec =  gsac$nzsec+msec/1000+ gsac$b 
 
      
      ##################################  the b is an offset
      
      
      thesta1= as.character(gsac$kstnm)

      sblank  = unlist( strsplit(thesta1, split="") )

      thesta=   paste(sblank[which(sblank!=" ")], collapse="")

     ###   thecomp= B4[18]
      aunits="volts"
     ###  print(paste(sep=' ', infile, thesta, thecomp, aunits, N, dt, sec))
      
      md = getmoday(jd, yr)

      t1 = 0
      t2 = dt*(N-1)
      

      tstart = list(yr=yr, jd=jd , mo=md$mo, dom=md$dom, hr=hr, mi=mi, sec=sec, msec=0, dt=dt, t1=t1,
        t2=t2, off=0)

      
####   you have to have a sample rate
      if(dt==0)
        {
          dt=0.025
        }
      
      DATIM =   c(yr, jd, hr, mi)

      aunits="unknown"
      coords = NULL

      y = ONEsac$amp
      
####  if(is.null(thesta))   thesta="XXX"
#### if(is.null(thecomp))  thecomp="X"
####  if(is.na(thesta))   thesta="XXX"
#### if(is.na(thecomp))  thecomp="X"
      
      GIVE[[i]] = list(fn=fn, sta=thesta,  comp=thecomp, dt=dt, DATTIM=tstart, N=N, units=aunits ,
            coords=coords ,  amp=y , HEAD=gsac, IO=list(kind=2, Iendian=Iendian,  BIGLONG=BIGLONG) )

      if(PLOT>=0)
        {
          tee = seq(from=0, by=dt, length=length(y))
          plot(tee, y, type='l', main=paste(sep=' ', thesta, thecomp), xlab="time, s", ylab=aunits )
          
          if(PLOT==0 | PLOT==TRUE)
            {
              print("left CLICK once in WINDOW for NEXT TRACE:")
              locator(1)
            }
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




