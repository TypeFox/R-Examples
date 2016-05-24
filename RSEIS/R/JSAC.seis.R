`JSAC.seis` <-function(fnames, Iendian=1, HEADONLY=FALSE , BIGLONG=FALSE, PLOT=-1, RAW=FALSE)
{
  ###  get a bunch of sac files from a directory and store in structure
  ####  
  
  if(missing(PLOT)) { PLOT=-1 }
  
  if(missing(Iendian)) { Iendian=1 }
  if(missing(HEADONLY)) {HEADONLY=FALSE }
  if(missing(BIGLONG)) { BIGLONG=FALSE}
  if(missing(RAW)) { RAW=FALSE }

isign = TRUE
   
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
 if(length(endianVEC)<length(fnames)) { endianVEC  = rep(endianVEC, times=length(fnames) ) }


  
  for(i in 1:length(fnames))
    {

      fn = fnames[i]
      infile = fn
      theENDIAN = endianVEC[i]
####   print(fn);
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
   sacheadnames = c("delta", "depmin", "depmax", "scale", "odelta", "b", 
     "e", "o", "a", "internal1", "t0", "t1", 
     "t2", "t3", "t4", "t5", "t6", "t7", 
     "t8", "t9", "f", "resp0", "resp1", "resp2", 
     "resp3", "resp4", "resp5", "resp6", "resp7", "resp8", 
     "resp9", "stla", "stlo", "stel", "stdp", "evla", 
     "evlo", "evel", "evdp", "unused1", "user0", "user1", 
     "user2", "user3", "user4", "user5", "user6", "user7", 
     "user8", "user9", "dist", "az", "baz", "gcarc", 
     "internal2", "internal3", "depmen", "cmpaz", "cmpinc", "unused2", 
     "unused3", "unused4", "unused5", "unused6", "unused7", "unused8", 
     "unused9", "unused10", "unused11", "unused12", "nzyear", "nzjday", 
     "nzhour", "nzmin", "nzsec", "nzmsec", "internal4", "internal5", 
     "internal6", "npts", "internal7", "internal8", "unused13", "unused14", 
     "unused15", "iftype", "idep", "iztype", "unused16", "iinst", 
     "istreg", "ievreg", "ievtyp", "iqual", "isynth", "unused17", 
     "unused18", "unused19", "unused20", "unused21", "unused22", "unused23", 
     "unused24", "unused25", "unused26", "leven", "lpspol", "lovrok", 
     "lcalda", "unused27", "kstnm", "kevnm[16]", "khole", "ko", 
     "ka", "kt0", "kt1", "kt2", "kt3", "kt4", 
     "kt5", "kt6", "kt7", "kt8", "kt9", "kf", 
     "kuser0", "kuser1", "kuser2", "kcmpnm", "knetwk", "kdatrd", 
     "kinst")
      
      zz <- file(fn, "rb")


      A1 = vector(length=70)
      for(j in 1:70)
        {
      A1[j] =  readBin(zz, numeric() , n = 1, size = ifloat, signed = isign,
        endian = theENDIAN)
    }


      A2 = vector(length=40)
      for(j in 1:40)
        {
      A2[j] =  readBin(zz, integer() , n = 1, size = ilong, signed = isign,
        endian = theENDIAN)
    }

   ## close(zz)

   ##  70+40+1+1+21
     ##  A4 =  getBINstr(zz, 8,  endian =theENDIAN)

     A4 = readChar(zz, 8, useBytes = FALSE)
      
      A5 = readChar(zz, 16, useBytes = FALSE)

      B4 = vector()
      for(k in 1:21) B4[k] = readChar(zz, 8, useBytes = FALSE)


      ALLHEAD = c(A1, A2, A4, A5, B4)


      SACH = data.frame(names=sacheadnames, values=ALLHEAD)

##############  this bizarre header is the beginning of the trace (like an offset)

      ###  b stands for beginning?????
      
      b = A1[6]
      
      N = A2[10]
      
        ##   close(zz)
      

     ## i1 = unlist( strsplit(split=" ", B4[4]) )
      
    thecomp1=  as.character(B4[18])

      sblank  = unlist( strsplit(thecomp1, split="") )

      thecomp=   paste(sblank[which(sblank!=" ")], collapse="")

      dt = as.numeric( A1[1] )
      
      DATIM =   as.numeric(A2[1:4] )

      ##################################  the b is an offset
      sec = A2[5]+ A2[6]/1000 + b 
      
      thesta1= as.character(A4)

      sblank  = unlist( strsplit(thesta1, split="") )

      thesta=   paste(sblank[which(sblank!=" ")], collapse="")

     ###   thecomp= B4[18]
      aunits="volts"
     ###  print(paste(sep=' ', infile, thesta, thecomp, aunits, N, dt, sec))
      
      md = getmoday(DATIM[2], DATIM[1])

      t1 = 0
      t2 = dt*(N-1)
      

      tstart = list(yr=DATIM[1], jd=DATIM[2] , mo=md$mo, dom=md$dom, hr=DATIM[3], mi=DATIM[4], sec=sec, msec=0, dt=dt, t1=t1,
        t2=t2, off=0)

  ###      if(is.null(thesta))   thesta="XXX"
  ###      if(is.null(thecomp))  thecomp="X"
   ###      if(is.na(thesta))   thesta="XXX"
   ###     if(is.na(thecomp))  thecomp="X"
           
      if(is.null(aunits))  aunits="volts"
     
      if(HEADONLY==TRUE)
        {
          
         ###  print(paste("headonly ", i))
          
          x = NULL
          aunits=NA
        }
      else
        {
         
          D1 = readBin(zz, numeric() , n = N , size =4,  endian = DATAendian[Iendian], signed = TRUE)
           N = length(D1)
          x = as.vector(D1)
        }

      
      close(zz)
      
      GIVE[[i]] = list(fn=fn, sta=thesta,  comp=thecomp, dt=dt, DATTIM=tstart, N=N, units=aunits , amp=x, HEAD=SACH, IO=list(kind=2, Iendian=Iendian,  BIGLONG=BIGLONG)   )

      
      if(PLOT>=0)
        {
          plot(x, type='l', main=paste(sep=' ', thesta, thecomp))
          print("left CLICK in WINDOW for NEXT TRACE:")
          if(PLOT==0) { locator(1) }
          else
            {
              Sys.sleep(PLOT);
            }
        }
    }
  invisible(GIVE)
}

