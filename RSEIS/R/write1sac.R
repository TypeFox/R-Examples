write1sac<-function(alist, fn=NULL, BIGLONG=FALSE  )
  {

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
      "unused9", "unused10", "unused11", "unused12",

      "nzyear", "nzjday", 
      "nzhour", "nzmin", "nzsec", "nzmsec", "internal4", "internal5", 
      "internal6", "npts", "internal7", "internal8", "unused13", "unused14", 
      "unused15",
      "iftype",
      "idep", "iztype", "unused16", "iinst", 
      "istreg", "ievreg", "ievtyp", "iqual", "isynth", "unused17", 
      "unused18", "unused19", "unused20", "unused21", "unused22", "unused23", 
      "unused24", "unused25", "unused26", "leven", "lpspol", "lovrok", 
      "lcalda", "unused27",

      "kstnm",

      "kevnm[16]",


      "khole","ko", 
      "ka", "kt0", "kt1", "kt2", "kt3", "kt4", 
      "kt5", "kt6", "kt7", "kt8", "kt9", "kf", 
      "kuser0", "kuser1", "kuser2", "kcmpnm", "knetwk", "kdatrd", 
      "kinst")


    

    
    theENDIAN =  .Platform$endian
    
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
##################   DEFAULT is for BIGLONG=FALSE ilong = ishort
        ishort = 2
        iint  = 4
        ilong = 4
        ifloat = 4
        idouble = 8
      }


    N = length(alist$amp)
    
    A1 = rep(-12345, length=70)
    NAM1 = sacheadnames[1:70]

    
    A2 = rep(-12345, length=40)
    NAM2  = sacheadnames[71:(70+40) ]


    
    A4 =format(alist$sta, width=8, justify="left" )
    NAM4 = sacheadnames[111]
    A5 = format("-12345", width=16, justify="left" )
    NAM5 = sacheadnames[112]
    
    B =  rep(format("-12345", width=8, justify="left" ), length=21)
    
    NAMB =  sacheadnames[113:length(sacheadnames) ]

    
    A1[1] = alist$dt
    
    A1[6] =alist$DATTIM$t1
    A1[7] = alist$DATTIM$t2


####  "stla", "stlo", "stel"

    
    A2[1] = alist$DATTIM$yr
    A2[2]=  alist$DATTIM$jd
    A2[3] = alist$DATTIM$hr
    A2[4] = alist$DATTIM$mi

    
    A2[5]  =trunc( alist$DATTIM$sec)
    A2[6] =  (alist$DATTIM$sec-A2[5])*1000
    A2[7]  =6
    
    A2[10]  =  N
    A2[16]  =  1    ##############  iftype = "Time Series File"

    A2[ which(NAM2=="leven") ] = 1

    B[ NAMB=="kcmpnm"] = format(alist$comp, width=8, justify="left" )
    
    sig = alist$amp

    if(is.null(fn))
      {
        staname =  alist$sta
        comp3 = alist$comp
        yr = alist$DATTIM$yr
        jd=  alist$DATTIM$jd
        hr= alist$DATTIM$hr
        mi = alist$DATTIM$mi
        trsec =trunc( alist$DATTIM$sec)
        
        sacfn = paste(sep=".",
           formatC(yr , width = 4, flag = "0") ,
            formatC(jd, width = 3, flag = "0"),
         formatC(hr, width = 2, flag = "0") ,
          formatC(mi, width = 2, flag = "0") ,
           formatC(trsec, width = 2, flag = "0")  ,
          staname,
          comp3,
          "SAC")

      }
    else
      {

        sacfn = fn
      }





    
#####################################   write the file
    zz <- file(sacfn, "wb")
    for(j in 1:70)
      {
        writeBin(A1[j] , zz, size = ifloat, 
                 endian = theENDIAN)
      }

    for(j in 1:40 )
      {
        writeBin(as.integer(A2[j]) , zz,  size = ilong, 
                 endian = theENDIAN)
      }

    writeChar(A4, zz, nchars = 8,eos=NULL)

    writeChar(A5, zz, nchars = 16,eos=NULL)

    for(k in 1:21)
      { writeChar(B[k] , zz, nchars = 8,eos=NULL) }


    ##  for(k in 1:N)

    writeBin(sig , zz,  size =ifloat,  endian =theENDIAN) 

    close(zz)

    invisible(sacfn)

  }
