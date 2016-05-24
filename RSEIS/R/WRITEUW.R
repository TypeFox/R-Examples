`writeUW.Acard` <-
function(LOC)
  {
    ## A 200406270357 30.03  0N4156  77W5069 12.00  4.1  0/000   0  0 0.00  0.0XX EC
    ############ writeUW.Acard(P$LOC)
    
      ID = paste(sep="",
        formatC(LOC$yr, format="d", width=4, flag="0"),
        formatC(LOC$mo, format="d", width=2, flag="0"), 
	formatC(LOC$dom, format="d", width=2, flag="0"), 
	formatC(LOC$hr, format="d", width=2,  flag="0"), 
	formatC(LOC$mi, format="d", width=2,flag="0"))

      L = abs(LOC$lat)
      LAT1 = floor(L)
      LAT2 = round((L - LAT1)*6000)
      if(LOC$lat<0){LATNS="S"}  else  {LATNS="N"}
          
      L = abs(LOC$lon)
      LON1 = floor(L)
      LON2 = round((L - LON1)*6000)
      if(LOC$lon<0) {LONEW="W"} else  {LONEW="E"}
          
      if(is.na(LOC$mag)) { LOC$mag=0 }
      if(is.na(LOC$z)) { LOC$z=0 }

### print(paste(sep=' ', LAT1, LATNS, LAT2, LON1, LONEW, LON2))
     
    AC = paste(sep='', "A ", ID, sprintf(fmt="%6.2f",LOC$se) , " ",
        sprintf(fmt="%02d", LAT1)   ,LATNS, sprintf(fmt="%04d", LAT2) ," ",
        sprintf(fmt="%03d",LON1),LONEW,  sprintf(fmt="%04d",LON2),
        sprintf(fmt="%6.2f",LOC$z), " ",
       sprintf(fmt="%4.1f", LOC$mag),
          "  0/000   0  0 0.00  0.0XX EC"
        )


      

      return(AC)

  }

`writeUW.Commentcard` <-
function(comments)
  {
    v = vector()
    for(i in 1:length(comments))
      {
        v[i] = paste(sep=" ", "C", comments[i])
      }
    return(v)

  }

`writeUW.DOTcard` <-
function(STAS)
  {

 
    STAS$pol[is.na(STAS$pol)] = "_"
    STAS$res[is.na(STAS$res)] = 0
    STAS$flg[is.na(STAS$flg)] = 0
    STAS$err[is.na(STAS$err)] = 0
    STAS$phase[is.na(STAS$phase)] = "P"
   

    
    
    tag = paste(sep='', ".", paste(sep=".", STAS$name, STAS$c3))

    pflag = STAS$phase=="P"

    
  
    pdots = paste(sep=" ", tag , rep(" (P",length=length(STAS$name)),     
          STAS$phase ,
          STAS$pol  ,
          STAS$sec ,
          STAS$flg ,
          STAS$err ,
          STAS$res   , rep(")",length=length(STAS$name)))

    wflag = !is.na(STAS$name) &   !is.na(STAS$c3)  & !STAS$name=="" & !STAS$c3==""
     
    alldots = c(pdots[wflag])
    
    
    return(alldots)
  }

`writeUW.Ecard` <-
function(E)
      {
       ecard = sprintf("E %2s%6.3f%6.3f%6.3f%6.3f%8.2f %3d%4s %5.2f%5.2f%5.2f%5.2f%5.2f%5.2f", E$LOC,  E$rms ,  E$meanres,  E$sdres, E$sdmean,  E$sswres,  E$ndf, E$fixflgs, E$sterrx, E$sterry, E$sterrz, E$sterrt, E$mag,  E$sterrmag )
       return(ecard)
      }

`writeUW.Fcard` <-
function(F)
  {
     fcard = sprintf("F %3.0f %2.0f %5.2f %3.0f %2.0f %5.2f %3.0f %2.0f %5.2f%7.2f%7.2f",
       F$azim1,  F$plunge1,  F$val1, F$azim2, F$plunge2, F$val2, F$azim3, F$plunge3, F$val3, F$herr, F$verr)
     
 return(fcard)
  }

`writeUW.Hcard` <-
function(H)
  {
   hcard =  paste(sep=" ", "H", H$yr, H$mo, H$dom, H$hr, H$mi, H$sec, H$lat, H$lon, H$z, H$mag) 
   return(hcard)
  }

`writeUW.Ncard` <-
function(N)
{
  ncard = paste(sep=" ", "N", N$name)
  
  return(ncard)
}

`writeUW.OSTAScard` <-
function(OSTAS)
  {
    
    v = vector()
    k = 0
    
    for(i in 0:(length(OSTAS)-1) )
      {
        #####print(paste(sep=' ', i, i %% 8))
        if( i %% 8 == 0 )
          {
            k = k + 1
            v[k] = paste(sep=" ", "O", OSTAS[i+1])
          }
        else
          {
            v[k] = paste(sep=" ", v[k], OSTAS[i+1])

          }

        
      }

        
    return(v)

  }

`writeUWpickfile` <-
function(A, output="")
  {
    if(missing(output)) output=""


############  need here to rectify the times in the pickfile
    ###  in UW format the times are all in seconds relative to the
    ####   minute on the Acard.
    ####   need to make sure these are reduced (i.e. less than 60 in general)
    



    
    
    cat(file=output  , writeUW.Acard(A$LOC), sep="\n")
    
    if(!is.null( A$E ))
      {
        if(!all(is.na( A$E )))
          {
            
            if(!is.na( A$E$sterrx ))
              {
                cat(file=output  , writeUW.Ecard(A$E), sep="\n", append=TRUE)
              }
          }
        
        
      }
    
    
    if(!is.null( A$F ))
      {
        
        if(!all(is.na( A$F )))
          {
            
            cat(file=output  , writeUW.Fcard(A$F), sep="\n", append=TRUE)
          }
      }

    if(!is.null( A$STAS ))
      {
        
        if(TRUE)
          {
            
            cat(file=output  , writeUW.DOTcard(A$STAS), sep="\n", append=TRUE)
          }
      }

    
    if(!is.null( A$N ))
      {
        if(!all(is.na( A$N )))
          {
            
            cat(file=output  ,   writeUW.Ncard(A$N) , sep="\n", append=TRUE)
          }
      }

    if(!is.null( A$H ))
      {
        if(!all(is.na( A$H )))
          {
            
            cat(file=output  ,   writeUW.Hcard(A$H) , sep="\n", append=TRUE)
          }
      }


    if(!is.null( A$OSTAS ))
      {
        if(!all(is.na( A$OSTAS )))
          {
            cat(file=output  ,   writeUW.OSTAScard(A$OSTAS) , sep="\n", append=TRUE)
          }
      }

    if(!is.null( A$comments ))
      {
        if(!all(is.na( A$comments )))
          {
            cat(file=output  ,   writeUW.Commentcard(A$comments) , sep="\n", append=TRUE)
          }
      }


    
  }

