`getpfile` <-
function(uwpickfile, stafile=NULL)
{
  if(missing(stafile))
    { stafile = NULL }
  else
    {
       sta =setstas(stafile)

    }

  if(is.null(stafile)) {  sta = NULL  }

  

  #########  if this already is in memory and is legit, return the pickfile as is
  if(is.list(uwpickfile))
    {
      if(length(uwpickfile$STAS)>1)
        {
          return(uwpickfile)

        }
      else
        {
          print("problem with pickfile")
          return(NULL)

        }
    }


  
  APF  = scan(file=uwpickfile, what="", sep="\n", quiet =TRUE)
  flg = substr(APF, 0,1)
  ACARD = APF[flg=="A"]
  LOC = unpackAcard(ACARD)
  PICS = APF[flg=="."]
  
  if(length(PICS)<1)
    {
      return(NULL)
    }
  
  MCARD = APF[flg=="M"]
  if(length(MCARD)>=1)
    {
      MC = list()
      for(  mnum in 1:length(MCARD)   )
        {
          F = list(az=0, dip=0)
          F$az = as.numeric(substr(MCARD, 5, 7))
          F$dip = as.numeric(substr(MCARD,8, 10))
          
          G = list(az=0, dip=0)
          G$az = as.numeric(substr(MCARD, 14, 16))
          G$dip = as.numeric(substr(MCARD,17, 19))
          
          U = list(az=0, dip=0)
          U$az = as.numeric(substr(MCARD, 22, 25))
          U$dip = as.numeric(substr(MCARD, 26, 28))
          
          V = list(az=0, dip=0)
          V$az = as.numeric(substr(MCARD, 32, 34))
          V$dip = as.numeric(substr(MCARD, 35, 37))
          
          P = list(az=0, dip=0)
          P$az = as.numeric(substr(MCARD, 41, 43))
          P$dip = as.numeric(substr(MCARD, 44, 46))
          
          T  = list(az=0, dip=0)
          T$az = as.numeric(substr(MCARD, 50, 52))
          T$dip = as.numeric(substr(MCARD, 53, 55))
          
          
          MC[[mnum]] = list(F=F, G=G, U=U, V=V, P=P, T=T, CNVRG=NA)
        }
    }
  else
    {
      MC = NULL
    }
  
  
  FCARD = APF[flg=="F"]
  if(length(FCARD)>=1)
    {

      EF = getFcard(FCARD)
      
    uf = unlist(strsplit(FCARD, " "))

    uf = uf[2:length(uf)]
    
    FF = as.numeric(uf)
    
    FF = FF[!is.na(FF)]
    phi = c(FF[2], FF[5], FF[8])
    lam = c(FF[1], FF[4], FF[7])
    valeig = c(FF[3], FF[6], FF[9])
    
    v = TOCART(phi, lam)

    LIP = v
    
    ## check vectors
   ##  sum(v[,3] * v[,2])
    ## sum(v[,2] * v[,1])
    ## sum(v[,3] * v[,1])

    
  }
  else
    {
      EF = NA
      LIP = NA
    }
  ECARD = APF[flg=="E"]
  if(length(ECARD)>=1)
    {
      E = getEcard(ECARD)
    }
  else
    {
      E = NA
    }
    
PICS = APF[flg=="."]

  if(length(PICS)<1)
    {
      return(NULL)
    }

##  fpics = substr(PICS, 0,1)

  STAS= NEW.getUWSTAS(PICS)


  
  if(!is.null(sta))
    {
      stmatch = match(STAS$name, sta$name)
      STAS$lat = sta$lat[stmatch]
      STAS$lon = sta$lon[stmatch]
      STAS$z = sta$z[stmatch]


    }

  UWFILEID = paste(sep="",
    formatC(LOC$yr, format="d", width=4, flag="0"),
    formatC(LOC$mo, format="d", width=2, flag="0"), 
    formatC(LOC$dom, format="d", width=2, flag="0"), 
    formatC(LOC$hr, format="d", width=2,  flag="0"), 
    formatC(LOC$mi, format="d", width=2,flag="0"),
    formatC(floor(LOC$sec), format="d", width=2,flag="0")
    )


  C1 =  APF[flg=="C"]
  if(length(C1)>0)
    {
      
      Comments = substr(C1, 3, 1000000)
    }
  else
    {
      Comments =NULL

    }

  
  OS1 = APF[flg=="O"]
  if(length(OS1)>0)
    {
      OSTAS = readUW.OSTAS(OS1)
    }
  else
    {
      OSTAS =NULL

    }

  
  ahh = APF[flg=="H"]
  
  
 if(length(ahh)>0)
    {
     AHH=getHcard(ahh)
    }
  else
    {
      AHH =NULL

    }




  
  nhh = APF[flg=="N"]
  
 
   if(length(nhh)>0)
    {
     Ncard = getNcard(nhh)
    }
  else
    {
       Ncard =NULL

    }


  return(list(PF=APF, AC=ACARD, LOC=LOC, MC=MC, STAS=STAS, LIP=LIP, E=E, F=EF, filename=uwpickfile, UWFILEID=UWFILEID, comments=Comments, OSTAS=OSTAS,H=AHH, N = Ncard ))
}

