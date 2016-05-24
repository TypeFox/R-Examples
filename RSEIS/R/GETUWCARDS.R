`getEcard` <-
function(ECARD)
  {
    E = list()
    
    E$LOC = as.character(substr(ECARD, 2+1, 2+2))
    E$rms = as.numeric(substr(ECARD, 4+1, 4+6))
    E$meanres = as.numeric(substr(ECARD, 10+1, 10+6))
    E$sdres = as.numeric(substr(ECARD, 16+1, 16+6))
    E$sdmean = as.numeric(substr(ECARD, 22+1, 22+6))
    E$sswres = as.numeric(substr(ECARD, 28+1, 28+8))
    E$ndf = as.numeric(substr(ECARD, 37+1, 37+3))
    E$fixflgs= as.character(substr(ECARD, 41+1, 41+4))

    E$sterrx= as.numeric(substr(ECARD, 45+1, 45+5))
    E$sterry= as.numeric(substr(ECARD, 50+1, 50+5))
    E$sterrz = as.numeric(substr(ECARD, 55+1,  55+5))     
    E$sterrt= as.numeric(substr(ECARD, 60+1, 60+5))

    ##  these might be non-numeric
    E$mag  = as.numeric(substr(ECARD, 65+1, 65+5))    
    E$sterrmag= as.numeric(substr(ECARD, 70+1, 70+5))
 
    return(E)

  }

`getFcard` <-
function(FCARD)
  {

    F  = list(
      azim1= as.numeric(substr(FCARD, 3, 3+3)),
      plunge1= as.numeric(substr(FCARD, 7, 7+2)),
      val1= as.numeric(substr(FCARD, 10, 10+5)),
      
      azim2= as.numeric(substr(FCARD, 16, 16+3)),
      plunge2= as.numeric(substr(FCARD, 20, 20+2)),
      val2= as.numeric(substr(FCARD, 23, 23+5)),
      
      
      azim3= as.numeric(substr(FCARD, 29, 29+3)),
      plunge3= as.numeric(substr(FCARD, 33, 33+3)),
      val3= as.numeric(substr(FCARD, 36, 36+5)),
      herr= as.numeric(substr(FCARD, 41, 41+7)),
      verr  = as.numeric(substr(FCARD, 48, 48+7)))
    return(F)

  }

`getHcard` <-
function(hcard)
{
  h = unlist(strsplit(split=" ", hcard))
  h = h[!(h=="")]

  
  H = list(yr=as.numeric(h[2]),
    mo=as.numeric(h[3]),
    dom=as.numeric(h[4]),
    hr=as.numeric(h[5]),
    mi=as.numeric(h[6]),
    sec=as.numeric(h[7]),
    lat=as.numeric(h[8]),
    lon=as.numeric(h[9]),
    z=as.numeric(h[10]),
    mag=as.numeric(h[11]))

  return(H)
}

`getNcard` <-
function(ncard)
{
  h = unlist(strsplit(split=" ", ncard))
  N = list(name=h[2])
  return(N)
}

