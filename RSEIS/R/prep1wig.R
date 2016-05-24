`prep1wig` <-
function(wig=vector(), dt=1, sta="STA", comp="CMP", units="UNITS", starttime=list(yr=0, jd=1,mo=1,dom=1,hr=1,mi=1,sec=0) )
{
  if(missing(dt)) {dt=1}
  if(missing(sta)) {sta="STA"}
  if(missing(comp)) {comp="CMP"}
  if(missing(units)) {units="UNITS"}
  if(missing(starttime)) {
    starttime=list(yr=1972, jd=1,mo=1,dom=1,hr=0,mi=0,sec=0)
  }

  
  KLIST = list()
  KLIST[[1]] = list(
         amp=wig,
         dt=dt,
         nzyear=2000,
         nzjday=1,
         nzhour=1,
         nzmin=1,
         nzsec=0,
         nzmsec=0,
         b=0,
         e=0,
         o=0,
         fn="inputFILEname",
         sta=sta,
         comp=comp,
         DATTIM=
         list(yr=starttime$yr,
              jd=starttime$jd,
              mo=starttime$mo,
              dom=starttime$dom,
              hr=starttime$hr,
              mi=starttime$mi,
              sec=starttime$sec,
              msec=0,
              dt=dt,
              t1=0,
              t2=0,
              off=0),
         N=length(wig),
         units=units)
  

  return(KLIST)
  

  
}

