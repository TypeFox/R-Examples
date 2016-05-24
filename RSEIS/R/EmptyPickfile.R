'EmptyPickfile'<-
  function(GH=list())
{
  #############   create an empty pick file for later processing

  PF =
    list(
         PF = "",   ## ascii version of the pickfile, character vector
         AC = "",   ## character string

         LOC=list(yr=0,
           jd=0,
           mo=0,
           dom=0,
           hr=0,
           mi=0,
           sec=0,
           lat=NA,
           lon=NA,
           z=NA,
           mag=0,
           gap=0,
           delta=0 ,
           rms=0,
           hozerr=0),
         

         MC=list(az1=0,
           dip1=0,
           az2=0,
           dip2=0,
           dir=0,
           rake1=0,
           dipaz1=0,
           rake2=0,
           dipaz2=0,
           F=list(az=0, dip=0),
           G=list(az=0, dip=0),
           U=list(az=0, dip=0),
           V=list(az=0, dip=0),
           P=list(az=0,dip=0),
           T=list(az=0, dip=0),
           sense=0,
           M=list( az1=0, d1=0,  az2=0, d2=0, uaz=0, ud=0, vaz=0, vd=0, paz=0, pd =0, taz=0, td=0),
           UP=TRUE,
           icol=1,
           ileg="",
           fcol='red',
           CNVRG="",
           LIM =c(0,0,0,0)),

         STAS=list(tag="",  ## character vector
           name="",
           comp="",
           c3="",
           phase="",
           sec=0,
           err=0,
           pol="",
           flg=0 ,
           res=0,
           lat=NA,
           lon=NA,
           z=NA,
           pdel=0,
           sdel=0
           ),

         LIP=vector(length=6),


         E=list(rms=NA,
           meanres=0,
           sdres=0,
           sdmean=0,
           sswres=0,
           ndf=0,
           fixflgs=0,
           sterrx=NA,
           sterry=0,
           sterrz=0,
           sterrt=0,
           mag=0,
           sterrmag=0),
         
         filename="", ## character string
         
         PICKER="", ## character string
         
         UWFILEID="", ## character string
         
         winID1="", ## character string
         
         comments="",## character vector
         
         OSTAS="",  ## character vector
         
         H=list(yr=0,
           mo=0,
           dom=0,
           hr=0,
           mi=0,
           sec=0,
           lat=0,
           lon=0,
           z=0,
           mag=0),
         
         N=list(name="")
         
         )


  if(missing(GH)) {
    return(PF)
  }
  else
    {
      
      PF$LOC=list(yr=GH$info$yr[1],
        jd=GH$info$jd[1],
        mo=GH$info$mo[1],
        dom=GH$info$dom[1],
        hr=GH$info$hr[1],
        mi=GH$info$mi[1],
        sec=GH$info$sec[1],
        lat=NA,
        lon=NA,
        z=NA,
        mag=0,
        gap=0,
        delta=0 ,
        rms=0,
        hozerr=0)
      
      PF$comments="FIRST TRACE"

      return(PF)
      
    }
  

  
}

  

