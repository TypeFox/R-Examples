ternfocgeo<-function(CMTSOL, PROJ=NULL, icut=5, ndivs=10,   bbox=c(0,1, 0, 1)  , PLOT=TRUE, add=FALSE, RECT=FALSE)
  {

    if(missing(PLOT)) PLOT=TRUE
    if(missing(add)) add=FALSE
      if(missing(RECT)) RECT=FALSE
    if(missing(PROJ))
      {

        PROJ = GEOmap::setPROJ(type=2, LAT0=mean(CMTSOL$lat) , LON0=mean(CMTSOL$lon) )   ##   utm


      }

    if (missing(icut)) {
        icut = 5
    }
    if (missing(ndivs)) {
        ndivs = 10
    }
   if (missing(bbox)) {
        bbox = par("usr")
    }
     
    ####   CMT are a list of CMT solutions downloaded from the web
    ####   they have 
 ####      lon lat str1 dip1 rake1 str2 dip2 rake2 sc iexp name Elat Elon jd yr mo dom
 ####   CMT= list(lon=0, lat=0, str1=0, dip1=0, rake1=0, str2=0, dip2=0, rake2=0, sc=0, iexp=0, name=0, Elat=0, Elon=0, jd=0, yr=0, mo=0, dom=0)
    ## for this program ELAT andELON are not important
### read these from the screen dump on the internet - these are in old GMT style

###### OLDCMT = scan(file="OLDmechsGMT.txt", skip=1, what=list(lon=0, lat=0, str1=0, dip1=0, rake1=0, str2=0, dip2=0, rake2=0, sc=0, iexp=0, name=""))

######  mechs on this site have two forms for the name

###  see getCMT for info
    PHAIT = prepFOCS(CMTSOL)
    

XY = GEOmap::GLOB.XY(CMTSOL$lat, CMTSOL$lon, PROJ)

if(PLOT)
  {

    if(add==FALSE) plot(XY$x, XY$y, asp=1, type="p", xlab="East, km", ylab="North, km" )


    
    u = c(min(XY$x, na.rm=TRUE), min(XY$y, na.rm=TRUE), max(XY$x, na.rm=TRUE), max(XY$y, na.rm=TRUE))

    RI = RectDense( XY$x, XY$y, icut=icut, u=bbox, ndivs=ndivs)

    i = 1
    sizy = RI$icorns[i,4]-RI$icorns[i,2]
    sizx = RI$icorns[i,3]-RI$icorns[i,1]
    siz = .5*min(c(sizy, sizx))
    
    
    
   if(RECT) { rect(RI$icorns[,1],RI$icorns[,2],RI$icorns[,3],RI$icorns[,4], col=NA, border='blue') }
    for(i in 1:length(RI$ipass))
      {
        flag = XY$x>RI$icorns[i,1]& XY$y>RI$icorns[i,2] & XY$x<RI$icorns[i,3] & XY$y<RI$icorns[i,4]
        jh =PHAIT$h[flag]
        jv= PHAIT$v[flag]
        PlotTernfoc(jh,jv,x=mean(RI$icorns[i,c(1,3)]), y=mean(RI$icorns[i,c(2,4)]), siz=siz, fcols=PHAIT$fcols[flag], add=TRUE)
      }


    invisible(PHAIT)

  }






    
  }


