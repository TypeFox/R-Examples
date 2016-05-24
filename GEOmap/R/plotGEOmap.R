`plotGEOmap` <-
function(MAP, LIM=c(-180, -90, 180, 90) , shiftlon=0, add=TRUE , NUMB=FALSE , SEL=NULL, MAPcol=NULL, MAPstyle=NULL, border=NA,  PLOT=TRUE, PRINT=FALSE, BB=FALSE, ...)
{
  if(missing(add)) { add=FALSE }
  if(missing(NUMB)) { NUMB=FALSE }
  if(missing(shiftlon)) {  shiftlon=0 }
  if(missing(SEL)) {  SEL=NULL }
  if(missing(PLOT)) {  PLOT=TRUE }
  if(missing(MAPcol)) { MAPcol=NULL  }
  if(missing(MAPstyle)) { MAPstyle=NULL  }
  if(missing(border)) { border=NA  }
   if(missing(PRINT)) {  PRINT=FALSE }
   if(missing(BB)) {  BB=FALSE }
  
  if(missing(LIM))
    {
      lon = RPMG::fmod(MAP$POINTS$lon-shiftlon, 360)
      
      LIMP=c( min(lon), min(MAP$POINTS$lat), max(lon), max(MAP$POINTS$lat))
      LIM=c( min(RPMG::fmod(MAP$POINTS$lon, 360)), min(MAP$POINTS$lat), max(fmod(MAP$POINTS$lon, 360)), max(MAP$POINTS$lat))
    }
  else
    {

      if(is.null(LIM))
        {
          lon = RPMG::fmod(MAP$POINTS$lon-shiftlon, 360)
          LIMP=c( min(lon), min(MAP$POINTS$lat), max(lon), max(MAP$POINTS$lat))
          LIM=c( min(RPMG::fmod(MAP$POINTS$lon, 360)), min(MAP$POINTS$lat), max(fmod(MAP$POINTS$lon, 360)), max(MAP$POINTS$lat))
        }

      
      if(is.list(LIM))
        {
          
          lon = RPMG::fmod(LIM$lon-shiftlon, 360)
          lat = LIM$lat
          LIMP=c( min(lon), min(lat), max(lon), max(lat))
          LIM=c( min(RPMG::fmod(lon, 360)), min(lat), max(fmod(lon, 360)), max(lat))
          
        }
      else
        {
          LIMP=LIM
        }
      
    }
  
if(!is.null(MAPcol)) {  MAP$STROKES$col = rep(MAPcol, length(MAP$STROKES$col))  }
if(!is.null(MAPstyle)) {  MAP$STROKES$style = rep(MAPstyle, length(MAP$STROKES$style))  }


###   cat(paste(sep=" ", "###### plotGEOmap", paste(collapse=" ",LIMP) ) , sep="\n" )
  
###  determine stroke inclusion

###  (x2>=x3)&&(x4>=x1)&& (y2>=y3)&&(y4>=y1)

if( is.null(MAP$STROKES$LAT1)) { MAP = boundGEOmap(MAP)  }



  
  y1 = MAP$STROKES$LAT1
  y2 = MAP$STROKES$LAT2
  x1 = RPMG::fmod(MAP$STROKES$LON1-shiftlon, 360)
  x2 = RPMG::fmod(MAP$STROKES$LON2-shiftlon, 360)

  y3 = LIM[2]
  y4 = LIM[4]
  
  x3 = RPMG::fmod(LIM[1], 360)
  x4 = RPMG::fmod(LIM[3], 360)
  

  
  OUT = y1>=y4 | x1>=x4 | y2 <= y3 | x2 <= x3

  IN = which(!OUT)
  
 ###   MAP$STROKES$LAT1>=LIM[1] 
     
##########    for selecting parts of the map
  if(!is.null(SEL))
    {

      slin = which(IN %in%SEL)
      if(length(slin)>=1)
        { IN = IN[slin] }
    }
  
 ####     Kstroke = length(MAP$STROKES$num)
  ####   for(i in 1:Kstroke)


    if(add==FALSE)
    {
      ##  plot(MAP$POINTS$lon, MAP$POINTS$lat, type='n')
      ##  xlab="Lon", ylab="Lat",
      ## if(is.null(xlab)) { xlab="Lon" }
      
      plot(RPMG::fmod(MAP$POINTS$lon-shiftlon, 360), MAP$POINTS$lat, xlim=c(LIMP[1], LIMP[3])  , ylim=c(LIM[2], LIM[4]),
           type='n', xlab="lon", ylab='lat', axes=FALSE,  ...)

      if(PLOT==FALSE) { return(0) }

      
      axis(2)
      
      pp = axTicks(1)

      xlabs = RPMG::fmod(pp, 360)
      axis(1, at=pp, labels = xlabs)

      xlabs[xlabs>180 & xlabs<=359.99] = xlabs[xlabs>180 & xlabs<=359.99]-360
      axis(3, at=pp, labels =xlabs )

      
      box()
      
      
    }

  if(length(IN)<1) return(0)

  ###

  if(PRINT)
    {
     ### print(IN)
      print(paste(collapse=",", IN))
 
     ### print(sep="", "C(",  paste(collapse=",", IN), ")")

    }

     for(i in IN)
    {

      
      j1 = MAP$STROKES$index[i]+1
      j2 = j1+MAP$STROKES$num[i]-1

      LONS = RPMG::fmod(MAP$POINTS$lon[j1:j2]-shiftlon, 360)
      LATS = MAP$POINTS$lat[j1:j2]
       if(NUMB==TRUE)
        {

          
          numblats = length(LATS)
          points(LONS[1], LATS[1], pch=1, col="red")
          points(LONS[numblats], LATS[numblats], pch=6, col="blue")
          
          text( LONS[numblats], LATS[numblats] , labels=i, pos=3)
          text( LONS[1], LATS[1], labels=i, pos=4)
          ##  cat(paste(sep=" ",i,",") )
        }

      


      if(MAP$STROKES$style[i]==1)
        {
          points(LONS, LATS, col=MAP$STROKES$col[i], pch=".")
        }
       if(BB==TRUE)
        {
          rect(x1[i], y1[i], x2[i], y2[i], col=grey(.7) )

        }
      
      if(MAP$STROKES$style[i]==2)
        {
          mline = LONS
          dline = c(0, abs(diff(mline)))
          mline[dline>100] = NA
          lines( mline, LATS, col=MAP$STROKES$col[i])        
        }

      if(MAP$STROKES$style[i]==3)
        {

          mline = LONS
           dline = abs(diff(mline))

          ww  = which(dline>100)
         ##  mline[dline>100] = NA
          if(any(dline>100))
            {
            ##   print(paste('plotGEOmap', i, "wrapping", length(ww)  ))
              
         ## MM = fixCoastwrap(list(x=LONS, y=LATS), 100)
              if(MAP$STROKES$nam[i]=="ANTARCMAP")
                {
                  u = par("usr")
                  MM = fixCoastwrap(list(x=LONS, y=LATS), 100)
                  polygon(MM$x, MM$y, border=border, col=MAP$STROKES$col[i])
                  
                }
              else
                {

              
              xshift = LONS[1]
              v = RPMG::fmod(LONS-xshift, 360)
              x = v+xshift
              polygon(x, LATS,border=border, col=MAP$STROKES$col[i] )
              wv = which.max(x)
              xshift = (LONS[wv])
              v = x-x[wv]+LONS[wv]
              polygon(v,LATS,border=border, col=MAP$STROKES$col[i] )
              
            }
         
            }
          else
            {
              MM = list(x=LONS, y=LATS)
              polygon(MM$x, MM$y, border=border, col=MAP$STROKES$col[i])
            }
          
          
          
          
          
        }



      
      ##  locator()

    }

}

