`plotworldmap` <-
function(MAP, LIM=c(-180, -90, 180, 90) , shiftlon=0, add=TRUE , NUMB=FALSE , ...)
{
  if(missing(add)) { add=FALSE }
  if(missing(NUMB)) { NUMB=FALSE }
  if(missing(shiftlon)) {  shiftlon=0 }

  
  if(missing(LIM))
    {
      lon = RPMG::fmod(MAP$POINTS$lon-shiftlon, 360)
      
      LIMP=c( min(lon), min(MAP$POINTS$lat), max(lon), max(MAP$POINTS$lat))
      LIM=c( min(RPMG::fmod(MAP$POINTS$lon, 360)), min(MAP$POINTS$lat), max(fmod(MAP$POINTS$lon, 360)), max(MAP$POINTS$lat))
    }
  else
    {

      if(LIM[1]==(-180) & LIM[3]==(180))
        {
          LIM[1] =  min(RPMG::fmod(MAP$POINTS$lon, 360))
          LIM[3] =  max(RPMG::fmod(MAP$POINTS$lon, 360))
        }

      if(LIM[1]==(0) & LIM[3]==(360))
        {
          LIM[1] =  min(RPMG::fmod(MAP$POINTS$lon, 360))
          LIM[3] =  max(RPMG::fmod(MAP$POINTS$lon, 360))
        }

      

      LIM = c(  RPMG::fmod(LIM[1], 360)  , LIM[2] , fmod(LIM[3], 360) , LIM[4]) 
      
      LIMP=LIM
    }
  

###  determine stroke inclusion

###  (x2>=x3)&&(x4>=x1)&& (y2>=y3)&&(y4>=y1)

  y1 = MAP$STROKES$LAT1
  y2 = MAP$STROKES$LAT2
  x1 = RPMG::fmod(MAP$STROKES$LON1, 360)
  x2 = RPMG::fmod(MAP$STROKES$LON2, 360)

  y3 = LIM[2]
  y4 = LIM[4]
  
  x3 = RPMG::fmod(LIM[1], 360)
  x4 = RPMG::fmod(LIM[3], 360)
  

  
  OUT = y1>=y4 | x1>=x4 | y2 <= y3 | x2 <= x3

  IN = which(!OUT)
  
 ###   MAP$STROKES$LAT1>=LIM[1] 
###  print(IN)
  
  Kstroke = length(MAP$STROKES$num)
  ####   for(i in 1:Kstroke)


    if(add==FALSE)
    {
      ##  plot(MAP$POINTS$lon, MAP$POINTS$lat, type='n')
      ##  xlab="Lon", ylab="Lat",
      ## if(is.null(xlab)) { xlab="Lon" }
      
      plot(RPMG::fmod(MAP$POINTS$lon-shiftlon, 360), MAP$POINTS$lat, xlim=c(LIMP[1], LIMP[3])  , ylim=c(LIM[2], LIM[4]), type='n', xlab="lon", ylab='lat', axes=FALSE,  ...)
      u = par('usr')
      axis(2)
      
      pp = axTicks(1)

      xlabs = RPMG::fmod(pp+1, 360)
      axis(1, at=pp, labels = xlabs)

      xlabs[xlabs>180 & xlabs<=359.99] = xlabs[xlabs>180 & xlabs<=359.99]-360
      axis(3, at=pp, labels =xlabs )

      j = RPMG::fmod(seq(from=(-180), to=180, by=6)-shiftlon, 360)
      jhalf = RPMG::fmod(seq(from=(-177), to=177, by=6)-shiftlon, 360)
      text(jhalf , u[3], labels=1:length(jhalf), pos=3)
      
      segments(j,  u[3]    , j   ,      u[4]      ,     lty=2, col=rgb(.9, .9, .9) )

      
       abline(v=RPMG::fmod(0-shiftlon, 360), col=rgb(.8, .8, .8) )
      abline(h=0, col=rgb(.8, .8, .8))


       htags =  LETTERS[6:23]
       htags = htags[-c(4, 10)]

      k = seq(from =-56, to =72, by=8)
      
      segments(u[1],  k    , u[2]   ,  k      ,     lty=2, col=rgb(.9, .9, .9) )
       ty = k+4
      text(u[2],ty[1:length(htags)],     labels=htags , xpd=TRUE, pos=4 )

      
    box()
     

    }



     for(i in IN)
    {

      
      j1 = MAP$STROKES$index[i]+1
      j2 = j1+MAP$STROKES$num[i]-1

      LONS = RPMG::fmod(MAP$POINTS$lon[j1:j2]-shiftlon, 360)

       if(NUMB==TRUE)
        {
          points(LONS[1], MAP$POINTS$lat[j1])
          text( LONS[1], MAP$POINTS$lat[j1], labels=i, pos=3)
        }


      if(MAP$STROKES$style[i]==1)
        {
          points(LONS, MAP$POINTS$lat[j1:j2], col=MAP$STROKES$col[i])
        }

      
      if(MAP$STROKES$style[i]==2)
        {
          mline = LONS
          dline = c(0, abs(diff(mline)))
          mline[dline>100] = NA
          lines( mline, MAP$POINTS$lat[j1:j2], col=MAP$STROKES$col[i])
        }

      if(MAP$STROKES$style[i]==3)
        {
          polygon(LONS, MAP$POINTS$lat[j1:j2], border=FALSE, col=MAP$STROKES$col[i])
        }



      
      ##  locator()

    }

}

