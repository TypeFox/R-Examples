`OverTurned` <-
function(x,y, syn=TRUE, spacing=NULL, N=1, r1= 1, r2= 1.2, h1= .5, h2= .5, endtol=.1,
                    REV=FALSE, col='black', ...)
{
  if(missing(spacing))  spacing=NULL
  if(missing(REV)) REV=FALSE
  if(missing(r1)) { r1 = 1 }
  if(missing(r2)) { r2 = 1.2 }
  if(missing(h1)) { h1 = .5 }
  if(missing(h2)) { h2 = .5 }

  if(missing(col)) { col='black' }
  if(missing(N)) { N = 1 }
    if(missing(syn)) { syn=TRUE }

  
  if(REV){ x= rev(x); y = rev(y) }
     if(missing(endtol)) { endtol=.1 }


  n = length(x)
  
 g = PointsAlong(x, y, N=N,  endtol=endtol)

  lines(x,y, col=col, ...)

 ##  arrows(x[n-1], y[n-1], x[n], y[n], col=col)
  
 ## g$rot$sn = -g$rot$sn

  
  HOR= horseshoe(g$x  , g$y , r1=r1, r2=r2, h2=h2, h1=h1, rot=g$rot, col=col)
 ##  fin = par("fin")
  pin = par("pin")
  u = par("usr")
   umm =   (u[4]-u[3])/pin[2]/25.4
  lenarr =  umm*h1


##  phi = 160*pi/180

##  vleg1 = list(x= (cos(phi)*g$rot$cs-sin(phi)*g$rot$sn), y=sin(phi)*g$rot$cs+cos(phi)*g$rot$sn   )

##  vlen = sqrt(vleg1$x^2+vleg1$y^2)

##  vleg = list(x= vleg1$x/vlen, y=vleg1$y/vlen)

##  p5 = list(x=g$x+lenarr*vleg$x  , y = g$y+lenarr*vleg$y  )
 
##  p6 = list(x= p3$x-shoff*m*vleg$x  , y = p3$y-shoff*m*vleg$y  )
##  segments( p2$x, p2$y,p5$x, p5$y,  col=col)
##  segments( p3$x, p3$y,p6$x, p6$y, col=col)
  
  
  if(syn==FALSE)
    {
      for(i in 1:length(HOR))
        {
          x = HOR[[i]]$x
          y = HOR[[i]]$y
          m = length(x)


         ## p5 = list(x=x[1]+lenarr*vleg$x[i]  , y = y[1]+lenarr*vleg$y[i]  )
          
         ## segments( x[1], y[1],p5$x, p5$y,  col=col)

          arrows(x[1],y[1], x[2], y[2], col=col, length=lenarr)
          arrows( x[m], y[m], x[m-1], y[m-1], col=col, length=lenarr)

          
        }
      
    }
  else
    {

      for(i in 1:length(HOR))
        {
          x = HOR[[i]]$x
          y = HOR[[i]]$y
          m = length(x)
          arrows( x[2], y[2], x[1],y[1],col=col, length=lenarr)
          arrows(  x[m-1], y[m-1], x[m], y[m],col=col, length=lenarr)
        }
      
      
    }

}

