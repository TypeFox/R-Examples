`thrust` <-
function(x,y, h=1, N=1, REV=FALSE,  endtol=.1, col='black', ...)
{
  if(missing(N))  N = 1
  if(missing(REV)) REV=FALSE
  if(missing(h))   { h=1  }

  if(missing(col)) { col='black' }
  if(missing(endtol)) { endtol=.1 }
  
  if(REV){ x= rev(x); y = rev(y) }


  pin = par("pin")
  u = par('usr')
##  h = 1 mm
##  r =0.875 mm
## uinch =   (u[4]-u[3])/pin[2]
###  2.54 cn per inch
  umm =   (u[4]-u[3])/pin[2]/25.4
###  2.54 cn per inch
###25.4 mm/inch
  stickh = h*umm
   
  
  g = PointsAlong(x, y, N=N,  endtol=endtol)

### vperp = list(x=-rot$sn, y=rot$cs)

###  plot(c(-1,1), c(-1,1), type='n', asp=TRUE)
  lines(x,y, col=col, ...)
  
  rot=list(cs = g$rot$cs , sn = g$rot$sn)


###  rot=list(cs = g$rot$cs , sn = g$rot$sn )

  teeth(g$x, g$y, stickh, g$rot, col=col)


  
}

