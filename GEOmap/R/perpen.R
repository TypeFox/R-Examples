`perpen` <-
function(x,y,h=1, rot=list(cs=1, sn=0) , col='black', lwd=1 )
  {
    if(missing(col)) col='black'
    if(missing(h)) h = 1
    if(missing(lwd))lwd=1

#########  plot(c(-1,1), c(-1,1), type='n', asp=1)
######   polygon(px,py)
######   polygon(ax,ay)
if(missing(col)) { col='black' }
###  if(missing(border)) { border='black' }
 if(missing(lwd)) { lwd=1 }

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
    
#######  zim = cos(60*pi/180)

    cs = -rot$sn
    sn = rot$cs

    for(i in 1:length(x))
      {
        
        ax = x[i]+(stickh*cs[i])
        ay = y[i]+(stickh*sn[i])
        segments(x[i], y[i], ax,ay, col=col)
      }


  }

