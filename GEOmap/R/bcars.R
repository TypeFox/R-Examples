`bcars` <-
function(x,y, h1=1, h2=.3, rot, col='black', border='black' )
  {

#########  plot(c(-1,1), c(-1,1), type='n', asp=1)
######   polygon(px,py)
######   polygon(ax,ay)
if(missing(col)) { col='black' }
if(missing(border)) {border =col }

###  if(missing(border)) { border='black' }

 ###     border=col 
#######  zim = cos(60*pi/180)
 
pin = par("pin")
u = par('usr')
##  h = 1 mm
##  r =0.875 mm
## uinch =   (u[4]-u[3])/pin[2]
###  2.54 cn per inch
umm =   (u[4]-u[3])/pin[2]/25.4
###  2.54 cn per inch
###25.4 mm/inch
stickh = h1*umm
   stickh2 = h2*umm


    px = c(-stickh, stickh, stickh, -stickh  )/2
    py=c(0, 0, stickh2, stickh2)
#####    polygon(px,py)


    cs = rot$cs
    sn = rot$sn

    for(i in 1:length(x))
      {
        
        ax = x[i]+(px*cs[i]-py*sn[i])
        ay = y[i]+(px*sn[i]+py*cs[i])
        polygon(ax,ay, col=col, border=border)
      }


  }

