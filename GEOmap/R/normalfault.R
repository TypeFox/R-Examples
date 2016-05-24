normalfault<-function(x,y,h=1,hoff=1, rot=list(cs=1, sn=0) , col='black')
{

if(missing(col)) {  col='black' }
if(missing(h)) { h=1}
if(missing(hoff)) { hoff=1}

pin = par("pin")
u = par('usr')
##  h = 1 mm
##  r =0.875 mm
uinch =   (u[4]-u[3])/pin[2]
###  2.54 cn per inch
umm =   (u[4]-u[3])/pin[2]/25.4
###  2.54 cn per inch
###25.4 mm/inch
stickh = h*umm

rd = hoff*0.875*umm/2

ci = RPMG::circle(n=36)
rdot = list(x=rd*ci$x, y=rd*ci$y)
    cs = -rot$sn
    sn = rot$cs

 for(i in 1:length(x))
      {
        
        ax = x[i]+(stickh*cs[i])
        ay = y[i]+(stickh*sn[i])
        segments(x[i], y[i], ax,ay, col=col)
        polygon(ax+rdot$x,ay+rdot$y, col=col, border=col)
        
      }

}
