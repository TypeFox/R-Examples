SSfault<-function(x,y,h=1,hoff=.15, rot=list(cs=1, sn=0) , col='black', dextral=TRUE, lwd=1)
{
if(missing(dextral)) {  dextral=TRUE }
if(missing(col)) {  col='black' }
if(missing(h)) {h=1}
if(missing(hoff)) {hoff=1}
 if(missing(lwd))lwd=1

pin = par("pin")
u = par('usr')
##  h = 1 mm
##  r =0.875 mm
## uinch =   (u[4]-u[3])/pin[2]
###  2.54 cn per inch
umm =   (u[4]-u[3])/pin[2]/25.4
###  2.54 cn per inch
###25.4 mm/inch
stickh = h*umm*5.25
shoff = hoff*umm*.875



vperp = list(x=-rot$sn, y=rot$cs)

phi = 160*pi/180

m = 2


p1 = list(x= x-0.5*stickh*rot$cs+shoff*vperp$x  , y = y-0.5*stickh*rot$sn+shoff*vperp$y  )

p2 = list(x= x+0.5*stickh*rot$cs+shoff*vperp$x  , y = y+0.5*stickh*rot$sn+shoff*vperp$y   )

p3 = list(x= x-0.5*stickh*rot$cs-shoff*vperp$x  , y = y-0.5*stickh*rot$sn-shoff*vperp$y  )

p4 = list(x= x+0.5*stickh*rot$cs-shoff*vperp$x  , y = y+0.5*stickh*rot$sn-shoff*vperp$y   )


segments(p1$x, p1$y, p2$x, p2$y, col=col)
segments(p3$x, p3$y, p4$x, p4$y, col=col)


if(dextral)
{
phi = 160*pi/180

vleg1 = list(x= (cos(phi)*rot$cs-sin(phi)*rot$sn), y=sin(phi)*rot$cs+cos(phi)*rot$sn   )

vlen = sqrt(vleg1$x^2+vleg1$y^2)

vleg = list(x= vleg1$x/vlen, y=vleg1$y/vlen)

p5 = list(x=p2$x+shoff*m*vleg$x  , y = p2$y+shoff*m*vleg$y  )
p6 = list(x= p3$x-shoff*m*vleg$x  , y = p3$y-shoff*m*vleg$y  )
segments( p2$x, p2$y,p5$x, p5$y,  col=col)
segments( p3$x, p3$y,p6$x, p6$y, col=col)
}
else
{
phi = -160*pi/180
vleg1 = list(x= (cos(phi)*rot$cs-sin(phi)*rot$sn), y=sin(phi)*rot$cs+cos(phi)*rot$sn   )
vlen = sqrt(vleg1$x^2+vleg1$y^2)
vleg = list(x= vleg1$x/vlen, y=vleg1$y/vlen)
p5 = list(x=p1$x-shoff*m*vleg$x  , y = p1$y-shoff*m*vleg$y  )
p6 = list(x= p4$x+shoff*m*vleg$x  , y = p4$y+shoff*m*vleg$y  )
segments( p1$x, p1$y,p5$x, p5$y,  col=col)
segments( p4$x, p4$y,p6$x, p6$y, col=col)
}

}


