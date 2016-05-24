set.template <-
function(void=TRUE) 
{

#### generates visual field template
#### void == TRUE produces a hashed area at 0 deg and 180 deg


par(pty='m')
par(xpd=TRUE)


eqscplot(c(-90,90), c(-70,70), ylim=c(-70, 70), axes=F,
type="n",xlab="",ylab="")
rad.aux<- 92.5

if(void)
{
polygon(c(-rad.aux, -70, -65, -70, -rad.aux),c( 5,5,0,-5,-5), density=25, col='lightgray')
polygon(c(rad.aux, 70, 65, 70, rad.aux),c( 5,5,0,-5,-5), density=25, col='lightgray')
}

draw.circle( 0,0, radius=30,border='gray', lwd=2 ,nv=500)
seq.angles<- seq(0, 180, by=15)
ind.angles1<- c(1:3, 9:11 )### complete lines
for ( i in ind.angles1) 
{
segments(0,0, 90*cos(2*pi*i/24), 90*sin(2*pi*i/24),
               lwd=1/2, col='gray')
segments(0,0, -90*cos(2*pi*i/24), -90*sin(2*pi*i/24),
               lwd=1/2, col='gray')
}

### obtain corners
upper.corner<- sqrt(c(80,90)^2 - 70^2)
segments(-upper.corner, 70, upper.corner, 70,col='gray', lwd=1)

lower.corner<- sqrt( c(80,90)^2 - 74^2)
segments(-lower.corner, -74, lower.corner, -74,col='gray', lwd=1)

draw.arc(0,0, 90, angle1=asin(70/90), angle2=asin(-74/90), col='gray',lwd=1)
draw.arc(0,0, 80, angle1=asin(70/80), angle2=asin(-74/80), col='gray',lwd=1)


draw.arc(0,0, 90, deg1= 180+asin(-70/90)*180/pi, 
                  deg2= 180+asin(74/90)*180/pi,
                  col='gray',lwd=1)

draw.arc(0,0, 80, deg1= 180+asin(-70/80)*180/pi, 
                  deg2= 180+asin(74/80)*180/pi,
                  col='gray',lwd=1)

### remove upper and lower caps

draw.arc(0,0,90, deg1=180+asin(-70/90)*180/pi,
                 deg2=asin(70/90)*180/pi, col='white', 
                       lwd=1)
### draw major axes

segments(0, -74, 0, 70, col='gray')
segments(-90, 0, 90, 0, col='gray')                      

segments(rep(0,4), rep(0,4),
         70/tan(rad(c(60,75,105,120))), rep(70,4), 
            col='gray', lwd=1)


segments(rep(0,4), rep(0,4),
         -74/tan(rad(c(240,255,285,300))), rep(-74,4), 
            col='gray', lwd=1)


for ( i in seq(10,70,by=10)) 
{
draw.circle( 0,0, radius=i,border='gray', nv=500)
text( 0, -i,  i, cex=2/3, col='slategray')
text( 0, i, i, cex=2/3, col='slategray')
text( -i,2, i, cex=2/3, col='slategray')
text(i,2,i,cex=2/3, col='slategray')
}


text(c(-80,-90,80,90),rep(2,4), rep(c(80,90),2),
cex=2/3, col='slategray')

text(par()$usr[1] -4, 0, "180", col='darkslategray',cex=3/5, family='sans')
text(par()$usr[2] + 2, 0, "0", col='darkslategray',cex=3/5, family='sans')


for (i in seq.angles[-c(1, 5:13, length(seq.angles))])  
   text(  rad.aux*cos(rad(i)), rad.aux*sin(rad(i)), i, 
   col='darkslategray',cex=3/5, family='sans',pos=4, 
   offset=1/4)

for (i in seq.angles[5:9])  
   text(  rad.aux*cos(rad(i)), rep(70,5), i, 
   col='darkslategray',cex=3/5, family='sans',pos=3, 
   offset=1/4)
   
for (i in seq.angles[-c(1, 2:9, length(seq.angles))])  
   text(  rad.aux*cos(rad(i)), rad.aux*sin(rad(i)), i, 
   col='darkslategray',cex=3/5, family='sans',pos=2, 
   offset=1/4)
   
for (i in 180+seq.angles[-c(1,5:13,length(seq.angles))]) 
   text( rad.aux*cos(rad(i)), rad.aux*sin(rad(i)), i, 
   col='darkslategray',cex=3/5, family='sans', pos=2,
   offset=1/4)

for (i in 180+seq.angles[5:9])  
   text(  rad.aux*cos(rad(i)), -rep(74,5), i, 
   col='darkslategray',cex=3/5, family='sans',pos=1, 
   offset=1/4)
   
 for (i in 180+seq.angles[-c(1, 2:9, length(seq.angles))])  
   text(  rad.aux*cos(rad(i)), rad.aux*sin(rad(i)), i, 
   col='darkslategray',cex=3/5, family='sans',pos=4, 
   offset=1/4)

}
