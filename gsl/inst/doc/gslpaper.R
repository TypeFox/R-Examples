### R code from vignette source 'gslpaper.Rnw'

###################################################
### code chunk number 1: gslpaper.Rnw:119-119
###################################################



###################################################
### code chunk number 2: gslpaper.Rnw:120-121
###################################################
library(gsl)


###################################################
### code chunk number 3: gslpaper.Rnw:126-127
###################################################
airy_Ai(1:3)


###################################################
### code chunk number 4: gslpaper.Rnw:140-154
###################################################
x <- seq(from=0,to=10,len=100)
plot(c(0,11),c(-1,1),type="n",main="Fig 10.6, p446",xlab="",ylab="",yaxt="n",xaxt="n",frame=FALSE)
axis(1,pos=0,at=c(0,2,4,6,8,10),labels=c("","2","4","6","8","10"))
axis(2,pos=0)
lines(x,airy_Ai       ( x),type="l",lty=1)
lines(x,airy_Ai       (-x),type="l",lty=2)
lines(x,airy_Ai_deriv ( x),type="l",lty=3)
lines(x,airy_Ai_deriv (-x),type="l",lty=4)
text(1,0.6     ,"Ai(-x)" )
text(0.85,0.33 ,"Ai(x)"  )
text(1.08,-0.26,"Ai'(x)" )
text(10.5,0.4  ,"Ai'(-x)")
arrows(10, 0, 11, 0,angle=11)
text(11,-0.1,"x")


###################################################
### code chunk number 5: gslpaper.Rnw:165-179
###################################################
x <- seq(from=0,to=10,len=100)
plot(c(0,10),c(-1,2.2),type="n",main="Fig 10.7, p446",xlab="",ylab="",yaxt="n",xaxt="n",frame=FALSE)
axis(1,pos=0,at=c(0,1:9),labels=c("","1","2","3","4","5","6","7","8","9"))
axis(2,pos=0)
lines(x,airy_Bi       ( x),type="l",lty=1)
lines(x,airy_Bi       (-x),type="l",lty=2)
lines(x,airy_Bi_deriv ( x),type="l",lty=3)
lines(x,airy_Bi_deriv (-x),type="l",lty=4)
text(0.15,1.44     ,"Bi(x)",pos=4)
text(1,0.90 ,"Bi'(x)",pos=4)
text(2.25,0.56,"Bi'(-x)")
text(0.7,-0.55,"Bi'(-x)",pos=4)
arrows(9, 0, 10, 0, angle=11)
text(10,-0.1,"x")


###################################################
### code chunk number 6: gslpaper.Rnw:277-290
###################################################
f <- function(r,n){ 
-airy_Ai(r+airy_zero_Ai(n+1))/airy_zero_Ai_deriv(n+1)}
plot(c(0,10),c(0,10),type="l",yaxt="n",xaxt="n",frame=FALSE,xlab="r",ylab="V(r)")
axis(1,pos=0)
axis(2,pos=0)

x <- seq(from=0,to=10,len=400)
for(i in 0:5){
  jj <- -airy_zero_Ai(i+1)
  lines(x=c(0,jj),y=c(jj,jj))
  lines(x=c(jj,10),y=c(jj,jj),col="gray",lty=2)
  points(x,(i+1)*(-1)^i*f(x,i)+jj,type="l")
}


