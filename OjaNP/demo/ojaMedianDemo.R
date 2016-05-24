# naive function to calculate the intersection point of two lines
intersection <- function(a,b,c,d){
	fab <-function(a,b,x){
		slope <- (b[2]-a[2])/(b[1]-a[1])
		if(slope==Inf) slope <- 0
		achse <- a[2]-a[1]*slope
		slope*x+achse
	}
	fcd <-function(c,d,x){
		slope <- (d[2]-c[2])/(d[1]-c[1])
		if(slope==Inf) slope <- 0
		achse <- c[2]-c[1]*slope
		slope*x+achse
		return(slope)
	}
	xValue<-uniroot(function(x) fab(a,b,x)-fab(c,d,x),lower=-9999999999999999999999999999999999,upper=999999999999999999999999999999999999)$root
	yValue<-fab(a,b,xValue) 
	result<-c(xValue,yValue) 
	return(result)
}

# The settings for the layout of the plot

#op<-par(mfrow=c(2,3),pty="s",mex=0.001)

# All graphics are plotted in the same way. The axis are drawn with arrows, the connection lines with lines and the datapoints with points.
# The intersection function is needed to determin the intersection of the connection lines.

# The n=2 case
#plot(c(7.5),c(2),,col="white",xlim=c(0,10),ylim=c(0,7.5),xaxt="n",yaxt="n",xlab="",ylab="",lty=0,mar=c(0,0,0,0),mfg=c(1,1),bty="n")
plot(c(7.5),c(2),,col="white",xlim=c(0,10),ylim=c(0,7.5),xaxt="n",yaxt="n",xlab="",ylab="",lty=0,mar=c(0,0,0,0),bty="n")
arrows(0.5,0.5,0.5,7,angle=30,lwd=1)
arrows(0.5,0.5,10,0.5,angle=30,lwd=1)
lines(c(2,9),c(1.2,6.8),lwd=1,col="red")
lines(c(9,10),c(6.8,7.6),lwd=1,col="red",lty=3)
lines(c(1,2),c(0.4,1.2),lwd=1,col="red",lty=3)
points(c(3,8),c(2,6),pch=20,lwd=4)
#text(9.5,6,paste("n=2"),cex=1.85)
title(main = "Oja median in the bivariate case for n=2",sub="red area corresponds to the Oja median")


# The n=3 case
#plot(c(7.5),c(2),xlim=c(0,10),ylim=c(0,7.5),xaxt="n",yaxt="n",xlab="",ylab="",lty=0,mar=c(0,0,0,0),mfg=c(1,2),bty="n")
plot(c(7.5),c(2),xlim=c(0,10),ylim=c(0,7.5),xaxt="n",yaxt="n",xlab="",ylab="",lty=0,mar=c(0,0,0,0),bty="n")
arrows(0.5,0.5,0.5,7,angle=30,lwd=1)
arrows(0.5,0.5,10,0.5,angle=30,lwd=1)
#text(9.5,6,paste("n=3"),cex=1.85)
polygon(c(7.5,1.5,4),c(2,4.5,7),col="red",lwd=1.5)
points(7.5,2,pch=20,lwd=4)
points(1.5,4.5,pch=20,lwd=4)
points(4,7,pch=20,lwd=4)
title(main = "Oja median in the bivariate case for n=3",sub="red area corresponds to the Oja median")


# The n=4 case
#plot(c(7.5),c(2),xlim=c(0,10),ylim=c(0,7.5),xaxt="n",yaxt="n",xlab="",ylab="",lty=0,mar=c(0,0,0,0),mfg=c(1,3),bty="n")
plot(c(7.5),c(2),xlim=c(0,10),ylim=c(0,7.5),xaxt="n",yaxt="n",xlab="",ylab="",lty=0,mar=c(0,0,0,0),bty="n")
arrows(0.5,0.5,0.5,7,angle=30,lwd=1)
arrows(0.5,0.5,10,0.5,angle=30,lwd=1)
#text(9.5,6,paste("n=4"),cex=1.85)
sch6 <- intersection(c(2.5,1.2),c(4,7),c(1.4,4.5),c(7.5,2))
polygon(c(7.5,1.4,4,2.5),c(2,4.5,7,1.2),lwd=1)
lines(c(1.4,2.5),c(4.5,1.2),lwd=1)
lines(c(4,7.5),c(7,2),lwd=1)
points(7.5,2,pch=20,lwd=4)
points(1.4,4.5,pch=20,lwd=4)
points(4,7,pch=20,lwd=4)
points(2.5,1.2,pch=20,lwd=4)
points(sch6[1],sch6[2],pch=4,lwd=3,cex=1.5,col="red")
title(main = "Oja median in the bivariate case for n=4",sub="red area corresponds to the Oja median")


# The n=5 case
#plot(c(7.5),c(2),xlim=c(0,10),ylim=c(0,7.5),xaxt="n",yaxt="n",xlab="",ylab="",lty=0,mar=c(0,0,0,0),mfg=c(2,1),bty="n")
plot(c(7.5),c(2),xlim=c(0,10),ylim=c(0,7.5),xaxt="n",yaxt="n",xlab="",ylab="",lty=0,mar=c(0,0,0,0),bty="n")
arrows(0.5,0.5,0.5,7,angle=30,lwd=1)
arrows(0.5,0.5,10,0.5,angle=30,lwd=1)
#text(9.5,6,paste("n=5"),cex=1.85)
sch1 <- intersection(c(2.5,1.2),c(4,7),c(1.4,4.5),c(6.7,6.4))
sch2 <- intersection(c(2.5,1.2),c(4,7),c(1.4,4.5),c(7.5,2))
sch3 <- intersection(c(2.5,1.2),c(6.7,6.4),c(1.4,4.5),c(7.5,2))
sch4 <- intersection(c(2.5,1.2),c(6.7,6.4),c(4,7),c(7.5,2))
sch5 <- intersection(c(1.4,4.5),c(6.7,6.4),c(4,7),c(7.5,2))
polygon(c(sch1[1],sch2[1],sch3[1],sch4[1],sch5[1]),c(sch1[2],sch2[2],sch3[2],sch4[2],sch5[2]),col="red",lwd=1.5)
polygon(c(6.7,1.4,7.5,4,2.5,6.7,4,1.4,2.5,7.5,6.7),c(6.4,4.5,2,7,1.2,6.4,7,4.5,1.2,2,6.4),lwd=1)
points(7.5,2,pch=20,lwd=4)
points(1.4,4.5,pch=20,lwd=4)
points(4,7,pch=20,lwd=4)
points(2.5,1.2,pch=20,lwd=4)
points(6.7,6.4,pch=20,lwd=4)
title(main = "Oja median in the bivariate case for n=5",sub="red area corresponds to the Oja median")


# The n=6 case
#plot(c(7.5),c(2),xlim=c(0,10),ylim=c(0,7.5),xaxt="n",yaxt="n",xlab="",ylab="",lty=0,mar=c(0,0,0,0),mfg=c(2,2),bty="n")
plot(c(7.5),c(2),xlim=c(0,10),ylim=c(0,7.5),xaxt="n",yaxt="n",xlab="",ylab="",lty=0,mar=c(0,0,0,0),bty="n")
arrows(0.5,0.5,0.5,7,angle=30,lwd=1)
arrows(0.5,0.5,10,0.5,angle=30,lwd=1)
#text(9.5,6,paste("n=6"),cex=1.85)
sch7 <- intersection(c(2.5,1.2),c(6.7,6.4),c(1.4,4.5),c(8,4))
polygon(c(6.7,1.4,7.5,4,2.5,6.7,4,1.4,2.5,7.5,6.7),c(6.4,4.5,2,7,1.2,6.4,7,4.5,1.2,2,6.4),lwd=1)
points(7.5,2,pch=20,lwd=4)
points(1.4,4.5,pch=20,lwd=4)
points(4,7,pch=20,lwd=4)
points(2.5,1.2,pch=20,lwd=4)
points(6.7,6.4,pch=20,lwd=4)
points(8.5,4,pch=20,lwd=4)
points(sch7[1],sch7[2],pch=4,lwd=3,cex=1.5,col="red")
lines(c(8.5,7.5),c(4,2),lwd=1)
lines(c(8.5,1.4),c(4,4.5),lwd=1)
lines(c(8.5,4),c(4,7),lwd=1)
lines(c(8.5,2.5),c(4,1.2),lwd=1)
lines(c(8.5,6.7),c(4,6.4),lwd=1)
title(main = "Oja median in the bivariate case for n=6",sub="red area corresponds to the Oja median")


# The n=7 case
#plot(c(7.5),c(2),col="white",xlim=c(0,10),ylim=c(0,7.5),xaxt="n",yaxt="n",xlab="",ylab="",lty=0,mar=c(0,0,0,0),mfg=c(2,3),bty="n")
plot(c(7.5),c(2),col="white",xlim=c(0,10),ylim=c(0,7.5),xaxt="n",yaxt="n",xlab="",ylab="",lty=0,mar=c(0,0,0,0),bty="n")
arrows(0.5,0.5,0.5,7,angle=30,lwd=1)
arrows(0.5,0.5,10,0.5,angle=30,lwd=1)
#text(9.5,6,paste("n=7"),cex=1.85)
points(c(1,8.5,2.5,5,8.5,3.5,7.5),c(1,2.5,5,6.5,5,3.5,1.5),pch=20,lwd=4)
lines(c(1,8.5),c(1,2.5),lwd=1)
lines(c(1,2.5),c(1,5),lwd=1)
lines(c(1,5),c(1,6.5),lwd=1)
lines(c(1,8.5),c(1,5),lwd=1)
lines(c(1,3.5),c(1,3.5),lwd=1)
lines(c(1,7.5),c(1,1.5),lwd=1)
lines(c(8.5,2.5),c(2.5,5),lwd=1)
lines(c(8.5,5),c(2.5,6.5),lwd=1)
lines(c(8.5,8.5),c(2.5,5),lwd=1)
lines(c(8.5,3.5),c(2.5,3.5),lwd=1)
lines(c(8.5,7.5),c(2.5,1.5),lwd=1)
lines(c(2.5,5),c(5,6.5),lwd=1)
lines(c(2.5,8.5),c(5,5),lwd=1)
lines(c(2.5,3.5),c(5,3.5),lwd=1)
lines(c(2.5,7.5),c(5,1.5),lwd=1)
lines(c(5,8.5),c(6.5,5),lwd=1)
lines(c(5,3.5),c(6.5,3.5),lwd=1)
lines(c(5,7.5),c(6.5,1.5),lwd=1)
lines(c(8.5,3.5),c(5,3.5),lwd=1)
lines(c(8.5,7.5),c(5,1.5),lwd=1)
lines(c(3.5,7.5),c(3.5,1.5),lwd=1)
sch1 <- intersection(c(3.5,3.5),c(8.5,5),c(2.5,5),c(8.5,2.5))
sch2 <- intersection(c(2.5,5),c(8.5,2.5),c(1,1),c(8.5,5))
sch3 <- intersection(c(1,1),c(8.5,5),c(2.5,5),c(7.5,1.5))
sch4 <- intersection(c(3.5,3.5),c(8.5,5),c(2.5,5),c(7.5,1.5))
polygon(c(sch1[1],sch2[1],sch3[1],sch4[1]),c(sch1[2],sch2[2],sch3[2],sch4[2]),col="red",lwd=1.5)
title(main = "Oja median in the bivariate case for n=7",sub="red area corresponds to the Oja median")
