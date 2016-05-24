perspec.dyna<-function(x,y,z,col="black",phi=10,theta=0)
{
persp(x=x,y=y,z=z,col=col,
xlab="level",ylab="h",zlab="",ticktype="detailed",
phi=phi,theta=theta)

loc<-locator(1)
ycor<-loc$y 

alaraja<--0.4
while (loc$y>=alaraja){

     if (loc$x>=0) theta<-theta+10 else theta<-theta-10
     if (loc$y>=0) phi<-phi+10 else phi<-phi-10

     persp(x=x,y=y,z=z,col=col,
     xlab="level",ylab="h",zlab="",ticktype="detailed",
     phi=phi,theta=theta)

     loc<-locator(1)
}
dev.off()
}

