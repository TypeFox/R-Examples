require(aws)
if(exists("X11")) X11(,12,6)
y <- rnorm(10000)

y1 <- rep(0,10000)
y1[999:1001]<-1
y1[1998:2002]<-sqrt(.5)
y1[2996:3004]<-sqrt(.25)
y1[3992:4008]<-sqrt(.125)
y1[4984:5016]<-sqrt(.0625)
y1[5968:6032]<-sqrt(.03125)
y1[6936:7064]<-sqrt(.03125/2)
y1[7872:8128]<-sqrt(.03125/4)
y1[8744:9256]<-sqrt(.03125/8)

delta <- as.numeric(readline("Width of central interval (default 0)"))
if(is.na(delta)||delta < 0) delta <- 0
delta <- delta/2

size <- as.numeric(readline("Signal size (default 5+intervalwidth/2)"))
if(is.na(size)||size <= 0) size <- 5+delta

plot(y+size*y1,col=3)
lines(size*y1,col=1)
lines(c(1,100000),c(delta,delta),col=2)
lines(c(1,100000),-c(delta,delta),col=2)
title("Data, true function (black) and interval range (red)" )
readline("Enter for next step")

z<-aws.segment(y+size*y1,level=0,delta=delta,hmax=1000,graph=TRUE)

ergs<-c(sum(z@segment>0&size*y1>delta),sum(z@segment>0&size*y1<=delta),sum(z@segment<=0&size*y1>delta),sum(z@segment<=0&size*y1<=delta))
names(ergs) <-c("Correct alternative","False positive","False negative","Correct hypothesis")
print(ergs)

readline("Enter for next step")

plot(size*y1>delta,col=1,type="l")
lines(z@segment,col=2)
title("True segments (black) and segmentation results (red)" )

if(! readline("keep files and device (N/Y) :") %in% c("y","Y")){ 
rm(y,y1,size,delta,z,ergs)
dev.off()
}
