diag2edge<-function(cordseg){

if (is.data.frame(cordseg)) 
    cordseg <- as.matrix(cordseg)
    cordseg <- matrix(cordseg, ncol = 2)
    d <- dim(cordseg)
    if (d[2] != 2 | d[1] != 2) 
        stop("Segment coordinates must be a 2 x 2 matrix")
    cordseg <- cordseg[order(cordseg[, 1]), ]
    if (cordseg[1,1]==cordseg[2,1]) cordseg <- cordseg[rev(order(cordseg[, 2])), ]
    Di<-sqrt((cordseg[1,1]-cordseg[2,1])^2+(cordseg[1,2]-cordseg[2,2])^2)
    d<-Di/sqrt(2)

    flag=TRUE; if (cordseg[1,2]>cordseg[2,2]) flag=FALSE

if (flag){
x3<-cordseg[1,1]+d*cos(pi/4+acos(abs(cordseg[1,1]-cordseg[2,1])/Di))
y3<-cordseg[2,2]-d*cos(pi/4+acos(abs(cordseg[1,2]-cordseg[2,2])/Di))
}
else
{
x3<-cordseg[2,1]-d*cos(pi/4+acos(abs(cordseg[1,1]-cordseg[2,1])/Di))
y3<-cordseg[1,2]-d*cos(pi/4+acos(abs(cordseg[1,2]-cordseg[2,2])/Di))
}
x<-rbind(cordseg[1,],cbind(x3,y3),deparse.level=0)
colnames(x)<-NULL
x
}
