"asyvar" <-
function(xgrid=seq(0,1,length=21),ygrid=seq(0,1,length=21)){

asyvar<-matrix(0,length(xgrid),length(ygrid))

for (i in 1:length(xgrid)){
for (j in 1:length(ygrid)){

p1<-xgrid[i]
p2<-ygrid[j]

q1 <- 1-p1
q2 <- 1-p2

asyvar[i,j] <- (p1*q1 + p2*q2)/((p1+p2)*(2-p1-p2)/2)

}
}

asyvar[1,1]<-0
asyvar[length(xgrid),length(xgrid)]<-0

asyvar

}

