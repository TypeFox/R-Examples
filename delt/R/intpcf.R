intpcf<-function(pcf)
{
value<-pcf$value
down<-pcf$down
high<-pcf$high
support<-pcf$support
N<-pcf$N

d<-length(N)  #dim(down)[2]
step<-stepcalc(support,N)
recnum<-length(value)
int<-0
rr<-1
while (rr<=recnum){
     recint<-matrix(0,2*d,1)
     for (dd in 1:d){
          recint[2*dd-1]<-down[rr,dd]
          recint[2*dd]<-high[rr,dd]
     }
     volint<-massone(recint)
     vol<-prod(step)*volint
     int<-int+value[rr]*vol
     rr<-rr+1
}

return(int)
}




