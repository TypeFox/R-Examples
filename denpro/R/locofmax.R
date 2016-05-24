locofmax<-function(pcf)
{
d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d){
   step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]
}

nod<-which.max(pcf$value)

lowi<-matrix(0,d,1)
uppi<-matrix(0,d,1)
for (jj in 1:d){
    lowi[jj]<-pcf$support[2*jj-1]+step[jj]*(pcf$down[nod,jj])
    uppi[jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
}
baryc<-lowi+(uppi-lowi)/2  

return(baryc)
}

