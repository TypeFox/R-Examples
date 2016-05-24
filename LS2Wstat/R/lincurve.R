lincurve <-function(x,start=1,end=2,a=0.25){

block<-a
lx<-length(x)

si<-floor(block*lx)

y<-x

y[1:si]<-start
y[(lx-si+1):lx]<-end

lpart<-lx-2*si

y[(si+1):(lx-si)]<-start+(1:lpart)*(end-start)/lpart

return(y)

}

