histo2data<-function(pcf){

d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d) step[i]=(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i];
xmin<-pcf$support[1]
xmax<-pcf$support[2]
ymin<-pcf$support[3]
ymax<-pcf$support[4]
zmin<-pcf$support[5]
zmax<-pcf$support[6]

nnew<-length(pcf$value)
desdat<-matrix(0,nnew,3)
for (i in 1:nnew){
     x1<-pcf$support[1]+step[1]*pcf$down[i,1]
     x2<-pcf$support[1]+step[1]*pcf$high[i,1] 
     y1<-pcf$support[3]+step[2]*pcf$down[i,2]
     y2<-pcf$support[3]+step[2]*pcf$high[i,2] 
     z1<-pcf$support[5]+step[3]*pcf$down[i,3]
     z2<-pcf$support[5]+step[3]*pcf$high[i,3] 
     desdat[i,]<-c((x1+x2)/2,(y1+y2)/2,(z1+z2)/2)
}

f0<-sqrt(pcf$value)
colo<-1-(f0-min(f0)+0.5)/(max(f0)-min(f0)+0.5)
col<-gray(colo)

return(list(dendat=desdat,col=col))
}


