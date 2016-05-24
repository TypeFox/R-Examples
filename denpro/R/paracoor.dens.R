paracoor.dens<-function(dendat,type="classical",h=1,b=0.25,k=100,m=100,alpha=1)
{
# k<-1000  # grid lkm vaakatasossa
# m<-1000  # grid lkm pystytasossa

n<-dim(dendat)[1]

if (type=="new"){

 vals<-matrix(0,n,1)
 for (i in 1:n){
    arg<-dendat[i,]
    vals[i]<-kernesti.dens(arg,dendat,h=h)
 }
 w<-(vals-min(vals))/(max(vals)-min(vals))
 or<-order(w)
 w2<-(1-w)^b
 paletti<-grey(w2)[or]    
 x<-dendat[or,]
 paracoor(x,paletti=paletti) 

}

if (type=="classical"){

 d<-dim(dendat)[2]
 maks<-matrix(0,d,1)
 mini<-matrix(0,d,1)
 for (i in 1:d){
    maks[i]<-max(dendat[,i])
    mini[i]<-min(dendat[,i])
 }
 dendat2<-dendat
 for (i in 1:d) dendat2[,i]<-(dendat[,i]-mini[i])/(maks[i]-mini[i])
 pc<-matrix(0,m,k*(d-1))
 for (dd in 1:(d-1)){
   for (kk in 1:k){
      x1<-dendat2[,dd]
      x2<-dendat2[,dd+1]
      t<-kk/(k+1)
      datai<-(1-t)*x1+t*x2
      ind<-(dd-1)*k+kk
      for (mm in 1:m){
          arg<-mm/m
          pc[mm,ind]<-kernesti.dens(arg,datai,h=h)
      }
   }
 }
 pc2<-t(pc)^b
 colo<-grey(seq(0,1,0.1),alpha=alpha)
 image(pc2,col=colo)  #image(pc2,col=topo.colors(120))
 #image(pc2,col=terrain.colors(50))
 #heatmap(pc2)
 #contour(pc2)

}

}



