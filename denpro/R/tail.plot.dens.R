tail.plot.dens<-function(denmat,h=1,k=100,b=0.25,alpha=1,type="left.tail",
minx=-0.2,plot=TRUE)
{
# log="y",cex.axis=1,pch=20,pchs=rep(20,1000))
lkm<-dim(denmat)[2]
n<-dim(denmat)[1]

if (type=="left.tail"){

m<-floor(n/2)
detmat<-matrix(0,m,lkm)
for (i in 1:lkm){
       dencur<-denmat[,i]
       ordi<-order(dencur)
       dendat.ord<-dencur[ordi]
       detmat[,i]<-dendat.ord[1:m]
       #split<-median(dencur)
       #redu.ind<-(dencur<split) 
       #dendat.redu<-dencur[redu.ind]
       #ordi<-order(dendat.redu)
       #dendat.ord<-dendat.redu[ordi]  #nredu<-length(dendat.redu)
       #detmat[,i]<-dendat.ord[1:m]
}
minu<-min(detmat,na.rm=TRUE)
maki<-max(detmat,na.rm=TRUE)

x<-matrix(0,k,1)
pc<-matrix(0,k,m)
for (mm in 1:m){
     datai<-detmat[mm,]
     if (is.null(h)){
        expon<-1/(1+4)
        sdev<-sd(datai,na.rm=TRUE)
        h<-(4/(1+2))^expon*sdev*n^(-expon)
     }  
     ini<-!is.na(datai)
     dataj<-datai[ini]   
     for (kk in 1:k){
       arg<-minx+(maki-minx)*kk/(k+1) 
       x[kk]<-arg
       pc[kk,mm]<-kernesti.dens(arg,dataj,h=h)
   }
}
y<-log(seq(1,m))
pc2<-(pc)^b
colo<-grey(seq(1,0,-0.01),alpha=alpha)

if (plot) image(x,y,pc2,col=colo)  #image(pc2,col=topo.colors(120))
else return(list(x=x,y=y,pc=pc,colo=colo))
}


}


