graph.matrix<-function(dendat,type="level",
tt=NULL,permu=seq(1:dim(dendat)[1]),col=seq(1:2000),
config="new",shift=0.1,segme=TRUE,poin=FALSE,epsi=0,ystart=0.5, 
pch=21, cex=1, cex.axis=1, yaxt="s",
# profile:
ylen=100,profcol=rep("black",n),texto=TRUE)
{
if (type=="level") 
 graph.matrix.level(dendat, tt=tt, permu=permu, col=col,
 config=config, shift=shift, segme=segme, poin=poin, epsi=epsi, ystart=ystart,
 pch=pch, cex=cex, yaxt=yaxt, cex.axis=cex.axis, texto=texto)

else{  # type="profile"
 n<-dim(dendat)[1]
 d<-dim(dendat)[2]
 x<-seq(1:n)
 y<-seq(1:ylen)
 z<-matrix(0,length(x),length(y))
 varit<-matrix("",length(x),length(y))
 ala<-matrix(0,d,1)
 for (i in 1:d) ala[i]<-min(dendat[,i])
 yla<-matrix(0,d,1)
 for (i in 1:d) yla[i]<-max(dendat[,i])
 range<-yla-ala
 alaind<-matrix(0,d,1)
 alaind[1]<-1
 for (i in 2:d) 
     alaind[i]<-min(alaind[i-1]+round(ylen*range[i]/sum(range))+1,ylen)
 ylaind<-matrix(0,d,1)
 ylaind[d]<-ylen
 for (i in 1:(d-1)) ylaind[i]<-alaind[i+1]+1
 plot(x="",y="",xlim=c(0,n),ylim=c(0,ylen),xlab="",ylab="",yaxt="n",xaxt="n")
 for (i in 1:n){
     for (j in 1:d){
         suht<-(dendat[i,j]-ala[j])/range[j]
         korkeus<-round(suht*(ylaind[j]-alaind[j]))
         alku<-alaind[j]
         loppu<-min(max(alku+korkeus,1),ylen)
         if (alku>=loppu) loppu<-loppu+1
         polygon(x=c(i-1,i-1,i,i),y=c(alku,loppu,loppu,alku),col=profcol[i],
                 lty="blank")
         #z[i,alku:loppu]<-1
     }
 }
 #image(x,y,z,col=c("white","black"),xlab="",ylab="",xaxt="n",yaxt="n")
}

}




