leafsfirst.visu<-function(tt,pcf,lev=NULL,refe=NULL,type="lst",
levmet="radius",ordmet="etaisrec",
lkmbound=NULL,radius=NULL,
orde="furthest",suppo=T,propor=NULL,lty=NULL,numbers=TRUE,
sigcol="lightblue",cex.axis=1,cex=1)
{

if ((!is.null(lev)) || (!is.null(propor))){
    type<-"shape"
    if (is.null(refe)) refe<-locofmax(pcf)
    if (!is.null(propor)) lev<-propor*max(pcf$value)
}
if (is.null(refe)) refe<-locofmax(pcf)

pp<-plotprof(tt,plot=FALSE,data=TRUE)
vecs<-pp$vecs

d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]

# order the atoms for the level set with level "lev"

lenni<-length(pcf$value)
distat<-matrix(0,lenni,1)
infopointer<-matrix(0,lenni,1)

if (type=="lst"){
  lkm<-lenni
  distat<-pcf$value
  infopointer<-seq(1,lkm)
}
else{

lkm<-0
for (i in 1:lenni){
  if (pcf$value[i]>=lev){
     lkm<-lkm+1
     nod<-i  #nod<-pcf$nodefinder[i]
     if (ordmet=="etaisrec"){
         recci<-matrix(0,2*d,1)
         for (jj in 1:d){
            recci[2*jj-1]<-pcf$support[2*jj-1]+step[jj]*pcf$down[nod,jj]
            recci[2*jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
         }
         distat[lkm]<-etaisrec(refe,recci)
     }
     else{
         lowi<-matrix(0,d,1)
         uppi<-matrix(0,d,1)
         for (jj in 1:d){
            lowi[jj]<-pcf$support[2*jj-1]+step[jj]*pcf$down[nod,jj]
            uppi[jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
         }
         baryc<-lowi+(uppi-lowi)/2
         distat[lkm]<-etais(baryc,refe)  #etais(baryc[lk m,],baryind)
     }
     infopointer[lkm]<-i
  }
}

}  #else

distat<-distat[1:lkm]
infopointer<-infopointer[1:lkm]   #pointe->pcf$value,pcf$nodefinder

ord<-order(distat)
infopointer<-infopointer[ord]

if (suppo){
  xmin<-pcf$support[1]
  xmax<-pcf$support[2]
  ymin<-pcf$support[3]
  ymax<-pcf$support[4]
}
else{
  xmin<-tt$refe[1]-tt$maxdis  #pcf$support[1]
  xmax<-tt$refe[1]+tt$maxdis  #pcf$support[2]
  ymin<-tt$refe[1]-tt$maxdis  #pcf$support[3]
  ymax<-tt$refe[2]+tt$maxdis  #pcf$support[4]
}

plot(x=refe[1],y=refe[2],xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
pch=20,cex.axis=cex.axis) #,col="red")

i<-1
while (i<=lkm){

     if (orde=="furthest") node<-lkm-i+1 else node<-i
     ip<-infopointer[node]   #ip<-pcf$nodefinder[infopointer[node]]

     x1<-pcf$support[1]+step[1]*pcf$down[ip,1]
     x2<-pcf$support[1]+step[1]*pcf$high[ip,1] 
     y1<-pcf$support[3]+step[2]*pcf$down[ip,2]
     y2<-pcf$support[3]+step[2]*pcf$high[ip,2] 
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col="gray",lty=lty)

     i<-i+1
}

if (!is.null(lkmbound)){
  i<-1
  while (i<=lkmbound){

     if (orde=="furthest") node<-lkm-i+1 else node<-i
     ip<-infopointer[node]  #ip<-pcf$nodefinder[infopointer[node]]

     x1<-pcf$support[1]+step[1]*pcf$down[ip,1]
     x2<-pcf$support[1]+step[1]*pcf$high[ip,1] 
     y1<-pcf$support[3]+step[2]*pcf$down[ip,2]
     y2<-pcf$support[3]+step[2]*pcf$high[ip,2] 
     dev.set(which = dev.next())
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col=sigcol,lty=lty)
     #points(x=refe[1],y=refe[2],pch=20,col="red")
     if (numbers) text(x=x1+(x2-x1)/2,y=y1+(y2-y1)/2,paste(i),cex=cex)

     i<-i+1
  }
}
else{
  i<-1
  radu<-tt$level[lkm]  #tt$madxdis
  while (radu>=radius){

     if (orde=="furthest") node<-lkm-i+1 else node<-i
     ip<-infopointer[node]  #ip<-pcf$nodefinder[infopointer[node]]

     x1<-pcf$support[1]+step[1]*pcf$down[ip,1]
     x2<-pcf$support[1]+step[1]*pcf$high[ip,1] 
     y1<-pcf$support[3]+step[2]*pcf$down[ip,2]
     y2<-pcf$support[3]+step[2]*pcf$high[ip,2] 
     dev.set(which = dev.next())
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col="blue",lty=lty)
     points(x=refe[1],y=refe[2],pch=20,col="red")

     i<-i+1
     radu<-tt$level[node]
  }
}

}

