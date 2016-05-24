plottwin<-function(tt,et,lev,bary,orde="furthest",ordmet="etaisrec")
{

#if (is.null(et$low)){
   d<-length(et$N)
   step<-matrix(0,d,1)
   for (i in 1:d) step[i]=(et$support[2*i]-et$support[2*i-1])/et$N[i];
   et$step<-step
   et$low<-et$down
   et$upp<-et$high
#}

pp<-plotprof(tt,plot=FALSE,data=TRUE)
vecs<-pp$vecs

d<-length(et$step)

# order the atoms for the level set with level "lev"

lenni<-length(et$value)
distat<-matrix(0,lenni,1)
infopointer<-matrix(0,lenni,1)
lkm<-0
for (i in 1:lenni){
  if (et$value[i]>=lev){
     lkm<-lkm+1
     nod<-i  #nod<-et$nodefinder[i]
     if (ordmet=="etaisrec"){
         recci<-matrix(0,2*d,1)
         for (jj in 1:d){
            recci[2*jj-1]<-et$support[2*jj-1]+et$step[jj]*et$low[nod,jj]
            recci[2*jj]<-et$support[2*jj-1]+et$step[jj]*et$upp[nod,jj]
         }
         distat[lkm]<-etaisrec(bary,recci)
     }
     else{
         lowi<-matrix(0,d,1)
         uppi<-matrix(0,d,1)
         for (jj in 1:d){
            lowi[jj]<-et$support[2*jj-1]+et$step[jj]*et$low[nod,jj]
            uppi[jj]<-et$support[2*jj-1]+et$step[jj]*et$upp[nod,jj]
         }
         baryc<-lowi+(uppi-lowi)/2
         distat[lkm]<-etais(baryc,bary)  #etais(baryc[lk m,],baryind)
     }
     infopointer[lkm]<-i
  }
}
distat<-distat[1:lkm]
infopointer<-infopointer[1:lkm]   #pointe->et$value,et$nodefinder

ord<-order(distat)
infopointer<-infopointer[ord]

xmin<-et$support[1]
xmax<-et$support[2]
ymin<-et$support[3]
ymax<-et$support[4]
plot(x=bary[1],y=bary[2],xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
pch=20,col="red")

i<-1
while (i<=lkm){

     if (orde=="furthest") node<-lkm-i+1 else node<-i
     ip<-infopointer[node]   #ip<-et$nodefinder[infopointer[node]]

     x1<-et$support[1]+et$step[1]*et$low[ip,1]
     x2<-et$support[1]+et$step[1]*et$upp[ip,1] 
     y1<-et$support[3]+et$step[2]*et$low[ip,2]
     y2<-et$support[3]+et$step[2]*et$upp[ip,2] 
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col="lightblue")

     i<-i+1
}

xmin2<-min(vecs[,1])
xmax2<-max(vecs[,3])
ymin2<-0
ymax2<-omamax(vecs[,2])
dev.new()
plot("","",xlab="",ylab="",xlim=c(xmin2,xmax2),ylim=c(ymin2,ymax2))

ycor<-ymax
i<-1
while ((i<=lkm) && (ycor>ymin2)){

     if (orde=="furthest") node<-lkm-i+1 else node<-i
     ip<-infopointer[node]  #ip<-et$nodefinder[infopointer[node]]

     x1<-et$support[1]+et$step[1]*et$low[ip,1]
     x2<-et$support[1]+et$step[1]*et$upp[ip,1] 
     y1<-et$support[3]+et$step[2]*et$low[ip,2]
     y2<-et$support[3]+et$step[2]*et$upp[ip,2] 
     dev.set(which = dev.next())
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col="blue")
     points(x=bary[1],y=bary[2],pch=20,col="red")

     ttnode<-node     
     vecci<-vecs[ttnode,]
     x0<-vecci[1]
     y0<-vecci[2]
     x1<-vecci[3]
     y1<-vecci[4] 
     dev.set(which = dev.next())
     segments(x0, y0, x1, y1)

     loc<-locator(1)
     ycor<-loc$y 

     i<-i+1
}



}

