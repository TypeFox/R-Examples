cluster.lst<-function(dendat,h,N=NULL,cut=NULL,lambda=NULL,complete=FALSE,
type="grid",labels="number",nodes=NULL,minobs=1)
{
# cut is in (0,1)
n<-dim(dendat)[1]

if (type=="grid"){ 
   pcf<-pcf.kern(dendat,h,N,radi=0)
   lst<-leafsfirst(pcf)
}
else{
    pcf<-pcf.greedy.kernel(dendat,h,minobs=minobs)
    lst<-leafsfirst.adagrid(pcf)
}
if (!is.null(cut)) clusterlevel<-cut*(max(lst$level)-min(lst$level))
if (!is.null(lambda)) clusterlevel<-lambda
if (!is.null(nodes)) clusterlevel<-NULL

if (type=="grid")
cd<-colors2data(dendat,pcf,lst,clusterlevel=clusterlevel,nodes=nodes,type="regular")
else
cd<-colors2data(dendat,pcf,lst,clusterlevel=clusterlevel,nodes=nodes,type="ada")

cls0<-cd$datacolo
mita0<-unique(cd$datacolo)
ind<-(mita0=="grey")
mita<-mita0
mita[length(mita0)]<-"grey"
mita[ind]<-mita0[length(mita0)]

clnum<-length(mita)
cnum<-clnum-1
nums<-matrix(0,cnum,1)
cls<-matrix(0,n,1)
for (i in 1:n){
    if (cls0[i]=="grey") cls[i]<-0
    else{ 
       for (j in 1:cnum) 
       if (mita[j]==cls0[i]){ 
          cls[i]<-j
          nums[j]<-nums[j]+1 
       }
    }
}

if (complete){
   indi<-(cls==0)
   data<-dendat[indi,]
   n0<-dim(data)[1]
   part<-matrix(0,n0,1)
   part0<-matrix("",n0,1)
   mito<-mita[1:cnum]
   for (i in 1:n0){
       arg<-data[i,]
       arvot<-matrix(0,cnum,1)
       for (j in 1:cnum){
           ota<-(cls==j)
           x<-dendat[ota,]
           arvot[j]<-(nums[j]/(n-n0))*kernesti.dens(arg,x,h=h)
       }
       makind<-(arvot==max(arvot))
       part[i]<-seq(1,cnum)[makind]
       part0[i]<-mito[makind]
   }
   cls[indi]<-part
   cls0[indi]<-part0
}

if (labels=="number") labels<-cls else labels<-cls0
return(labels)
}

