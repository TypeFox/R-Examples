paraclus<-function(dendat,algo="kmeans",k=2,method="complete",
scatter=FALSE,coordi1=1,coordi2=2,levelmethod="center",
startind=c(1:k),range="global",terminal=TRUE,coordi=1,
paletti=NULL,xaxt="s",yaxt="s",cex.axis=1,pch.paletti=NULL)
{
if (is.null(paletti)) paletti<-seq(1,2000)
if (is.null(pch.paletti)) pch.paletti<-rep(21,2000)
if (algo!="kmeans"){ 
      method<-algo
      algo<-"hclust"
}

n<-dim(dendat)[1]
d<-dim(dendat)[2]
colot<-c(colors()[2],colors()[3])

if (algo=="kmeans"){
    starters<-dendat[startind,]
    cl<-kmeans(dendat,k,centers=starters)
    ct<-cl$cluster
    centers<-cl$centers
}
else if (algo=="hclust"){
       dis<-dist(dendat)
       hc <- hclust(dis, method=method)
       ct<-cutree(hc,k=k)
       centers<-matrix(0,k,d)
       for (ij in 1:k) centers[ij,]<-mean(data.frame(dendat[(ct==ij),]))
}

# calculate innerlevel
innerlevel<-matrix(0,n,1)
maxlevel<-matrix(0,k,1)
for (i in 1:n){
  classlabel<-ct[i]
  cente<-centers[classlabel,]
  if (levelmethod=="random"){ 
        set.seed(i)
        luku<-runif(1)
  }
  else{ 
        luku<-sqrt(sum((dendat[i,]-cente)^2))
  }
  innerlevel[i]<-luku  
  maxlevel[classlabel]<-max(maxlevel[classlabel],luku)
}
# calculate classlevel
classlevel<-matrix(0,k,1)
i<-2
while (i<=k){
   if (levelmethod=="random"){ 
        classlevel[i]<-classlevel[i-1]+1
   }
   else{
        classlevel[i]<-classlevel[i-1]+maxlevel[i-1]
   }
   i<-i+1
}
# calculate level
level<-matrix(0,n,1)
for (i in 1:n){
    classlabel<-ct[i]
    level[i]<-innerlevel[i]+classlevel[classlabel]
}

if (d<=5){ 
   lkm<-d  
   times<-0
   reminder<-d
} 
else{
   lkm<-5
   times<-floor(d/lkm)
   reminder<-d-lkm*times
}
curcolo<-1
ymin<-0  #min(level)
ymax<-max(level)

if (!terminal){
       coordinate<-coordi
       x<-dendat[,coordinate]
       if (range=="global"){
          xmin<-min(dendat) 
          xmax<-max(dendat)
       }
       else{
          xmin<-min(x)
          xmax<-max(x)
       }
       plot(x="",y="",xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
            xaxt=xaxt,yaxt=yaxt,cex.axis=cex.axis)
       for (j in 1:k){
         if (curcolo==1) curcolo<-2 else curcolo<-1
         polygon(c(xmin,xmax,xmax,xmin),
                 c(classlevel[j],classlevel[j],
                   classlevel[j]+maxlevel[j],classlevel[j]+maxlevel[j]),
                 col=colot[curcolo]) 
       }
       points(x,level,col=paletti[ct],pch=pch.paletti[ct])
       if (scatter) plot(dendat[,coordi1],dendat[,coordi2], col = paletti[ct],
                         xaxt=xaxt,yaxt=yaxt,pch=pch.paletti[ct])

}
########################################################
else{

t<-1
while (t<=times){
   mat<-matrix(c(1:lkm),1,lkm)
   dev.new()
   layout(mat)
   for (i in 1:lkm){
       coordinate<-(times-1)*lkm+i
       x<-dendat[,coordinate]
       if (range=="global"){
          xmin<-min(dendat) 
          xmax<-max(dendat)
       }
       else{
          xmin<-min(x)
          xmax<-max(x)
       }
       plot(x="",y="",xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
            xaxt=xaxt,yaxt=yaxt,cex.axis=cex.axis)
       for (j in 1:k){
         if (curcolo==1) curcolo<-2 else curcolo<-1
         polygon(c(xmin,xmax,xmax,xmin),
                 c(classlevel[j],classlevel[j],
                   classlevel[j]+maxlevel[j],classlevel[j]+maxlevel[j]),
                 col=colot[curcolo]) 
       }
       points(x,level,col=ct)
   }
   t<-t+1
}
if (reminder>0){
   lkm<-reminder
   mat<-matrix(c(1:lkm),1,lkm)
   dev.new()
   layout(mat)
   for (i in 1:lkm){
       coordinate<-i
       x<-dendat[,coordinate]
       if (range=="global"){
          xmin<-min(dendat) 
          xmax<-max(dendat)
       }
       else{
          xmin<-min(x)
          xmax<-max(x)
       }
       plot(x="",y="",xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
            xaxt=xaxt,yaxt=yaxt)
       for (j in 1:k){
         if (curcolo==1) curcolo<-2 else curcolo<-1
         polygon(c(xmin,xmax,xmax,xmin),
                 c(classlevel[j],classlevel[j],
                   classlevel[j]+maxlevel[j],classlevel[j]+maxlevel[j]),
                 col=colot[curcolo]) 
       }
       points(x,level,col=ct)
   }
}

# scatter plot
if (scatter){
   dev.new()
   plot(dendat[,coordi1],dendat[,coordi2], col = ct, xaxt=xaxt, yaxt=yaxt)
}

} # if terminal

}


