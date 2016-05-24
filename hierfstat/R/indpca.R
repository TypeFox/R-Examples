
indpca<-function(dat,ind.labels=NULL){
#requires ade4
#given a genotype data set dat
#and individual labels lab
#performs a PCA on individuals
  if (is.genind(dat)) dat<-genind2hierfstat(dat)
  
indp<-pop.freq(cbind(1:dim(dat)[1],dat[,-1]))
mati<-NULL
for (i in 1:length(indp)) mati<-rbind(mati,indp[[i]])
mati<-t(mati) #ind rows, nbal col
mp<-apply(mati,2,mean,na.rm=T) 
matic<-sweep(mati,2,FUN="-",mp)
matic[is.na(matic)]<-0.0 #replace NA with 0
pca.matic<-dudi.pca(matic,scannf=FALSE,nf=min(min(dim(matic)),50))
if (is.null(ind.labels)) pca.matic$rownames<-as.character(dat[,1]) 
else pca.matic$rownames<-ind.labels
res<-list(call=match.call(),ipca=pca.matic,ifreq=mati)
class(res)<-"indpca"
res
}

print.indpca<-function(x,...){
print(names(x))
print("$call:")
print(x$call)
print("Dimension of individidual frequency matrix is: ")
print(dim(x$ifreq))
print("$ipca:")
print(x$ipca)
invisible(x)
}

plot.indpca<-function(x,eigen=FALSE,ax1=1,ax2=2,...){
if(eigen){
par(mfrow=c(2,1))
plot(x$ipca$eig/sum(x$ipca$eig),type="h",xlab="axes",ylab="Prop. variance")
}
plot(x$ipca$li[,ax1],x$ipca$li[,ax2],xlab=paste("Axis: ",ax1,sep=""),ylab=paste("Axis: ",ax2,sep=""),type="n")
text(x$ipca$li[,ax1],x$ipca$li[,ax2],labels=x$ipca$rownames,...)
}
