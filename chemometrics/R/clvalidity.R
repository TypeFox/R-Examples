clvalidity <-
function(x,clnumb=c(2:10)){
#
# Compute and plot cluster validity measure
#
#require(e1071)
#require(mclust)
vali=matrix(NA,nrow=3,ncol=length(clnumb))
for (j in 1:length(clnumb)){
	k=clnumb[j]
	#kmeans:
	set.seed(100)
	res=kmeans(x,k)
	wss=rep(NA,k)
	for (i in 1:k){
	  wss[i]=sum(apply(t(t(x[res$clu==i,])-res$cent[i,])^2,1,sum))
	}
	bss=dist(res$cent)
	vali[1,j]=sum(wss)/sum(bss)

        #fuzzy:
        set.seed(100)
        res=cmeans(x,k)
        wss=rep(NA,k)
        for (i in 1:k){
          wss[i]=sum(apply(t(t(x[res$clu==i,])-res$cent[i,])^2,1,sum))
        }
        bss=dist(res$cent)
        vali[2,j]=sum(wss)/sum(bss)

        #mclust
        set.seed(100)
        res=Mclust(x,k)
        wss=rep(NA,k)
	centr=matrix(NA,nrow=k,ncol=ncol(x))
        for (i in 1:k){
	  if (sum(res$class==i)==1)
		centr[i,]=x[res$class==i,]
	  else 
		centr[i,]=apply(x[res$class==i,],2,mean)
	}	
        for (i in 1:k){
          wss[i]=sum(apply(t(t(x[res$class==i,])-centr[i,])^2,1,sum))
        }
        bss=dist(centr)
        vali[3,j]=sum(wss)/sum(bss)
}
matplot(t(vali),xlab="Number of clusters",ylab="Cluster validity",cex.lab=1.2,type="b",
   ylim=c(0,max(vali)),xaxt="n",col=1,lty=1,pch=c(1,2,3))
axis(1,at=1:length(clnumb),clnumb)
legend("topright",legend=c("k-means","fuzzy","model-based"),pch=c(1,2,3))
dimnames(vali)=list(c("kmeans","fuzzy","mclust"),clnumb)
list(validity=vali)
}

