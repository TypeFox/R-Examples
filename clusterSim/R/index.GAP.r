index.Gap<-function(x,clall,reference.distribution="unif",B=10,method="pam",d=NULL,centrotypes="centroids") {
GAP<-function(X,cl,referenceDistribution,B,method,d,centrotypes){
simgap<-function(Xvec) {
	ma<-max(Xvec)
	mi<-min(Xvec)
	Xout<-runif(length(Xvec),min=mi,max=ma)
	return(Xout) 
}



# Thanks to Chuen Seng for this fix
diagvar<-function(X){
X<-as.matrix(X)
D<-matrix(1,nrow=nrow(X),ncol=1)
res0 = X - (D %*% solve(t(D) %*% D) %*% t(D)) %*% X
return(colSums(res0^2)/(nrow(X)-1))
}



pcsim<-function(X,d,centrotypes) {
  if(centrotypes=="centroids"){
    Xmm<-apply(X,2,mean)
  }
  else{
    Xmm<-.medoid(x,d)
  }
	for (k in (1:dim(X)[2])) 
	{
		X[,k]<-X[,k]-Xmm[k] 
	}
	ss<-svd(X)
	Xs<-X%*%ss$v
	Xnew<-apply(Xs,2,simgap)
	Xt<-Xnew%*%t(ss$v)
	for (k in (1:dim(X)[2])) 
	{
		Xt[,k]<-Xt[,k]+Xmm[k] 
	}
	return(Xt) 
}
  if(is.null(dim(x))){
    dim(x)<-c(length(x),1)
  }
	ClassNr<-max(cl)
	#print(ClassNr)
	Wk0<-0
	WkB<-matrix(0,1,B)
	for (bb in (1:B)) {
		if (reference.distribution=="unif")
			Xnew<-apply(X,2,simgap) 
		else if (reference.distribution=="pc") 
			Xnew<-pcsim(X,d,centrotypes) 
		else
			stop("Wrong reference distribution type")	
		if (bb==1) {
			pp<-cl
			if (ClassNr==length(cl))
			pp2<-1:ClassNr
			else if (method=="pam")
			pp2<-pam(Xnew,ClassNr)$cluster 
			else if (method=="k-means")
			pp2<-kmeans(Xnew,ClassNr,100)$cluster 
			else if(method=="diana")
			pp2<-cutree(as.hclust(diana(dist(Xnew))),k=ClassNr)
			else if (method=="single" || method=="complete"
			|| method=="average" || method=="ward"
			|| method=="mcquitty" || method=="median" || method=="centroid")
				pp2<-cutree(hclust(dist(Xnew), method = method),ClassNr)
			else	
				stop ("Wrong clustering method")	
			if (ClassNr>1) {
				for (zz in (1:ClassNr)) {
				Xuse<-X[pp==zz,]
				Wk0<-Wk0+sum(diagvar(Xuse))*(length(pp[pp==zz])-1)/(dim(X)[1]-ClassNr)
				Xuse2<-Xnew[pp2==zz,]			
				WkB[1,bb]<-WkB[1,bb]+sum(diagvar(Xuse2))*(length(pp2[pp2==zz])-1)/(dim(X)[1]-ClassNr)  
				} 
			}
			if (ClassNr==1) {
				Wk0<-sum(diagvar(X))
				WkB[1,bb]<-sum(diagvar(Xnew)) 
			} 
		}
		if (bb>1) { 
			if (ClassNr==length(cl))
			pp2<-1:ClassNr
			else if (method=="pam")
				pp2<-pam(Xnew,ClassNr)$cluster 
			else if (method=="k-means")
				pp2<-kmeans(Xnew,ClassNr,100)$cluster
      else if (method=="diana")
        pp2<-cutree(as.hclust(diana(dist(Xnew))),k=ClassNr)
			else if (method=="single" || method=="complete"
				|| method=="average" || method=="ward"
				|| method=="mcquitty" || method=="median" || method=="centroid")
				pp2<-cutree(hclust(dist(Xnew), method = method),ClassNr)
			else	
				stop ("Wrong clustering method")	
			if (ClassNr>1) {
				for (zz in (1:ClassNr)) {
					Xuse2<-Xnew[pp2==zz,]
					WkB[1,bb]<-WkB[1,bb]+sum(diagvar(Xuse2))*length(pp2[pp2==zz])/(dim(X)[1]-ClassNr) 
				} 
			}
			if (ClassNr==1) {
				WkB[1,bb]<-sum(diagvar(Xnew)) 
			} 
		} 
	}

	Sgap<-mean(log(WkB[1,]))-log(Wk0)
Sdgap<-sqrt(1+1/B)*sqrt(var(log(WkB[1,])))*sqrt((B-1)/B) 
	resul<-list(Sgap=Sgap,Sdgap=Sdgap) 
	resul
}

  if(sum(c("centroids","medoids")==centrotypes)==0)
  stop("Wrong centrotypes argument")
  if("medoids"==centrotypes && is.null(d))
  stop("For argument centrotypes = 'medoids' d cannot be null")
  if(!is.null(d)){
  if(!is.matrix(d)){
    d<-as.matrix(d)
  }
  row.names(d)<-row.names(x)
  }
	#print(x)
	X<-as.matrix(x)
	gap1<-GAP(X,clall[,1],reference.distribution,B,method,d,centrotypes)
	gap<-gap1$Sgap
	gap2<-GAP(X,clall[,2],reference.distribution,B,method,d,centrotypes)
	diffu<-gap-(gap2$Sgap-gap2$Sdgap)
	resul<-list(gap=gap,diffu=diffu)
	resul
}
