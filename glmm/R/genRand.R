#genRandOld <-
#function(sigma.star,s.star,z.list,m,distrib,gamm){
#	nrand<-lapply(z.list,ncol)
#	nrandom<-unlist(nrand)
#	q<-sum(nrandom)
#	if(q!=length(s.star)) stop("Number of random effects inconsistent between s.star and mod.mcml$z")
#	
#	eek<-getEk(z.list)
#	Aks<-Map("*",eek,sigma.star)
#	A.star<-addVecs(Aks) #at this point still a vector
#	#A.star<-diag(A.star)
#	u.star<-A.star*s.star #this is the mean vector if normal, location if t
#	
#	u<-matrix(data=NA,nrow=m,ncol=q)
#	if(distrib=="normal"){	
#		for(k in 1:m){
#			u.swoop<-rnorm(q)
#			u[k,]<-u.swoop*A.star+u.star
#		}
#	}

#	if(distrib=="tee"){
#		for(k in 1:m){ 
#		u[k,]<-	rmvt(1,diag(A.star),delta=u.star,df=gamm,type="shifted")
#		}
#	}
#	list(u=u,u.star=u.star,distrib=distrib)
#}


genRand <-function(mew,Sigma,m){
	q<-length(mew)
	if(length(mew)!=nrow(Sigma)) stop("problem in genRand: mew and Sigma dims disagree")
	u<-matrix(data=NA,nrow=m,ncol=q)

	#calculate Sigma^.5
	eigenstuff<-eigen(Sigma,symmetric=TRUE)
	lambda.half<-sqrt(eigenstuff$values)
	Sigma.half<-eigenstuff$vectors%*%diag(lambda.half)%*%t(eigenstuff$vectors)

	for(k in 1:m){
		u.swoop<-rnorm(q)
		u[k,]<-u.swoop%*%Sigma.half+mew
		}
	u
}
