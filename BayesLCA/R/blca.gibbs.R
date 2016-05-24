blca.gibbs <-
function(X,G, alpha=1, beta=1, delta=1, start.vals= c("prior","single","across"), counts.n=NULL, iter=5000, thin=1, accept=thin, burn.in=100, relabel=TRUE, verbose=TRUE, verbose.update=1000)
{
	if(is.null(counts.n))
	{
		if(inherits(X,"data.blca")){
			counts.n<- X$counts.n
			X<- X$data
		}else{
			Xdat<- data.blca(X)
			X<- Xdat$data
			counts.n<- Xdat$counts.n
			}# else class(X)
		} else{ 
			X<- as.matrix(X)
			if(any(X[X>0]!=1))
			stop("If vector of counts is supplied separately, then data must be binary.")
			}
	
	M<-ncol(X); N<-nrow(X)
	if(2^M <= (M+1)*G) warning(paste("Model may be improperly specified. Maximum number of classes that should be run is", floor(2^M/(M+1)), "."))
	
		if(length(delta)==1) delta<-rep(delta,G)
		if(length(delta)!=G) stop("delta prior is wrong length (i.e., !=1 or G)")
			if(!is.matrix(alpha)){
			if(any(length(alpha)==c(1,G)) ){
				alpha<-matrix(alpha,G,M)
			} else {
				if(length(alpha)==M){
					alpha<- matrix(alpha,G,M, byrow=TRUE)
				} else stop("Item probability prior improperly specified.")
			}
	} #else {
	#	if(!is.matrix(alpha)) stop("Item probability prior improperly specified.")
	#}	
	
	if(!is.matrix(beta)){
		if(any(length(beta)==c(1,G)) ){
			beta<-matrix(beta,G,M)
			} else {
				if(length(beta)==M){
					beta<- matrix(beta,G,M, byrow=TRUE)
				} else stop("Item probability prior improperly specified.")
			}
		} #else {
			#if(!is.matrix(beta)) stop("Item probability prior improperly specified.")
	#}	
	if(any(dim(alpha)!=c(G,M), dim(beta)!=c(G,M))) stop("Item probability prior improperly specified.")


	
	#Prior
	tau<- rep(0,G)
	theta<- matrix(0,G,M)
	if(is.character(start.vals)){
	  if(start.vals[1]=="prior"){
	    tau<-rdirichlet(1,delta)
	    theta<-matrix(rbeta(G*M,alpha,beta),nrow=G,ncol=M)
	    }else{
	  if(start.vals[1]=="single"){
	    Z<-unMAP(sample(1:G,size=N,replace=TRUE))
	    if(ncol(Z)<G) Z<-cbind(Z, matrix(0,nrow=N, ncol=(G-ncol(Z))))
	    tau<-rdirichlet(1,delta+colSums(Z))
	    for(g in 1:G) theta[g,]<-rbeta(M,alpha+colSums(Z[,g]*X), beta+colSums(Z[,g]*(1-X)))
	  }else{
	    if(start.vals[1]=="across"){
		Z<- matrix(runif(N*G), N,G)
		Z<- Z/rowSums(Z)
		tau<-rdirichlet(1,delta+colSums(Z))
		for(g in 1:G) theta[g,]<-rbeta(M,alpha+colSums(Z[,g]*X), beta+colSums(Z[,g]*(1-X)))
		} else stop("start.vals improperly specified. See help files for more details.")
	      }
	     }
	    } else{
	      if(is.matrix(start.vals) & all(dim(as.matrix(start.vals)) == c(N,G))){
		Z<- start.vals 
		tau<-rdirichlet(1,delta+colSums(Z))
		for(g in 1:G) theta[g,]<-rbeta(M,alpha+colSums(Z[,g]*X), beta+colSums(Z[,g]*(1-X)))
	      }else{
		  if(is.numeric(start.vals) & length(as.numeric(start.vals))==N){ 
		    Z<- unMAP(start.vals)
		    tau<-rdirichlet(1,delta+colSums(Z))
		    for(g in 1:G) theta[g,]<-rbeta(M,alpha+colSums(Z[,g]*X), beta+colSums(Z[,g]*(1-X)))
		    } else stop("start.vals improperly specified. See help files for more details.")
		 }
	  }	

	#Store Values
	maxiter<- iter
	if(accept!=thin) thin<- accept
	Kstore<-maxiter*thin
	taustore<-matrix(NA,Kstore,G)
	thetastore<-array(NA,c(Kstore,G,M))
	labelstore<-matrix(NA, Kstore,G)
	logpost.store<- rep(NA, Kstore)
	counter<-1
	#Dummy matrix
	W<-matrix(nrow=N, ncol=G)
	Zstore<-matrix(0, N, G)
	label.swap<- FALSE

	if(verbose==TRUE) cat("Initialising sampler...starting burn-in.\n")
	#Gibbs Sampler
	for (iter in 1:(maxiter+burn.in))
	{
		for(g in 1:G)	W[,g]<-tau[g]*apply(theta[g,]^t(X)*(1-theta[g,])^t(1-X),2,prod) # g
	
		Z<-Zsamp(W, counts.n) 

		if(iter==burn.in+1){
		  if(verbose==TRUE) cat("Burn-in completed...\n")
		  Zmatch<- counts.n*Z
		  }
		
		tau<-rdirichlet(1,delta+colSums(Z))

		for(g in 1:G) theta[g,]<-rbeta(M,alpha+colSums(Z[,g]*X), beta+colSums(Z[,g]*(1-X)))

		if((iter>burn.in)&&(iter%%round(1/thin)==0))
		{
			if(verbose==TRUE & (iter-burn.in)%%verbose.update == 0) cat(iter - burn.in, "of", maxiter, "samples completed...\n")
			if(G > 1) { 
				match1<- matchClasses(t(Z)%*%Zmatch, method="exact", verbose=FALSE) } else match1 <- 1:G
			if(any(match1!=1:G)) label.swap<- TRUE
			labelstore[counter,]<- match1
			if(relabel){
				taustore[counter,labelstore[counter,]]<-tau
				thetastore[counter, labelstore[counter,], ]<-theta
				Zstore<- ((iter-burn.in-1)*Zstore + Z[, labelstore[counter,]])/(iter - burn.in)				
			}else{
				taustore[counter, ]<-tau
				thetastore[counter, , ]<-theta
				Zstore<- ((iter-burn.in-1)*Zstore + Z)/(iter - burn.in)	
			}
			
			logpost.store[counter]<- sum(log(rowSums(W))*counts.n) + sum(xlogy(alpha-1, theta) + xlogy(beta-1, 1-theta) + sum(xlogy(delta-1, tau)))
			
			counter<-counter+1
			if(counter>Kstore) break
			}
		} #iter
	if(verbose==TRUE) cat("Sampling run completed.\n")
	tau<- apply(taustore, 2, mean)
	o<- order(tau, decreasing=TRUE)

	x<-NULL
	x$call<- match.call()
	x$classprob<- tau[o]
	x$itemprob<- apply(thetastore, c(2,3), mean)[o,]
	x$classprob.sd<- x$classprob.se<- apply(taustore, 2, sd)[o]
	x$itemprob.sd<- x$itemprob.se<- apply(thetastore, c(2,3), sd)[o, ]
	
	if(G == 1){
		 x$itemprob <- matrix(x$itemprob, G, M)
		 x$itemprob.sd <- x$itemprob.se <- matrix(x$itemprob.sd, G, M)
		 }
	
	if(G > 1)	dum<-array(apply(x$itemprob,1,dbinom, size=1, x=t(X)), dim=c(M,N,G)) else dum<-array( dbinom(t(X), 1, x$itemprob), dim=c(M,N,G))

	Z1<-t(x$classprob*t(apply(dum, c(2,3), prod)))

	x$logpost<- sum(log(rowSums(Z1))*counts.n) + sum(xlogy(alpha-1, x$itemprob) + xlogy(beta-1, 1- x$itemprob) + sum(xlogy(delta-1, x$classprob)))
	
	if(G > 1) x$Z<- (Zstore/counts.n)[, o] else x$Z<- matrix( (Zstore/counts.n)[, o], nrow = N, ncol = G)
	
	rownames(x$Z)<- names(counts.n)
	colnames(x$Z)<- paste("Group", 1:G)
	
	x$samples<-NULL
	if(G > 1){
	x$samples$classprob<-taustore[, o]
	x$samples$itemprob<-thetastore[, o, ]
	} else {
	x$samples$classprob<- taustore
	x$samples$itemprob<- thetastore
	}
	x$samples$logpost<- logpost.store
	
	if(!is.null(colnames(X))){
		colnames(x$itemprob)<- colnames(x$itemprob.sd)<- colnames(x$itemprob.se)<- colnames(X)
		dimnames(x$samples$itemprob)<- list(NULL, NULL, colnames(X)) 
	}
	
	Dbar<- mean(logpost.store)
	S2<- var(logpost.store)
	
	x$DIC<- 2*(2*Dbar - x$logpost)
	x$BICM<- 2*(Dbar - S2*(log(sum(counts.n))-1))
	x$AICM<- 2*(Dbar - S2)
	
	x$counts<- counts.n
	
	x$prior<-NULL
	x$prior$alpha<- alpha[o, ]
	x$prior$beta<- beta[o, ]
	x$prior$delta<- delta[o]

	x$thin<-thin
	x$burn.in<-burn.in
	x$relabel<- relabel
	x$labelstore<-labelstore

	class(x)<-c("blca.gibbs", "blca")
	#if(relabel && matchClasses(t(Z)%*%Z1, method="exact", verbose=FALSE)) warning("Label-switching (provisionally) corrected for - proceed with caution")
	if(relabel && label.swap) warning("Label-switching (provisionally) corrected for - diagnostic plots are recommended. Use '?plot.blca' for details.")
	if(!relabel && label.swap) warning("Label-switching may have occurred - diagnostic plots are recommended. Use '?plot.blca' for details.")
	
x
}
