blca.em <-
function(X,G,alpha=1, beta=1, delta=1, start.vals= c("single","across"), counts.n=NULL, iter=500, restarts=5, verbose=TRUE, sd=FALSE, se=sd, conv=1e-6, small=1e-100)
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

	N<-nrow(X); M<-ncol(X) 
	
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

	if(2^M <= (M+1)*G) warning(paste("Model may be improperly specified. Maximum number of classes that should be run is", floor(2^M/(M+1)), "."))
	
	N1<- sum(counts.n)
	
	conv<-N1*M*conv
	eps<-N1*M
		
	counter<-0
	llstore<-0
	llcheck<- -Inf

	if(is.numeric(restarts)){
	  if(length(restarts)>1){
	    restarts<- restarts[1]
	    warning("restarts improperly specified - first value will be used, other values will be ignored")
	    }# else {stop("restarts improperly specified. Must be an integer of length 1.")}
	} else {stop("restarts improperly specified. Must be an integer of length 1.")}

	multistart.lp.store<- rep(0, restarts)
	if(sd!=se) se<- sd
	#Set Parameters
	for(r in 1:restarts){

	if(is.character(start.vals)){
	  
	  if(start.vals[1]=="single"){
	    Z<-unMAP(sample(1:G,size=N,replace=TRUE))
	    if(ncol(Z)<G) Z<-cbind(Z, matrix(0,nrow=N, ncol=(G-ncol(Z))))
	  }else{
	    if(start.vals[1]=="across"){
	      Z<- matrix(runif(N*G), N,G)
	      Z<- Z/rowSums(Z)
	      } else stop("start.vals improperly specified. See help files for more details.")
	    }
	 } else{
	      if(is.matrix(start.vals) & all(dim(as.matrix(start.vals)) == c(N,G)))  Z<- start.vals else{
		  if(is.numeric(start.vals) & length(as.numeric(start.vals))==N){ 
		    Z<- unMAP(start.vals)
		    } else stop("start.vals improperly specified. See help files for more details.")
		 }
	  }	

	while(abs(eps)>conv || counter< 20)
	{	
		Z.sum<-colSums(Z*counts.n)
		
		#M-step
		Taut<-(Z.sum+delta-1)/(N1 + sum(delta)-G)
		Thetat<-(t(Z)%*%(X*counts.n) + alpha-1)/(Z.sum + alpha + beta -2 + small)
		
		if(any(Taut<0)){Taut[Taut<0]<-0; Taut<-Taut/sum(Taut)}
		if(any(Thetat<0)) Thetat[Thetat<0]<-0
		if(any(Thetat>1)) Thetat[Thetat>1]<- 1 ##This should only be due to precision error
				
		#E-Step
		dum<-array(apply(Thetat,1,dbinom, size=1, x=t(X)), dim=c(M,N,G))
		Z1<-t(Taut*t(apply(dum, c(2,3), prod))) + small
		Z<-Z1/rowSums(Z1)

		#Log-Posterior
		l<-sum(log(rowSums(Z1))*counts.n)+sum(xlogy(alpha-1,Thetat)+xlogy(beta-1,1-Thetat))+sum(xlogy(delta-1, Taut)) + lgamma(sum(delta)) - sum(lgamma(delta)) + sum(lgamma(alpha + beta)) - sum(lgamma(alpha) + lgamma(beta))

		llstore[counter]<-l
		if(counter>2)
		{
			ll.inf<-(llstore[counter-1] - llstore[counter-2])/(llstore[counter] - llstore[counter-1])*(llstore[counter] - llstore[counter-2]) + llstore[counter-1]
		
			if(llstore[counter] == llstore[counter-1]){ ll.inf<- llstore[counter]}
		
			eps<- ll.inf - llstore[counter]
		}			
				
		counter<-counter+1		
		if(counter>iter) break 
		}#while

		if(l > llcheck){
		if(r>1 & verbose==TRUE)  cat("New maximum found... ")
		  rstore<- list(Thetat=Thetat, Taut=Taut,Z=Z, l=l, llstore=llstore, counter=counter, eps=eps)
		  
		  llcheck<- l
		  }
		
		if(verbose==TRUE){ cat(paste("Restart number ", r, ", logpost = ", round(l, 2), "... \n", sep=""))}
		multistart.lp.store[r]<- l
		
		eps<-N1*M  ## Very important to reset these!!!
		counter<-0
		llstore<-0
		if(r==1 & (is.matrix(start.vals) | is.numeric(start.vals)))  start.vals<- "single"
		}#r
		l<- llcheck 
		o<- order(rstore$Taut, decreasing=TRUE)

		x<-NULL
		x$call<- match.call()
		#Z.return<-matrix(0, N1, G)
		#for(g in 1:G) Z.return[,g]<-rep(Z[,g], counts.n)
 		if(G>1) x$itemprob<-rstore$Thetat[o, ] else x$itemprob<-rstore$Thetat
 		if(!is.null(colnames(X)))if (G>1) colnames(x$itemprob)<- colnames(X)# else names(as.numeric(x$itemprob))<- colnames(X)
		
		x$classprob<- rstore$Taut[o]
		x$Z<- rstore$Z[,o];
 		if(G>1) rownames(x$Z)<- names(counts.n) else names(Z)<- names(counts.n)
 		if(G>1) colnames(x$Z)<- paste("Group", 1:G)

		x$logpost<- l
  
		likl<- l - sum(xlogy(alpha-1,Thetat)+xlogy(beta-1,1-Thetat))+sum(xlogy(delta-1, Taut)) + lgamma(sum(delta)) - sum(lgamma(delta)) + sum(lgamma(alpha + beta)) - sum(lgamma(alpha) + lgamma(beta))
		x$BIC<- 2*likl-(G*M + G-1)*log(N1)
		x$AIC<- 2*likl - 2*(G*M + G-1)
		x$iter<- length(rstore$llstore)
		x$poststore<- rstore$llstore
		x$eps<- rstore$eps
		x$counts<- counts.n
		x$lpstarts<- multistart.lp.store
		
		x$prior<-NULL
		x$prior$alpha<- alpha[o,]
		x$prior$beta<- beta[o,]
		x$prior$delta<- delta[o]

		if(G>1){
		if(sd){
#			if(any(x$itemprob==0)){ warning("some item probability estimates are exactly zero. standard errors in this case are undefined.")}
#			if(any(x$classprob==0)){ warning("some class probability estimates are exactly zero. standard errors in this case are undefined.")}
			s.e.<- blca.em.sd(x,X,counts.n)
			x$itemprob.sd<- x$itemprob.se<- s.e.$itemprob
			x$classprob.sd<- x$classprob.se<- s.e.$classprob
			convergence<- s.e.$convergence
		} else convergence<- 0
		
		if(counter>iter){ 
		  convergence<- 3 
		  warning("maximum iteration reached - algorithm not deemed to have converged. rerunning the function with 'iter' set to a higher value is recommended.")
		} #else{ convergence<- s.e.$convergence}
		if(convergence==2) {warning("some point estimates likely converged at saddle-point. at least some points will not be at a local maximum. \n rerunning the function with a larger number of restarts is recommended.")}
		if(convergence==4){ warning("some point estimates located at boundary (i.e., are 1 or 0). posterior standard deviations will be 0 for these values.")}
		x$convergence<- convergence
		} else{ 
		  x$convergence<- 1
		  x$classprob.sd<- x$classprob.se<- 0
		  x$itemprob.sd<- x$itemprob.se<- sqrt( ((Thetat*N1 + alpha)*( (1 - Thetat)*N1 + beta))/( (N1 + alpha + beta + small)^2 * (N1 + alpha + beta + 1 + small) ) )
		}
		x$small<- small
		if((se==TRUE)&&(is.null(s.e.$classprob))) se<- FALSE
		x$sd<- x$se<- se
		class(x)<-c("blca.em", "blca")

		x
		}
