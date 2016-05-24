blca.vb <-
function(X,G, alpha=1, beta=1, delta=1, start.vals= c("single","across"), counts.n=NULL, iter=500, restarts=1, verbose=TRUE, conv=1e-6, small=1e-100)
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
	if(2^M <= (M+1)*G) warning(paste("Model may be improperly specified. Maximum number of classes that should be fitted is", floor(2^M/(M+1)), "."))
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

	
	conv<-N*M*conv
	eps<-N*M

	counter<-0
	llstore<-0
	llcheck<- -Inf

	if(is.numeric(restarts)){
	  if(length(restarts)>1){
	    restarts<- restarts[1]
	    warning("restarts improperly specified - first value will be used, other values will be ignored")
	    } #else {stop("restarts improperly specified. Must be an integer of length 1.")}
	} else {stop("restarts improperly specified. Must be an integer of length 1.")}

	multistart.lp.store<- rep(0, restarts)

	for(r in 1:restarts){
	
	#Set Parameters
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
		gamma<-colSums(Z*counts.n) + delta
		
		zeta1<- (t(Z)%*%(X*counts.n) + alpha)
		zeta2<- (t(Z)%*%((1-X)*counts.n) + beta)
				
		Elog.gamma<-digamma(gamma) - digamma(sum(gamma))
		
		Elog.zeta1<-digamma(zeta1)-digamma(zeta1 + zeta2)
		
		Elog.zeta2<-digamma(zeta2)-digamma(zeta1 + zeta2)
		
		Z<- t(exp(Elog.gamma+ t(X%*%t(Elog.zeta1) + (1-X)%*%t(Elog.zeta2))))
		
		Z<-Z+ small
		
		Z<- Z/rowSums(Z)
		
		l1<- sum(Z*(counts.n*(X%*%t(Elog.zeta1)) + counts.n*((1-X)%*%t(Elog.zeta2))))
		
		l2<- sum(t(counts.n*Z)*(Elog.gamma))
		
		l3<- sum((delta - 1)*Elog.gamma)+lgamma(sum(delta))-sum(lgamma(delta))
		
		l4<- sum((alpha-1)*Elog.zeta1 +(beta-1)*Elog.zeta2)+sum(lgamma(alpha+beta))-sum(lgamma(alpha+beta))
		
		l5<- sum(zeta1*Elog.zeta1 + zeta2*Elog.zeta2)+sum(lgamma(zeta1+zeta2+2))-sum(lgamma(zeta1+1),lgamma(zeta2+1))
	
		l6<- sum(gamma*Elog.gamma)+lgamma(sum(gamma+1))-sum(lgamma(gamma+1))
		
		l7<- sum(xlogy(counts.n*Z,counts.n*Z))
		
		lnew<- l1+l2+l3+l4-l5-l6-l7
		
		llstore[counter]<-lnew
		ll<-lnew

		if(counter>2)
		{
			ll.inf<-(llstore[counter-1] - llstore[counter-2])/(llstore[counter] - llstore[counter-1])*(llstore[counter] - llstore[counter-2]) + llstore[counter-1]
			
			if(llstore[counter] == llstore[counter-1]){ ll.inf<- llstore[counter]}
			
			eps<- ll.inf - llstore[counter]
			}
		
		counter<-counter+1
	
		if(counter==iter) {print("Maximum iteration reached."); break}

		}#while eps>conv

		if(ll > llcheck){
		if(r>1 & verbose==TRUE)  cat("New maximum found... ")
		  rstore<- list(gamma=gamma, zeta1=zeta1, zeta2=zeta2, Z=Z, l=ll, llstore=llstore, counter=counter, eps=eps)
		  
		  llcheck<- ll
		  }
		
		if(verbose==TRUE){ cat(paste("Restart number ", r, ", logpost = ", round(ll, 2), "... \n", sep=""))}
		multistart.lp.store[r]<- ll
		
		eps<-N*M  ## Very important to reset these!!!
		counter<-0
		llstore<-0
		if(r==1 & (is.matrix(start.vals) | is.numeric(start.vals)))  start.vals<- "single"
		}#r
	ll<- llcheck
	o<- order(rstore$gamma, decreasing=TRUE)
	gamma<- rstore$gamma[o]
	zeta1<- rstore$zeta1[o,]
	zeta2<- rstore$zeta2[o,]
	Z<- rstore$Z[,o]
	ll<- rstore$l
	llstore<- rstore$llstore
	counter<- rstore$counter
	eps<- rstore$eps
 
	x<-NULL
	x$call<- match.call()
	x$itemprob<- (zeta1 - 1)/(zeta1+zeta2 - 2)
	x$itemprob[x$itemprob<0]<- 0
	x$classprob<- (gamma-1)/(sum(gamma)-G)
	
	if(any(x$classprob<0)){
	x$classprob[x$classprob<0]<- 0
	x$classprob<- x$classprob/sum(x$classprob)
	}
	
	x$itemprob.sd<- x$itemprob.se<- sqrt((zeta1*zeta2)/((zeta1+zeta2)^2*(zeta1+zeta2+1)) )
	x$classprob.sd<- x$classprob.se<- sqrt( (gamma*(sum(gamma)-gamma))/(sum(gamma)^2*(sum(gamma)+1)))
	
	x$parameters$itemprob<-array(cbind(zeta1, zeta2), dim=c(G,M,2))
	x$parameters$classprob<-gamma
	
	x$Z<- Z
	rownames(x$Z)<- names(counts.n)
	colnames(x$Z)<- paste("Group", 1:G)
	x$LB<-ll 
	if(!is.null(colnames(X))){
		colnames(x$itemprob) <- colnames(x$itemprob.se) <-  colnames(X)
		dimnames(x$parameters$itemprob)<- list(NULL, colnames(X), NULL) 
	}
		
	x$lbstore<-llstore
	x$iter<-counter
	x$eps<-eps
	x$counts<- counts.n
	x$prior<-NULL
	x$prior$alpha<- alpha[o,]
	x$prior$beta<- beta[o, ]
	x$prior$delta<- delta[o]
	class(x)<-c("blca.vb", "blca")
	x
		
	}
