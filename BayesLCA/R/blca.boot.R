blca.boot <-
function(X,G=NULL, alpha=1,beta=1, delta=rep(1,G), start.vals= c("single","across"), counts.n=NULL, fit=NULL, iter=50,  B=100, relabel=FALSE, verbose=TRUE, verbose.update=10, small=1e-100)
{			
	if(is.null(fit)){
		if(verbose==TRUE) cat("Object 'fit' not supplied. Obtaining starting values via blca.em...\n")
		if(is.null(G)) warning("Number of groups must be specified.")
		x<-blca.em(X, G, iter=500, conv=1e-10, alpha=alpha, beta=beta, delta=delta, start.vals= start.vals)
		conv<-x$eps
		if(verbose==TRUE) cat("Starting values obtained...\n")
	} else {
		x<- fit
		G<- length(fit$classprob)
		conv<- x$eps
		alpha<- x$prior$alpha
		beta<- x$prior$alpha
		delta<- x$prior$delta
	}
	Tn.Theta<-x$itemprob 
	Tn.Tau<-x$classprob
	Zorig<- x$Z
	
	if(is.null(counts.n))
	{
		if(class(X)=="data.blca"){
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
	
	N1<-sum(counts.n) 
	countsorig.n<-counts.n 
	Zcheck<- sqrt(countsorig.n)*Zorig
	Nstar<-N
	swapcheck<- FALSE

	##Store Values
	Ti.Theta<-array(0,dim=c(B,G,M))
	Ti.Tau<-matrix(0,B,G)
	Z0<- matrix(0, N,G)
	rownames(Z0)<- names(countsorig.n)
	colnames(Z0)<- paste("Group", 1:G)
	if(relabel) labelstore<- matrix(0,B,G)
	if(verbose==TRUE) cat("Beginning bootstrapping run...\n")
	####
	##Bootstrap Estimates
	####
	for(b in 1:B)
	{
	if(verbose==TRUE & b%%verbose.update == 0) cat(b, "of", B, "samples completed...\n")
			
	counts.n<-table(sample(1:N,N1, replace=TRUE, prob=countsorig.n))
	lab1<- as.numeric(names(counts.n))
	if(length(counts.n)<N)
	{
		Xstar<-X 
		X<-X[lab1,]
		N<-nrow(X)
		}
	counts.n<-as.numeric(counts.n)
	
	#Set Parameters
	Thetat<-x$itemprob; Taut<-x$classprob
	
	eps<-N*M
	llstore<-0; counter<-0
	while(abs(eps)>conv)
	{
			#E-Step
			dum<-array(apply(Thetat,1,dbinom, size=1, x=t(X)), dim=c(M,N,G))
			Z1<-t(Taut*t(apply(dum, c(2,3), prod)))
			Z<-Z1/apply(Z1,1,sum)
			
			##M-Step
			Z.sum<-colSums(Z*counts.n)
			
			Taut<-(Z.sum+delta-1)/(sum(counts.n) + sum(delta)-G)
			Thetat<-(t(Z)%*%(X*counts.n) + alpha-1)/(Z.sum + alpha + beta -2 + small)
			if(any(Thetat>1)) Thetat[Thetat>1]<-1 ##Rounding Error
			if(any(Thetat<0)) Thetat[Thetat<0]<-0 ##Rounding Error
			
			#Log-Posterior
			l<-sum(log(rowSums(Z1))*counts.n)+sum(xlogy(alpha-1,Thetat)+xlogy(beta-1,1-Thetat))+sum(xlogy(delta-1,Taut))
			
			llstore[counter]<-l
			if(counter>2)
			{
				ll.inf<-(llstore[counter-1] - llstore[counter-2])/(llstore[counter] - llstore[counter-1])*(llstore[counter] - llstore[counter-2]) + llstore[counter-1]
		
				if(llstore[counter] == llstore[counter-1]){ ll.inf<- llstore[counter]}
		
				eps<- ll.inf - llstore[counter]
					} #Convergence
					
			counter<-counter+1		
					
			if(counter>iter) break			
			}#while
		  Z1<- matrix(0, Nstar, G)
		  Z1[lab1, ]<- Z

		if(relabel){
		  Z1<- matrix(0, Nstar, G)
		  Z1[lab1, ]<- Z
		  c1<- rep(0,Nstar)
		  c1[lab1]<- counts.n
		#if(N<Nstar){ c1<- rep(0,Nstar)
		  #Z1[-lab1,]<- Zscore.internal(Thetat, Taut, Xstar[-lab1,])
		 # counts
		 # }
		  match1<- matchClasses(t(Zcheck)%*%(sqrt(c1)*Z1), method='exact', verbose=FALSE)
		  labelstore[b, ]<- match1
		  if(swapcheck == FALSE & any(match1 != 1:G)) swapcheck<- TRUE
		  } else match1<- 1:G
		#if(any(is.nan(Z1))) match1<- matchClasses(t(Zcheck[lab1, ])%*%Z1[lab1,], method='exact', verbose=FALSE) else match1<- matchClasses(t(Zcheck)%*%Z1, method='exact', verbose=FALSE)

		#print(match1)
		
		Ti.Theta[b, , ]<-Thetat[match1, ]
		Ti.Tau[b, match1]<-Taut[match1]
		Z0<- ((b-1)*Z0 + Z1[, match1])/b
		
		if(N<Nstar){ 
			X<-Xstar
			N<-nrow(X)
			}
		}#b
	if(verbose==TRUE) cat("Bootstrap sampling run completed.\n")
	if(swapcheck==TRUE) warning("Some samples re-labelled to prevent label switching occurring. Some checking of density plots is recommended. Use '?plot.blca' for more details. ")
		
	######
	##Return Values
	######
	tau<- apply(Ti.Tau,2,mean)
	o<- order(tau, decreasing=TRUE)

	boot<-NULL
	boot$call<- match.call()
	boot$itemprob<-apply(Ti.Theta, c(2,3),mean)[o,]
	boot$classprob<- tau[o]
	boot$Z<- Z0[,o]
	boot$itemprob.sd<- boot$itemprob.se<- sqrt(apply(Ti.Theta, c(2,3),var))[o,]
	boot$classprob.sd<- boot$classprob.se<- sqrt(apply(Ti.Tau,2,var))[o]
	
	boot$classprob.initial<-x$classprob
	boot$itemprob.initial<-x$itemprob
	
	boot$samples<-NULL
	
	boot$samples$classprob<-Ti.Tau[, o]
	boot$samples$itemprob<-Ti.Theta[, o, ]
	
	dum<-array(apply(boot$itemprob, 1, dbinom, size=1, x=t(X)), dim=c(M,N,G))
	Z1<-t(boot$classprob*t(apply(dum, c(2,3), prod)))
	
	boot$logpost<-sum(log(rowSums(Z1))*countsorig.n)+sum(xlogy(alpha[o, ]-1,boot$itemprob)+xlogy(beta[o, ]-1,1-boot$itemprob))+sum(xlogy(delta[o]-1, boot$classprob))
	
	boot$BIC<- 2*boot$logpost - (G*M + G-1)*log(N1)
	boot$AIC<- 2*boot$logpost - 2*(G*M + G-1)
	
	if(!is.null(colnames(X))){
		colnames(boot$itemprob) <- colnames(boot$itemprob.se) <- colnames(boot$itemprob.initial) <-  colnames(X)
		dimnames(boot$samples$itemprob)<- list(NULL, NULL, colnames(X)) 
	}
	if(relabel) boot$label<- labelstore else boot$label<- FALSE
	boot$counts<- countsorig.n
	
	boot$prior<-NULL
	boot$prior$alpha<- alpha[o, ]
	boot$prior$beta<- beta[o, ]
	boot$prior$delta<- delta[o]
	
	class(boot)<-c("blca.boot", "blca")
	boot
	}
