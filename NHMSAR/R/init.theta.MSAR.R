init.theta.MSAR <-
function(data,...,M,order,regime_names=NULL,nh.emissions=NULL,nh.transitions=NULL,label=NULL,ncov.emis = 0,ncov.trans=0,cl.init="mean"
) {
	if (missing(M) || is.null(M) || M==0) {print("Need at least one regime : M=1"); M <- 1 }
	if (missing(order)) { order <- 0 } # AR order
    if (missing(label)) {label = 'HH'}
    T = dim(data)[1]
    N.samples = dim(data)[2]
    d = dim(data)[3]
    if (is.null(d) ) {d=1}
    else if ( is.na(d)) {d=1}
    #if (d==1) {data <- as.matrix(data)}

     if (length(dim(data))<3) {d = 1}
     else {d = dim(data)[3]} 
     
     
     if (cl.init == "mean") {
		 d.data = matrix(0,(T-1)*N.samples,d)
		 data.mat = matrix(0,T*N.samples,d)
		 for (id in 1:d) { 
			 d.data[,id] = c(data[2:T,,id]-data[1:(T-1),,id])
			 data.mat[,id] = c(data[,,id])
		 }
		 if (M>1) {
		 	centers=d.data[sample(1:dim(d.data)[1],M),] 
		 	if (any(duplicated(centers))) {
		 		cnt=0
		 		while (any(duplicated(centers)) & cnt<5){ 
		 			cnt = cnt+1
		 			centers = d.data[sample(1:dim(d.data)[1],M),]
		 		}
		 		if (any(duplicated(centers))) {print("Is there a lot of equal data?")}
		 	}
			mc = kmeans(d.data,centers=centers)
			class = mc$cluster
		 }
		 else {class = 1}
	 }
     else {    
     	d.var = array(0,c((T-1),N.samples,d))
		ht = 3
    	for (id in 1:d) { 
			for (t in 1:(T-1)) {
        		ii = max(t-ht,1):min(T-1,t+ht)
     			d.var[t,,id] = apply(data[ii+1,,id]-data[ii,,id],2,var)    
     		}
     	}
		qv = quantile(d.var,probs=(1:M)/M)
		d.var = matrix(d.var,(T-1)*N.samples,d)    
		if (M>1) {
			mc = kmeans(d.var,centers=d.var[sample(1:dim(d.var)[1],M),])
			class = mc$cluster
		}
	 }
	 res = Mstep.classif(data,array(class,c(T,N.samples,1)),order=order)
     A = res$A
     A0 = res$A0
     sigma = res$sigma
     
     
      
     #TRANSITION PROBABILITY MATRIX / INITIAL DISTRIBUTION
     if (M>1) {
     	prior <- normalise(runif(M)) 
     	transmat <- mk_stochastic(diag(1,M)+matrix(runif(M*M),M)) # homogeneous Markov chain
     } else {
     	prior = 1
     	transmat = 1
     } 
     prior=matrix(prior,M,1)
     transmat=matrix(transmat,M,M)
    
     n_par= M*(M-1)+(order*d^2+d+d*(d+1)/2)*M 
     
     #EMISSION PARAMETERS
     emis.linear = FALSE
     if (substr(label,2,2)=="N") {
     	par.emis <- list()
      	if (missing(nh.emissions)) {nh.emissions = 'linear'}
      	if (!is.function(nh.emissions)){
     	if (nh.emissions == 'linear') {
     		emis.linear = TRUE
     		if (ncov.emis<1) {ncov.emis=1}
     		nh.emissions <- function(covar,par.emis){
     			d <- dim(par.emis)[1]
     			if(is.null(d) || is.na(d)){d <- 1}
     			f <- matrix(0,d,dim(covar)[1]) # if covar is a vector?
     			for(i in 1:d){
     				f[i,] <- par.emis[i,1:dim(par.emis)[2]]%*%t(covar)
     			}
     			return(f)
     		}
     	}
     	}
     	for(i in 1:M){
     		par.emis[[i]] <- matrix(0,d,ncov.emis)
     	}

        n_par <- n_par+ncov.emis
     }
      
      #TRANSITION PARAMETERS
      if (substr(label,1,1)=="N") { 
           par.trans <- array(1,c(M,max(2,ncov.trans+1)))
     	   if (missing(nh.transitions)) { nh.transition = 'gauss'}
     	   if (!is.function(nh.transitions)){if (nh.transitions=='gauss') {     	   		nh.transitions <- function(covar,par.trans,transmat){		
  					T <- dim(covar)[1]
  					N.samples = dim(covar)[2]
  					ncov = dim(covar)[3]
  					M = dim(transmat)[1]
  					f <- array(0,c(M,M,T))
  					p.t = array(0,c(T,N.samples,ncov))
  					for (j in 1:M) {	
  						for (nc in 1:ncov) {p.t[,,nc] = par.trans[j,nc]}
  						xx = (covar-p.t)^2 
  						sxx = apply(xx,1,sum)
  						temp=exp( -.5*sxx/par.trans[j,ncov+1]^2)/par.trans[j,ncov+1]
  						for (i in 1:M) {f[i,j,] = transmat[i,j]*temp }  
  					}	
  					f.sum = apply(f , c(1,3), sum)
  					for (i in 1:M) {f[i,,] = f[i,,]/t(matrix(f.sum[i,],T,M))}
 					return(f)
				}
	  }
     	   else if (nh.transitions=='logistic') {
     	   		nh.transitions <- function(covar,par.trans,transmat){
  					eps = 1e-10 
  					T <- dim(covar)[1]
  					nc <- dim(covar)[3]
  					covar = matrix(covar,T,nc)
  					M = dim(transmat)[1]
  					f <- array(0,c(M,M,T))
  					for (m in 1:M) {
  						f[m,m,] = eps+(1-2*eps)/(1 + exp(par.trans[m,1] + par.trans[m,2:(nc+1)] %*% t(covar)))
  					}
  					f[1,2,] = 1-f[1,1,]
  					f[2,1,] = 1-f[2,2,]
  					return(f)
  				}
  			}
			   else if (nh.transitions=='VM') {
				   nh.transitions <- function(covar,par.trans,transmat){		
					   T <- dim(covar)[1]
					   N.samples = dim(covar)[2]
					   ncov = dim(covar)[3]
					   M = dim(transmat)[1]
					   par.trans = repmat(t(par.trans[,2])*exp(1i*par.trans[,1]),M,1)%*%(matrix(1,M,M)-diag(1,nrow=M))

					   f <- array(0,c(M,M,T))
					   for (j in 1:M) {	
						   for (i in 1:M) {
							   f[i,j,] = transmat[i,j]*abs(exp(par.trans[i,j]*exp(-1i*covar))) ;
						   }  
					   }	
					   f.sum = apply(f , c(1,3), sum)
					   for (i in 1:M) {f[i,,] = f[i,,]/t(matrix(f.sum[i,],T,M))}
					   return(f)
				   }
			   } else if (nh.transitions=='probit') {
     	   		nh.transitions <- function(covar,par.trans,transmat){		
  					T <- dim(covar)[1]
  					N.samples = dim(covar)[2]
  					ncov = dim(covar)[3]
  					M = dim(transmat)[1]
  					f <- array(0,c(M,M,T))
  					for (j in 1:M) {	
  						temp=pnorm(par.trans[j,1]+par.trans[j,2:(ncov+1)]%*%t(as.matrix(covar)))
  						 for (i in 1:M) {
  							f[i,j,] = transmat[i,j]*temp ; # On perd la monotonie...
  							#f[i,j,] = temp ;

  						 }  
  					}	
  					f.sum = apply(f , c(1,3), sum)
  					for (i in 1:M) {f[i,,] = f[i,,]/t(matrix(f.sum[i,],T,M))}
 					return(f)
				}
  			}
  			}

            n_par <- n_par+max(2,ncov.trans+1)
      }
  
    if (order==0 && label=='HH'){theta=list(A0,sigma,prior,transmat)}
    if (order==0 && label=='HN'){theta=list(A0,sigma,prior,transmat,par.emis)}
    if (order==0 && label=='NH'){theta=list(A0,sigma,prior,transmat,par.trans)}
    if (order==0 && label=='NN'){theta=list(A0,sigma,prior,transmat,par.trans,par.emis)}
    else if (order>0 && label=='HH'){theta=list(A,A0,sigma,prior,transmat)}
    else if (order>0 && label=='HN'){theta=list(A,A0,sigma,prior,transmat,par.emis)}
    else if (order>0 && label=='NH'){theta=list(A,A0,sigma,prior,transmat,par.trans)}
   # else if (order>0 && label=='NN'){theta=list(A,A0,sigma,prior,transmat,par.emis,par.trans)}
    else if (order>0 && label=='NN'){theta=list(A,A0,sigma,prior,transmat,par.trans,par.emis)}
    
    attr(theta,'NbComp') <- d
    attr(theta,'NbRegimes') <- M
    attr(theta,'order') <- order
    attr(theta,'label') <- label  
    attr(theta,'nh.emissions') <- nh.emissions
    attr(theta,'nh.transitions') <- nh.transitions 
    attr(theta,'n_par') <- n_par 
    attr(theta,'emis.linear') <- emis.linear
    
    theta=as.thetaMSAR(theta,label=label,ncov.emis=ncov.emis,ncov.trans=ncov.trans)
	#class(theta) <- "MSAR"
	#theta$call <- match.call()
    return(theta)
}
