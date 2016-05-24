init.theta.MSAR.VM <-
function(data,...,M,order,regime_names=NULL,nh.emissions=NULL,nh.transitions=NULL,label=NULL,ncov.emis = 0,ncov.trans=0
) {
	if (missing(M) || is.null(M) || M==0) {print("Need at least one regime : M=1"); M <- 1 }
	if (missing(order)) { order <- 0 } # AR order
    if (missing(label)) {label = 'HH'}
    d = dim(data)[3]
    if (is.null(d) ) {d=1}
    else if ( is.na(d)) {d=1}
    #if (d==1) {data <- as.matrix(data)}

     if (length(dim(data))<3) {d = 1}
     else {d = dim(data)[3]} 
     
     #MEAN
     mu <- matrix(6*runif(M*d),M,d) ;
      
     #kappa (autoregressive coefficients)
     kappa = matrix(1+2*runif(M*(order+1)),M,order+1)
     
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
    
     n_par=order*M+M+M^2
     
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
     			f <- matrix(0,d,dim(covar)[1]) # if covar is a vecteur?
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
     	   if (missing(nh.transitions)) { nh.transition = 'VM'} 
     	   if (!is.function(nh.transitions)){
     	   		if (nh.transitions=='VM') {
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
			} else if (nh.transitions=='gauss') {
     	   		nh.transitions <- function(covar,par.trans,transmat){		
  					T <- dim(covar)[1]
  					N.samples = dim(covar)[2]
  					ncov = dim(covar)[3]
  					M = dim(transmat)[1]
  					f <- array(0,c(M,M,T))
  					for (j in 1:M) {	
  						xx = (covar-array(par.trans[j,1:ncov],c(T,N.samples,ncov)))^2 
  						sxx = apply(xx,1,sum)
  						temp=exp( -sxx/par.trans[j,ncov+1]^2/2)
  						for (i in 1:M) {
  							f[i,j,] = transmat[i,j]*temp ;
  						}  
  					}	
  					f.sum = apply(f , c(1,3), sum)
  					for (i in 1:M) {f[i,,] = f[i,,]/t(matrix(f.sum[i,],T,M))}
 					return(f)
				}
	  } else if (nh.transitions=='logistic') {
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
  			}}

            n_par <- n_par+max(2,ncov.trans+1)
      }
  
    if (label=='HH'){theta=list(mu,kappa,prior,transmat)}
    if (label=='HN'){theta=list(mu,kappa,prior,transmat,par.emis)}
    if (label=='NH'){theta=list(mu,kappa,prior,transmat,par.trans)}
   # if (order==0 && label=='NN'){theta=list(mu,sigma,prior,transmat,par.emis,par.trans)}
    if (label=='NN'){theta=list(mu,kappa,prior,transmat,par.trans,par.emis)}
    
    attr(theta,'NbComp') <- d
    attr(theta,'NbRegimes') <- M
    attr(theta,'order') <- order
    attr(theta,'label') <- label  
    attr(theta,'nh.emissions') <- nh.emissions
    attr(theta,'nh.transitions') <- nh.transitions 
    attr(theta,'n_par') <- n_par 
    attr(theta,'emis.linear') <- emis.linear
    
    theta=as.thetaMSAR.VM(theta,label=label,ncov.emis=ncov.emis,ncov.trans=ncov.trans)
	class(theta) <- "MSAR.VM"
	#theta$call <- match.call()
    return(theta)
}
