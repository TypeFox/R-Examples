hybridMC <-
function(y.start,n.samp=1,logDens,dLogDens,epsilon,LFsteps,compWeights=NULL,MPwidth=1,MPweights=NULL,progress=0,...){

## Sanity checks and definitions
	if(mode(logDens) != "function") stop("logDens is not a function.")
	if(mode(dLogDens) != "function") stop("dLogDens is not a function.")
	if(!is.numeric(y.start) | length(y.start)<1 | length(dim(y.start))>1) stop("Invalid starting values.")
	if(n.samp<1) stop("Invalid number of samples.")
	if(LFsteps<1) stop("Number of leapfrog steps must be at least 1.")
	if(any(compWeights<0)) stop("Component weights must be nonnegative")
	if(MPwidth>LFsteps) stop("Multipoint window cannot be larger than number of leapfrog steps.")
	if(MPwidth<1) stop("Invalid multipoint window length.")

	MPwidth= as.integer(MPwidth)
	n.samp = as.integer(n.samp)
	LFsteps= as.integer(LFsteps)
	Ex     = function(x) -logDens(x,...)
	dEx    = function(x) -dLogDens(x,...)
	if(is.null(compWeights) | length(compWeights)==1) compWeights=rep(1,length(y.start))	
	if(is.null(MPweights) | length(MPweights)==1) MPweights=rep(1,MPwidth)
	if(length(MPweights)!=MPwidth) stop("Invalid MPweights vector.")
	if(length(compWeights)!=length(y.start)) stop("Invalid compWeights vector.")

	if(length(MPweights)<MPwidth) stop("Multipoint weight vector must be at least the length of the Multipoint width.")	
	if(any(MPweights<0)) stop("Multipoint weights must be nonnegative.")

	if(any(is.na(Ex(y.start)))) stop("Function could not be evaluated at initial values.")
	if(any(is.na(dEx(y.start))) | length(dEx(y.start))!=length(y.start)) stop("Gradient function returned NA or wrong length.")

	if(length(epsilon)>1){ 
		epsLo=epsilon[1]
		epsRng=epsilon[2]-epsilon[1]
	}else if(length(epsilon)==1){
		epsLo=epsilon
		epsRng=0
	}else{
		stop("Invalid epsilon specified.")	
	}
	
	if(epsLo<=0 | epsRng<0) stop("Invalid epsilon specified.")

	progress = as.integer(abs(progress))

## Create progress bar
	if(progress){ pb = txtProgressBar(min = 0, max = n.samp, style = 3) }else{ pb=NULL }
	pbFun = function(samps){ if(progress) setTxtProgressBar(pb, samps)}
## Call the C code
	
	if(MPwidth==1){
		samples = .Call("hybridMC",y.start,n.samp,Ex,dEx,epsLo,epsRng,LFsteps,compWeights,progress,pbFun,new.env(),PACKAGE="HybridMC")
	}else{
		samples = .Call("MPhybridMC",y.start,n.samp,Ex,dEx,epsLo,epsRng,LFsteps,compWeights,MPwidth,MPweights,progress,pbFun,new.env(),PACKAGE="HybridMC")	
	}
	
	if(progress & is.list(pb)) close(pb)
	dim(samples)=c(length(y.start),n.samp)
	return(mcmc(t(samples)))
}
