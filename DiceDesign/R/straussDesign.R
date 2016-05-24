	straussDesign <- function(n,dimension,RND,alpha=0.5,repulsion=0.001,NMC=1000,constraints1D=0,repulsion1D=0.0001,seed=NULL){
			
    # if no seed is provided in argument, choice of the seed for 'runif'
	if (is.null(seed)){
		seed <- as.numeric(Sys.time())
	}
	set.seed(seed)	
	
	# Case with potential different to zero
	# the C code is computed for beta = -ln(gamma)
	if (alpha != 0){
		repulsion <- -log(repulsion)
	} # End of the case with no constraint
	
	# Case with no constraints
	if(constraints1D!=0){
		R1D <- 0.75/n 
	}

	# initial design
	init <- runif(n*dimension)
		
 	out <- .C("C_StraussDesign", as.double(init), as.integer(n), 
	as.integer(dimension), as.integer(constraints1D), as.integer(NMC), 
	as.double(RND), as.double(alpha), as.double(repulsion), 
	as.double(repulsion1D),as.integer(seed),ans = double(n * dimension), 
	PACKAGE="DiceDesign")
		
	if (alpha != 0){
		repulsion <- exp(-repulsion)
	}
		
	# Outputs
	if (constraints1D==0){
		return(list(n=n, dimension=dimension, 
			design_init=t(matrix(init,ncol=n,nrow=dimension,byrow=FALSE)), 
			radius=RND, alpha=alpha, repulsion=repulsion, NMC=NMC,    
			seed=seed, constraints1D=constraints1D,
			design=matrix(out$ans,nrow=n,ncol=dimension,byrow=TRUE)))
	} else {
		return(list(n=n,dimension=dimension,
			design_init=matrix(init,ncol=dimension,nrow=n,byrow=TRUE), 	
			radius=RND, alpha=alpha, repulsion=repulsion,NMC=NMC,  	
			seed=seed, constraints1D=constraints1D, repulsion1D=repulsion1D, 
			design=matrix(out$ans,nrow=n, ncol=dimension,byrow=TRUE)))
	}

}