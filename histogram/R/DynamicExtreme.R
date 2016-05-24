`DynamicExtreme` <-
function(FctWeight,n,Dmax,mini=TRUE,msg=TRUE) {
	
	###################################################################
    # use only for vectorizing the computation
    OptimizePath <- function(D) {
    	ancestor0 = matrix(nrow=n-D+1)
    	cumWeight0 = matrix(nrow=n-D+1)

        for(i in D:n) {
            tmp <- cumWeight[(D-1):(i-1),D-1]+weight[D:i,i+1]
            
            if (mini) ancestor0[i-D+1] <- which.min(tmp)
            else ancestor0[i-D+1] <- which.max(tmp)
            cumWeight0[i-D+1] <- tmp[ancestor0[i-D+1]]
            }
        ancestor[D:n,D] <<- ancestor0 + (D-2)
        cumWeight[D:n,D] <<- cumWeight0
        }
 	###################################################################
       
	if (n==1) return(list(extreme=FctWeight(1,2),ancestor=cbind(0,0) ))

    # matrix of extreme values reached in d=1:Dmax steps
    cumWeight <- matrix(nrow=n,ncol=Dmax)
    
    # matrix of the ancestors
    ancestor <- matrix(0,nrow=n,ncol=Dmax)
    
    # weight matrix for combinatorial
    weight <- matrix(nrow=n+1,ncol=n+1)
    
    # first compute all n*(n-1)/2 weights
    if (msg) message("- Computing weights for dynamic programming algorithm.")
    
    for(i in 1:n) 
    	weight[i,(i+1):(n+1)] = mapply(FctWeight,rep(i,n-i+1),(i+1):(n+1));
    	    	    
    cumWeight[,1] <- weight[1,2:(n+1)] 
    
    if (msg) message("- Now performing dynamic optimization.")
    
    mapply(OptimizePath,2:Dmax)    
    extreme <- cumWeight[n,]


    list(extreme=extreme,ancestor=ancestor)
}

