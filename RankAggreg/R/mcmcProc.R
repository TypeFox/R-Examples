`mcmcProc` <-
function(v, N, thin, burn.in, comp.list, verbose=TRUE)
{
    if(is.null(thin))
	   thin = 20 # default value for thinning
    if(is.null(burn.in))
	   burn.in=N # default value for burn.in

    n <- nrow(v)
    k <- ncol(v)
    M <- (N+burn.in)*thin

    samples <- matrix(0,N,k)
    ind <- 1
	
	# start with the best candidate
	last <- rep(0,n)
	excl <- n+1
    for(j in 1:k){
         last[j] <- rmult(1, (1:n)[-excl], v[,j][-excl])
         excl <- last[1:j]
	}
	last[(k+1):n] <- sample((1:n)[-last[1:k]], n-k)

	#last <- comp.list[sample(1:n,n)]  
    
    vs <- sample(1:k, M, replace=TRUE)
    hs <- sample(1:n, M, replace=TRUE)
    unifs <- runif(M)   
    
    if(verbose)
        cat("MCMC: ")   
    for (i in 1:M)
    {      
        x <- last
        vert <- vs[i]
        horz <- hs[i]

        temp <- last[vert]      
        x[vert] <- last[horz]
        x[horz] <- temp
    
        if(horz > k)        
            A <- v[x[vert],vert]/v[last[vert],vert]
        else
            A <- v[x[vert],vert]*v[x[horz],horz]/
                (v[last[vert],vert]*v[last[horz],horz])
    
        if (A > unifs[i]){
            if(i%%thin==0 && i>=burn.in*thin)      
                samples[i/thin-burn.in,] <- x[1:k]
            last <- x
        }
        else
            if(i%%thin==0 && i>=burn.in*thin)
                samples[i/thin-burn.in,] <- last[1:k]

        if(i%%(M/10)==0 && verbose)
            cat(paste(i*100/M, "% ",sep=""))
    }
    samples
}

