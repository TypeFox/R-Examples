genperms <-
function(   Z,
                        blockvar = NULL,
                        clustvar = NULL,
                        maxiter=10000) {
                        	
	if (is.null(clustvar)) clustvar <- c(1:length(Z))

    desmat.out <- desmat.sanitize(Z,blockvar,clustvar)
    desmat <- desmat.out$desmat
    B <- nrow(desmat)

    rands <- prod(desmat$rand)
    if(rands > maxiter){
        cat(paste("Too many permutations to use exact method.\nDefaulting to approximate method.\nIncrease maxiter to at least ", rands, " to perform exact estimation.\n",sep=""))
        permclus <- replicate(maxiter,do.call(randfun.default,list(desmat.out)))
    }

    if(rands <= maxiter) {
        perms <- as.matrix(1)
        
        unitind <- rep(NA,length(Z))
        for(b in 1:B){
                Z.b <- desmat.out$Z[desmat.out$blockvar==b]
                N.b <- length(Z.b)
                m.b <- sum(Z.b)
                perms.b <- combn(N.b, m.b, tabulate, nbins = N.b)
                perms <- rbind( t(rep(1,ncol(perms.b))%x%t(perms)),
                                t(t(perms.b)%x%rep(1,ncol(perms))))
            }
            
        perms <- perms[-1,]

        rownames(perms) <- c(1:nrow(perms))

        permclus <- matrix(NA, nrow=length(desmat.out$clustvar), ncol=ncol(perms))
        permclus <- perms[match(desmat.out$clustvar, rownames(perms)),]
    }
    return(permclus)
}
