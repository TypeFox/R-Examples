revise.ispd <- function(tpm=NULL,gamma=NULL,lns=NULL,cis=TRUE) {
# Function revise.ispd.  To revise the initial state probability
# distribution

if(cis) {
    if(is.null(tpm) + (is.null(gamma) & is.null(lns)) != 1)
    	stop(paste("When \"cis\" is TRUE either \"tpm\" should be given\n",
                   "OR \"gamma\" and \"lns\" should be given.\n"))
    
# If ispd is taken to be the steady state distribution:
    if(!is.null(tpm)) {
    	eee   <- eigen(t(tpm))
    	k <- match(1,round(eee$values,6))
    	if(length(k) != 1) {
    		cat('Problems with eigenvalues:\n')
    		print(eee$values)
    		stop()
    	}
    	v <- Re(eee$vectors[,k])
        ispd <- v/sum(v)
    } else {
    
# Steady state not assumed:
        v <- 0
        j <- 1
        nseq <- length(lns)
        for(i in 1:nseq) {
            v <- v + gamma[,j]
            j <- j + lns[i]
        }
        ispd <- v/nseq
    }
} else {
    if(is.null(gamma) | is.null(lns))
    	stop(paste("When \"cis\" is FALSE \"gamma\" and \"lns\"\n",
                   "must be given.\n"))
    jg   <- 1
    nseq <- length(lns)
    K    <- nrow(gamma)
    ispd <- matrix(0,K,nseq)
    for(j in 1:nseq) {
        ispd[,j] <- gamma[,jg]
        #ispd[which.max(gamma[,jg]),j] <- 1
        jg <- jg + lns[j]
    }
}
return(ispd)
}
