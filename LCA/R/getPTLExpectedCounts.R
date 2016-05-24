getPTLExpectedCounts <- function(alpha,beta,gamma,bin_limits,ntrials){
        counts <- rep(0,length(bin_limits)-1)
        if(beta>0){
                for(i in 1:length(counts)){
                        counts[i] <- pPTL(q=bin_limits[i+1],alpha=alpha,beta=beta,gamma=gamma)-pPTL(q=bin_limits[i],alpha=alpha,beta=beta,gamma=gamma)
                }
        }
        counts*ntrials
}

