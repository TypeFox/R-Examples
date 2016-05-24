
BPSpost2mcmc<-function(sampost){
    chain<-sampost$eta
    colnames(chain)<-paste(rep("eta",ncol(sampost$eta)),seq(1,ncol(sampost$eta)),sep="")
    if(require(coda)){ # convert to MCMC
        chain<-mcmc(chain,start=sampost$burnin+1,thin=sampost$thin)
    }else warning("Package 'coda' is not available: a 'matrix' is returned...")
    return(chain)
} # end BPSpostSample
