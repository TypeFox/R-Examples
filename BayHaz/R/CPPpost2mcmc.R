
CPPpost2mcmc<-function(sampost){
    chain<-cbind(sampost$sgm,sampost$xi0,sampost$csi)
    colnames(chain)<-paste(c(rep("si",sampost$hyp$F),rep("xi",sampost$hyp$F+1)),
                           c(1:sampost$hyp$F,0,1:sampost$hyp$F),sep="")
    if(require(coda)){ # convert to MCMC
        chain<-mcmc(chain,start=sampost$burnin+1,thin=sampost$thin)
    }else warning("Package 'coda' is not available: a 'matrix' is returned...")
    return(chain)
} # end CPPpostSample
