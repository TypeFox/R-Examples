bayesmbc<-function(G,Z,prior,Z.K.ref=NULL,sample.size=2000,interval=10,burnin=500,verbose=FALSE){
  start<-mbc.VII.EM(G,Z)
  state<-start<-with(start,list(Z=Z,
                                Z.mean=Z.mean,
                                Z.var=Z.var,
                                Z.K=Z.K,
                                Z.pK=Z.pK)
                     )
  if(verbose>1) cat("Running MBC MCMC... ")
  state<-bayesmbc.MCMC.C(G,start,prior,
                         sample.size=1,interval=burnin)[["sample"]][[1]]
  state[["Z"]]<-Z
  mcmc.out<-bayesmbc.MCMC.C(G,state,prior,sample.size,interval)
  if(verbose>1) cat("Finished.\n")

  if(verbose>1) cat("Running label switching... ")
  Q.start <- switch.Q.K(if(is.null(Z.K.ref)) start[["Z.K"]] else Z.K.ref,G)
  mcmc.out[["sample"]]<-klswitch.C(Q.start,mcmc.out[["sample"]],Z,verbose=verbose)
  if(verbose>1) cat("Finished.\n")
  
  mcmc.out[["pmean"]]<-with(mcmc.out[["sample"]],
                       list(Z.mean=apply(Z.mean,2:3,mean),
                            Z.var=apply(Z.var,2,mean),
                            Z.K=apply(Z.K,2,function(x)which.max(tabulate(x,G))),
                            Z.pZK=t(apply(Z.K,2,function(x)tabulate(x,G)/length(x))),
                            Z.pK=tabulate(c(Z.K),G)/length(Z.K)
                            )
                       )
  mcmc.out
}

bayesmbc.snowFT<-function(threads,G,Z,prior,Z.K.ref=NULL,sample.size=2000,interval=10,burnin=500,verbose=FALSE){
  if(!requireNamespace("snowFT",quietly=TRUE)) stop("Package 'snowFT' required for multithreaded model based clustering MCMC.")
  start<-mbc.VII.EM(G,Z)
  state<-start<-with(start,list(Z=Z,
                                Z.mean=Z.mean,
                                Z.var=Z.var,
                                Z.K=Z.K,
                                Z.pK=Z.pK)
                     )
  if(verbose>1) cat("Running MBC MCMC... ")
  state<-bayesmbc.MCMC.C(G,start,prior,
                         sample.size=1,interval=burnin)[["sample"]][[1]]
  state[["Z"]]<-Z
  mcmc.out<-bayesmbc.MCMC.C(G,state,prior,sample.size,interval)
  if(verbose>1) cat("Finished.\n")

  if(verbose>1) cat("Running label switching... ")
  Q.start <- switch.Q.K(if(is.null(Z.K.ref)) start[["Z.K"]] else Z.K.ref,G)
  mcmc.out[["sample"]]<-klswitch.snowFT(threads,Q.start,mcmc.out[["sample"]],Z,verbose=verbose)
  if(verbose>1) cat("Finished.\n")
  
  mcmc.out[["pmean"]]<-with(mcmc.out[["sample"]],
                       list(Z.mean=apply(Z.mean,2:3,mean),
                            Z.var=apply(Z.var,2,mean),
                            Z.K=apply(Z.K,2,function(x)which.max(tabulate(x,G))),
                            Z.pZK=t(apply(Z.K,2,function(x)tabulate(x,G)/length(x))),
                            Z.pK=tabulate(c(Z.K),G)/length(Z.K)
                            )
                       )
  mcmc.out
}
