bayesmbc.MCMC.C<-function(G, start, prior, sample.size=NULL, interval=NULL){
  ## Note that passing NULL as a parameter will cause the corresponding parameter in
  ## the C function be set to NULL when NULL is coerced to double.
  ## (as.double(NULL)==double(0))

  Z<-start[["Z"]]
  n<-dim(Z)[1]
  d<-dim(Z)[2]

  ## Sanity checks: the following block of code checks that all dimensionalities and
  ## dimensions are correct, and those optional parameters that are required by the presence
  ## of other optional parameters are present.

  if(!all(dim(start[["Z"]])==c(n,d))) stop("Incorrect size for the starting latent positions.")
  if(G > 0){
    if(length(start[["Z.K"]])!=n) stop("Incorrect length for the vector of starting cluster assignments.")
    if(length(start[["Z.pK"]])!=G) stop("Incorrect length for the vector of starting cluster probabilities.")
    if(!all(dim(start[["Z.mean"]])==c(G,d))) stop("Incorrect size for the starting cluster means.")
    if(length(start[["Z.var"]])!=G) stop("Incorrect size for the starting cluster variances.")
  }

  ## End Sanity checks.

  RESERVED<-1

#  cat("Entering C routine... ")
  Cret <- .C("MBC_MCMC_wrapper",
             sample.size=as.integer(sample.size),
             interval=as.integer(interval),
             
             n=as.integer(n),
             d=as.integer(d),
             G=as.integer(G), 
             
             lpZ.mcmc=double(sample.size+RESERVED),
             lpLV.mcmc=double(sample.size+RESERVED),
             
             Z=as.double(Z),
             
             Z.pK=if(G > 0) as.double(start[["Z.pK"]]) else double(0),
             Z.mean=if(G > 0) as.double(start[["Z.mean"]]) else double(0),
             Z.var=as.double(start[["Z.var"]]),
             Z.K=if(G > 0) as.integer(start[["Z.K"]]) else integer(0),
             
             prior.Z.var=as.double(prior[["Z.var"]]),
             prior.Z.mean.var=if(G > 0) as.double(prior[["Z.mean.var"]]) else double(0),
             prior.Z.pK=if(G > 0) as.double(prior[["Z.pK"]]) else double(0),
             prior.Z.var.df=as.double(prior[["Z.var.df"]]),
             
             K.mcmc = if(G > 0) integer(n*(sample.size+RESERVED)) else integer(0),
             Z.pK.mcmc = if(G > 0) double(G*(sample.size+RESERVED)) else double(0),
             mu.mcmc = double(d*G*(sample.size+RESERVED)),
             Z.var.mcmc = double(max(G,1)*(sample.size+RESERVED)),

             PACKAGE="latentnet")
#  cat("Finished C routine.\n")
  
  sample<-list(## MCMC Sample
                lpZ=Cret[["llike.mcmc"]],
                Z.K = if(G>0) matrix(Cret[["K.mcmc"]],ncol=n),
                Z.mean = if(G>0) array(Cret[["mu.mcmc"]],dim=c((sample.size+RESERVED),G,d)),
                Z.var = if(d>0) matrix(Cret[["Z.var.mcmc"]],ncol=max(G,1)),
                Z.pK = if(G>0) matrix(Cret[["Z.pK.mcmc"]],ncol=G)
                )
  class(sample)<-"ergmm.par.list"
  
  
  mcmc.mle<-sample[[1]]
  
  sample<-del.iteration(sample,1)
  
  ## Construct the list (of lists) for return.
  out<-list(sample=sample,
            mcmc.mle=mcmc.mle)

  out
}
