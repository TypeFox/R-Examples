### ergmm.MCMC.C: This is a pretty minimal R function that prepares the R data to be
### passed into the C function, calls the C function to estimate the latent space model,
### and then puts the C data back into readable R storage. The hope is to separate the
### part specific to interfacing with C from the rest of the program.
### It does NOT come up with initial values --- those must be passed to it.

### If present, parameters 'sample.size' and 'interval' override those in 'control'.

### Also note that it does NOT perform the burnin. To perform the burnin,
### pass sample.size=1 and interval=burnin, and then pass the last (and only) iteration
### to the actual run.


ergmm.MCMC.C<-function(model, start, prior, control, sample.size=NULL, interval=NULL){
  Ym.noNA<- Ym <-model[["Ym"]]
  Ym.noNA[is.na(Ym.noNA)]<-0
    
  n<-network.size(model[["Yg"]])
  d<-model[["d"]]
  G<-model[["G"]]
  p<-model[["p"]]

  if(is.null(sample.size)) sample.size<-control[["sample.size"]]
  if(is.null(interval)) interval<-control[["interval"]]
  
  if(length(prior[["beta.mean"]])==1) prior[["beta.mean"]]<-rep(prior[["beta.mean"]],p)
  if(length(prior[["beta.var"]])==1) prior[["beta.var"]]<-rep(prior[["beta.var"]],p)

  ## Figure out the design matrix.
  observed<-observed.dyads(model[["Yg"]])
  
  if((all(observed==(diag(n)==0)) && is.directed(model[["Yg"]])) ||
     (all(observed==lower.tri(diag(n))) && !is.directed(model[["Yg"]])))
    observed<-NULL

  ## Sanity checks: the following block of code checks that all dimensionalities and
  ## dimensions are correct, and those optional parameters that are required by the presence
  ## of other optional parameters are present.

  if(sample.size!=round(sample.size)) stop("Non-integer MCMC sample size.")
  if(interval!=round(interval)) stop("Non-integer MCMC interval.")
  
  for(i in 1:p)
    if(!all(dim(model[["X"]][[i]])==c(n,n))) stop("Incorrect size for covariate matrices.")

  if(!is.null(start[["Z"]])){
    if(!all(dim(start[["Z"]])==c(n,d))) stop("Incorrect size for the starting latent positions.")
    if(is.null(control[["Z.delta"]])) stop("Need Z-move proposal standard deviation (control[[\"Z.delta\"]]).")
    if(G > 0){
      if(length(start[["Z.K"]])!=n) stop("Incorrect length for the vector of starting cluster assignments.")
      if(length(start[["Z.pK"]])!=G) stop("Incorrect length for the vector of starting cluster probabilities.")
      if(!all(dim(start[["Z.mean"]])==c(G,d))) stop("Incorrect size for the starting cluster means.")
      if(length(start[["Z.var"]])!=G) stop("Incorrect size for the starting cluster variances.")
    } else{
       if(length(start[["Z.var"]])!=1) stop("Missing starting latent space variance.")
    }
  }  
  if(length(start[["beta"]])!=p) stop("Incorrect length for the starting beta vector.")
  if(length(prior[["beta.mean"]])!=p) stop("Incorrect length for the prior beta mean vector.")
  if(length(prior[["beta.var"]])!=p) stop("Incorrect length for the prior beta standard deviation vector.")

  if(model[["sender"]] || model[["receiver"]] || model[["sociality"]]){
    if(is.null(control[["RE.delta"]])) stop("Need random effect proposal sd.")
  }

  if(model[["sender"]]){
    if(length(start[["sender"]])!=n) stop("Incorrect length for the vector of starting sender effects.")
    if(is.null(start[["sender.var"]])) stop("Need starting sender variance.")
    if(is.null(prior[["sender.var"]])) stop("Need sender prior variance.")
    if(is.null(prior[["sender.var.df"]])) stop("Need sender prior variance df.")
  }
  if(model[["receiver"]]){
    if(length(start[["receiver"]])!=n) stop("Incorrect length for the vector of starting receiver effects.")
    if(is.null(start[["receiver.var"]])) stop("Need starting receiver variance.")
    if(is.null(prior[["receiver.var"]])) stop("Need receiver prior variance.")
    if(is.null(prior[["receiver.var.df"]])) stop("Need receiver prior variance df.")
  }

  if(model[["sociality"]]){
    if(length(start[["sociality"]])!=n) stop("Incorrect length for the vector of starting sociality effects.")
    if(is.null(start[["sociality.var"]])) stop("Need starting sociality variance.")
    if(is.null(prior[["sociality.var"]])) stop("Need sociality prior variance.")
    if(is.null(prior[["sociality.var.df"]])) stop("Need sociality prior variance df.")
  }

  if(model[["dispersion"]]){
    if(is.null(start[["dispersion"]])) stop("Need starting dispersion.")
    if(is.null(prior[["dispersion"]])) stop("Need dispersion prior.")
    if(is.null(prior[["dispersion.df"]])) stop("Need dispersion prior df.")
  }

  ## End Sanity checks.

  RESERVED<-2

#  cat("Entering C routine... ")
  Cret <- .C("ERGMM_MCMC_wrapper",
             # 1:
             sample.size=as.integer(sample.size),
             interval=as.integer(interval),
             # 3:
             n=as.integer(n),
             p=as.integer(p),
             d=as.integer(d),
             G=as.integer(NVL(G,0)),
             latent=as.integer(NVL(model[["latentID"]],0)),
             family=as.integer(NVL(model[["familyID"]],0)),
             res=as.integer(with(model,c(sender,receiver,sociality,dispersion))),
             # 10:
             dir=as.integer(is.directed(model[["Yg"]])),
             viY=as.integer(Ym.noNA),
             vdY=as.double(Ym.noNA),
             iconsts=as.integer(model[["iconsts"]]),
             dconsts=as.double(model[["dconsts"]]),
             # 15:
             vX=as.double(unlist(model[["X"]])),  
             # 16:
             lpY.mcmc=double(sample.size+RESERVED),
             lpZ.mcmc=if(!is.null(start[["Z"]])) double(sample.size+RESERVED) else double(0),
             lpbeta.mcmc=if(p>0) double(sample.size+RESERVED) else double(0),
             lpRE.mcmc=if(model[["sender"]] || model[["sociality"]] || model[["receiver"]]) double(sample.size+RESERVED) else double(0),
             lpLV.mcmc=if(!is.null(start[["Z"]])) double(sample.size+RESERVED) else double(0),
             lpREV.mcmc=if(model[["sender"]] || model[["sociality"]] || model[["receiver"]]) double(sample.size+RESERVED) else double(0),
             lpdispersion.mcmc=if(model[["dispersion"]]) double(sample.size+RESERVED) else double(0),
             # 23:
             Z=as.double(start[["Z"]]),
             # 24:
             Z.pK=if(G > 0) as.double(start[["Z.pK"]]) else double(0),
             Z.mean=if(G > 0) as.double(start[["Z.mean"]]) else double(0),
             Z.var=as.double(start[["Z.var"]]),
             Z.K=if(G > 0) as.integer(start[["Z.K"]]) else integer(0),
             # 28:
             prior.Z.var=as.double(prior[["Z.var"]]),
             prior.Z.mean.var=if(G > 0) as.double(prior[["Z.mean.var"]]) else double(0),
             prior.Z.pK=if(G > 0) as.double(prior[["Z.pK"]]) else double(0),
             prior.Z.var.df=as.double(prior[["Z.var.df"]]),
             # 32:
             Z.mcmc = double((sample.size+RESERVED)*n*d),
             Z.rate = if(d > 0 || model[["sender"]] || model[["sociality"]] || model[["receiver"]]) double((sample.size+RESERVED)) else double(0),
             # 34:
             K.mcmc = if(G > 0) integer(n*(sample.size+RESERVED)) else integer(0),
             Z.pK.mcmc = double(G*(sample.size+RESERVED)),
             mu.mcmc = double(d*G*(sample.size+RESERVED)),
             Z.var.mcmc = double(max(G,d>0)*(sample.size+RESERVED)),
             # 38:
             start.beta=as.double(start[["beta"]]),
             prior.beta.mean=as.double(prior[["beta.mean"]]),
             prior.beta.var=as.double(prior[["beta.var"]]),
             beta.mcmc=double((sample.size+RESERVED)*p),
             beta.rate=double((sample.size+RESERVED)),
             # 43:
             start.sender=if(model[["sociality"]]) as.double(start[["sociality"]]) else as.double(start[["sender"]]),
             start.receiver=as.double(start[["receiver"]]),
             sender.var=if(model[["sociality"]]) as.double(start[["sociality.var"]]) else as.double(start[["sender.var"]]),
             receiver.var=as.double(start[["receiver.var"]]),
             # 47:
             prior.sender.var=if(model[["sociality"]]) as.double(prior[["sociality.var"]]) else as.double(prior[["sender.var"]]),
             prior.sender.var.df=if(model[["sociality"]]) as.double(prior[["sociality.var.df"]]) else as.double(prior[["sender.var.df"]]),
             prior.receiver.var=as.double(prior[["receiver.var"]]),
             prior.receiver.var.df=as.double(prior[["receiver.var.df"]]),
             # 51:
             sender.mcmc=if(model[["sender"]] || model[["sociality"]]) double(n*(sample.size+RESERVED)) else double(0),
             receiver.mcmc=if(model[["receiver"]]) double(n*(sample.size+RESERVED)) else double(0),
             sender.var.mcmc=if(model[["sender"]] || model[["sociality"]]) double((sample.size+RESERVED)) else double(0),
             receiver.var.mcmc=if(model[["receiver"]]) double((sample.size+RESERVED)) else double(0),
             # 55:
             start.dispersion=if(model[["dispersion"]]) as.double(start[["dispersion"]]) else as.double(0),
             prior.dispersion=if(model[["dispersion"]]) as.double(prior[["dispersion"]]) else as.double(0),
             prior.dispersion.df=if(model[["dispersion"]]) as.double(prior[["dispersion.df"]]) else as.double(0),
             dispersion.mcmc=if(model[["dispersion"]]) double(sample.size+RESERVED) else double(0),
             # 59:
             observed=as.integer(NVL(observed,-1)),
             # 60:
             deltas=with(control,as.numeric(c(Z.delta,RE.delta,group.deltas))),
             beta.eff.sender=as.double(c(model[["beta.eff.sender"]],model[["beta.eff.sociality"]])),
             beta.eff.sender.size=as.integer(NROW(rbind(model[["beta.eff.sender"]],model[["beta.eff.sociality"]]))),
             beta.eff.receiver=as.double(model[["beta.eff.receiver"]]),
             beta.eff.receiver.size=as.integer(NROW(model[["beta.eff.receiver"]])),
             # 65:
             accept.all=control[["accept.all"]],

             PACKAGE="latentnet")
#  cat("Finished C routine.\n")
  
  sample<-list(## MCMC Sample
                lpY=Cret[["lpY.mcmc"]],
                lpZ=Cret[["lpZ.mcmc"]],
                lpbeta=Cret[["lpbeta.mcmc"]],
                lpRE=Cret[["lpRE.mcmc"]],
                lpLV=Cret[["lpLV.mcmc"]],
                lpREV=Cret[["lpREV.mcmc"]],
                lpdispersion=Cret[["lpdispersion.mcmc"]],
                beta=matrix(Cret[["beta.mcmc"]],ncol=p),
                beta.rate=Cret[["beta.rate"]],
                Z.K = if(G>0) matrix(Cret[["K.mcmc"]],ncol=n),
                Z.mean = if(G>0) array(Cret[["mu.mcmc"]],dim=c((sample.size+RESERVED),G,d)),
                Z.var = if(d>0) matrix(Cret[["Z.var.mcmc"]],ncol=max(G,1)),
                Z.pK = if(G>0) matrix(Cret[["Z.pK.mcmc"]],ncol=G),
                Z=if(d>0)array(Cret[["Z.mcmc"]],dim=c((sample.size+RESERVED),n,d)),
                Z.rate=if(d>0 || model[["sender"]] || model[["sociality"]] || model[["receiver"]]) Cret[["Z.rate"]],
                sender=if(model[["sender"]] && !model[["sociality"]]) matrix(Cret[["sender.mcmc"]],ncol=n),
                receiver=if(model[["receiver"]]) matrix(Cret[["receiver.mcmc"]],ncol=n),
                sociality=if(model[["sociality"]]) matrix(Cret[["sender.mcmc"]],ncol=n),
                sender.var=if(model[["sender"]] && !model[["sociality"]]) Cret[["sender.var.mcmc"]],
                receiver.var=if(model[["receiver"]]) Cret[["receiver.var.mcmc"]],
                sociality.var=if(model[["sociality"]]) Cret[["sender.var.mcmc"]],
                dispersion=if(model[["dispersion"]]) Cret[["dispersion.mcmc"]]
                )
  class(sample)<-"ergmm.par.list"
  
  
  mcmc.mle<-sample[[1]]
  mcmc.pmode<-sample[[2]]
  sample<-del.iteration(sample,1:2)
  
  
  ## Construct the list (of lists) for return.
  out<-list(sample=sample,
            mcmc.mle=mcmc.mle,
            mcmc.pmode=mcmc.pmode)

  out
}



ergmm.MCMC.snowFT<-function(threads, reps, model.l, start.l, prior.l, control.l, sample.size.l=NULL, interval.l=NULL){
  l.sizes<-c(length(model.l),
             length(start.l),
             length(prior.l),
             length(control.l),
             length(sample.size.l),
             length(interval.l))
  l.sizes<-l.sizes[l.sizes>0]
  param.sets<-max(l.sizes)
  if(any(l.sizes!=param.sets & l.sizes!=1)) stop("Length of each input list must be either 1 or a single other number.")

  if(!requireNamespace("snowFT",quietly=TRUE)) stop("Package 'snowFT' is required for multithreaded MCMC.")
  mcmc.out.l<-snowFT::performParallel(threads,rep(1:param.sets,reps),
                              ergmm.MCMC.snowFT.slave,
                              lib=.latentnetEnv$path.to.me,
                              model.l=model.l,
                              start.l=start.l,
                              prior.l=prior.l,
                              control.l=control.l,
                              sample.size.l=sample.size.l,
                              interval.l=interval.l,
                              seed=floor(runif(6,0,.Machine[["integer.max"]])))
  mcmc.mle<-mcmc.out.l[[which.max(sapply(1:length(mcmc.out.l),
                                         function(i) mcmc.out.l[[i]][["mcmc.mle"]][["lpY"]]))]][["mcmc.mle"]]
  mcmc.pmode<-mcmc.out.l[[which.max(sapply(1:length(mcmc.out.l),
                                         function(i) lpsum(mcmc.out.l[[i]][["mcmc.pmode"]])))]][["mcmc.pmode"]]
  result.list<-list(sample=list(),mcmc.mle=mcmc.mle,mcmc.pmode=mcmc.pmode)
  for(i in 1:length(mcmc.out.l)) result.list[["sample"]][[i]]<-mcmc.out.l[[i]][["sample"]]
  class(result.list$sample) <- "ergmm.mcmc.list.list"
  result.list
}

ergmm.MCMC.snowFT.slave<-function(i, lib, model.l, start.l, prior.l, control.l, sample.size.l=NULL, interval.l=NULL){
  library(latentnet,lib.loc=lib)
  ergmm.MCMC.C(model.l[[min(length(model.l),i)]],
               start.l[[min(length(start.l),i)]],
               prior.l[[min(length(prior.l),i)]],
               control.l[[min(length(control.l),i)]],
               sample.size.l[[min(length(sample.size.l),i)]],
               interval.l[[min(length(interval.l),i)]])
}

