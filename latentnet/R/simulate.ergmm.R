simulate.ergmm<-function(object, nsim=1, seed=NULL,...){
  extraneous.argcheck(...)
  
  ## If the random seed has been specified, save the old seed, to
  ## pick up where it left off. If not, don't.
  if(!is.null(seed)){
    old.seed<-.Random.seed
    .Random.seed<-seed
  }else runif(1) # This is needed to initialize .Random.seed if it isn't already.
  start.seed<-.Random.seed
  
  l<-list()
  for(i in 1:nsim){
    iter<-floor(runif(1,1,object[["control"]][["sample.size"]]+1))
    l[[i]]<-sim.1.ergmm(object[["model"]],object[["sample"]][[iter]],object[["prior"]])
  }
  
  if(!is.null(seed)) .Random.seed<-old.seed

  if(nsim > 1){
    l <- list(formula = object[["model"]][["formula"]], networks = l,
                     stats = NULL, coef=NULL)
    attr(l,"class")<-"network.series"
  }else{
    l <- l[[1]]
  }
  return(l)
}

simulate.ergmm.model<-function(object,nsim=1,seed=NULL,par,prior=list(),...){
  extraneous.argcheck(...)
  
  ## If the random seed has been specified, save the old seed, to
  ## pick up where it left off. If not, don't.
  if(!is.null(seed)){
    old.seed<-.Random.seed
    .Random.seed<-seed
  }else runif(1) # This is needed to initialize .Random.seed if it isn't already.
  start.seed<-.Random.seed

  l<-list()
  for(i in 1:nsim){
    l[[i]]<-sim.1.ergmm(object,par,prior)
  }
  
  if(!is.null(seed)) .Random.seed<-old.seed

  if(nsim==1) return(l[[1]])
  else{
    attr(l,"class")<-"network.series"
    return(l)
  }
}
  
sim.1.ergmm<-function(model,par,prior=list()){
  nv<-network.size(model[["Yg"]])
  mypar<-par
  
  if(length(model[["X"]])>0 && is.null(mypar[["beta"]]))
    mypar[["beta"]]<-rnorm(length(model[["X"]]),prior[["beta.mean"]],sqrt(prior[["beta.var"]]))

  if(model[["d"]]>0 && is.null(mypar[["Z"]])){
    if(model[["G"]]>0){
      if(is.null(mypar[["Z.mean"]]))
        mypar[["Z.mean"]]<-matrix(rnorm(model[["G"]]*model[["d"]],0,sqrt(prior[["Z.mean.var"]])),nrow=model[["G"]])
      if(is.null(mypar[["Z.K"]]))
        mypar[["Z.K"]]<-sample(1:model[["G"]],nv,replace=TRUE)
    }
    
    if(is.null(mypar[["Z.var"]]))
      mypar[["Z.var"]]<-prior[["Z.var"]]*prior[["Z.var.df"]]/rchisq(max(model[["G"]],1),prior[["Z.var.df"]])

    mypar[["Z"]]<-matrix(rnorm(nv*model[["d"]],
                            if(model[["G"]]>0) mypar[["Z.mean"]][mypar[["Z.K"]],] else 0,
                            if(model[["G"]]>0) mypar[["Z.var"]][mypar[["Z.K"]]] else mypar[["Z.var"]]
                            ),nrow=nv)
  }

  if(model[["sociality"]] && is.null(mypar[["sociality"]])){
    if(is.null(mypar[["sociality.var"]]))
      mypar[["sociality.var"]]<-with(prior,sociality.var*sociality.var.df/rchisq(1,sociality.var.df))
    model[["sociality"]]<-rnorm(nv,0,sqrt(mypar[["sociality.var"]]))
  }

  if(model[["sender"]] && is.null(mypar[["sender"]])){
    if(is.null(mypar[["sender.var"]]))
      mypar[["sender.var"]]<-with(prior,sender.var*sender.var.df/rchisq(1,sender.var.df))
    mypar[["sender"]]<-rnorm(nv,0,sqrt(mypar[["sender.var"]]))
  }
  
  if(model[["receiver"]] && is.null(mypar[["receiver"]])){
    if(is.null(mypar[["receiver.var"]]))
      mypar[["receiver.var"]]<-with(prior,receiver.var*receiver.var.df/rchisq(1,receiver.var.df))
    mypar[["receiver"]]<-rnorm(nv,0,sqrt(mypar[["receiver.var"]]))
  }

  if(model[["dispersion"]] && is.null(mypar[["dispersion"]])){
    mypar[["dispersion"]]<-with(prior,dispersion*dispersion.df/rchisq(1,dispersion.df))
  }

  eta<-ergmm.eta(model,mypar)

  sm<-rsm.fs[[model[["familyID"]]]](eta,dispersion=mypar[["dispersion"]],fam.par=model[["fam.par"]])
  if(!has.loops(model[["Yg"]]))
    diag(sm)<-0
  net<-with(model,as.network.matrix(sm,matrix.type="adjacency",
                                    directed=is.directed(Yg),
                                    loops=has.loops(Yg)))
  if(!is.null(model[["response"]]))
    net<-set.edge.value(net,model[["response"]],sm)

  attr(net,"ergmm.par")<-mypar

  net
}
