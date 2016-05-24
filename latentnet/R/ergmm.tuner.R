nterms.model<-function(model){
  (model[["p"]]
   +(if(model[["d"]]) 1 else 0)
   +(if(model[["sender"]]) nrow(model[["beta.eff.sender"]]) else 0)
   +(if(model[["receiver"]]) nrow(model[["beta.eff.receiver"]]) else 0)
   +(if(model[["sociality"]]) nrow(model[["beta.eff.sociality"]]) else 0)
   +(if(model[["dispersion"]]) 1 else 0)
   )
}

get.init.deltas<-function(model, control){
  nterms<-nterms.model(model)

  ## If proposal coefficient matrix is given, leave it alone.
  if(is.matrix(control[["group.deltas"]])){
    if(any(dim(control[["group.deltas"]])!=nterms))
      stop(paste("Incorrect proposal coefficient matrix size: (",
                 paste(dim(control[["group.deltas"]]),collapse=", "),"), ",
                 "while the model calls for ",nterms,".",sep=""))
    return(control)
  }
  
  ## If a vector of appropriate length is given, use a diagonal matrix
  if(length(control[["group.deltas"]])==nterms){
    control[["group.deltas"]]<-diag(control[["group.deltas"]],nrow=nterms)
    return(control)
  }
  
  ## If a scalar is given, construct a diagonal matrix that's in the ballpark.
  if(length(control[["group.deltas"]])==1){
    group.deltas.scale<-control[["group.deltas"]]
    control[["group.deltas"]]<-1/sapply(1:model[["p"]],function(i) sqrt(mean((model[["X"]][[i]][observed.dyads(model[["Yg"]])])^2)))
    if(model[["d"]]) control[["group.deltas"]]<-c(control[["group.deltas"]], 0.05)
    control[["group.deltas"]]<-c(control[["group.deltas"]],rep(1/model[["p"]],nterms-length(control[["group.deltas"]])))
    control[["group.deltas"]]<-diag(group.deltas.scale*control[["group.deltas"]]*2/(1+nterms),nrow=nterms)
  }

  control[["group.acf.adjust"]]<-rep(1,nrow(control[["group.deltas"]]))
  
  control
}

get.sample.deltas<-function(model,sample,control){
  ## Convert the "stacked" list of draws into a list of threads:
  samples<-if(control[["threads"]]>1) unstack.ergmm.par.list(sample) else list(sample)

  Z.rate<-0
  beta.rate<-0
  cov.beta.ext<-0

  acf.adjust<-rep(1,length(control[["group.acf.adjust"]]))
  
  for(thread in 1:control[["threads"]]){
    sample<-samples[[thread]]
    use.draws<-ceiling(length(sample)*control[["pilot.discard.first"]]):(length(sample))
    Z.rate<-if(!is.null(sample[["Z.rate"]]))Z.rate+mean(sample[["Z.rate"]][use.draws])/control[["threads"]] else 0.5
    beta.rate<-beta.rate+mean(sample[["beta.rate"]][use.draws])/control[["threads"]]

    beta.ext<-get.beta.ext(model,sample)
    beta.ext.ar.eff<-1/(1-apply(beta.ext,2,function(b) ar(b,order.max=1,aic=FALSE)$ar))
    acf.adjust<-acf.adjust*beta.ext.ar.eff/exp(mean(log(beta.ext.ar.eff)))

    ## Note that this only measures the covariances _within_ the thread.    
    cov.beta.ext<-cov.beta.ext+cov(beta.ext)/control[["threads"]]
  }

  control[["group.acf.adjust"]]<-control[["group.acf.adjust"]]*acf.adjust^(1/control[["threads"]])
  
  cov.beta.ext<-t(cov.beta.ext*control[["group.acf.adjust"]])*control[["group.acf.adjust"]]

  if(model[["d"]]) control[["Z.delta"]]<-control[["Z.delta"]]*Z.rate/control[["target.acc.rate"]]
  if(model[["sender"]] || model[["receiver"]] || model[["sociality"]]) control[["RE.delta"]]<-control[["RE.delta"]]*Z.rate/control[["target.acc.rate"]]
  control[["pilot.factor"]]<-control[["pilot.factor"]]*beta.rate/control[["target.acc.rate"]]
  
  ## Take the Choletsky decomposition of the empirical covariance matrix.
  control[["group.deltas"]]<-try(chol(cov.beta.ext)*control[["pilot.factor"]])
  if(inherits(control[["group.deltas"]],"try-error")) stop("Looks like a pilot run did not mix at all (practically no proposals were accepted). Try using a smaller starting proposal variance.")
  control
}

## Compute the empirical covariance of coefficients, latent scale, and random effect means.
get.beta.ext<-function(model,sample){
  n<-network.size(model[["Yg"]])
  ## Construct the "extended" beta: not just the coefficients, but also the scale of latent space and all random effects.
  betas<-cbind(if(model[["p"]]) sample[["beta"]], # covariate coefs
               if(model[["d"]]) log(apply(sqrt(apply(sample[["Z"]]^2,1:2,sum)),1,mean)), # scale of Z
               if(model[["sender"]]) tcrossprod(sample[["sender"]],model[["beta.eff.sender"]])/n, # sender eff.
               if(model[["receiver"]]) tcrossprod(sample[["receiver"]],model[["beta.eff.receiver"]])/n, # receiver eff.
               if(model[["sociality"]]) tcrossprod(sample[["sociality"]],model[["beta.eff.sociality"]]/n), # sociality eff.
               if(model[["dispersion"]]) log(sample[["dispersion"]]) # dispersion scale.
        )
  sweep(betas,2,apply(betas,2,mean))
}

backoff.check<-function(model,sample,control){
  bt <-control[["backoff.threshold"]]
  if(bt>0.5) bt<-1-bt
  
  backoff<-0
  if((model[["d"]] || model[["sender"]] || model[["receiver"]] || model[["sociality"]])){
    if(mean(sample[["Z.rate"]] < bt)){
      control[["Z.delta"]]<-control[["Z.delta"]]*control[["backoff.factor"]]
      control[["RE.delta"]]<-control[["RE.delta"]]*control[["backoff.factor"]]
      backoff<--1
    }
    if(mean(sample[["Z.rate"]] > 1-bt)){
      control[["Z.delta"]]<-control[["Z.delta"]]/control[["backoff.factor"]]
      control[["RE.delta"]]<-control[["RE.delta"]]/control[["backoff.factor"]]
      backoff<-+1
    }
  }
  
  if(mean(sample[["beta.rate"]]) < bt) {
    control[["group.deltas"]]<-control[["group.deltas"]]*control[["backoff.factor"]]
    control[["pilot.factor"]]<-control[["pilot.factor"]]*control[["backoff.factor"]]
    backoff<--1
  }
  if(mean(sample[["beta.rate"]]) > 1-bt) {
    control[["group.deltas"]]<-control[["group.deltas"]]/control[["backoff.factor"]]
    control[["pilot.factor"]]<-control[["pilot.factor"]]/control[["backoff.factor"]]
    backoff<-+1
  }
  if(backoff){
    backoff.str<-paste("Backing off: too",if(backoff<0) "few" else "many", "acceptances. If you see this message several times in a row, use a longer burnin.",sep=" ")
    if(control[["verbose"]]) cat(paste(backoff.str,"\n",sep=""))
    else warning(backoff.str)
  }
  control[["backedoff"]]<-backoff
  control
}
