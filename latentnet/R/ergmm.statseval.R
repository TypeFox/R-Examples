ergmm.statseval <- function (mcmc.out, model, start, prior, control,Z.ref=NULL,Z.K.ref=NULL){
  if(control[["verbose"]]) cat("Post-processing the MCMC output:\n")
   
  Yg<- model[["Yg"]]
  G <- model[["G"]]
  d <- model[["d"]]
  n <- network.size(Yg)
  p <- model[["p"]]
  family <- model[["family"]]
  sample.size <- control[["sample.size"]]

  if(!is.null(Z.ref)){
    control[["Z.ref"]]<-.scale.ergmm.model(model, list(Z=Z.ref))[["Z"]]
  }
  if(!is.null(Z.K.ref)){
    control[["Z.K.ref"]]<-Z.K.ref
  }
  
  ## Start putting together the output structure.
  if(is.null(mcmc.out)) mcmc.out<-list()
  mcmc.out[["model"]]<-model
  mcmc.out[["prior"]]<-prior
  mcmc.out[["control"]]<-control
  mcmc.out[["start"]]<-start
  class(mcmc.out)<-"ergmm"


  if(control[["tofit"]][["klswitch"]]) mcmc.out<-labelswitch.sample.ergmm(mcmc.out)
  if(control[["tofit"]][["pmode"]]) mcmc.out<-add.mcmc.pmode.pmode.ergmm(mcmc.out)
  if(control[["tofit"]][["mle"]]) mcmc.out<-add.mcmc.mle.mle.ergmm(mcmc.out)
  if(control[["tofit"]][["mkl"]]) mcmc.out<-add.mkl.pos.ergmm(mcmc.out)
  if(control[["tofit"]][["mkl.mbc"]]) mcmc.out<-add.mkl.mbc.ergmm(mcmc.out)
  if(control[["tofit"]][["procrustes"]]) mcmc.out<-.procr.sample.ergmm(mcmc.out)

  class(mcmc.out)<-"ergmm"
  return(mcmc.out)
}

statsreeval.ergmm<-function(x,Z.ref=NULL,Z.K.ref=NULL,rerun=FALSE){
  if(!is.null(Z.ref)){
    control[["Z.ref"]]<-.scale.ergmm.model(x[["model"]], list(Z=Z.ref))
  }
  
  if(!is.null(Z.K.ref)){
    x[["control"]][["Z.K.ref"]]<-Z.K.ref
  }
  x<-labelswitch.sample.ergmm(x)
  if(rerun) x<-add.mcmc.pmode.pmode.ergmm(x)
  if(rerun) x<-add.mcmc.mle.mle.ergmm(x)
  if(rerun) x<-add.mkl.pos.ergmm(x)
  x<-add.mkl.mbc.ergmm(x)
  x<-.procr.sample.ergmm(x)
  x
}

find.mkl<-function(model,sample,control){
  if(control[["verbose"]]>1) cat("Evaluating matrix of predicted dyad values and finding initial value... ")
  EY<-post.predict.C(model,sample,control,TRUE)
  EY[!observed.dyads(model[["Yg"]])]<-NA
  if(control[["verbose"]]>1) cat("Finished.\nMaximizing...")

  mkl<-sample[[attr(EY,"s.MKL")]]
  model[["Ym"]]<-EY
  for(i in 1:control[["mle.maxit"]]){
    if(control[["verbose"]]>1) cat(i,"")
    mkl.old<-mkl
    mkl<-find.mle(model,start=mkl,control=control,mllk=FALSE)
    if(is.null(mkl)) stop("MKL failed!")
    if(isTRUE(all.equal(mkl.old,mkl))) break
  }
  mkl
}

nullapply<-function(X,margin,FUN,...){
  if(is.null(X)) return(X)
  else apply(X,margin,FUN,...)
}

add.mkl.mbc.ergmm<-function(x,Z.K.ref=best.avail.Z.K.ref.ergmm(x)){
  if(is.null(x[["mkl"]]) || x[["model"]][["d"]]<=0 || x[["model"]][["G"]]<=0){
    if(x[["control"]][["verbose"]]) cat("MKL MBC is not available or non-latent-cluster model.\n")
    return(x)
  }else Z<-x[["mkl"]][["Z"]]

  if(x[["control"]][["verbose"]]) cat("Fitting MBC conditional on MKL locations... ")
  x[["mkl"]][["mbc"]] <- {
    if(x[["control"]][["kl.threads"]]==1)
      bayesmbc(x[["model"]][["G"]],Z,x[["prior"]],Z.K.ref,verbose=x[["control"]][["verbose"]])[["pmean"]]
    else
      bayesmbc.snowFT(x[["control"]][["kl.threads"]],x[["model"]][["G"]],Z,x[["prior"]],Z.K.ref,verbose=x[["control"]][["verbose"]])[["pmean"]]
  }
  if(x[["control"]][["verbose"]]) cat("Finished.\n")
  x
}

add.mcmc.mle.mle.ergmm<-function(x,Z.ref=best.avail.Z.ref.ergmm(x)){
  if(!is.null(x[["start"]])){
    if(x[["control"]][["verbose"]]) cat("Using the conditional posterior mode to seed an MLE fit... ")
    mle1<-find.mle.loop(x[["model"]],x[["start"]],control=x[["control"]])
    if(x[["control"]][["verbose"]]) cat(" Finished.\n")
  }else mle1<-list(lpY=-Inf)
  
  if(!is.null(x[["mcmc.mle"]])){
    ## Use the iteration with the highest probability to seed another
    ## shot at the MLE.
    
    if(x[["control"]][["verbose"]]) cat("Using the highest-likelihood iteration to seed another MLE fit... ")
    mle2 <- find.mle.loop(x[["model"]],x[["mcmc.mle"]],control=x[["control"]])
    mle2<-.scale.ergmm.model(x[["model"]],mle2)
    if(x[["control"]][["verbose"]]) cat("Finished.\n")
  }
  else mle2<-list(lpY=-Inf)
  if(is.null(mle2)) mle2<-list(lpY=-Inf)
  
  
  if(mle2[["lpY"]]>mle1[["lpY"]]) x[["mle"]]<-mle2 else x[["mle"]]<-mle1
  
  x[["mle"]]<-.scale.ergmm.model(x[["model"]],x[["mle"]])

  if(x[["model"]][["d"]] && "rotation" %in% latent.effect.invariances[[x[["model"]][["familyID"]]]]){
    x[["mle"]][["Z"]]<-x[["mle"]][["Z"]]%*%.procr(x[["model"]],Z.ref,x[["mle"]][["Z"]])
  }
  
  x
}

find.mle.loop<-function(model,start,control){
  mle<-start
  for(i in 1:control[["mle.maxit"]]){
    mle.old<-mle
    mle<-find.mle(model,mle,control=control)
    if(all.equal(mle.old,mle)[1]==TRUE) break
  }
  mle
}


add.mcmc.pmode.pmode.ergmm<-function(x,Z.ref=best.avail.Z.ref.ergmm(x)){
  if(!is.null(x[["start"]])){
    if(x[["control"]][["verbose"]]) cat("Double-checking conditional posterior mode estimate... ")
    pmode2<-find.pmode.loop(x[["model"]],x[["start"]],prior=x[["prior"]],control=x[["control"]])
    if(x[["control"]][["verbose"]]) cat("Finished.\n")
  }else pmode2<-list(mlp=-Inf)
  
  if(!is.null(x[["mcmc.pmode"]])){
    if(x[["control"]][["verbose"]]) cat("Using MCMC posterior mode to seed another conditional posterior mode fit... ")
    x[["pmode"]]<-find.pmode.loop(x[["model"]],x[["mcmc.pmode"]],prior=x[["prior"]],control=x[["control"]])
    if(x[["control"]][["verbose"]]) cat("Finished.\n")
    
  }
  if(is.null(x[["pmode"]]) || pmode2[["mlp"]]>x[["pmode"]][["mlp"]]) x[["pmode"]]<-pmode2
    
  if(x[["model"]][["d"]]>0 && !is.null(x[["pmode"]])){
    x[["pmode"]]<-.scale.ergmm.model(x[["model"]], x[["pmode"]])
    P<-.procr(x[["model"]],x[["pmode"]][["Z"]],Z.ref)
    x[["pmode"]][["Z"]]<-x[["pmode"]][["Z"]]%*%P
    if(!is.null(x[["pmode"]][["Z.mean"]]))
      x[["pmode"]][["Z.mean"]]<-x[["pmode"]][["Z.mean"]]%*%P
  }
  x
}

find.pmode.loop<-function(model,start,prior,control){
  pmode<-start
  for(i in 1:control[["mle.maxit"]]){
    if(control[["verbose"]]>1) cat(i,"")
    pmode.old<-pmode
    pmode<-find.mpe(model,pmode,prior=prior,given=list(Z.K=pmode[["Z.K"]]),control=control)
    if(model[["G"]]>1) pmode[["Z.K"]]<-mbc.VII.EM(model[["G"]],pmode[["Z"]],resume=list(Z.mean=pmode[["Z.mean"]],Z.var=pmode[["Z.var"]],Z.pK=pmode[["Z.pK"]]))[["Z.K"]]
    if(all.equal(pmode.old,pmode)[1]==TRUE) break
  }
  pmode
}

add.mkl.pos.ergmm<-function(x, Z.ref=best.avail.Z.ref.ergmm(x)){
  if(!is.null(x[["sample"]])){
    if(x[["control"]][["verbose"]]) cat("Fitting the MKL locations... ")
    x[["mkl"]]<-find.mkl(x[["model"]],x[["sample"]],x[["control"]])
  }
  if(!is.null(x[["mkl"]][["Z"]])) x[["mkl"]]<-.scale.ergmm.model(x[["model"]],x[["mkl"]])
  if(!is.null(x[["mkl"]][["Z"]])) x[["mkl"]][["Z"]]<-x[["mkl"]][["Z"]]%*%.procr(x[["model"]],x[["mkl"]][["Z"]],Z.ref)
  if(x[["control"]][["verbose"]]) cat("Finished.\n")
  x
}

.procr.sample.ergmm<-function(x,Z.ref=best.avail.Z.ref.ergmm(x),...){
  if(!is.null(x[["sample"]]) && x[["model"]][["d"]]>0 && "rotation" %in% latent.effect.invariances[[x[["model"]][["latentID"]]]] && "reflection" %in% latent.effect.invariances[[x[["model"]][["latentID"]]]]){
    if(x[["control"]][["verbose"]]) cat("Performing Procrustes transformation... ")
    x[["sample"]]<-.procrustes.Z.mean.C(x[["sample"]],Z.ref,verbose=x[["control"]][["verbose"]])
    if(x[["control"]][["verbose"]]) cat("Finished.\n")
  }
  x
}

labelswitch.sample.ergmm<-function(x,Z.K.ref=best.avail.Z.K.ref.ergmm(x)){
  if(!is.null(x[["sample"]]) && x[["model"]][["G"]]>1){
    if(x[["control"]][["verbose"]]) cat("Performing label-switching... ")
    Q.start<-switch.Q.K(Z.K.ref,x[["model"]][["G"]])
    x[["sample"]] <- {
      if(x[["control"]][["kl.threads"]]==1)
        klswitch.C(Q.start,x[["sample"]],verbose=x[["control"]][["verbose"]])
      else
        klswitch.snowFT(x[["control"]][["kl.threads"]],Q.start,x[["sample"]],verbose=x[["control"]][["verbose"]])
    }
    if(x[["control"]][["verbose"]]) cat("Finished.\n")
  }
  x
}

best.avail.Z.ref.ergmm<-function(x){
  if(x[["model"]][["d"]]==0) return(NULL)
  
  if(!is.null(x[["control"]][["Z.ref"]])) return(x[["control"]][["Z.ref"]])
  if(!is.null(x[["mkl"]][["Z"]])) return(x[["mkl"]][["Z"]])
  if(!is.null(x[["pmode"]][["Z"]])) return(x[["pmode"]][["Z"]])
  if(!is.null(x[["mcmc.pmode"]][["Z"]])) return(x[["mcmc.pmode"]][["Z"]])
  if(!is.null(x[["mle"]][["Z"]])) return(x[["mle"]][["Z"]])
  if(!is.null(x[["mcmc.mle"]][["Z"]])) return(x[["mcmc.mle"]][["Z"]])
  if(!is.null(x[["start"]][["Z"]])) return(x[["start"]][["Z"]])
}

best.avail.Z.K.ref.ergmm<-function(x){
  if(x[["model"]][["G"]]==0) return(NULL)
  
  if(!is.null(x[["control"]][["Z.K.ref"]])) return(x[["control"]][["Z.K.ref"]])
  if(!is.null(attr(x[["sample"]],"Q"))) return(apply(attr(x[["sample"]],"Q"),1,which.max))
  if(!is.null(x[["pmode"]][["Z.K"]])) return(x[["pmode"]][["Z.K"]])
  if(!is.null(x[["mkl"]][["mbc"]][["Z.K"]])) return(x[["mkl"]][["mbc"]][["Z.K"]])
  if(!is.null(x[["mcmc.pmode"]][["Z.K"]])) return(x[["mcmc.pmode"]][["Z.K"]])
  if(!is.null(x[["start"]][["Z.K"]])) return(x[["start"]][["Z.K"]])
  
  return(mbc.VII.EM(x[["model"]][["G"]],best.avail.Z.ref.ergmm(x))[["Z.K"]])
}

.scale.ergmm.model<-function(x,theta,...){
  extraneous.argcheck(...)
  if(!is.null(theta[["Z"]])){
    if("translation" %in% latent.effect.invariances[[x[["latentID"]]]]){
      if(!is.null(theta[["Z.mean"]])){
        Z.center<-colMeans(theta[["Z.mean"]])
        theta[["Z.mean"]]<-sweep(theta[["Z.mean"]],2,Z.center,check.margin=FALSE)
      }else{
        Z.center<-colMeans(theta[["Z"]])
      }
      theta[["Z"]]<-sweep(theta[["Z"]],2,Z.center,check.margin=FALSE)
    }
    if("scaling" %in% latent.effect.invariances[[x[["latentID"]]]]){
      theta[["Z"]]<-scale(theta[["Z"]],center=FALSE)
      theta[["Z.mean"]]<-scale(theta[["Z.mean"]],center=FALSE)
      ## Fixme: Should Z.var be scaled as well? How?
    }
  }
  if(x[["intercept"]]){
    shift<-0

    if(!is.null(theta[["sociality"]])){
      shift.mul<-1+is.directed(x[["Yg"]])
      shift<-shift+mean(theta[["sociality"]])*shift.mul
      theta[["sociality"]]<-theta[["sociality"]]-mean(theta[["sociality"]])
    }
    
    for(eff in c("sender","receiver")){
      if(!is.null(theta[[eff]])){
        shift<-shift+mean(theta[[eff]])
        theta[[eff]]<-theta[[eff]]-mean(theta[[eff]])
      }
    }

    theta[["beta"]][1]<-theta[["beta"]][1]+shift
  }
  theta
}
