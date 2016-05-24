predict.ergmm<-function(object,...,type="post"){
  if(class(type)=="list"){
    type<-type
  }else if(type=="start"){
    type<-object[["start"]]
  }else if(type=="mle"){
    type<-object[["mle"]]
  }else if(type=="pmean"){
    type<-summary(object,point.est=c("pmean"))[["pmean"]]
  }else if(type=="mkl"){
    type<-object[["mkl"]]
  }else if(type=="pmode"){
    type<-object[["pmode"]]
  }else if(type=="post"){
    return(with(object,post.predict.C(model,sample,control)))
  }else if(is.numeric(type) && round(type)==type){
    type<-object[["sample"]][[type]]
  }else stop("Invalid parameter structure.")

  ergmm.EY(object[["model"]],type,NA.unobserved=FALSE)
}


post.predict.C<-function(model,sample,control,MKL=FALSE){
  n<-network.size(model[["Yg"]])
  d<-model[["d"]]
  p<-model[["p"]]

  ## Figure out the design matrix.
  observed<-observed.dyads(model[["Yg"]])

  if((observed==(diag(n)==0) && is.directed(model[["Yg"]])) ||
     (observed==lower.tri(diag(n)) && !is.directed(model[["Yg"]])))
    observed<-NULL

  ret<-.C("post_pred_wrapper",
          S = as.integer(control[["sample.size"]]),
          
          n = as.integer(n),
          p = as.integer(p),
          d = as.integer(d),
          latent=as.integer(NVL(model[["latentID"]],0)),
          family=as.integer(NVL(model[["familyID"]],0)),
          res=as.integer(with(model,c(sender,receiver,sociality,dispersion))),
          
          dir=as.integer(is.directed(model[["Yg"]])),
          iconsts=as.integer(model[["iconsts"]]),
          dconsts=as.double(model[["dconsts"]]),

          
          X=as.double(unlist(model[["X"]])),
          
          Z = as.double(sample[["Z"]]),
          beta = as.double(sample[["beta"]]), # coef
          sender = if(model[["sociality"]]) as.double(sample[["sociality"]]) else as.double(sample[["sender"]]),
          receiver = as.double(sample[["receiver"]]),
          dispersion = as.double(sample[["dispersion"]]),
          
          observed=as.integer(NVL(observed,-1)),
          
          EY=double(n*n),
          s.MKL=if(MKL) as.integer(TRUE) else as.integer(FALSE),
          verbose=as.integer(control[["verbose"]]),

          PACKAGE="latentnet")
  EY<-array(ret[["EY"]],dim=c(1,n,n))[1,,] 
  if(MKL) attr(EY,"s.MKL")<-ret[["s.MKL"]]+1 # C counts from 0; R counts from 1
  EY
}

post.predict.R<-function(model,sample,control,MKL=FALSE){
  EY.f<-EY.fs[[model[["familyID"]]]]
  EY<-matrix(0,network.size(model[["Yg"]]),network.size(model[["Yg"]]))
  for(i in 1:control[["sample.size"]]){
    state<-sample[[i]]
    eta<-ergmm.eta(model,state)
    EY<-EY+EY.f(eta,model[["fam.par"]])
  }
  EY<-EY/control[["sample.size"]]

  if(MKL){
    min.MKL<-NA
    min.dev<-Inf
    model[["Ym"]]<-EY
    model[["Ym"]][!observed.dyads(model[["Yg"]])]<-NA
    for(i in 1:control[["sample.size"]]){
      state<-sample[[i]]
      dev<--ergmm.lpY(model,state,up.to.const=TRUE)
      if(dev<min.dev){
        min.dev<-dev
        min.MKL<-i
      }
    }
    attr(EY,"s.MKL")<-min.MKL
  }
  
  EY
}
