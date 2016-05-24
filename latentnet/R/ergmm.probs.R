not.given<-function(name,theta,given){
  is.null(given[[name]]) && !is.null(theta[[name]])
}

PRIOR_NAMES<-list(sender=c("sender.var"),
                  receiver=c("receiver.var"),
                  sociality=c("sociality.var"),
                  Z=c("Z.var"))

lp.works<-function(name,theta,given){
  not.given(name, theta, given) && all(PRIOR_NAMES[[name]]%in%names(.merge.lists(theta,given)))
}

bipartite.augment<-function(m){
  actors<-dim(m)[1]
  events<-dim(m)[2]
  rbind(matrix(NA,actors,actors+events),
        cbind(t(m),matrix(NA,events,events)))
}

getYm<-function(Yg,response=NULL){
  m <-
    if(is.null(response)){
      as.matrix.network(Yg, matrix.type="adjacency")
    }else{
      if(is.matrix(response)) response
      else as.matrix.network(Yg, response, matrix.type="adjacency")
    }

  ## If bipartite, augment into a matrix of the form
  ##  N    m
  ## t(m)  N
  ## where N is a matrix of NAs of appropriate dimension.
  if(is.bipartite(Yg)) m<-bipartite.augment(m)

  m[!observed.dyads(Yg)]<-NA
  m
}

ergmm.eta<-function(model,theta){
  n<-network.size(model[["Yg"]])
  dir<-is.directed(model[["Yg"]])

  eta<-matrix(0,n,n)
  
  if(!is.null(theta[["Z"]]))
    eta<-eta+latent.effect.fs[[model[["latentID"]]]](theta[["Z"]])

  if(!is.null(theta[["beta"]]))
    for(k in 1:length(theta[["beta"]]))
      eta<-eta+theta[["beta"]][k]*model[["X"]][[k]]

  if(!is.null(theta[["sociality"]])){
    eta<-eta+theta[["sociality"]]
    eta<-t(t(eta)+theta[["sociality"]])
  }
  
  if(!is.null(theta[["sender"]]))
    eta<-eta+theta[["sender"]]

  if(!is.null(theta[["receiver"]]))
    eta<-t(t(eta)+theta[["receiver"]])

  return(eta)
}

ergmm.EY<-function(model,theta,NA.unobserved=TRUE){
  eta<-ergmm.eta(model,theta)
  if(NA.unobserved) eta[!observed.dyads(model[["Yg"]])]<-NA
  EY.fs[[model[["familyID"]]]](eta,dispersion=theta[["dispersion"]],fam.par=model[["fam.par"]])
}

ergmm.lpY<-function(model,theta,given=list(),up.to.const=FALSE){
  theta<-.merge.lists(theta,given)
  Yg<-model[["Yg"]]
  Ym<-model[["Ym"]]
  n<-network.size(Yg)
  eta<-ergmm.eta(model,theta)
  obs<-observed.dyads(Yg)
  lpY<-if(up.to.const) lpYc.fs[[model[["familyID"]]]](Ym[obs],eta[obs],dispersion=theta[["dispersion"]],model[["fam.par"]]) else lpY.fs[[model[["familyID"]]]](Ym[obs],eta[obs],dispersion=theta[["dispersion"]],model[["fam.par"]])
  return(sum(lpY))
}

ergmm.lpY.grad<-function(model,theta,given=list()){
  theta<-.merge.lists(theta,given)
  n<-network.size(model[["Yg"]])
  obs<-observed.dyads(model[["Yg"]])
  eta<-ergmm.eta(model,theta)
  
  dlpY.deta <- dlpY.deta.fs[[model[["familyID"]]]](model[["Ym"]],eta,dispersion=theta[["dispersion"]],model[["fam.par"]])
  dlpY.deta[!obs] <- 0

  grad<-list()
  
  if(not.given("beta",theta,given)) grad[["beta"]] <- sapply(1:length(theta[["beta"]]),function(k) sum((dlpY.deta*model[["X"]][[k]])[obs]))

  if(not.given("Z",theta,given)){
    grad[["Z"]]<-dlpY.dZ.fs[[model[["latentID"]]]](theta[["Z"]],dlpY.deta)
  }

  if(not.given("sociality",theta,given))
    grad[["sociality"]] <- sapply(1:n,function(i) sum(dlpY.deta[i,][obs[i,]])+sum(dlpY.deta[,i][obs[,i]]))
  else{
    if(not.given("sender",theta,given)) grad[["sender"]] <- sapply(1:n,function(i) sum(dlpY.deta[i,][obs[i,]]))
    if(not.given("receiver",theta,given)) grad[["receiver"]] <-  sapply(1:n,function(i) sum(dlpY.deta[,i][obs[,i]]))
  }
  
  if(not.given("dispersion",theta,given)){
    grad[["dispersion"]] <- dlpY.ddispersion.fs[[model[["familyID"]]]](model[["Ym"]],eta,dispersion=theta[["dispersion"]],model[["fam.par"]])
  }
  
  grad
}
  
ergmm.lpY.C<-function(model,theta){
  Y <- model[["Ym"]]
  Y[is.na(Y)] <- 0
  n <- network.size(model[["Yg"]])

  ## Figure out the design matrix.
  observed<-observed.dyads(model[["Yg"]])

  if((observed==(diag(n)==0) && is.directed(model[["Yg"]])) ||
     (observed==lower.tri(diag(n)) && !is.directed(model[["Yg"]])))
    observed<-NULL
  
  ## Sanity checks: the following block of code checks that all dimensionalities and
  ## dimensions are correct, and those optional parameters that are required by the presence
  ## of other optional parameters are present.
  
  for(i in 1:model[["p"]])
    if(!all(dim(model[["X"]][[i]])==c(n,n))) stop("Incorrect size for covariate matrices.")

  if(!is.null(theta[["Z"]])){
    if(!all(dim(theta[["Z"]])==c(n,model[["d"]]))) stop("Incorrect size for the latent positions.")
  }  
  if(length(theta[["beta"]])!=model[["p"]]) stop("Incorrect length for the beta vector.")

  if(!is.null(theta[["sociality"]])){
    if(length(theta[["sociality"]])!=n) stop("Incorrect length for the vector of sociality effects.")
  }
  if(!is.null(theta[["sender"]])){
    if(length(theta[["sender"]])!=n) stop("Incorrect length for the vector of sender effects.")
  }
  if(!is.null(theta[["receiver"]])){
    if(length(theta[["receiver"]])!=n) stop("Incorrect length for the vector of receiver effects.")
  }
  ## End Sanity checks.
  
  ret <- .C("ERGMM_lp_Y_wrapper",
            n=as.integer(n), p=as.integer(model[["p"]]),
            d=as.integer(model[["d"]]), latent=as.integer(NVL(model[["latentID"]],0)), family=as.integer(NVL(model[["familyID"]],0)), res=as.integer(with(model,c(sender,receiver,sociality))),
            
            dir=as.integer(is.directed(model[["Yg"]])),
            viY=as.integer(Y),
            vdY=as.double(Y),
            
            iconsts=as.integer(model[["iconsts"]]), dconsts=as.integer(model[["dconsts"]]),
            
            vX=as.double(unlist(model[["X"]])),
            Z=as.double(theta[["Z"]]),
            
            beta=as.double(theta[["beta"]]),

            sender=if(is.null(theta[["sociality"]])) as.double(theta[["sender"]]) else as.double(theta[["sociality"]]), receiver=as.double(theta[["receiver"]]), lock.RE=as.integer(!is.null(theta[["sociality"]])),
            dispersion=NVL(theta[["dispersion"]],0),
            
            observed=as.integer(NVL(observed,-1)),

            lpY=double(1),

            PACKAGE="latentnet")
  

  ret[["lpY"]]
}

observed.dyads<-function(Yg){
  observed.dyads<-get.network.attribute(Yg,"design")
  if(is.null(observed.dyads)){
    if(!is.bipartite(Yg))
      observed.dyads<-!is.na(as.matrix.network(Yg,matrix.type="adjacency"))
    else
      observed.dyads<-bipartite.augment(!is.na(as.matrix.network(Yg,matrix.type="adjacency")))
  }
  else{
    observed.dyads<-bipartite.augment(as.matrix.network(observed.dyads,matrix.type="adjacency")==0)
  }

  observed.dyads[is.na(observed.dyads)]<-FALSE
  if(!is.directed(Yg)) observed.dyads[upper.tri(observed.dyads)]<-FALSE
  if(!is.bipartite(Yg) && !has.loops(Yg)) diag(observed.dyads)<-FALSE
  
  observed.dyads
}

pack.optim<-function(theta,fit.vars=NULL){
  if(is.null(fit.vars))
    return(c(theta[["beta"]],theta[["Z"]],
             theta[["sender"]],theta[["receiver"]],theta[["sociality"]],
             theta[["Z.var"]],theta[["Z.mean"]],
             theta[["sender.var"]],theta[["receiver.var"]],theta[["sociality.var"]],
             theta[["dispersion"]]))
  else
    return(c(if(fit.vars[["beta"]])theta[["beta"]],
             if(fit.vars[["Z"]])theta[["Z"]],
             
             if(fit.vars[["sender"]])theta[["sender"]],
             if(fit.vars[["receiver"]])theta[["receiver"]],
             if(fit.vars[["sociality"]])theta[["sociality"]],

             if(fit.vars[["Z.var"]])theta[["Z.var"]],
             if(fit.vars[["Z.mean"]])theta[["Z.mean"]],
             
             if(fit.vars[["sender.var"]])theta[["sender.var"]],
             if(fit.vars[["receiver.var"]])theta[["receiver.var"]],
             if(fit.vars[["sociality.var"]])theta[["sociality.var"]],
             if(fit.vars[["dispersion"]])theta[["dispersion"]]))
}

reg.fit.vars<-function(fit.vars){
  for(name in ERGMM.PAR_VAR_NAMES)
    if(!(name %in% names(fit.vars)))fit.vars[[name]]<-FALSE
  fit.vars
}

inv.fit.vars<-function(fit.vars){
  for(name in ERGMM.PAR_VAR_NAMES)
    if(name %in% names(fit.vars))
      fit.vars[[name]]<-!fit.vars[[name]]
  fit.vars
}

FIT_ALL<-list(beta=TRUE,Z=TRUE,sender=TRUE,receiver=TRUE,sociality=TRUE,
              Z.var=TRUE,Z.mean=TRUE,
              sender.var=TRUE,receiver.var=TRUE,sociality.var=TRUE,dispersion=TRUE)

FIT_MLE<-list(beta=TRUE,Z=TRUE,sender=TRUE,receiver=TRUE,sociality=TRUE,dispersion=TRUE)

unpack.optim<-function(v,fit.vars,model){
  p<-model[["p"]]
  n<-network.size(model[["Yg"]])
  G<-model[["G"]]
  d<-model[["d"]]
  v.must.be<-with(fit.vars,
                  beta*p +
                  Z*n*d +
                  sender*n +
                  receiver*n +
                  sociality*n +
                  Z.var*(d>0)*max(1,G) +
                  Z.mean*G*d +
                  sender.var +
                  receiver.var +
                  sociality.var +
                  dispersion)
  if(length(v)!=v.must.be){
    stop(paste("Input vector wrong length: ", length(v),
               " but should be ",v.must.be,".",
               sep=""))
  }
  pos<-0
  ret<-list()
  if(fit.vars[["beta"]] && p>0){
    ret[["beta"]]<-v[pos+1:p]
    pos<-pos+p
  }

  if(fit.vars[["Z"]] && d>0){
    ret[["Z"]]<-matrix(v[pos+1:(n*d)],nrow=n,ncol=d)
    pos<-pos+n*d
  }

  if(fit.vars[["sender"]]){
    ret[["sender"]]<-v[pos+1:n]
    pos<-pos+n
  }

  if(fit.vars[["receiver"]]){
    ret[["receiver"]]<-v[pos+1:n]
    pos<-pos+n
  }

  if(fit.vars[["sociality"]]){
    ret[["sociality"]]<-v[pos+1:n]
    pos<-pos+n
  }

  if(fit.vars[["Z.var"]] && d>0){
    ret[["Z.var"]]<-v[pos+1:max(1,G)]
    pos<-pos+max(1,G)
  }
  
  if(fit.vars[["Z.mean"]] && d>0 && G>0){
    ret[["Z.mean"]]<-matrix(v[pos+1:(G*d)],nrow=G,ncol=d)
    pos<-pos+G*d
  }
  
  if(fit.vars[["sender.var"]]){
    ret[["sender.var"]]<-v[pos+1]
    pos<-pos+1
  }
  
  if(fit.vars[["receiver.var"]]){
    ret[["receiver.var"]]<-v[pos+1]
    pos<-pos+1
  }

  if(fit.vars[["sociality.var"]]){
    ret[["sociality.var"]]<-v[pos+1]
    pos<-pos+1
  }

  if(fit.vars[["dispersion"]]){
    ret[["dispersion"]]<-v[pos+1]
    pos<-pos+1
  }

  ret
}

mk.lp.optim.fs<-function(fit.vars,model,prior,given=list(),opt=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV","lpdispersion")){
  fit.vars<-reg.fit.vars(fit.vars)
  return(list(
              f=function(v){
                theta<-unpack.optim(v,fit.vars,model)
                ergmm.lp(model,theta,prior=prior,given=given,
                         opt=opt,
                         up.to.const=TRUE)
              },
              grad.f=function(v){
                theta<-unpack.optim(v,fit.vars,model)
                # note that ergmm.lp.grad doesn't take the up.to.const parameter
                gr<-ergmm.lp.grad(model,theta,prior=prior,given=given,
                                  opt=opt)
                pack.optim(gr,fit.vars)
              }
              )
         )
}

find.mle<-function(model,start,given=list(),control,
                     hessian=FALSE,mllk=TRUE){
  fit.vars<-list()
  for(name in ERGMM.PAR_LLK_NAMES)
    fit.vars[[name]]<-not.given(name,start,given)
  mpe<-find.mpe(model,start,given=given,control=control,
                hessian=hessian,mlp=mllk,opt="lpY",fit.vars=fit.vars)
  if(mllk) mpe[["lpY"]]<-mpe[["mlp"]]
  mpe
}



find.mpe<-function(model,start,given=list(),prior=list(),control,fit.vars=NULL,opt=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV","lpdispersion"),
                   hessian=FALSE,mlp=TRUE){
  if(is.null(fit.vars)){
    fit.vars<-list()
    for(name in names(start))
      fit.vars[[name]]<-not.given(name,start,given)
  }else{
    fit.vars<-reg.fit.vars(fit.vars)
    for(name in names(fit.vars))
      if(!not.given(name,start,given)) fit.vars[[name]]<-FALSE
    
  }
  
  fit.vars<-reg.fit.vars(fit.vars)

  optim.control<-list(fnscale=-1,
                      maxit=control[["mle.maxit"]],
                      trace=max(0,control[["verbose"]]-2))
  
  optim.fs<-mk.lp.optim.fs(fit.vars,model,prior=prior,given=given,opt=opt)
  
  start.vals<-pack.optim(start,fit.vars)

  p<-model[["p"]]
  n<-network.size(model[["Yg"]])
  G<-model[["G"]]
  d<-model[["d"]]
  
  vmpe <- ##try(
              optim(par=start.vals,fn=optim.fs[["f"]],gr=optim.fs[["grad.f"]],
                    method="L-BFGS-B",
                    lower=pack.optim(list(
                      beta=rep(-Inf,p),
                      Z=rep(-Inf,n*d),
                      sender=rep(-Inf,n),
                      receiver=rep(-Inf,n),
                      sociality=rep(-Inf,n),
                      Z.var=rep(sqrt(.Machine[["double.eps"]]),(d>0)*max(G,1)),
                      Z.mean=rep(-Inf,d*G),
                      sender.var=sqrt(.Machine[["double.eps"]]),
                      receiver.var=sqrt(.Machine[["double.eps"]]),
                      sociality.var=sqrt(.Machine[["double.eps"]]),
                      dispersion=sqrt(.Machine[["double.eps"]])),
                      fit.vars=fit.vars),
                    control=optim.control,hessian=hessian)
            ##  )

  if(inherits(vmpe,"try-error")) return(NULL)
  mpe<-unpack.optim(vmpe[["par"]],fit.vars,model)

  mpe<-.merge.lists(mpe,given)
  mpe[["Z.K"]]<-.merge.lists(start,given)[["Z.K"]]
  mpe[["Z.pK"]]<-if(!is.null(mpe[["Z.K"]])) tabulate(mpe[["Z.K"]])/n
  
  if(mlp)
    mpe[["mlp"]]<-ergmm.lp(model,mpe,prior=prior,given=given,opt=opt)
  
  if(hessian) mpe[["hessian"]]<-vmpe[["hessian"]]
  mpe
}

ergmm.lp<-function(model,theta,prior,given=list(),opt=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV"),up.to.const=FALSE){

  lpY<-if("lpY" %in% opt) ergmm.lpY(model,theta,
                                        given=given,up.to.const=up.to.const) else 0
  
  lpZ<-if("lpZ" %in% opt) ergmm.lpZ(theta,given=given) else 0
  lpRE<-if("lpRE" %in% opt) ergmm.lpRE(theta,given=given) else 0
  lpBeta<-if("lpBeta" %in% opt) ergmm.lpBeta(theta,prior,given=given) else 0
  lpREV<-if("lpREV" %in% opt) ergmm.lpREV(theta,prior,given=given) else 0
  lpLV<-if("lpLV" %in% opt) ergmm.lpLV(theta,prior,given=given) else 0
  lpdispersion<-if("lpdispersion" %in% opt) ergmm.lpdispersion(theta,prior,given=given) else 0
  
  lpAll<-lpY+lpZ+lpRE+lpBeta+lpREV+lpLV+lpdispersion
                                                      
  return(lpAll)
}

.merge.lists<-function(...){
  out<-list()
  for(l in list(...)){
    for(name in names(l))
      out[[name]]<-l[[name]]
    if(class(l)!="list")
      class(out)<-class(l)
  }

  out
}

.sum.lists<-function(...){
  out<-list()
  for(l in list(...)){
    for(name in names(l)){
      if(name %in% names(out))
        out[[name]]<-out[[name]]+l[[name]]
      else
        out[[name]]<-l[[name]]
    }
  }
  class(out)<-class(list(l))
  out
}

zero.list<-function(x){
  out<-x
  for(name in names(x)){
    d<-dim(x[[name]])
    if(is.null(d)) out[[name]]<-rep(0,length(x[[name]]))
    else out[[name]]<-array(0,dim=d)
  }
}

filter.list<-function(x,y){
  out<-x
  for(name in names(y)){
    if(!isTRUE(y)) out[[name]]<-NULL
  }
  out
}


cmp.lists<-function(x,y){
  out<-list()
  for(name in names(x))
    if(name %in% names(y))
      out[[name]]<-mean(abs(x[[name]]-y[[name]]))/mean(abs(x[[name]]))
  out
}

ergmm.lp.grad<-function(model,theta,prior,given=list(),opt=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV","lpdispersion")){
  
  grad<-.sum.lists(if("lpY" %in% opt) if(not.given("beta",theta,given)||
                                        not.given("Z",theta,given)||
                                        not.given("sender",theta,given)||
                                        not.given("receiver",theta,given)||
                                        not.given("sociality",theta,given)||
                                        not.given("dispersion",theta,given)) ergmm.lpY.grad(model,theta,given=given),
                  if("lpZ" %in% opt) ergmm.lpZ.grad(theta,given=given),
                  if("lpRE" %in% opt) ergmm.lpRE.grad(theta,given=given),
                  if("lpBeta" %in% opt) ergmm.lpBeta.grad(theta,prior,given=given),
                  if("lpREV" %in% opt) ergmm.lpREV.grad(theta,prior,given=given),
                  if("lpLV" %in% opt) ergmm.lpLV.grad(theta,prior,given=given),
                  if("lpdispersion" %in% opt) ergmm.lpdispersion.grad(theta,prior,given=given))
  
  grad
}

ergmm.lp.grad.approx<-function(which.vars,model,theta,prior,delta,given=list(),opt=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV","lpdispersion")){
  which.vars[["Z.K"]]<-FALSE
  which.vars<-reg.fit.vars(which.vars)
  for(var in names(which.vars)) if(!(var %in% names(theta))) which.vars[[var]]<-FALSE

  v<-pack.optim(theta,which.vars)
  dlpdv<-numeric(length(v))
  
  for(i in 1:length(v)){
    v.m<-v.p<-v
    v.m[i]<-v.m[i]-delta
    theta.m<-.merge.lists(theta,unpack.optim(v.m,which.vars,model))
    lp.m<-ergmm.lp(model,theta.m,prior,given=given,opt=opt)

    
    v.p[i]<-v.p[i]+delta
    theta.p<-.merge.lists(theta,unpack.optim(v.p,which.vars,model))
    lp.p<-ergmm.lp(model,theta.p,prior,given=given,opt=opt)

    dlpdv[i]<-(lp.p-lp.m)/(2*delta)
  }
  
  return(unpack.optim(dlpdv,which.vars,model))
}

ergmm.lpZ<-function(theta,given=list()){
  theta<-.merge.lists(theta,given)
  if(lp.works("Z",theta,given)){
    n<-dim(theta[["Z"]])[1]
    d<-dim(theta[["Z"]])[2]
    if(is.null(theta[["Z.K"]])){
      if(!is.null(theta[["Z.mean"]])) stop("Given cluster means without cluster assignments!")
      theta[["Z.K"]]<-rep(1,n)
      theta[["Z.mean"]]<-matrix(0,nrow=1,ncol=d)
    }
    sum(dnorm(theta[["Z"]],theta[["Z.mean"]][theta[["Z.K"]],],matrix(sqrt(theta[["Z.var"]][theta[["Z.K"]]]),nrow=n,ncol=d,byrow=FALSE),TRUE))
  }
  else 0
}

ergmm.lpZ.grad<-function(theta,given=list()){
  theta<-.merge.lists(theta,given)
  deriv<-list()
  if(lp.works("Z",theta,given)){
    n<-dim(theta[["Z"]])[1]
    d<-dim(theta[["Z"]])[2]
    G<-if(is.null(theta[["Z.K"]])) 0 else dim(theta[["Z.mean"]])[1]

    if(is.null(theta[["Z.K"]])){
      if(!is.null(theta[["Z.mean"]])) stop("Given cluster means without cluster assignments!")
      theta[["Z.K"]]<-rep(1,n)
      theta[["Z.mean"]]<-matrix(0,nrow=1,ncol=d)
    }
    
    Z.dev<-(theta[["Z"]]-theta[["Z.mean"]][theta[["Z.K"]],])/matrix(theta[["Z.var"]][theta[["Z.K"]]],nrow=n,ncol=d,byrow=FALSE)
    deriv[["Z"]]<--Z.dev
    if(not.given("Z.var",theta,given)) deriv[["Z.var"]]<-sapply(1:max(G,1),
                                                           function(g)
                                                           (sum(Z.dev[theta[["Z.K"]]==g,,drop=FALSE]^2)-d*sum(theta[["Z.K"]]==g)/theta[["Z.var"]][g])/2)
    
    if(not.given("Z.mean",theta,given) && G) deriv[["Z.mean"]]<-do.call(rbind,lapply(1:G,function(g) apply(Z.dev[theta[["Z.K"]]==g,,drop=FALSE],2,sum)))
  }
  deriv
}

ergmm.lpRE<-function(theta, given=list()){
  theta<-.merge.lists(theta,given)
  ({if(lp.works("sender",theta,given)) sum(dnorm(theta[["sender"]],0,sqrt(theta[["sender.var"]]),TRUE)) else 0}+
   {if(lp.works("receiver",theta,given)) sum(dnorm(theta[["receiver"]],0,sqrt(theta[["receiver.var"]]),TRUE)) else 0}+
   {if(lp.works("sociality",theta,given)) sum(dnorm(theta[["sociality"]],0,sqrt(theta[["sociality.var"]]),TRUE)) else 0})
}

ergmm.lpRE.grad<-function(theta,given=list()){
  theta<-.merge.lists(theta,given)
  deriv<-list()
  if(lp.works("sender",theta,given)) deriv[["sender"]]<--theta[["sender"]]/theta[["sender.var"]]
  if(not.given("sender.var",theta,given)) deriv[["sender.var"]]<-(sum(theta[["sender"]]^2)/theta[["sender.var"]]-length(theta[["sender"]]))/theta[["sender.var"]]/2

  if(lp.works("receiver",theta,given)) deriv[["receiver"]]<--theta[["receiver"]]/theta[["receiver.var"]]
  if(not.given("receiver.var",theta,given)) deriv[["receiver.var"]]<-(sum(theta[["receiver"]]^2)/theta[["receiver.var"]]-length(theta[["receiver"]]))/theta[["receiver.var"]]/2

  if(lp.works("sociality",theta,given)) deriv[["sociality"]]<--theta[["sociality"]]/theta[["sociality.var"]]
  if(not.given("sociality.var",theta,given)) deriv[["sociality.var"]]<-(sum(theta[["sociality"]]^2)/theta[["sociality.var"]]-length(theta[["sociality"]]))/theta[["sociality.var"]]/2
  
  deriv
}

ergmm.lpBeta<-function(theta,prior,given=list()){
  theta<-.merge.lists(theta,given)
  return({if(not.given("beta",theta,given)) sum(dnorm(theta[["beta"]],prior[["beta.mean"]],sqrt(prior[["beta.var"]]),TRUE)) else 0})
}

ergmm.lpBeta.grad<-function(theta,prior,given=list()){
  theta<-.merge.lists(theta,given)
  deriv<-list()
  if(not.given("beta",theta,given)) deriv[["beta"]]<--(theta[["beta"]]-prior[["beta.mean"]])/prior[["beta.var"]]

  deriv
}

dsclinvchisq<-function(x,df,scale=1,log=FALSE){
  if(log) dchisq(df*scale/x,df,log=TRUE)+log(df)+log(scale)-2*log(x)
  else dchisq(df*scale/x,df,log=FALSE)*df*scale/x/x
}

ergmm.lpREV<-function(theta,prior,given=list()){
  theta<-.merge.lists(theta,given)
  ({if(not.given("sender.var",theta,given)) dsclinvchisq(theta[["sender.var"]],prior[["sender.var.df"]],prior[["sender.var"]],TRUE) else 0}+
   {if(not.given("receiver.var",theta,given)) dsclinvchisq(theta[["receiver.var"]],prior[["receiver.var.df"]],prior[["receiver.var"]],TRUE) else 0}+
   {if(not.given("sociality.var",theta,given)) dsclinvchisq(theta[["sociality.var"]],prior[["sociality.var.df"]],prior[["sociality.var"]],TRUE) else 0})
}

ergmm.lpREV.grad<-function(theta,prior,given=list()){
  theta<-.merge.lists(theta,given)
  deriv<-list()
  if(not.given("sender.var",theta,given)) deriv[["sender.var"]]<-prior[["sender.var.df"]]*prior[["sender.var"]]/theta[["sender.var"]]^2/2-(prior[["sender.var.df"]]/2+1)/theta[["sender.var"]]
  if(not.given("receiver.var",theta,given)) deriv[["receiver.var"]]<-prior[["receiver.var.df"]]*prior[["receiver.var"]]/theta[["receiver.var"]]^2/2-(prior[["receiver.var.df"]]/2+1)/theta[["receiver.var"]]
  if(not.given("sociality.var",theta,given)) deriv[["sociality.var"]]<-prior[["sociality.var.df"]]*prior[["sociality.var"]]/theta[["sociality.var"]]^2/2-(prior[["sociality.var.df"]]/2+1)/theta[["sociality.var"]]
  
  deriv
}

ergmm.lpLV<-function(theta,prior,given=list()){
  theta<-.merge.lists(theta,given)
  ({if(not.given("Z.var",theta,given)) sum(dsclinvchisq(theta[["Z.var"]],prior[["Z.var.df"]],prior[["Z.var"]],log=TRUE)) else 0}+
   {if(not.given("Z.mean",theta,given)) sum(dnorm(theta[["Z.mean"]],0,sqrt(prior[["Z.mean.var"]]),log=TRUE)) else 0})
}

ergmm.lpLV.grad<-function(theta,prior,given=list()){
  theta<-.merge.lists(theta,given)
  deriv<-list()
  if(not.given("Z.var",theta,given)) deriv[["Z.var"]]<-prior[["Z.var.df"]]*prior[["Z.var"]]/theta[["Z.var"]]^2/2-(prior[["Z.var.df"]]/2+1)/theta[["Z.var"]]
  if(not.given("Z.mean",theta,given)) deriv[["Z.mean"]]<--theta[["Z.mean"]]/prior[["Z.mean.var"]]
  deriv
}

ergmm.lpdispersion<-function(theta,prior,given=list()){
  theta<-.merge.lists(theta,given)
  (if(not.given("dispersion",theta,given)) dsclinvchisq(theta[["dispersion"]],prior[["dispersion.df"]],prior[["dispersion"]],TRUE) else 0)
}

ergmm.lpdispersion.grad<-function(theta,prior,given=list()){
  theta<-.merge.lists(theta,given)
  deriv<-list()
  if(not.given("dispersion",theta,given)) deriv[["dispersion"]]<-prior[["dispersion.df"]]*prior[["dispersion"]]/theta[["dispersion"]]^2/2-(prior[["dispersion.df"]]/2+1)/theta[["dispersion"]]
  deriv
}


lpsum<-function(theta,which=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV","lpdispersion")){
  sum(sapply(which,function(lp) if(is.null(theta[[lp]])) 0 else theta[[lp]]))
}
