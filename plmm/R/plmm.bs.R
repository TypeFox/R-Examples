plmm.bs <-
function(object, B, data, h0,...){
  arg<-list(...)
  plmm<-new.env()
  if(any(class(object)=="wplmm")){# wplmm
    call0<-object$plmm.call
    if(missing(data)){
      plmm$hetMat<-eval(object$call$heteroX, parent.frame())
    }else{ 
      plmm$hetMat<-eval(object$call$heteroX, envir=data, enclos=parent.frame())
    }

    if(is.atomic(plmm$hetMat)){
      if(!is.matrix(plmm$hetMat)){# vector or list
        plmm$hetMat<-as.matrix(unlist(plmm$hetMat))
      }
    }else{#2 dim list 
      plmm$hetMat<-matrix(unlist(plmm$hetMat), ncol=length(plmm$hetMat))
    }
    
  }else{ call0<-object$call }# plmm
  
  Form<-Formula(object$formula)
  match_<-match(c("formula", "data"), names(call0), 0L)
  mframe<-call0[c(1L, match_)]
  mframe[[1]]<-as.name("model.frame")
  mframe[[2]]<-Form
  dat<-eval(mframe, parent.frame())
  
  plmm$y<-model.response(dat)
  plmm$X<-as.matrix(model.matrix(Form, data=dat, rhs=1)[,-1])
  plmm$T_mat<-as.matrix(model.matrix(Form, data=dat, rhs=2)[,-1])
  
  if(missing(data)){
    plmm$cls<-as.factor(eval(call0$random, parent.frame()))
  }else{
    plmm$cls<-as.factor(eval(call0$random, envir=data, enclos=parent.frame()))
  }  
  plmm$clsLev<-levels(plmm$cls)
  plmm$ni<-table(plmm$cls)
  dat[[deparse(call0$random)]]<-plmm$cls
  
  N<-length(plmm$cls)
  m<-nlevels(plmm$cls)
  p<-ncol(plmm$X)
  
  if(match("hetero.prop", names(call0), 0L)){
    if(missing(data)){
      plmm$alpha<-eval(call0$hetero.prop, parent.frame())
    }else{
      plmm$alpha<-eval(call$hetero.prop, envir=data, enclos=parent.frame())
    }
    
    if(length(plmm$alpha)==m){
      plmm$alpha<-rep(plmm$alpha, times=plmm$ni)
      plmm$alpha<-plmm$alpha[order(order(plmm$cls))]
    }
    plmm$aInv0<-split(1/plmm$alpha, plmm$cls)
  }
  
  #### Preparation for iteration() ####
  vc.method<-ifelse(is.null(call0$vc.method), "FC", eval(call0$vc.method))
  
  plmm$X_tilde0<-split(plmm$X, plmm$cls)
  PX<-matrix(NA, nrow=N, ncol=p)
  for(i in 1:p){ PX[,i]<-ave(plmm$X[,i], plmm$cls) } 
  
  QX<-plmm$X-PX # plmm$X = X_tilde in the iteration
  plmm$invarIND<-apply(QX, 2, function(x){identical(all.equal(0, sum(x%*%x)), T)})# invarIND can be used for hetero case
  if(any(plmm$invarIND)){ QX<-QX[,!plmm$invarIND] }
  plmm$ncolQX<-ncol(QX)
  
  if(vc.method=="SA" || vc.method=="FC"){
    if(any(plmm$invarIND)){# invariant variables
      plmm$invXQX<-solve(crossprod(plmm$X[,!plmm$invarIND], QX))
    }else{ plmm$invXQX<-solve(crossprod(plmm$X, QX)) }
    
    if(vc.method=="SA"){
      plmm$invXPX<-solve(crossprod(plmm$X, PX))
      
      sum_xixi<-sapply(as.list(plmm$clsLev), SA_xixi, p=p, plmm=plmm)
      sum_xixi<-matrix(unlist(sum_xixi), nrow=p*p)
      plmm$sum_xixi<- matrix(rowSums(sum_xixi), ncol=p)
      
    }else if(vc.method=="FC"){      
      plmm$XX<-crossprod(plmm$X) # (p by p)
      
      sum_nixx<-sapply(as.list(plmm$clsLev), FC_xixi, p=p, plmm=plmm)
      sum_nixx<-matrix(unlist(sum_nixx), nrow=p*p)
      plmm$sum_nixx<-matrix(rowSums(sum_nixx), ncol=p)
    }
  }
  #####################################  
  
  order_cls<-order(plmm$cls)
  u_hat<-ranef(object, data=dat)
  #e_hat=object$residuals[order_cls] - rep(u_hat, plmm$ni)
  e_hat<-object$residuals # original order
  #y_hat=object$fitted.values[order_cls]
  y_hat<-plmm$X%*%object$coefficients+object$nonpar.values# original order

  #h0.call=NULL
  if(missing(h0)){
    h0.call<-object$h0.call
    #h0.call$BSenv=plmm
    h0<-eval(h0.call)
    h0<-h0[[1]]
  #}else if(!is.null(h0$h0.call)){# h0 from select.h0
  }else if(is.atomic(h0)){# user-specified vector or matrix
    h0.call<-object$h0.call
    h0<-matrix(h0, nrow=1)
  }else{# h0 from select.h0
    h0.call<-h0$h0.call    
    h0<-h0[[1]]
  }  
  
  #}else if(!is.null(h0$h0.call)){# h0 from select.h0
    #h0.call=h0$h0.call    
    #h0=h0[[1]]
  #}else{# user-specified vector or matrix
    #h0.call=object$h0.call
    #h0=matrix(h0, nrow=1)
  #}
  
  h0.call$BSenv<-plmm
  h0.call$h_y<-T # to update only h_y in select.h0()
  beta_BS<-matrix(NA, nrow=B, ncol=object$rank)
  VC_BS<-numeric(B)
  
  if(any(class(object)=="wplmm")){
    plmm.call<-object$plmm.call
    pos<-match(c("formula", "h0"), names(plmm.call))
    plmm.call[[pos[1]]]<-Form
    if(is.na(pos[2])){
      plmm.call$h0<-as.name("h0")
    }else{ plmm.call[[pos]]<-as.name("h0") }    
    plmm.call$BSenv<-plmm

    LS.call<-object$call
    #pos=match("plmm.obj", names(LS.call))
    pos<-match("object", names(LS.call))
    LS.call[[pos]]<-as.name("plmm.obj")
    LS.call$BSenv<-plmm
    
    for(b in 1:B){
      #u_BS=rep(rnorm(n=m)*u_hat, ni); e_BS=rnorm(n=N)*e_hat
      #y_BS=object$fitted.values+u_BS+e_BS
      rn<-rnorm(n=N+m)
      e_BS<-rn[1:N]*e_hat
      u_BS<-rep(rn[(N+1):(N+m)]*u_hat, plmm$ni) # sorted order
      #y_BS=y_hat+u_BS+e_BS
      #y_BS=y_hat+u_BS[order(order_cls)]+e_BS
      #y_BS=y_BS[order(order_cls)]
      plmm$y<-y_hat+u_BS[order(order_cls)]+e_BS
   
      ### plmm estimation
      h0[,1]<-eval(h0.call) # update of h_y
      plmm.obj<-eval(plmm.call) # update of plmm.obj

      #dummy$call=object$h0.call# plmm obj h0.call      
      #h_y=update(dummy, BSenv=plmm)
      #dummy$call=object$plmm.call# plmm obj call
      #mod=update(dummy, h_y=h_y, BSenv=plmm)
      
      ### LSplmm estimation
      #LS_mod=update(object, plmm.obj=mod, BSenv=plmm)
      LS_mod<-eval(LS.call)
      beta_BS[b,]<-LS_mod$coefficients          
      VC_BS[b]<-LS_mod$var.comp[1]
    }
  }else{# class(ojbect) = plmm
    plmm.call<-object$call
    pos<-match(c("formula", "h0"), names(plmm.call))
    plmm.call[[pos[1]]]<-Form
    if(is.na(pos[2])){
      plmm.call$h0<-as.name("h0")
    }else{ plmm.call[[pos]]<-as.name("h0") }
    plmm.call$BSenv<-plmm
    
    for(b in 1:B){
      rn<-rnorm(n=N+m)
      e_BS<-rn[1:N]*e_hat
      u_BS<-rep(rn[(N+1):(N+m)]*u_hat, plmm$ni)# sorted order
      plmm$y<-y_hat+u_BS[order(order_cls)]+e_BS

      h0[,1]<-eval(h0.call)
      mod<-eval(plmm.call)
      beta_BS[b,]<-mod$coefficients
      VC_BS[b]<-mod$var.comp[1]
   
      #dummy$call=object$h0.call      
      #h_y=update(dummy, BSenv=plmm)
      #mod=update(object, h_y=h_y, BSenv=plmm)
      #beta_BS[b,]=mod$coefficients
      #VC_BS[b]=mod$var.comp[1]
    }   
  }
  
  #colnames(beta_BS)=attr(terms(Formula(object$formula), rhs=1), "term.labels")
  res<-cbind(beta_BS, VC_BS)
  colnames(res)<-c(attr(terms(Formula(object$formula), rhs=1), "term.labels"), "var.comp")
  #result=list(coef=beta_BS, var.comp=VC_BS)
  res<-list(res)
  class(res)<-"bs.plmm"    
  return(res)
}
