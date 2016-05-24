ranef.plmm <-
function(object, data, ...){
  
  if(any(class(object)=="wplmm")){# wplmm
    call0<-object$plmm.call            
  }else{ call0<-object$call }# plmm
  
  if(is.null(call0$vc.method)){ vc.method<-"FC"
  }else{ vc.method<-eval(call0$vc.method) }
  # eval() is necessary in a case like VC="SA"; plmm(..., vc.method=VC)
  
  plmm<-new.env()
  Form<-Formula(object$formula)
  match_<-match(c("formula", "data"), names(call0), 0L)
  mframe<-call0[c(1L, match_)]
  mframe[[1]]<-as.name("model.frame")
  mframe[[2]]<-Form
  dat<-eval(mframe, parent.frame())
  
  y<-model.response(dat)                         
  X<-as.matrix(model.matrix(Form, data=dat, rhs=1)[,-1])

  if(missing(data)){
    cls<-as.factor(eval(call0$random, parent.frame()))
  }else{
    cls<-as.factor(eval(call0$random, envir=data, enclos=parent.frame()))
  }  
  plmm$clsLev<-levels(cls)
  plmm$ni<-table(cls)  
  VC<-object$var.comp

  if(!is.null(call0$hetero.prop)){
    if(missing(data)){
      #alpha=eval(parse(text=call0$hetero.prop), parent.frame())   
      alpha<-eval(call0$hetero.prop, parent.frame())
    }else{
      alpha<-eval(call0$hetero.prop, envir=data, enclos=parent.frame())
      #alpha=ifelse(any(deparse(call0$hetero.prop)==names(data)), eval(parse(text=call0$hetero.prop), data), eval(parse(text=call0$hetero.prop), parent.frame()))
    }
    if(length(alpha)==nlevels(cls)){
      alpha<-rep(alpha, times=plmm$ni)
      alpha<-alpha[order(order(plmm$cls))]
    }
    plmm$aInv0<-split(1/alpha, cls)
    plmm$ViInv<-sapply(as.list(plmm$clsLev), inv_Vi_hetero, var_u=VC[1], var_e=VC[2], alpha=alpha, plmm=plmm)
  }
  
  v<-y-X%*%object$coefficients-object$nonpar.values
  plmm$v0<-split(v, cls)

  if(any(class(object)=="wplmm")){# wplmm object
    plmm$nu2Inv0<-split(1/object$var.fun.values, cls)
    u_hat<-sapply(as.list(plmm$clsLev), cal_u_LS, var_u=VC[1], LS=T, plmm=plmm)
  }else{# plmm object
    if(vc.method=="FC" || vc.method=="SA"){
      u_hat<-sapply(as.list(plmm$clsLev), cal_u, var_u=VC[1], var_e=VC[2], plmm=plmm)
    }else if(vc.method=="FChetero" || vc.method=="SAhetero"){
      #plmm$aInv0=split(1/alpha, cls)
      #plmm$ViInv=sapply(plmm$clsLev, inv_Vi_hetero, var_u=VC[1], var_e=VC[2], alpha=alpha, plmm=plmm)
      u_hat<-sapply(as.list(plmm$clsLev), cal_u_LS, var_u=VC[1], LS=F, plmm=plmm)
    }
  }
  
  names(u_hat)<-plmm$clsLev
  return(u_hat)
}
