wplmm <-
function(object, heteroX, data, nonpar.bws="h.select", poly.index=1, var.fun.bws="ROT", var.fun.poly.index=0, scale.h=1, trim=0.01, lim.binning=100, ...){
  arg<-list(...)
  call<-match.call()
  if(is.null(arg$BSenv)){
    plmm<-new.env()
    Form<-Formula(object$formula)
    call0<-object$call
  
    match_<-match(c("formula", "data"), names(call0), 0L)
    mframe<-call0[c(1L, match_)]
    mframe[[1]]<-as.name("model.frame")
    mframe[[2]]<-Form
    dat<-eval(mframe, parent.frame()) 
  
    plmm$y<-model.response(dat)
    plmm$X<-as.matrix(model.matrix(Form, data=dat, rhs=1)[,-1])
    xName<-attr(terms(Form, rhs=1), "term.labels")
    plmm$T_mat<-as.matrix(model.matrix(Form, data=dat, rhs=2)[,-1])

    if(missing(data)){
      plmm$cls<-as.factor(eval(call0$random, parent.frame()))
      hetMat<-eval(call$heteroX, parent.frame())
    }else{
      plmm$cls<-as.factor(eval(call0$random, envir=data, enclos=parent.frame()))
      hetMat<-eval(call$heteroX, envir=data, enclos=parent.frame())
    }
    plmm$clsLev<-levels(plmm$cls)
    plmm$ni<-table(plmm$cls)

    plmm$X_tilde0<-split(plmm$X, plmm$cls)
  
  #### To obtain the colnames of hetMat ####
    if(is.atomic(hetMat)){# var_name in data; vec or mat from parent.frame  
      if(!is.null(colnames(hetMat))){
        name_<-colnames(hetMat)
      }else{
        name_<-deparse(call$heteroX)
        if(is.matrix(hetMat)){name_<-NULL}# mat without colnames
      }
    }else{# named or unnamed list in data or from parent.frame
      name_<-rep(NA, length(hetMat))
      if(!is.null(names(hetMat))){# named
        name_<-names(hetMat)
      }else{# unnamed
        if(length(call$heteroX)==1){
          name_<-deparse(call$heteroX)
# name_ = heteroX name even if it is for a multi-dim list in parent.frame
        }else{
          for(i in 1:length(name_)){
            name_[i]<-deparse(call$heteroX[1+i][[1]])
          }  
        }    
      }
      if(length(name_)<length(hetMat)){name_<-NULL}
# to deal with the above problem
    }
  ##########################################

    if(is.atomic(hetMat)){
      if(!is.matrix(hetMat)){# vector or list
        hetMat<-as.matrix(unlist(hetMat))
      }
    }else{#2 dim list 
      hetMat<-matrix(unlist(hetMat), ncol=length(hetMat))
    }
    plmm$hetMat<-hetMat

  }else{ plmm<-arg$BSenv }
  
  N<-length(plmm$y)
  p<-ncol(plmm$X)
  d<-ncol(plmm$T_mat)
  m<-nlevels(plmm$cls)
  
  if(!is.null(arg$nbins)){nbins<-arg$nbins
  }else{ nbins<-round(8*log(length(plmm$y))/d) }
# nbins for cal_df_gamma1 and cal_df_gamma2; nbins in ... is used for h.select
  
  plmm$y_tilde0<-split(plmm$y-object$nonpar.values, plmm$cls)
  var_u<-object$var.comp[1]
  
  v<-plmm$y-plmm$X%*%object$coefficients-object$nonpar.values
  
  #res=LiStengos(v=v, var_u=var_u, hetMat=hetMat, LS.bws=var.fun.bws, poly.index=poly.index, nonpar.bws=nonpar.bws, N=N, p=p,d=d, plmm=plmm, nbins=nbins, LS.poly.index=var.fun.poly.index, trim=trim, scale.h=scale.h, lim.binning=lim.binning, ...)
  res<-LiStengos(v=v, var_u=var_u, LS.bws=var.fun.bws, poly.index=poly.index, nonpar.bws=nonpar.bws, N=N, p=p,d=d, nbins_=nbins, LS.poly.index=var.fun.poly.index, trim=trim, scale.h=scale.h, lim.binning=lim.binning, plmm=plmm, ...)

  VC<-c(var_u, NA)
  names(VC)[2]<-"var_e"
  
  res<-list(
    coefficients=res$beta,
    residuals=NULL,
    fitted.values=NULL,
    var.comp=VC,
    
    nonpar.values=res$gamma,
    h.nonpar=res$h,
    var.fun.values=res$nu2,
    h.var.fun=res$LS_h,

    rank=p,
    df.residual=N-p-res$df_gamma,
    nbins=nbins,
    formula=object$formula,

    call=call,     
    h0.call=NULL,
    plmm.call=object$call,
    xlevels=NULL,
    heteroX=NULL
  )
  class(res)<-c("wplmm", "plmm")  
    
  if(is.null(arg$BSenv)){
    dat[[deparse(res$plmm.call$random)]]<-plmm$cls
    u_hat<-ranef(res, data=dat)
    order_cls<-order(plmm$cls)
    fitted.values<-(plmm$X%*%res$coefficients+res$nonpar.values)[order_cls]+rep(u_hat, times=plmm$ni)
    res$fitted.values<-c(fitted.values[order(order_cls)])
    res$residuals<-c(plmm$y-fitted.values) 
    
    colnames(res$coefficients)<-"Estimate"
    rownames(res$coefficients)<-xName
   
    res$h0.call<-object$h0.call
    res$xlevels<-object$xlevels
    res$heteroX<-name_
  }
  
  return(res)
}
