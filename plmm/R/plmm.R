plmm <-
function(formula, random, h0, data, vc.method="FC", nonpar.bws="h.select", poly.index=1, iter=20, scale.h=1, epsilon=0.003, lim.binning=100, hetero.prop=NULL, ...){
  arg<-list(...)
  call<-match.call()
  if(is.null(arg$BSenv)){
    BS<-F
    plmm<-new.env()
    Form<-Formula(formula)
    match_<-match(c("formula", "data"), names(call), nomatch=0L)
    mframe<-call[c(1L, match_)]
    mframe[[1]]<-as.name("model.frame")
    mframe[[2]]<-Form
    dat<-eval(mframe, parent.frame())# if data is given, parent.frame falls into eclosure argument
  
    plmm$y<-model.response(dat)
    plmm$X<-as.matrix(model.matrix(Form, data=dat, rhs=1)[,-1])
    xName<-attr(terms(Form, rhs=1), "term.labels")
    plmm$T_mat<-as.matrix(model.matrix(Form, data=dat, rhs=2)[,-1])
    N<-length(plmm$y)
    p<-ncol(plmm$X)
    d<-ncol(plmm$T_mat)
  
    if(missing(data)){
      plmm$cls<-as.factor(eval(call$random, parent.frame()))
    }else{
      plmm$cls<-as.factor(eval(call$random, envir=data, enclos=parent.frame()))
    }  
    plmm$clsLev<-levels(plmm$cls)
    plmm$ni<-table(plmm$cls) 
    m<-nlevels(plmm$cls)
  
    mterms<-attr(dat,"terms")
  #if(any("factor"==attr(terms(dat), "dataClasses"))){
    if(any("factor"==attr(mterms, "dataClasses"))){
      ind<-match("factor", attr(mterms, "dataClasses")) 
      xlevels<-list(NULL)
      for(i in 1:length(ind)){
        xlevels[[i]]<-levels(dat[[ind[i]]])
        names(xlevels)[[i]]<-names(attr(terms(dat), "dataClasses"))[ind]
      }
    }else{xlevels<-NULL}
# xlevels is used in prediction. c.f. lm()
  
    if(!missing(hetero.prop)){
      if(missing(data)){
        plmm$alpha<-eval(call$hetero.prop, parent.frame())
      }else{
      plmm$alpha<-eval(call$hetero.prop, envir=data, enclos=parent.frame())
      }
    
      if(length(plmm$alpha)==m){
        plmm$alpha<-rep(plmm$alpha, times=plmm$ni)
        plmm$alpha<-plmm$alpha[order(order(plmm$cls))]
      }
      plmm$aInv0<-split(1/plmm$alpha, plmm$cls)
    }

    if(missing(h0)){#without h0 input
      h0<-select.h0(formula=formula, data=dat)
      h0$h0.call[[2]]<-Form
      h0.call<-h0$h0.call
      h0<-h0[[1]]
    }else if(!is.null(h0$h0.call)){# from select.h0
      h0.call<-h0$h0.call    
      h0<-h0[[1]]
    }else{# h0 is user-specified vector or matrix
      if(is.vector(h0)){ h0<-matrix(h0, nrow=1) }
      h0.call<-NULL
    }
  
  }else{# plmm.bs
    BS<-T
    plmm<-arg$BSenv
    N<-length(plmm$y)
    p<-ncol(plmm$X)
    d<-ncol(plmm$T_mat)
    m<-nlevels(plmm$cls)# clustering variable
  }  

  if(!is.null(arg$nbins)){
    nbins<-arg$nbins
  }else{ nbins<-round(8*log(length(plmm$y))/d) }
# nbins for cal_df_gamma1 and cal_df_gamma2; nbins in ... is used for h.select

  res<-iteration0(N=N, p=p, m=m, d=d, h0=h0, vc.method=vc.method, poly.index=poly.index, scale.h=scale.h, lim.binning=lim.binning, nbins_=nbins, plmm=plmm, ...)

### Iteratve Process
  res<-iteration(
      beta=res$beta, 
      VC=res$VC, 
      gamma=res$gamma, 
      h=res$h,
      df_gamma=res$df_gamma,
      vc.method=vc.method, 
      iter=iter, 
      nonpar.bws=nonpar.bws, 
      poly.index=poly.index, 
      N=N, m=m, p=p, d=d, 
      scale.h=scale.h, epsilon=epsilon, lim.binning=lim.binning, 
      nbins_=nbins, plmm=plmm, BS=BS, ...)
  
  if(res$iter>0){
      var.comp<-res$VC_iter[res$iter+1,]
  }else{ var.comp<-res$VC_iter }
  
  res<-list(
          coefficients=res$beta, 
		  fitted.values=NULL,
          residuals=NULL,
          var.comp=var.comp,
          nonpar.values=res$gamma,
          h.nonpar=res$h,
          rank=p,
          df.residual=N-p-res$df_gamma,
          nbins=nbins,
          iter=res$iter,
          coef.iter=res$beta_iter,
          vc.iter=res$VC_iter,
          formula=formula,
          call=call, 
          h0.call=NULL,
          xlevels=NULL
          )
  class(res)<-"plmm"

  if(!BS){
    dat[[deparse(call$random)]]<-plmm$cls
    #u_hat<-ranef(res, data=dat, vc.method=vc.method)
# in case vc.method is not specified in call
    u_hat<-ranef(res, data=dat)
    order_cls<-order(plmm$cls)
    fitted.values<-(plmm$X%*%res$coefficients+res$nonpar.values)[order_cls]+rep(u_hat, times=plmm$ni)
    res$fitted.values<-as.vector(fitted.values[order(order_cls)])
    res$residuals<-as.vector(plmm$y-res$fitted.values) 
    
    colnames(res$coefficients)<-"Estimate"
    rownames(res$coefficients)<-xName
    if(iter>0){
      colnames(res$coef.iter)<-xName
    }else{attr(res$coef.iter,"names")<-xName} 
    
    res$h0.call<-h0.call
    res$xlevels<-xlevels
  }
  
  return(res)
}
