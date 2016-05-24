var.plot <-
function(object, heteroX, data, var.fun.bws="ROT", var.fun.poly.index=0,...){
  if(class(object)[1]!="plmm"){stop("object must be of plmm class")}
  if(missing(heteroX)){stop("heteroX is missing")}
  
  plmm<-new.env()  
  call<-match.call()
  arg<-list(...)
  Form<-Formula(object$formula)
  call0<-object$call

  match_<-match(c("formula", "data"), names(call0), 0L)
  mframe<-call0[c(1L, match_)]
  mframe[[1]]<-as.name("model.frame")
  mframe[[2]]<-Form
  dat<-eval(mframe, parent.frame())  
  
  y<-model.response(dat)
  X<-as.matrix(model.matrix(Form, data=dat, rhs=1)[,-1])
  T_mat<-as.matrix(model.matrix(Form, data=dat, rhs=2)[,-1])
  N<-length(y)
  d<-ncol(T_mat)
  
  v<-y-X%*%object$coefficients-object$nonpar.values
  var_u <- object$var.comp[1]

  if(missing(data)){
    hetMat<-eval(call$heteroX, parent.frame())
  }else{
    hetMat<-eval(call$heteroX, envir=data, enclos=parent.frame())
  }
  
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
# name_=var_name if var_name is for a multi-dim list in parent.frame
      }else{
        for(i in 1:length(name_)){
          name_[i]<-deparse(call$heteroX[1+i][[1]])
        }  
      }    
    }
    if(length(name_)<length(hetMat)){ name_<-NULL }
# to deal with the above problem
  }
  
  if(is.atomic(hetMat)){
    if(!is.matrix(hetMat)){#var_name in data; vec from parent.frame   
      hetMat<-as.matrix(unlist(hetMat))
    }
  }else{ hetMat<-matrix(unlist(hetMat), ncol=length(hetMat)) }# list
  
  if(!is.null(name_)){colnames(hetMat)<-name_}
  
  plmm$hetMat<-hetMat
  
  #LS_h_var.ij=.var.ij(v=object$residuals, var_u=var_u, hetMat=hetMat, LS.bws=var.fun.bws, LS.poly.index=var.fun.poly.index, d=length(object$nonpar), plmm=plmm,...)
  LS_h_var.ij<-.var.ij(v=v, LS.bws=var.fun.bws, LS.poly.index=var.fun.poly.index, N=N, d=d, plmm=plmm,...)
  
  LS_h<-LS_h_var.ij$LS_h
  nu2<-LS_h_var.ij$var.ij-var_u # nu2 not trimmed
   
  #v2=object$residuals^2   
  v2<-v^2

### plot ###
  if(is.null(arg$xlab)){ xlab<-colnames(hetMat) }

  if(ncol(hetMat)==1){
    if(is.null(arg$ylab)){ ylab<-"squared residuals" }
    if(!is.null(arg$ann) && arg$ann==F ){
      sm.regression(y=v2, x=hetMat, h=LS_h, poly.index=var.fun.poly.index, xlab="", ylab="", ...)
    }else if(is.null(arg$ylab)){ # ylab not given
      if(is.null(arg$xlab)){ # xlab not given
        sm.regression(y=v2, x=hetMat, h=LS_h, poly.index=var.fun.poly.index, ylab=ylab, xlab=xlab, ...)  
      }else{  # xlab given
        sm.regression(y=v2, x=hetMat, h=LS_h, poly.index=var.fun.poly.index, ylab=ylab, ...)        
      }
    }else{# ylab given
      if(is.null(arg$xlab)){ # xlab not given
        sm.regression(y=v2, x=hetMat, h=LS_h, poly.index=var.fun.poly.index, xlab=xlab, ...)
      }else{ # xlab given
        sm.regression(y=v2, x=hetMat, h=LS_h, poly.index=var.fun.poly.index, ...)
      }
    }      
  
  }else if(ncol(hetMat)==2){
    if(is.null(arg$zlab)){ zlab<-"squared residuals" }
    if(!is.null(arg$ann) && arg$ann==F){
      sm.regression(y=v2, x=hetMat, h=LS_h, poly.index=var.fun.poly.index, ylab="", xlab="", zlab="", ...)
    }else if(is.null(arg$zlab)){# zlab missing
      if(is.null(arg$xlab)){ # xlab missing
        sm.regression(y=v2, x=hetMat, h=LS_h, poly.index=var.fun.poly.index, zlab=zlab, xlab=xlab[1], ylab=xlab[2],...)
      }else{ # xlab given
        sm.regression(y=v2, x=hetMat, h=LS_h, poly.index=var.fun.poly.index, zlab=zlab, ...)      
      }
    }else{ # zlab given
      if(is.null(arg$xlab)){ # xlab missing 
        sm.regression(y=v2, x=hetMat, h=LS_h, poly.index=var.fun.poly.index, xlab=xlab[1], ylab=xlab[2], ...)
      }else{
        sm.regression(y=v2, x=hetMat, h=LS_h, poly.index=var.fun.poly.index, ...)      
      }            
    }
  }
  
  return(invisible(list(var.fun.values=nu2, var.comp=as.numeric(var_u), h.var.fun=LS_h)))
}
