predict.plmm <-
function(object, newdata, data, cond=TRUE, ...){
  Form<-Formula(object$formula)
  if(any(class(object)=="wplmm")){#wplmm
    call0<-object$plmm.call            
  }else{call0<-object$call} #plmm
  match_<-match(c("formula", "data"), names(call0), 0L)
  mframe<-call0[c(1L, match_)]
  mframe[[1]]<-as.name("model.frame")
  mframe[[2]]<-Form
  dat<-eval(mframe, parent.frame())

  if(cond){
    if(missing(data)){ cls<-as.factor(eval(call0$random, parent.frame()))
    }else{
      cls<-as.factor(eval(call0$random, envir=data, enclos=parent.frame()))
    }
    dat[[deparse(call0$random)]]<-cls
    u<-ranef(object, data=dat)
  }  

  if(missing(newdata)){# prediction of the response 
    if(cond){
      pred<-object$fitted.values
    }else{# unconditional
      X<-as.matrix(model.matrix(Form, data=dat, rhs=1)[,-1])
      pred<-X%*%object$coefficients+object$nonpar.values     
#      order_cls=order(cls)
#      u_hat=rep(u, times=table(cls))
#      pred=object$fitted.values[order_cls]+u_hat
#      pred=pred[order(order_cls)]
    }
    
  }else{# newdata is given
    y<-model.response(dat)
    T_mat<-as.matrix(model.matrix(Form, dat, rhs=2)[,-1])
    new_t<-as.matrix(model.matrix(Form, data=newdata, rhs=2)[,-1])
    gamma<-sm.regression(y=y, x=T_mat, h=object$h.nonpar, eval.grid=F, eval.points=new_t, display="none")$estimate
    new_X<-as.matrix(model.matrix(Form, data=newdata, rhs=1, xlev=object$xlevels)[,-1])

    Xb<-ifelse(dim(new_X)[2]==1, t(new_X)%*%object$coefficients, new_X%*%object$coefficients)

    if(cond){ 
## check if the condictioning variable is given
      pos<-match(deparse(call0$random), names(newdata))
      #if(is.na(match(deparse(call0$random), names(newdata)))){ stop("random factor is not found in newdata") }
      if(is.na(pos)){ stop("random factor is not found in newdata") }
      new_cls<-as.factor(newdata[[pos]])
      new_cls<-match(new_cls, levels(cls))
#      u=ranef(object, data=dat)[new_cls_ind]
      u_hat<-u[new_cls]
      pred<-Xb+gamma+u_hat
    }else{ pred<-Xb+gamma }
  }
  return(as.numeric(pred))
#  return(prediction=as.numeric(pred))
}
