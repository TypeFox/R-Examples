plot.plmm <-
function(x, data, ...){
  arg<-list(...)
  Form<-Formula(x$formula)

  if(any(class(x)=="wplmm")){#wplmm
    call0<-x$plmm.call            
  }else{ call0<-x$call }# ifelse() can't be used for this assignment

  match_<-match(c("formula", "data"), names(call0), 0L)
  mframe<-call0[c(1L, match_)]
  mframe[[1]]<-as.name("model.frame")
  mframe[[2]]<-Form
  dat<-eval(mframe, parent.frame())   
  
  y<-model.response(dat)                         
  X<-as.matrix(model.matrix(Form, data=dat, rhs=1)[,-1])  
  T_mat<-as.matrix(model.matrix(Form, data=dat, rhs=2)[,-1])
  d<-ncol(T_mat)

  yXb<-y-X%*%x$coefficients
  h<-x$h.nonpar
  
  #ifelse(is.null(object$call$poly.index), poly.index<-1, poly.index<-object$call$poly.index)
  if(is.null(x$call$poly.index)){poly.index<-1    
  }else{ poly.index<-x$call$poly.index }

#  grwindow=par()$mfrow
#  oldpar=par(no.readonly=T)
#  oldpar$mfrow=NULL           
#  par(mfrow=c(grwindow))
#  on.exit(par(oldpar))
#browser()
  if(is.null(arg$xlab)){ xlab<-attr(terms(Form, rhs=2), "term.labels") }
  if(is.null(arg$ylab)){ylab<-""} # no label for vertical axis in default

  if(d==1){
    if(!is.null(arg$ann) && arg$ann==F){
      sm.regression(y=yXb, x=T_mat, h=h, ylab="", xlab="", poly.index=poly.index, ...)
    }else if(is.null(arg$ylab)){ # ylab missing
      if(is.null(arg$xlab)){ # xlab missing
        sm.regression(y=yXb, x=T_mat, h=h, poly.index=poly.index, ylab=ylab, xlab=xlab,...)
      }else{ # xlab given
        sm.regression(y=yXb, x=T_mat, h=h, poly.index=poly.index, ylab=ylab, ...)
      }
    }else{ #ylab given
      if(is.null(arg$xlab)){ # xlab missing
        sm.regression(y=yXb, x=T_mat, h=h, poly.index=poly.index, xlab=xlab, ...)
      }else{ # xlab given
        sm.regression(y=yXb, x=T_mat, h=h, poly.index=poly.index, ...)}
    }
  
  }else if(d==2){
    if(!is.null(arg$ann) && arg$ann==F){
      sm.regression(y=yXb, x=T_mat, h=h, poly.index=poly.index, ylab="", xlab="", zlab="", ...)
    }else if(is.null(arg$zlab)){ #zlab missing
      zlab<-""
      if(is.null(arg$xlab)){ # xlab missing
        sm.regression(y=yXb, x=T_mat, h=h, poly.index=poly.index, zlab=zlab, xlab=xlab[1], ylab=xlab[2],...)
      }else{# xlab given
        sm.regression(y=yXb, x=T_mat, h=h, poly.index=poly.index, zlab=zlab, ...)
      }
    }else{ # zlab given
      if(is.null(arg$xlab)){ # xlab missing
        sm.regression(y=yXb, x=T_mat, h=h, poly.index=poly.index, xlab=xlab[1], ylab=xlab[2],...)      
      }else{
        sm.regression(y=yXb, x=T_mat, h=h, poly.index=poly.index, ...)
      } 
    }
  }
}
