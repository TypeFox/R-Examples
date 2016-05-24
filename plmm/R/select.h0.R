select.h0 <-
function(formula, data, nonpar.bws="h.select", poly.index=1,...){
  arg<-list(...)

  if(is.null(arg$BSenv)){
    plmm<-new.env()
    call<-match.call()
    Form<-Formula(formula)
    match_ <- match(c("formula", "data"), names(call), 0L)
    mframe<-call[c(1L, match_)]
    mframe[[1]]<-as.name("model.frame")
    mframe[[2]]<-Form
    dat<-eval(mframe, parent.frame())  
  
    plmm$y<-model.response(dat)
    yName<-row.names(attr(terms(Form), "factors"))[1]
    plmm$X<-as.matrix(model.matrix(Form, data=dat, rhs=1)[,-1])
    xName<-attr(terms(Form, rhs=1), "term.labels")
    plmm$T_mat<-as.matrix(model.matrix(Form, data=dat, rhs=2)[,-1])
    tName<-attr(terms(Form, rhs=2), "term.labels")

  }else{ # from plmm.bs
    plmm<-arg$BSenv
    tName<-NULL
    yName<-NULL
    xName<-NULL
  }   
  
  p<-ncol(plmm$X)
  d<-ncol(plmm$T_mat)
    
### E[y|T] ###
  if(nonpar.bws=="h.select"){
    h<-h.select(y=plmm$y, x=plmm$T_mat, method="cv", poly.index=poly.index, ...)
  }else if(nonpar.bws=="hcv"){
    if(d==1){
      h<-hcv(y=plmm$y, x=as.vector(plmm$T_mat), poly.index=poly.index, ...)
    }else if(d==2){
      h<-hcv(y=plmm$y, x=plmm$T_mat, poly.index=poly.index, ...)
    }
  }
  if(!is.null(arg$h_y)){return(h)}
  
  h0Mat<-matrix(NA, nrow=d, ncol=p+1)
  rownames(h0Mat)<-tName
  colnames(h0Mat)<-c(yName, xName)
  h0Mat[,1]<-h

### E[x|T] ###
  for(i in 1:p){
    if(nonpar.bws=="h.select"){
      h<-h.select(y=plmm$X[,i], x=plmm$T_mat, method="cv", poly.index=poly.index, ...)
    }else if(nonpar.bws=="hcv"){
      if(d==1){
        h<-hcv(y=plmm$X[,i], x=as.vector(plmm$T_mat), poly.index=poly.index, ...)
      }else{
        h<-hcv(y=plmm$X[,i], x=plmm$T_mat, poly.index=poly.index, ...)
      }
    }
    h0Mat[,i+1]<-h
  }
  
  if(!is.null(arg$nbins)){
    nbins<-arg$nbins
  }else{ nbins<-round(8*log(length(plmm$y))/d) }
  
  res<-list(h0=h0Mat, nbins=nbins, h0.call=call)
#  class(result)<-"plmm"
  return(res)
}
