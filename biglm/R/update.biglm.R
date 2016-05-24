
"update.biglm" <-
function(object,moredata,...){
    mf<-model.frame(object$terms, moredata)
    mm<-model.matrix(object$terms, mf)
    if (is.null(object$weights))
      w<-NULL
    else
      w<-model.frame(object$weights, moredata)[[1]]
    if (!identical(object$assign, attr(mm,"assign")))
        stop("model matrices incompatible")
    if(is.null(off<-model.offset(mf))) off<-0
    object$qr<-update.bigqr(object$qr, mm, model.response(mf)-off,w)
    object$n<-object$n+NROW(mm)
    if(!is.null(object$sandwich)){
      p<-ncol(mm)
      n<-nrow(mm)
      xx<-matrix(nrow=n,ncol=p*(p+1))
      xx[,1:p]<-mm*(model.response(mf)-off)
      for(i in 1:p)
        xx[,p*i+(1:p)]<-mm*mm[,i]
      xyqr<-update(object$sandwich$xy,xx,rep(0,n),w*w)
      object$sandwich<-list(xy=xyqr)
    }
    object
  }

