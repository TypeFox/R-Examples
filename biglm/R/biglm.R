"biglm" <-
function(formula, data, weights=NULL, sandwich=FALSE){
   tt<-terms(formula)
   if (!is.null(weights)){
     if (!inherits(weights, "formula"))
       stop("`weights' must be a formula")
     w<-model.frame(weights, data)[[1]]
   } else w<-NULL
   mf<-model.frame(tt,data)
   if (is.null(off<-model.offset(mf))) off<-0
   mm<-model.matrix(tt,mf)
   qr<-bigqr.init(NCOL(mm))
   qr<-update(qr,mm,model.response(mf)-off,w)
   rval<-list(call=sys.call(), qr=qr,assign=attr(mm,"assign"), terms=tt, 
              n=NROW(mm),names=colnames(mm), weights=weights)

   if (sandwich){
     p<-ncol(mm)
     n<-nrow(mm)
     xyqr<-bigqr.init(p*(p+1))
     xx<-matrix(nrow=n, ncol=p*(p+1))
     xx[,1:p]<-mm*(model.response(mf)-off)
     for(i in 1:p)
       xx[,p*i+(1:p)]<-mm*mm[,i]
     xyqr<-update(xyqr,xx,rep(0,n),w*w)
     rval$sandwich<-list(xy=xyqr)
   }
   rval$df.resid<-rval$n-length(qr$D)
   class(rval)<-"biglm"
   rval	
}



print.biglm<-function(x,...){
  cat("Large data regression model: ")
  print(x$call)
  cat("Sample size = ",x$n,"\n")
  invisible(x)
}

summary.biglm<-function(object,...){
   beta<-coef(object)
   se<-sqrt(diag(vcov(object)))
   mat<-cbind(`Coef`=beta, `(95%`=beta-2*se, `CI)`=beta+2*se, `SE`=se, 
               `p`=2*pnorm(abs(beta/se),lower.tail=FALSE))
   rownames(mat)<-object$names
   rval<-list(obj=object, mat=mat)
   if (attr(object$terms,"intercept"))
     rval$nullrss<-object$qr$ss+sum(object$qr$D[-1]*object$qr$thetab[-1]^2)
   else
     rval$nullrss<-object$qr$ss+sum(object$qr$D*object$qr$thetab^2)
   rval$rsq<-1-deviance(object)/rval$nullrss
   class(rval)<-"summary.biglm"
   rval
}

print.summary.biglm<-function(x,digits=getOption("digits")-3,...){
   print(x$obj)
   print(round(x$mat,digits))
   if(!is.null(x$obj$sandwich))
     cat("Sandwich (model-robust) standard errors\n")
   invisible(x)
}


