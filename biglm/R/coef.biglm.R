"coef.biglm" <-
function(object,...){
    if (!object$qr$checked)
       object$qr<-singcheck.bigqr(object$qr)
    rval<-coef(object$qr)
    rval[object$qr$D==0]<-NA
    names(rval)<-object$names
    rval
}

