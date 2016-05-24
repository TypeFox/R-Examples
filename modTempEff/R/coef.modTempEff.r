`coef.modTempEff`<-function(object, which=c("cold","heat"), L, ...){
      if(length(object$ToTheat)<=0 && length(object$fit.seas)>0) {
        stop("the model does not include csdl coefficients")}
      which<-match.arg(which, several.ok = TRUE)
      cold<-match("cold",which,nomatch = 0)>0
      heat<-match("heat",which,nomatch = 0)>0
      xf<-NULL
      xc<-NULL
      if(cold) xf<- -object$betaCold
      if(heat) xc<- object$betaHeat
      if(length(xf)!=length(xc)) {
        na<-rep(NA,abs(length(xf)-length(xc)))
        if(length(xf)<length(xc)) {
          xf<-c(xf,na)} else {xc<-c(xc,na)}
          }
      m<-cbind(xf,xc)
      rownames(m)<-paste("lag",0:(nrow(m)-1),sep="")
      colnames(m)<-c("cold","heat")
      if(missing(L)) L<-1:nrow(m) else L<-seq_len(L)
      round(m[L,c(cold,heat),drop=FALSE],5)
      }

