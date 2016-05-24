predict.GORH <-
function(object,...){
     arg<-list(...)
     P<-length(object$ParEst$Beta)

     if(is.null(arg$len)) arg$len<-100
     if(is.null(arg$new.x)) arg$new.x<-rep(0,P)
     new.x<-as.vector(arg$new.x)
     mdata<-object$mdata
     ti<-unique(c(0,na.omit(mdata$Li),na.omit(mdata$Ri)))
     if(is.null(arg$tp)) arg$tp<-seq(0,max(ti)+1e-5,length.out=arg$len)
     exb<-exp(sum(object$ParEst$Beta*new.x))
     surv<-(1+object$ParEst$r*t(Ispline(arg$tp,order=object$ParEst$order,knots=object$ParEst$knots))%*%object$ParEst$gl*exb)^(-1/object$ParEst$r)
     pred<-list(SurvTime=arg$tp,SurvProb=surv)
     class(pred)<-"GORH"
     class(pred)<-"predict.GORH"
return(pred)
}
