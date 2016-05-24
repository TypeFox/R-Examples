predict.GORMC <-
function(object,...){
     arg<-list(...)
     M<-length(object$ParEst$Eta)
     P<-length(object$ParEst$Beta)

     if(is.null(arg$len)) arg$len<-100
     if(is.null(arg$new.z)) arg$new.z<-c(1,rep(0,M-1))
     if(is.null(arg$new.x)) arg$new.x<-rep(0,P)
     new.z<-as.vector(arg$new.z)
     new.x<-as.vector(arg$new.x)
     mdata<-object$mdata
     ti<-unique(c(0,na.omit(mdata$Li),na.omit(mdata$Ri)))
     tp<-seq(0,max(ti)+1e-5,length.out=arg$len)
     pzi<-exp(sum(object$ParEst$Eta*new.z))/(1+exp(sum(object$ParEst$Eta*new.z)))
     exb<-exp(sum(object$ParEst$Beta*new.x))
     if(object$ParEst$r>0) surv<-1-pzi+pzi*(1+object$ParEst$r*t(Ispline(tp,order=object$ParEst$order,knots=object$ParEst$knots))%*%object$ParEst$gl*exb)^(-1/object$ParEst$r)
     if(object$ParEst$r==0)  surv<-1-pzi+pzi*exp(-t(Ispline(tp,order=object$ParEst$order,knots=object$ParEst$knots))%*%object$ParEst$gl*exb)
     pred<-list(CureRate=1-pzi,Survival=data.frame(Time=tp,SurvProb=surv))
     class(pred)<-"GORMC"
     class(pred)<-"predict.GORMC"
return(pred)
}
