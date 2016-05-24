PVR.confint<-function(obj, variable, skewbeta=FALSE){
   if(missing(obj)) stop("Please provide an object with a list of logistf fits for analysis.\n")
   fit.list<-obj
   
   if(missing(variable)) variable=names(fit.list[[1]]$coefficients)
   nimp<-length(fit.list)
   nvar<-length(variable)
   zalpha<-qnorm(1-fit.list[[1]]$alpha/2)
   pseudo.var<-function(x,tail){
      meanx<-mean(x)
      if(tail==1) signx<--1
      else signx<-1
      diff<-x-meanx
      if(skewbeta) thistail<-(sign(diff)==signx) 
      else thistail<-(diff==diff)
      pseudovar<-sum(diff[thistail]^2)/(sum(thistail)-1)
      return(pseudovar)
    }
   if(nvar>1){
     lower.pse<-(t(sapply(1:nimp,function(X) (fit.list[[X]]$coefficient[variable]-fit.list[[X]]$ci.lower[variable])/zalpha)))
    upper.pse<-(t(sapply(1:nimp,function(X) -(fit.list[[X]]$coefficient[variable]-fit.list[[X]]$ci.upper[variable])/zalpha)))
    lower.var<-lower.pse^2
    upper.var<-upper.pse^2
    coef<-t(sapply(1:nimp,function(X) fit.list[[X]]$coefficient[variable]))
    var.coef.lo<-apply(coef,2,pseudo.var,1)
    var.coef.up<-apply(coef,2,pseudo.var,2)
    mean.coef<-apply(coef,2,mean)
    RR.lower.var<-apply(lower.var,2,mean)+(1+1/nimp)*var.coef.lo
    RR.lower.ci<-mean.coef-zalpha*sqrt(RR.lower.var)
    RR.upper.var<-apply(upper.var,2,mean)+(1+1/nimp)*var.coef.up
    RR.upper.ci<-mean.coef+zalpha*sqrt(RR.upper.var)
    ret<-cbind(RR.lower.ci, RR.upper.ci)
    rownames(ret)<-variable
    colnames(ret)<-c("lower","upper")
   } else {
    lower.pse<-((sapply(1:nimp,function(X) (fit.list[[X]]$coefficient[variable]-fit.list[[X]]$ci.lower[variable])/zalpha)))
    upper.pse<-((sapply(1:nimp,function(X) -(fit.list[[X]]$coefficient[variable]-fit.list[[X]]$ci.upper[variable])/zalpha)))
    lower.var<-lower.pse^2
    upper.var<-upper.pse^2
    coef<-(sapply(1:nimp,function(X) fit.list[[X]]$coefficient[variable]))
    var.coef.lo<-pseudo.var(coef,1)
    var.coef.up<-pseudo.var(coef,2)
    mean.coef<-mean(coef)
    RR.lower.var<-mean(lower.var)+(1+1/nimp)*var.coef.lo
    RR.lower.ci<-mean.coef-zalpha*sqrt(RR.lower.var)
    RR.upper.var<-mean(upper.var)+(1+1/nimp)*var.coef.up
    RR.upper.ci<-mean.coef+zalpha*sqrt(RR.upper.var)
    ret<-cbind(RR.lower.ci, RR.upper.ci)
    colnames(ret)<-c("lower","upper")
    rownames(ret)<-variable
   }
 res<-list(estimate=mean.coef, ci=ret,lower.var=RR.lower.var, upper.var=RR.upper.var,conflev=1-fit.list[[1]]$alpha, call=match.call(), variable=variable)
 attr(res,"class")="PVR.confint"
 return(res)
 }


print.PVR.confint<-function(x,exp=FALSE,...){
  object<-x
  print(object$call)
  cat("Pseudo-variance modification of Rubins Rules\n")
  cat("Confidence level: ", object$conflev*100, "%\n")
  mat<-cbind(object$estimate, object$ci[,1], object$ci[,2], object$lower.var, object$upper.var)
  colnames(mat)<-c("Estimate", "Lower", "Upper", "Lower pseudo variance", "Upper pseudo variance")
  if(exp) {
    mat[,1:3]<-exp(mat[,1:3])
    colnames(mat)[1]<-"Odds ratio"
    }
  rownames(mat)<-object$variable
  print(mat)
 }
 

