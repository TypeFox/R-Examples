"summary.MCMCglmm"<-function(object, random=FALSE, ...){

  DIC<-object$DIC
  fixed.formula<-object$Fixed$formula
  nF<-object$Fixed$nfl
  nL<-object$Fixed$nll

  if(random){
    nF<-sum(rep(object$Random$nrl, object$Random$nfl))+nF
    if(nF!=dim(object$Sol)[2]){stop("random effects not saved and cannot be summarised")}    
  }

  solutions<-cbind(colMeans(object$Sol[,1:nF,drop=FALSE]), coda::HPDinterval(object$Sol[,1:nF,drop=FALSE]), effectiveSize(object$Sol[,1:nF,drop=FALSE]), 2*pmax(0.5/dim(object$Sol)[1], pmin(colSums(object$Sol[,1:nF,drop=FALSE]>0)/dim(object$Sol)[1], 1-colSums(object$Sol[,1:nF,drop=FALSE]>0)/dim(object$Sol)[1])))
  if(nL>0){
  solutions<-rbind(solutions, cbind(colMeans(object$Lambda), coda::HPDinterval(object$Lambda),effectiveSize(object$Lambda), 2*pmax(0.5/dim(object$Lambda)[1], pmin(colSums(object$Lambda>0)/dim(object$Lambda)[1], 1-colSums(object$Lambda>0)/dim(object$Sol)[1]))))

  }
  colnames(solutions)<-c("post.mean", "l-95% CI", "u-95% CI", "eff.samp", "pMCMC")

  random.formula=object$Random$formula
  residual.formula=object$Residual$formula

  gterms<-sum(object$Random$nfl^2)
  
  rterms<-sum(object$Residual$nfl^2)
  covariances<-cbind(colMeans(object$VCV), coda::HPDinterval(object$VCV), effectiveSize(object$VCV))
  colnames(covariances)<-c("post.mean", "l-95% CI", "u-95% CI","eff.samp")
  if(gterms>0){
   Gcovariances<-covariances[1:gterms,,drop=FALSE]
  }else{
   Gcovariances<-NULL
  }
  Rcovariances<-covariances[gterms+1:rterms,,drop=FALSE]
  cstats<-attr(object$VCV, "mcpar")
  cstats[4]<-dim(object$VCV)[1]
  if(is.null(object$CP)){
     cutpoints<-NULL
  }else{
    cutpoints<-cbind(colMeans(object$CP), coda::HPDinterval(object$CP), effectiveSize(object$CP))
    colnames(cutpoints)<-c("post.mean", "l-95% CI", "u-95% CI", "eff.samp")

  }
  if(is.null(object$Random$nrt)){
    Gterms<-NULL
  }else{
    Gterms<-rep(rep(1:length(object$Random$nrt), object$Random$nrt), object$Random$nfl^2) 
  }
  Rterms<-rep(rep(1:length(object$Residual$nrt), object$Residual$nrt), object$Residual$nfl^2) 

  output<-list(DIC=DIC, fixed.formula=fixed.formula, random.formula=random.formula,residual.formula=residual.formula, solutions=solutions, Gcovariances=Gcovariances, Gterms=Gterms,  Rcovariances=Rcovariances,Rterms=Rterms, cstats=cstats,cutpoints=cutpoints)
  attr(output, "class")<-c("summary.MCMCglmm", "list")
  output
}

"print.summary.MCMCglmm"<-function (x, digits = max(3, getOption("digits") - 3), has.Pvalue=TRUE, eps.Pvalue = 1/(x$cstats[4]-1), cstats=TRUE,  ...) 
{

 if(cstats){
   cat("\n Iterations =", paste(x$cstats[1], ":", x$cstats[2], sep=""))
   cat("\n Thinning interval  =" , x$cstats[3]) 
   cat("\n Sample size  =" , x$cstats[4], "\n") 
 }

 cat("\n DIC:", x$DIC, "\n")
 if(is.null(x$random.formula)==FALSE){
   rcomponents<-split.direct.sum(as.character(x$random.formula)[2])
   for(i in 1:length(rcomponents)){
     if(i==1){
     cat(paste("\n G-structure:  ~", rcomponents[i], "\n\n", sep=""))
     }else{
     cat(paste("\n               ~", rcomponents[i], "\n\n", sep=""))
     }
     if(i%in%x$Gterms){
       print(as.data.frame(x$Gcovariance[x$Gterms==i,,drop=FALSE]), digits=digits, ...)
     }else{
       cat(" G-R structure below\n")
     }
   }
 }
 rcomponents<-split.direct.sum(as.character(x$residual.formula)[2])
 for(i in 1:length(rcomponents)){
   if(i==1){
     cat(paste("\n R-structure:  ~", rcomponents[i], "\n\n", sep=""))
   }else{
     cat(paste("\n               ~", rcomponents[i], "\n\n", sep=""))
   }
   print(as.data.frame(x$Rcovariance[x$Rterms==i,,drop=FALSE]), digits=digits, ...)
 }

 cat("\n Location effects:", paste(as.expression(x$fixed.formula)), "\n\n")
 printCoefmat(as.data.frame(x$solutions), has.Pvalue=has.Pvalue, digits=digits, eps.Pvalue=eps.Pvalue, ...)

 if(is.null(x$cutpoints)==FALSE){
   cat("\n Cutpoints:", "\n")
   print(as.data.frame(x$cutpoints), digits=digits, ...)
 }
}

