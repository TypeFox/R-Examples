svyloglin<-function(formula,design,...) UseMethod("svyloglin",design)

withOptions<-function(optlist,expr){
	oldopt<-options(optlist)
	on.exit(options(oldopt))
	expr<-substitute(expr)
	eval.parent(expr)
	}

tr<-function(m)sum(diag(m))
tr2<-function(m) sum(m*m)


svyloglin.survey.design<-function(formula,design,...){
  if (length(formula)!=2) stop("needs a one-sided formula")
  mdata<-model.frame(design)[,all.vars(formula)]
  mf<-model.frame(formula,mdata,na.action=na.pass)
  n<-as.numeric(nrow(mf))
  hatp<-svymean(~I(do.call(interaction,mf)),design,na.rm=TRUE)
  dat<-do.call(expand.grid,lapply(mdata,function(x) sort(unique(x))))
  dat<-as.data.frame(lapply(dat,as.factor))
  dat$y<-coef(hatp)*n
  ff<-update(formula, y~.)
  m1<-withOptions(list(contrasts=c("contr.sum","contr.poly")), 
                  glm(ff, data=dat,family=quasipoisson)
                  )
  P1<-(diag(fitted(m1)/n)-tcrossprod(fitted(m1)/n))/n
  V<-vcov(hatp)
  
  XX<-model.matrix(m1)[,-1,drop=FALSE]
  XX<-sweep(XX,2,colMeans(XX))
  Vtheta<-solve(t(XX)%*%P1%*%XX)%*%(t(XX)%*%V%*%XX)%*%solve(t(XX)%*%P1%*%XX)/(n*n)
  
  rval<-list(model=m1, var=Vtheta, prob.table=hatp,df.null=degf(design),n=n)
  call<-sys.call()
  call[[1]]<-as.name(.Generic)
  rval$call<-call
  class(rval)<-"svyloglin"
  rval
}

print.svyloglin<-function(x,...)  {cat("Loglinear model: ");print(x$call)}

coef.svyloglin<-function(object,...,intercept=FALSE){
  if (intercept)
    coef(object$model)
  else
    coef(object$model)[-1]
}
vcov.svyloglin<-function(object,...) object$var
deviance.svyloglin<-function(object,...) deviance(object$model)
degf.svyloglin<-function(design,...) length(design$prob.table)-length(coef(design))-1
terms.svyloglin<-function(x,...) terms(x$model,...)
model.matrix.svyloglin<-function(object,...) model.matrix(object$model,...)

update.svyloglin<-function(object, formula,...){
  n<-object$n
  model<-withOptions(list(contrasts=c("contr.sum","contr.sum")),
                     update(object$model, formula,data=object$model$model))
  P1<-(diag(fitted(model)/n)-tcrossprod(fitted(model)/n))/n
  V<-vcov(object$prob.table)
  
  XX<-model.matrix(model)[,-1,drop=FALSE]
  XX<-sweep(XX,2,colMeans(XX))
  A<-solve(t(XX)%*%P1%*%XX)
  B<-t(XX)%*%V%*%XX
  Vtheta<-A%*%B%*%A/(n*n)
  
  rval<-list(model=model, var=Vtheta, prob.table=object$prob.table,
             df.null=object$df.null,n=n)
  call<-sys.call()
  call[[1]]<-as.name(.Generic)
  rval$call<-call
  class(rval)<-"svyloglin"
  rval
}

anova.svyloglin<-function(object,object1,...,integrate=FALSE){
  if(length(coef(object1))<length(coef(object))) {
    tmp<-object1
    object1<-object
    object<-tmp
  }
  n<-object$n
  m0<-object$model
  m1<-object1$model
  dfnull<-object$df.null
  pi1<-fitted(m1)/n
  pi0<-fitted(m0)/n
  X1<-model.matrix(m0)[,-1,drop=FALSE]
  X2<-model.matrix(m1)[,-1,drop=FALSE]
  X1<-sweep(X1,2,colMeans(X1))
  X2<-sweep(X2,2,colMeans(X2))
  P1<-(diag(fitted(m1)/n)-tcrossprod(fitted(m1)/n))/n
  P0<-(diag(fitted(m0)/n)-tcrossprod(fitted(m0)/n))/n
  Psat<-(diag(coef(object$prob.table))-tcrossprod(coef(object$prob.table)))/n
  
  if (!all.equal(object$prob.table,object1$prob.table))
    stop("models must be fitted to the same data.")
  V<-vcov(object$prob.table)
  
  wX2sing<-X2-X1%*%solve(t(X1)%*%P1%*%X1,t(X1)%*%P1%*%X2)
  qwX2<-qr(round(wX2sing,11))
  wX2<-wX2sing[,qwX2$pivot[1:qwX2$rank],drop=FALSE]
  Delta<-solve(t(wX2)%*%Psat%*%wX2,t(wX2)%*%V%*%wX2)
  
  an<-anova(m0,m1)
  dev<-an$Deviance[2]
  if (integrate){
    pdev<-pchisqsum(dev,rep(1,ncol(wX2)), a=eigen(Delta,only.values=TRUE,symmetric=TRUE)$values,
                    lower.tail=FALSE,method="integration")
  } else pdev<-NA
  pdev1<-pchisq(dev*ncol(wX2)/tr(Delta),df=ncol(wX2),lower.tail=FALSE)
  pdev2a<-pf(dev/tr(Delta), tr(Delta)^2/tr2(Delta),dfnull*tr(Delta)^2/tr2(Delta),
             lower.tail=FALSE)
  pdevsad<-pchisqsum(dev,rep(1,ncol(wX2)), a=eigen(Delta,only.values=TRUE)$values,
                    lower.tail=FALSE,method="saddlepoint")
  
  pearson<-n*sum( (pi1-pi0)^2/pi0 )
  if (integrate){
    pearsonp<-pchisqsum(pearson, rep(1,ncol(wX2)), a=eigen(Delta,only.values=TRUE,symmetric=TRUE)$values,
                        lower.tail=FALSE,method="integration")
  } else pearsonp<-NA
  prs1<-pchisq(pearson*ncol(wX2)/tr(Delta),df=ncol(wX2),lower.tail=FALSE)
  prs2<-pchisq(pearson*ncol(wX2)/tr(Delta),df=tr(Delta)^2/tr2(Delta),lower.tail=FALSE)
  prs2a<-pf(pearson/tr(Delta), tr(Delta)^2/tr2(Delta),dfnull*tr(Delta)^2/tr2(Delta),
            lower.tail=FALSE)
  pchisqsad<-pchisqsum(pearson, rep(1,ncol(wX2)), a=eigen(Delta,only.values=TRUE,symmetric=TRUE)$values,
                        lower.tail=FALSE,method="saddlepoint")

  rval<-list(an, dev=list(dev=dev, p=c(pdev,pdev1,pdev2a,pdevsad)),
             score=list(chisq=pearson,p=c(pearsonp,prs1,prs2a,pchisqsad)),
             integrate=integrate,a=eigen(Delta,only.values=TRUE,symmetric=TRUE)$values,p=ncol(wX2))
  class(rval)<-"anova.svyloglin"
  rval
}

print.anova.svyloglin<-function(x,pval=c("F","saddlepoint","lincom","chisq"),...){
  cat(attr(x[[1]],"heading"),"\n")
  pval<-match.arg(pval)
  if (pval=="lincom" && !x$integrate){
    x$dev$p[1]<-pchisqsum(x$dev$dev, rep(1,x$p), a=x$a,
                        lower.tail=FALSE,method="integration")
    x$score$p[1]<-pchisqsum(x$score$chisq, rep(1,x$p), a=x$a,
                        lower.tail=FALSE,method="integration")
  }
  cat("Deviance=",x$dev$dev,"p=",
      switch(pval,lincom=x$dev$p[1],
             saddlepoint=x$dev$p[4],
             chisq=x$dev$p[2],
             F=x$dev$p[3]),"\n")
  cat("Score=",x$score$chisq,"p=",
      switch(pval,lincom=x$score$p[1],
             saddlepoint=x$score$p[4],
             chisq=x$score$p[2],
             F=x$score$p[3]),"\n")
  invisible(x)
}

summary.svyloglin<-function(object,...){
	rval<-list(ll=object)
	class(rval)<-"summary.svyloglin"
	rval
	}

print.summary.svyloglin<-function(x,...){
	print(x$ll)
	print(cbind(coef=coef(x$ll),
                    se=SE(x$ll),
                    p=2*pnorm(abs(coef(x$ll)/SE(x$ll)),lower.tail=FALSE)))
	invisible(x)
	}

svyloglin.svyrep.design<-svyloglin.survey.design

svyloglin.DBIsvydesign<-function (formula, design, ...) 
{
    design$variables <- dropFactor(getvars(formula, design$db$connection, 
        design$db$tablename, updates = design$updates), weights(design))
    class(design)<-c("survey.design2","survey.design")    
    rval<-svyloglin(formula,design)
    rval$call<-sys.call()
    rval$call[[1]]<-as.name(.Generic)
    rval
}
svyloglin.ODBCsvydesign<-function (formula, design, ...) 
{
    design$variables <- dropFactor(getvars(formula, design$db$connection, 
        design$db$tablename, updates = design$updates), weights(design))
    class(design)<-c("survey.design2","survey.design")    
    rval<-svyloglin(formula,design)
    rval$call<-sys.call()
    rval$call[[1]]<-as.name(.Generic)
    rval
}
