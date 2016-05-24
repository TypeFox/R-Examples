svyttest<-function(formula, design,...) UseMethod("svyttest",design)

svyttest.default<-function(formula, design, ...){
  if (formula[[3]]==1 || formula[[3]]==0){
    ## one-sample
    tt <- eval(bquote(svymean(~.(formula[[2]]),design)))
    rval<-list(statistic=coef(tt)[1]/SE(tt)[1],
               parameter=degf(design),
               estimate=coef(tt)[1],
               null.value=0,
               alternative="two.sided",
               method="Design-based one-sample t-test",
               data.name=deparse(formula))
    rval$p.value<-2*pt(-abs(rval$statistic),df=rval$parameter)
    names(rval$statistic)<-"t"
    names(rval$parameter)<-"df"
    names(rval$estimate)<-"mean"
    names(rval$null.value)<-"mean"
    class(rval)<-"htest"
  } else {
    ## two-sample
    m <- eval(bquote(svyglm(formula,design, family=gaussian())))
    rval<-list(statistic=coef(m)[2]/SE(m)[2],
               parameter=degf(design)-1,
               estimate=coef(m)[2],
               null.value=0,
               alternative="two.sided",
               method="Design-based t-test",
               data.name=deparse(formula))
    rval$p.value<-2*pt(-abs(rval$statistic),df=rval$parameter)
    names(rval$statistic)<-"t"
    names(rval$parameter)<-"df"
    names(rval$estimate)<-"difference in mean"
    names(rval$null.value)<-"difference in mean"
    class(rval)<-"htest"
  }
  return(rval)
  
}

expit<-function(eta) exp(eta)/(1+exp(eta))
	
svyciprop<-function(formula, design, method=c("logit","likelihood","asin","beta","mean"),
                    level=0.95,df=degf(design),...)	{
  method<-match.arg(method)
  if (method=="mean"){
    m<-eval(bquote(svymean(~as.numeric(.(formula[[2]])),design,...)))
    ci<-as.vector(confint(m,1,level=level,df=df,...))
    rval<-coef(m)[1]
    attr(rval,"var")<-vcov(m)
  } else if (method=="asin"){
    m<-eval(bquote(svymean(~as.numeric(.(formula[[2]])),design,...)))
    names(m)<-1
    xform<-svycontrast(m,quote(asin(sqrt(`1`))))
    ci<-sin(as.vector(confint(xform,1,level=level,df=df,...)))^2
    rval<-coef(m)[1]
    attr(rval,"var")<-vcov(m)
  } else if (method=="beta"){
    m<-eval(bquote(svymean(~as.numeric(.(formula[[2]])),design,...)))
    n.eff <- coef(m)*(1-coef(m))/vcov(m)
    rval<-coef(m)[1]
    attr(rval,"var")<-vcov(m)
    alpha<-1-level
    n.eff<-n.eff*( qt(alpha/2, nrow(design)-1)/qt(alpha/2, degf(design)) )^2
    ci<-c(qbeta(alpha/2, n.eff*rval,n.eff*(1-rval)+1),
          qbeta(1-alpha/2, n.eff*rval+1, n.eff*(1-rval)))
  } else {
    m<-eval(bquote(svyglm(.(formula[[2]])~1,design, family=quasibinomial)))
    cimethod<-switch(method, logit="Wald",likelihood="likelihood")
    ci<-suppressMessages(as.numeric(expit(confint(m,1,level=level,method=cimethod,ddf=df))))
    rval<-expit(coef(m))[1]
    attr(rval,"var")<-vcov(eval(bquote(svymean(~as.numeric(.(formula[[2]])),design,...))))
  }
  halfalpha<-(1-level)/2
  names(ci)<-paste(round(c(halfalpha,(1-halfalpha))*100,1),"%",sep="")
  names(rval)<-deparse(formula[[2]])
  attr(rval,"ci")<-ci
  class(rval)<-"svyciprop"
  rval
}

confint.svyciprop<-function(object,parm,level=NULL,...){
  if (!is.null(level)) stop("need to re-run svyciprop to specify level")
  rval<-t(as.matrix(attr(object,"ci")))
  rownames(rval)<-names(object)
  rval
}	

coef.svyciprop<-function(object,...) object

vcov.svyciprop<-function(object,...) attr(object,"var")

print.svyciprop<-function(x,digits=max(3,getOption("digits")-4),...){
  m <- cbind(coef(x), confint(x))
  printCoefmat(m,digits=digits)
  invisible(x)
}

