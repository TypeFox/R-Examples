##deviance methods not exported, used by method="LRT"
deviance.svycoxph<-function(object,...) 2 * (object$ll[1] - object$ll[2])
deviance.coxph<-function(object,...) 2 * (object$loglik[1] - object$loglik[2])

regTermTest<-function(model, test.terms, null=NULL, df=NULL, method=c("Wald","LRT"), lrt.approximation="saddlepoint"){

  method<-match.arg(method)
  
  canonicalOrder<-function(term){
    tt<-strsplit(term,":")
    tt<-lapply(tt,sort)
    sapply(tt,paste,collapse=":")
  }
    
  
  if(inherits(test.terms,"formula"))
    test.terms<-attr(terms(test.terms),"term.labels")
  
  okbeta<-!is.na(coef(model,na.rm=FALSE)) ## na.rm for svyglm
  tt<-attr(terms(model),"term.labels")
  aa<-attr(model.matrix(model),"assign")[okbeta]
  if((inherits(model,"coxph")|| inherits(model,"svyloglin") || inherits(model,"svyolr"))  && attr(terms(model),"intercept"))
    aa<-aa[-1]
  index<-which(aa %in% match(canonicalOrder(test.terms),canonicalOrder(tt)))
  if (any(is.na(index)))
    stop("Terms didn't match:",canonicalOrder(test.terms),canonicalOrder(tt))
  
  beta<-coef(model)[index]

  if (!is.null(null))
    beta<-beta-null
  V<-vcov(model)[index,index]

  ## this should be rewritten as methods, but that's not happening any time soon.
  if (is.null(df)){
    if (inherits(model,"svyglm"))
      df<-model$df.residual
    else if (inherits(model, "svycoxph"))
      df<-model$degf.resid
    else if (inherits(model,"lm"))
      df<-model$df.residual
    else if (inherits(model,"coxph"))
      df<-model$n-length(coef(model))
    else if (inherits(model, "MIresult"))
      df<-min(model$df[index])
    else if (inherits(model,"svyloglin"))
      df<-model$df+1-length(index)
    else if (inherits(model, "svyolr"))
      df<-model$df.residual
    else
      df<-length(resid(model))-length(coef(model))
  }

  if (method=="LRT"){
    if (inherits(model,"svyglm"))
      V0<-model$naive.cov
    else if (inherits(model, "svycoxph"))
      V0<-model$inv.info
    else if (inherits(model,"lm"))
      V0<-vcov(model)
    else if (inherits(model,"coxph")){
      if (is.null(model$naive.var))
        V0<-model$var
      else
        V0<-model$naive.var
    } else if (inherits(model,"svyolr")) {
      V0<-solve(model$Hess)
    } else stop("method='LRT' not supported for this model")
    V0<-V0[index,index]
    test.formula<-make.formula(test.terms)[[2]]
    if (!("formula") %in% names(model$call))
      names(model$call)[[2]]<-"formula"
    
    model0<-eval(bquote(update(model, .~.-(.(test.formula)))))
    chisq<-deviance(model0)-deviance(model)
    misspec<-eigen(solve(V0)%*%V, only.values=TRUE)$values
    if (df==Inf)
      p<-pchisqsum(chisq,rep(1,length(misspec)),misspec,method=lrt.approximation,lower.tail=FALSE)
    else
      p<-pFsum(chisq,rep(1,length(misspec)),misspec,ddf=df,method=lrt.approximation,lower.tail=FALSE)
      
    rval<-list(call=sys.call(),mcall=model$call,chisq=chisq,
               df=length(index),test.terms=test.terms, 
               p=p,lambda=misspec,ddf=df)
    class(rval)<-"regTermTestLRT"
    return(rval)
  }
  
  
  chisq<-beta%*%solve(V)%*%beta
  if (df<Inf){
    Ftest<-chisq/length(index)
    rval<-list(call=sys.call(),mcall=model$call, Ftest=Ftest,
             df=length(index),ddf=df,test.terms=test.terms,
             p=pf(Ftest,length(index),df,lower.tail=FALSE))
  } else {
    rval<-list(call=sys.call(),mcall=model$call,chisq=chisq,
               df=length(index),test.terms=test.terms,
               p=pchisq(chisq,length(index),lower.tail=FALSE))
  }
  class(rval)<-"regTermTest"
  rval
}


print.regTermTest<-function(x,...){
  cat("Wald test for ")
  cat(x$test.terms)
  cat("\n in ")
  print(x$mcall)
  if(is.null(x$Ftest))
    cat("Chisq = ",x$chisq," on ",x$df," df: p=",format.pval(x$p),"\n")
  else
    cat("F = ",x$Ftest," on ",x$df," and ",x$ddf," df: p=",format.pval(x$p),"\n")
  invisible(x)
}

print.regTermTestLRT<-function(x,...){
  if (is.null(x$ddf) || x$ddf==Inf)
    cat("Working (Rao-Scott) LRT for ")
  else
    cat("Working (Rao-Scott+F) LRT for ")
  cat(x$test.terms)
  cat("\n in ")
  print(x$mcall)
  chisq<-x$chisq/mean(x$lambda)
  cat("Working 2logLR = ",chisq, 'p=',format.pval(x$p),"\n")
  if (length(x$lambda)>1)
    cat("(scale factors: ",signif(x$lambda/mean(x$lambda),2),")")
  else cat("df=1")
  if (!is.null(x$ddf) && is.finite(x$ddf))
    cat(";  denominator df=",x$ddf)
  cat("\n")
  invisible(x)
}

svycontrast<-function(stat, contrasts,...) UseMethod("svycontrast")

match.names <- function(nms,contrasts){
  l<-length(nms)
  ll<-sapply(contrasts,length)
  if (all(ll==l)) return(contrasts)

  if (l==0) stop("No names to match")
  if( !all( unlist(sapply(contrasts,names)) %in% nms))
    stop("names not matched")
  
  lapply(contrasts,
         function(con) {
           r<-numeric(l)
           names(r)<-nms
           r[names(con)]<-con
           r
         })
  
}

contrast<-function(coef,var,contrasts){
  nas<-is.na(var[,1])
  drop<-nas & apply(contrasts,2,function(v) all(v==0))
  if(any(drop)){
    contrasts<-contrasts[,!drop,drop=FALSE]
    coef<-coef[!drop]
    var<-var[!drop,!drop,drop=FALSE]
  }
  if (any(is.na(coef))){
    badin<-is.na(coef)
    bad<-((contrasts!=0)%*%is.na(coef))>0
    rval<-rep(NA,NROW(contrasts))
    rval[!bad]<-contrasts[!bad,!badin,drop=FALSE]%*%coef[!badin]
    v<-matrix(NA,length(rval),length(rval))
    v[!bad,!bad]<-contrasts[!bad,!badin,drop=FALSE]%*%var[!badin,!badin,drop=FALSE]%*%t(contrasts[!bad,!badin,drop=FALSE])
    dimnames(v)<-list(names(rval),names(rval))
    attr(rval, "var")<-v
  } else{
    rval<-contrasts%*%coef
    v<-contrasts%*%var%*%t(contrasts)
    dimnames(v)<-list(names(rval),names(rval))
    attr(rval,"var")<-v
  }
  rval
}

svycontrast.svystat<-function(stat, contrasts,...){
  if (!is.list(contrasts))
    contrasts<-list(contrast=contrasts)
  if (is.call(contrasts[[1]])){
    rval<-nlcon(contrasts,as.list(coef(stat)), vcov(stat))
    class(rval)<-"svrepstat"
    attr(rval,"statistic")<-"nlcon"
    return(rval)
  }
  contrasts<-match.names(names(coef(stat)),contrasts)
  contrasts<-do.call(rbind,contrasts)
  coef<-contrast(coef(stat),vcov(stat),contrasts)
  class(coef)<-"svystat"
  attr(coef,"statistic")<-"contrast"
  coef
}

svycontrast.svystat<-function(stat, contrasts,...){
  if (!is.list(contrasts))
    contrasts<-list(contrast=contrasts)
  if (is.call(contrasts[[1]])){
    rval<-nlcon(contrasts,as.list(coef(stat)), vcov(stat))
    class(rval)<-"svrepstat"
    attr(rval,"statistic")<-"nlcon"
    return(rval)
  }
  contrasts<-match.names(names(coef(stat)),contrasts)
  contrasts<-do.call(rbind,contrasts)
  coef<-contrast(coef(stat),vcov(stat),contrasts)
  class(coef)<-"svystat"
  attr(coef,"statistic")<-"contrast"
  coef
}

svycontrast.svyolr<-function(stat, contrasts,...){
  if (!is.list(contrasts))
    contrasts<-list(contrast=contrasts)
  if (is.call(contrasts[[1]])){
      rval<-nlcon(contrasts,as.list(c(coef(stat),stat$zeta)), vcov(stat))
      class(rval)<-"svystat"
      attr(rval,"statistic")<-"nlcon"
      return(rval)
  }
  contrasts <- match.names(names(coef(stat)), contrasts)
  contrasts<-do.call(rbind,contrasts)
  coef<-contrast(as.vector(as.matrix(coef(stat))),
                 vcov(stat),contrasts)
  class(coef)<-"svystat"
  attr(coef,"statistic")<-"contrast"
  coef
}


svycontrast.svyglm<-svycontrast.svystat
svycontrast.svycoxph<-svycontrast.svystat
svycontrast.svyby<-svycontrast.svystat
svycontrast.default<-svycontrast.svystat

svycontrast.svrepstat<-function(stat, contrasts,...){
  if (!is.list(contrasts))
    contrasts<-list(contrast=contrasts)
  if (is.call(contrasts[[1]])){
    if (is.list(stat)){ ##replicates
      rval<-list(nlcon=nlcon(contrasts,as.list(coef(stat)),vcov(stat)))
      colnames(stat$replicates)<-names(coef(stat))
      rval$replicates<-t(apply(stat$replicates,1,
                             function(repi) nlcon(datalist=as.list(repi),
                                                  exprlist=contrasts, varmat=NULL)))
      attr(rval$nlcon,"statistic")<-"nlcon"
    } else {
      rval<-nlcon(contrasts,as.list(coef(stat)), vcov(stat))
      attr(rval,"statistic")<-"nlcon"
    }
    class(rval)<-"svrepstat"
    return(rval)
  }
  contrasts<-match.names(names(coef(stat)), contrasts)
  contrasts<-do.call(rbind,contrasts)
  
  coef<-contrast(coef(stat), vcov(stat), contrasts)
  if (is.list(stat)){
    coef<-list(contrast=coef,
               replicates=crossprod(stat$replicates, contrasts))
  }
  class(coef)<-"svrepstat"
  attr(coef,"statistic")<-"contrast"
  coef
}



nlcon<-function(exprlist, datalist, varmat){
  if (!is.list(exprlist)) exprlist<-list(contrast=exprlist)
  dexprlist<-lapply(exprlist,
                    function(expr) deriv(expr, names(datalist))[[1]])
  values<-lapply(dexprlist,
                 function(dexpr) eval(do.call(substitute, list(dexpr,datalist))))
  if (is.null(varmat))
    return(do.call(c,values))
  jac<-do.call(rbind,lapply(values,
                            function(value) attr(value,"gradient")))
  var<-jac%*%varmat%*%t(jac)
  values<-do.call(c, values)
  dimnames(var)<-list(names(values),names(values))
  attr(values, "var")<-var
  values
}
