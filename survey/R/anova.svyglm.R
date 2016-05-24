
oneanova.svyglm<-function(object,test,method){
	tt<-terms(object)
	tlbls<-attr(tt,"term.labels")
	nt<-length(tlbls)
	if (nt<2) return(NULL)
	seqtests<-vector("list",nt)
	if(test=="F") ddf<-NULL else ddf<-Inf
	thismodel<-object
	if (!("formula") %in% names(thismodel$call)) 
            names(thismodel$call)[[2]] <- "formula"
	for(i in nt:1){
		thisterm<-tlbls[i]
		seqtests[[i]]<-regTermTest(thismodel,tlbls[i],method=method,df=ddf)
		thisformula<-make.formula(thisterm)[[2]]
        thismodel <- eval(bquote(update(thismodel, . ~ . - (.(thisformula)))))
		}
	class(seqtests)<-"seqanova.svyglm"
        attr(seqtests,"method")<-method
        attr(seqtests,"test")<-test
	seqtests
	}

print.seqanova.svyglm<-function(x,...){
  isWald<-attr(x,"method")=="Wald"
  isF<-attr(x,"test")=="F"

  cat("Anova table: ")
  if (isWald) cat("(Wald tests)\n") else cat(" (Rao-Scott LRT)\n")
  print(x[[1]]$mcall)

  terms<-sapply(x,"[[","test.terms")
  stats<-if(isF && isWald) sapply(x,"[[","Ftest") else sapply(x,"[[","chisq")
  if(!isWald) stats<-cbind(stats,DEff=sapply(x,function(xi) mean(xi$lambda)))
  df<-sapply(x,"[[","df")
  p<-sapply(x,"[[","p")
  if (!isF){
    rval<-cbind(stats,df,p)
  } else {
    ddf<-sapply(x,"[[","ddf")
    rval<-cbind(stats,df,ddf,p)
  }
  rownames(rval)<-terms
  printCoefmat(rval,tst.ind=1,zap.ind=2:3,has.Pvalue=TRUE)
  invisible(x)
}


SD<-function(x) if (NCOL(x)>1) apply(x,2,sd) else sd(x)
		
anova.svyglm<-function(object, object2=NULL,test=c("F","Chisq"),method=c("LRT","Wald"),tolerance=1e-5,...,force=FALSE){
  test<-match.arg(test)
  method<-match.arg(method)
  if(is.null(object2)) ## sequential tests
    return(oneanova.svyglm(object,test,method))
  
  t1<-attr(terms(object),"term.labels")
  t2<-attr(terms(object2),"term.labels")
  if ((all(t1 %in% t2) || all(t2 %in% t1)) && !force){
    ## symbolically nested, call regTermTest
    biggerobject<-if(all(t1 %in% t2)) object2 else object
    termdiff<-make.formula(if(all(t1 %in% t2)) setdiff(t2,t1) else setdiff(t1,t2))
    if(test=="F") ddf<-NULL else ddf<-Inf
    return(regTermTest(biggerobject,termdiff,df=ddf,method=method))   
  }
  
  ## not symbolically nested, need to project explicitly
  X<-model.matrix(object)
  Z<-model.matrix(object2)
  if (nrow(X)!=nrow(Z)) stop("models have different numbers of observations")
  if (ncol(X)>ncol(Z)) {
    tmp<-X
    X<-Z
    Z<-tmp
    bigger<-1
  } else bigger<-2
  if (any(sapply(suppressWarnings(summary(lm(X~Z))), "[[","sigma")/(tolerance+SD(X))>tolerance)) stop("models not nested")
  
  XX<-matrix(nrow=nrow(Z),ncol=ncol(Z))
  xform<-lm(Z[,1]~X+0)
  XX[,1]<-resid(xform)
  for(i in 2:ncol(Z)){
    XX[,i]<-resid(xform<-lm(Z[,i]~X+Z[,1:(i-1)]+0))
  }
  colkeep<-colMeans(abs(XX))/(tolerance+colMeans(abs(Z))) > tolerance	
  XX<-XX[,colkeep,drop=FALSE]	
  index<-ncol(X)+(1:ncol(XX)) 
  
  ## and now need to refit the model
  ## ugly, but svyglm demands that all variables are in the design argument.
  ## We do know the fitted values at convergence, so one iteration suffices.
  mu<-if(bigger==1) fitted(object) else fitted(object2)
  eta<-if(bigger==1) object$linear.predictors else object2$linear.predictors
  offset<-if(bigger==1) object$offset else object2$offset
  if (is.null(offset)) offset<-0
  pweights<-weights(object$survey.design,"sampling")
  y<-object$y
  if (length(pweights)!=length(y)){
    pweights<-pweights[pweights>0]
    if (length(pweights)!=length(y)) stop("number of observations does not match design")
  }
  pweights<-pweights/mean(pweights)
  ywork<-eta-offset+(y-mu)/object$family$mu.eta(eta)
  wwork<-((pweights * object$family$mu.eta(eta)^2)/object$family$variance(mu))
  wlm<-lm.wfit(cbind(X,XX),ywork,wwork)
    p1<-1:wlm$rank
  Ainv<-chol2inv(wlm$qr$qr[p1,p1,drop=FALSE])
  
  estfun<-cbind(X,XX)*wwork*((y-mu)/object$family$mu.eta(eta))
  design<-object$survey.design
  if (inherits(design, "survey.design2")) 
    V<-svyrecvar(estfun %*% Ainv, design$cluster, design$strata, 
                 design$fpc, postStrata = design$postStrata)
  else if (inherits(design, "twophase")) 
    V<-twophasevar(estfun %*% Ainv, design)
  else if (inherits(design, "twophase2")) 
    V<-twophase2var(estfun %*% Ainv, design)
  else if (inherits(design, "pps")) 
        V<-ppsvar(estfun %*% Ainv, design)
  
  V<-V[index,index]
  df<-min(object$df.residual, object2$df.residual)
  
  if(method=="LRT"){
    V0<-Ainv[index,index]
    chisq <- if(bigger==1) deviance(object2) - deviance(object) else deviance(object)-deviance(object2)
    misspec <- eigen(solve(V0) %*% V, only.values = TRUE)$values
    
    
    if (test=="Chisq") 
      p <- pchisqsum(chisq, rep(1, length(misspec)), misspec, 
                     method = "sad", lower.tail = FALSE)
    else p <- pFsum(chisq, rep(1, length(misspec)), misspec, 
                    ddf = df, method = "sad", lower.tail = FALSE)
    
    rval <- list(call = sys.call(),  chisq = chisq, 
            df = length(index),  p = p, 
                 lambda = misspec, ddf = df, mcall=if(bigger==1) object$call else object2$call,
                 test.terms=if(bigger==1) c(setdiff(t1,t2),"-",setdiff(t2,t1)) else  c(setdiff(t2,t1),"-",setdiff(t1,t2))
                 )
    
    class(rval)<-"regTermTestLRT"
  } else {
    ## method=Wald	
    beta<- wlm$coefficients[index]	
    chisq<-crossprod(beta,solve(V,beta))
    if (test=="Chisq"){
      p<-pchisq(chisq,df=length(index),lower.tail=FALSE)
    } else {
      p<-pf(chisq/length(index),df1=length(index),df2=df,lower.tail=FALSE)	
    }
    rval <- list(call = sys.call(),  Ftest = chisq/length(index), 
                 df = length(index),  p = p, 
                 ddf = df, mcall=if(bigger==1) object$call else object2$call,
                 test.terms=if(bigger==1) c(setdiff(t1,t2),"-",setdiff(t2,t1)) else  c(setdiff(t2,t1),"-",setdiff(t1,t2))
                 )
    class(rval)<-"regTermTest"
    
  }
  
  rval 
}	
