

svyranktest<-function(formula,design,test=c('wilcoxon','vanderWaerden','median',"KruskalWallis"),...){
        UseMethod("svyranktest", design)
}

svyranktest.survey.design<-svyranktest.svyrep.design<-function(formula, design, test=c('wilcoxon','vanderWaerden','median',"KruskalWallis"),...)
{
  mf<-model.frame(formula,model.frame(design),na.action=na.omit)
  if (!is.null(naa<-attr(mf,"na.action"))){
    design<-design[-naa,]
    mf<-model.frame(formula,model.frame(design),na.action=na.fail)
  }
  y<-mf[,1]
  g<-mf[,2]

  if (length(unique(g))!=2) {
      return(multiranktest(formula,design, test,...))
  }
  
  if (is.character(test)) {
    test<-match.arg(test)
    testf<-switch(test, wilcoxon=,KruskalWallis=function(r,N) r/N,
		  vanderWaerden=function(r,N) qnorm(r/N),
		  median=function(r,N) as.numeric(r>N/2))
  } else{
    testf<-test
  }	
  if (identical(test,"wilcoxon")) test<-"KruskalWallis"
  
  ii<-order(y)
  n<-length(y)
  rankhat<-numeric(n)
  w<-weights(design,"sampling")
  
  N<-sum(w)
  rankhat[ii]<-ave(cumsum(w[ii])-w[ii]/2,factor(y[ii]))
  rankscore<-testf(rankhat,N)
  m <- lm(rankscore~g, weights=w)
  delta<-coef(m)[2]
  xmat<-model.matrix(m)
  infn<- (xmat*(rankscore-fitted(m)))%*%summary(m)$cov.unscaled
  tot.infn<-svytotal(infn,design)
  if (is.character(test))
    method<-paste("Design-based",test,"test")
  else if (!is.null(attr(test,"name")))
    method<-paste("Design-based",attr(test,"name"),"test")
  else method<-"Design-based rank test"
  
  rval <- list(statistic = coef(m)[2]/SE(tot.infn)[2], parameter = degf(design) - 
               1, estimate = coef(m)[2], null.value = 0, alternative = "two.sided", 
               method = method, data.name = deparse(formula))
  rval$p.value <- 2 * pt(-abs(rval$statistic), df = rval$parameter)
  names(rval$statistic) <- "t"
  names(rval$parameter) <- "df"
  names(rval$estimate) <- "difference in mean rank score"
  names(rval$null.value) <- "difference in mean rank score"
  class(rval) <- "htest"
  rval
}

multiranktest<-function(formula,design,test=c('wilcoxon','vanderWaerden','median','KruskalWallis'),...){
    mf<-model.frame(formula,model.frame(design),na.action=na.omit)
    if (!is.null(naa<-attr(mf,"na.action"))){
        design<-design[-naa,]
        mf<-model.frame(formula,model.frame(design),na.action=na.fail)
    }
    y<-mf[,1]
    g<-mf[,2]

    if (is.character(test)) {
        test<-match.arg(test)
        testf<-switch(test, wilcoxon=,KruskalWallis=function(r,N) r/N,
                      vanderWaerden=function(r,N) qnorm(r/N),
                      median=function(r,N) as.numeric(r>N/2))
    } else{
        testf<-test
    }	
    if (identical(test,"wilcoxon")) test<-"KruskalWallis"

    ii<-order(y)
    n<-length(y)
    rankhat<-numeric(n)
    w<-weights(design,"sampling")
    
    N<-sum(w)
    rankhat[ii]<-ave(cumsum(w[ii])-w[ii]/2,factor(y[ii]))
    rankscore<-testf(rankhat,N)
    m <- glm(rankscore~factor(g),weights=w)
    V<-svy.varcoef(m,design)
    ndf<-length(unique(g))-1
    beta<-coef(m)[-1]
    V<-V[-1,-1]
    chisq<-beta%*%solve(V,beta)
    ddf<-degf(design)-ndf
    if (is.character(test))
        method<-paste("Design-based",test,"test")
    else if (!is.null(attr(test,"name")))
        method<-paste("Design-based",attr(test,"name"),"test")
    else method<-"Design-based rank test"
    names(chisq)<-"Chisq"
    names(ndf)<-"df"
    rval<-list(parameter=chisq,statistic=ndf,ddf=ddf,p.value=pf(chisq/ndf,ndf,ddf,lower.tail=FALSE),
               method=method, data.name = deparse(formula))
    class(rval)<-"htest"
    rval
}
