gsi.betterPvalue <- function(pval,digits) {
  pval <- round(pval,digits)
  if( pval == 0 )
    pval <- 10^-digits
  structure(pval,digits=digits)
}

Gauss.test <- function(x,y=NULL,mean=0,sd=1,alternative = c("two.sided", "less", "greater")) {
  parameter <-c(mean=mean,sd=sd)
  if(is.null(y)) {
    statistic <- c(T=mean(x))
    isd <- sqrt((sd^2)/length(x))
  } else {
    statistic <- c(T=mean(x)-mean(y))
    sdq <- c(sd^2,sd^2)
    isd <- sqrt(sdq[1]/length(x)+sdq[2]/length(y))
  }
  p1 <- 1-pnorm(statistic,mean=parameter["mean"],sd=isd)
  p2 <- pnorm(statistic,mean=parameter["mean"],sd=isd)
  alternative = match.arg(alternative)
  p.value <- gsi.betterPvalue(switch(alternative,
                    "two.sided"=2*min(p1,p2),
                    "less"=p2,
                    "greater"=p1
                    ),6)
  structure(list(
                 data.name=deparse(substitute(x)),
                 method="Ein Stichproben Gauss-Test",
                 alternative=alternative,
                 parameter=parameter,
                 statistic=statistic,
                 p.value=p.value
                 ),
            class="htest")
}


gsiStepper <- function(x,update) {
  y = x+update
  if( any(y <= 0) ) {
    lf <- x / -update
    lf <- lf[lf>0]
    lf <- min(lf)/2 
    x+lf*update
  } else y
  
}


fitDirichlet <- function(x,elog=mean(ult(x)),alpha0=rep(1,length(elog)),maxIter=20,n=nrow(x)) {
  alpha <- alpha0
  for( i in 1:maxIter ) {
    E <- digamma(alpha)-digamma(sum(alpha))
    V <- gsi.diagGenerate(trigamma(alpha))-trigamma(sum(alpha))
    delta=sqrt(sum((elog-E)^2))
    if( i==1 )
      delta1=delta
    update <- solve(V,(elog-E))
    alpha <- gsiStepper(alpha,update)
#    print(list(alpha=alpha,E=E,V=V,update=update,deltaR=delta/delta1))
  }
  list(alpha=alpha,
       loglikelihood=n*(sum(lgamma(alpha))-lgamma(sum(alpha))+sum(elog*(alpha-1))),
       df=n*(length(elog)-1)-length(elog)
       )
}




acompNormalGOF.test <- function(x,...,method="etest") {
  method <- match.arg(method)
  switch(method,
         "etest"={
           mvnorm.etest(ilr(x),...)
         }
         )
}



gsi.acompUniformityGOF.test<- function(x,
                                samplesize=nrow(x)*20,
                                R=999
                                ) {
  data.name<- paste(deparse(substitute(x)),collapse="",sep="")
  theSample <- runif.acomp(samplesize,ncol(x))
  erg <- acompTwoSampleGOF.test(x,theSample,R=R) 
  erg$data.name<-data.name
  erg$method<-"Compositional test for uniformity"
  erg
}

ccompMultinomialGOF.test<-function(x,
                                   simulate.p.value=TRUE,
                                   R=1999
                                   ){
  chisq.test(unclass(x),simulate.p.value=simulate.p.value,B=R)
}


ccompPoissonGOF.test<-function(x,
                               simulate.p.value=TRUE,
                               R=1999
                               ){
  x <- ccomp(x)
  M <- mean(ccomp(x))
  N <- totals(x)
  erg1 <-chisq.test(unclass(x),simulate.p.value=simulate.p.value,B=R)
  erg2 <-PoissonGOF.test(N,R=R)
  statistic <- min(erg1$p.value , erg2$p.value)
  p.value <- pbeta(statistic,1,2)
  structure(list(
                 data.name=deparse(substitute(x)),
                 method="Count Composition Poission Goodness of fit test",
                 alternative="Count Composition is not a Multi Poission distribution with constant mean",
                 parameter=c(shape1=1,shape2=2),
                 sample=sample,
                 statistic=statistic,
                 p.value=gsi.betterPvalue(p.value,4),
                 erg1=erg1,
                 erg2=erg2
                 ),
            class="htest")
}

gsi.sortedUniforms <- function(n) {
  n = as.integer(c(n,0)[1])
  .C("gsiKSsortedUniforms",
     n    = as.integer(n),
     data = numeric(n),
     getRng= as.integer(1)
     )$data
}
  
PoissonGOF.test <- function(x,lambda=mean(x),R=999,estimated=missing(lambda)) {
    x <- as.integer(x)  
    Max <- as.integer(max(x))
    ps = dpois(0:Max,lambda)
    statistic <- .C("gsiKSPoisson",
                    nd  =as.integer(c(length(x),1)),
                    data=as.integer(x),
                    nps =as.integer(length(ps)),
                    ps  =as.numeric(ps),
                    n   =integer(length(ps)),
                    statistic = numeric(1)
                    )$statistic
    if( estimated ) {
      N <- sum(x)
      xsample <- rmultinom(R,N,rep(1,length(x)))
      ksample <- .C("gsiKSPoisson",
                   nd  =as.integer(dim(xsample)),
                   data=as.integer(xsample),
                   nps =as.integer(length(ps)),
                   ps  =as.numeric(ps),
                   n   =integer(length(ps)),
                   statistic = numeric(R)
                   )$statistic
    } else {
      ksample <- if(R>0)
        .C("gsiKSPoissonSample",
           n=as.integer(length(x)),
           data=numeric(length(x)),
           Nps =as.integer(length(ps)),
           ps  =as.numeric(ps),
           nsamples=as.integer(R),
           statistics=numeric(R)
           )$statistics
    }
    p.value <- gsi.betterPvalue((sum( statistic <= c(ksample) )+1)/(R+1),floor(log(R,10)))
    structure(list(
                   data.name=deparse(substitute(x)),
                   method=if(estimated) "Poisson KS-GOF-Test (with unknown lambda)" else "Poisson KS-GOF-Test (with known lambda)",
                   alternative="Random variable is not Poisson distributed with the given parameter",
                   parameter=c(lambda=lambda),
                   sample=ksample,
                   statistic=c(D=statistic),
                   p.value=p.value
                   ),
              class="htest")
}
 


acompTwoSampleGOF.test <- function(x,y,...,method="etest",data=NULL) {
  acompGOF.test(list(x,y),...,method=method)
}

acompGOF.test <- function(x,...) UseMethod("acompGOF.test")

acompGOF.test.formula <- function(formula, data,...,method="etest") {
  
    if (missing(formula) || (length(formula) != 3) || (length(attr(terms(formula[-2]), "term.labels")) != 1)) 
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    DATA <- split(mf[[response]], g)
    names(DATA) <- levels(g)
    y <- do.call("acompGOF.test", c(DATA, list(...,method=method)))
    y$data.name <- DNAME
    y
  }

gsi.AcompGOFEtest <- function(x,distance=FALSE,R=999,...,dname="data") {
  X <- lapply(x,ilr)
  N <- sapply(x,nrow)
  D <- as.data.frame(do.call("rbind",x))
  erg <- eqdist.etest(D,N,R=R)
  erg$data.name <- dname
  erg
}

acompGOF.test.list <- function(x,...,method="etest") {
  method <- match.arg(method)
  data.name <- paste(deparse(substitute(x)),collapse=" ",sep=" ")
  switch(method,
         etest=gsi.AcompGOFEtest(x,...,dname=data.name)
         )
}




gsi.AddMatrices <- function(M) 
  structure(structure(unlist(M),dim=c(length(M[[1]]),length(M))) %*% rep(1,length(M)),dim=dim(M[[1]]))


fitSameMeanDifferentVarianceModel <- function(x) {
  N <- sapply(x,nrow)
  n <- sum(N)
  G <- length(x)
  m <- ncol(x[[1]])
  a1 <- function(SigmaInv,n) n*SigmaInv
  a2 <- function(SigmaInv,n,mean) n*SigmaInv%*%mean
  a3 <- function(mean,Sigma) {D<-mean-M;Sigma+D%o%D} 
  Sigma0 <- lapply(x,var)
  means  <-  lapply(x,function(x) mean(rmult(x)))
  Sigma  <-  Sigma0
  M <- mean(rmult(do.call(cbind,means)))
  
  for(i in 1:20) {
    Mold <- M
    SigmaInv <- lapply(Sigma,solve)
    M <- rmult(solve(gsi.AddMatrices(mapply(a1,SigmaInv,N,SIMPLIFY=FALSE)),
               gsi.AddMatrices(mapply(a2,SigmaInv,N,means,SIMPLIFY=FALSE))))
    Sigma <- mapply(a3,means,Sigma0,SIMPLIFY=FALSE)
 #   print(norm(M-Mold))
  }
  list(mean=M,vars=Sigma,N=N)
}


acompNormalLocation.test <- function(x,g=NULL,var.equal=FALSE,paired=FALSE,R=ifelse(var.equal,999,0) ) {
  if( paired ) {
    erg <- acompNormalLocation.test(acomp(x)-acomp(g))
    erg$method<-"Compositional paired normal location test"
    return(erg)
  }
  if( inherits(x,"formula") ) {
    # formula interface
    mf <- model.frame(x,g)
    x  <- acomp(mf[[1]])
    g  <- mf[[2]]
    data.name <- mf$names[1]
    Splitted <- split(x,g)
  } else if( !is.null(g)  ) {
    if( is.acomp(g) ) {
          data.name <- paste(deparse(substitute(x)),deparse(substitute(g)),sep=",")
          Splitted <- list(x=x,y=g)
        } else {
          data.name <- deparse(substitute(x))
          Splitted <- split(x,g)
        }
  } else if( !is.list(x) ) {
    # One Sample Test
    data.name <- deparse(substitute(x))
    v <- unclass(ilr(x-mean(x)))
    n <- nrow(x)
    m <- ncol(v)
    w <- unclass(ilr(x))
    tss <- t(w)%*%w
    rss <- t(v)%*%v
    logLik <- n*(log(det(tss/n))-log(det(rss/n)))
    df <- m
    if( R > 0 ) {
      lS1 <- function(x) {
        v <- unclass(rmult(x)-mean(rmult(x)))
        tss <- t(x)%*%x
        rss <- t(v)%*%v
        n*(log(det(tss/n))-log(det(rss/n)))
      }
      sample <- replicate(R,lS1(structure(rnorm(n*m),dim=c(n,m))))
      p.value <- gsi.betterPvalue(mean(logLik<=c(sample,logLik)),floor(log(R,10)))
    } else {
      p.value <- gsi.betterPvalue(pchisq(logLik,df,lower.tail=FALSE),3)
      sample<-NULL
    } 
    return(
    structure(list(
                   data.name=data.name,
                   method="Compositional one sample normal location test",
                   alternative="location is neutral composition",
                   parameter=c(df=df),
                   sample=sample,
                   statistic=c(logLik=logLik),
                   p.value=p.value
                   ),
              class="htest")
    )
  } else {
    data.name <- deparse(substitute(x))
    Splitted <- x
  }
  Splitted <- lapply(Splitted,ilr)
  N <- sapply(Splitted,nrow)
  n <- sum(N)
  G <- length(Splitted)
  m <- ncol(Splitted[[1]])
  css <- function(x) {
    x <- rmult(x)
    x <- unclass(x-mean(x))
    t(x) %*% x
  }
  if( var.equal ) {
    TSS <- css(do.call("rbind",Splitted))
    iRSS <- lapply(Splitted,css)
    tRSS <- gsi.AddMatrices(iRSS)
    logLik <- n*(log(det(TSS/n))-log(det(tRSS/n)))
    df <- (G-1)*m
    if( R>0 ) {
      lS2 <- function(u) {
        Splitted <- lapply(Splitted,function(x) structure(rnorm(length(x)),dim=dim(x)))
        TSS <- css(do.call("rbind",Splitted))
        iRSS <- lapply(Splitted,css)
        tRSS <- gsi.AddMatrices(iRSS)
        n*(log(det(TSS/n))-log(det(tRSS/n)))
      }
      sample <- replicate(R,lS2())
      p.value <- gsi.betterPvalue(mean(logLik<=c(sample,logLik)),floor(log(R,10)))
    } else {
      p.value <- gsi.betterPvalue(pchisq(logLik,df,lower.tail=FALSE),3)
      sample  <- NULL
    }
    structure(list(
                   data.name=data.name,
                   method="Compositional normal location test with equal variances",
                   alternative="locations of groups are not equal",
                   parameter=c(df=df),
                   sample=sample,
                   statistic=c(logLik=logLik),
                   p.value=p.value
                   ),
              class="htest")
  } else {
    smdvm <- fitSameMeanDifferentVarianceModel(Splitted)$vars
    dmdvm <- lapply(Splitted,function(x) css(x)/nrow(x))
    a0 <- function(n,V1,V2) n*(log(det(V1))-log(det(V2)))
    logLik <- sum(mapply(a0,N,smdvm,dmdvm))
    df <- (G-1)*m
    if( R>0 ) {
      lS3 <- function(u) {
        Splitted <- lapply(Splitted,function(x) structure(rnorm(length(x)),dim=dim(x)))
        smdvm <- fitSameMeanDifferentVarianceModel(Splitted)$vars
        dmdvm <- lapply(Splitted,function(x) css(x)/nrow(x))
        a0 <- function(n,V1,V2) n*(log(det(V1))-log(det(V2)))
        sum(mapply(a0,N,smdvm,dmdvm))
      }
      sample <- replicate(R,lS3())
      p.value <- gsi.betterPvalue(mean(logLik<=c(sample,logLik)),floor(log(R,10)))
    } else {
      p.value <- gsi.betterPvalue(pchisq(logLik,df,lower.tail=FALSE),3)
      sample  <- NULL
    }
    structure(list(
                   data.name=data.name,
                   method="Compositional Normal Location test with different variances",
                   alternative="locations of groups are not equal",
                   parameter=c(df=df),
                   sample=sample,
                   statistic=c(logLik=logLik),
                   p.value=p.value
                   ),
              class="htest")
  }
}




  
