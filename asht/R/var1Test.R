var1Test<-function(x,nullVar=1, 
    alternative=c("two.sided","less","greater"),conf.level=0.95){

    alternative<-match.arg(alternative)
    alpha<-1-conf.level
    if (alternative=="two.sided") alpha<-alpha/2
    n<-length(x)
    # test based on normality assumption
    s2<- var(x)
    ci<- ((n-1)*s2)/ qchisq(c(1-alpha,alpha),n-1)
    # one-sided p-value for alternative: sigma2 < nullVar
    pAL<-pchisq( (n-1)*s2/nullVar,n-1) 
    pAG<- 1- pchisq((n-1)*s2/nullVar,n-1)

    if (alternative=="two.sided"){
         p.value<-2*min(pAL,pAG)
    } else if (alternative=="less"){
         p.value<-pAL
         ci<-c(0,ci[2])
    } else if (alternative=="greater"){
         p.value<-pAG
         ci<-c(ci[1],Inf)
    }
    df<-n-1
    names(df)<-"degrees of freedom"
    dname<-deparse(substitute(x))
    names(s2)<-"sample variance"
    ss<-n
    names(nullVar)<-"variance"
    names(ss)<-"sample size"
    attr(ci,"conf.level")<-conf.level
    out<-list(statistic=ss,
          parameter=NULL,
          p.value=p.value,
          conf.int=ci,
          estimate=s2,
          null.value=nullVar,
          alternative=alternative,
          method="One-Sample Normal Test on Variance",
          data.name=dname)
   class(out)<-"htest"
   return(out)
}
  
#var1Test(rnorm(25,sd=.3),alternative="less")


