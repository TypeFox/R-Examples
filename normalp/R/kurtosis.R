kurtosis<-function(x=NULL,p,value=c("estimate","parameter")){
if (missing(p) && is.null(x)) stop("no arguments inserted")
if (!missing(p) && !is.numeric(p) || !is.null(x) && !is.numeric(x)) stop(" Non-numeric argument to mathematical function")
if (!missing(p) && p<1) stop("p must be at least 1") 
value <- match.arg(value)
 if (is.null(x)){
        vi<-sqrt(gamma(1/p)*gamma(3/p))/gamma(2/p)
        b2<-(gamma(1/p)*gamma(5/p))/(gamma(3/p))^2
        bp<-p+1
 }
 else {
 if (!is.numeric(x)) stop ("x must be a numerical vector")
    n<-length(x)
    if (value=="estimate") {
       cmp<-paramp(x)
       p<-cmp$p
       mp<-cmp$mp
     }
    else  mp<-paramp(x,p=p)$mp
    vi<-(sqrt(n*(sum((x-mp)^2))))/(sum(abs(x-mp)))
    b2<-(n*sum((x-mp)^4))/(sum((x-mp)^2))^2
  Sp<-((sum((abs(x-mp))^p))/n)^(1/p)
  S<-Sp^(2*p)
  bp<-sum((abs(x-mp))^(2*p))/(n*S)
  }
  RIS<-c(VI=vi,B2=b2,Bp=bp)
  RIS
  }

