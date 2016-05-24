pFsum<-function(x,df,a,ddf=Inf,lower.tail=TRUE,method=c("saddlepoint","integration","satterthwaite"),...){
  if (ddf==Inf) return(pchisqsum(x,df=df,a=a,lower.tail=lower.tail,...))
  
  method<-match.arg(method)
  if (method=="integration" && !suppressWarnings(require("CompQuadForm"))){
      warning("Package 'CompQuadForm' not found, using saddlepoint approximation")
      method<-"saddlepoint"
    }
  
  
  if (method=="integration"){
    
    int<-davies(0,lambda=c(a,-x/ddf), h=c(df,ddf),acc=1e-7)
    if ( (int$ifault %in% c(0,2))){
      rval<-int$Qq
    } else { 
      rval<-davies(0,lambda=c(a,-x/ddf), h=c(df,ddf),acc=1e-5)$Qq
    }
    if(lower.tail)
      return(1-rval)
    else
      return(rval)
  } else if (method %in% c("satterthwaite","saddlepoint")){
    if(any(df>1)){
      a<-rep(a,df)
    }
    tr<-mean(a)
    tr2<-mean(a^2)/(tr^2)   
    scale=tr*tr2
    ndf=length(a)/tr2
    rval<-pf(x/ndf/scale, ndf,ddf,lower.tail=lower.tail)
    
    if (method=="saddlepoint"){
      a<-c(a,-x/ddf)
      df<-c(df,ddf)
      if(any(df>1))
        a<-rep(a,df)
      s<-saddle(0,a)
      if (!is.na(s)) {
        if (lower.tail)
          rval<-1-s
        else
          rval<-s
      }
    }
    rval
  }
}
