`assoc` <-
function(y,x,test="lrt",quantitative)
 {
 lrt<-function(m) 
  {
   if (m$family$family=="gaussian") {
    df1<-m$df.null
    df2<-m$df.residual
    df<-df1-df2 
    ans<-1-pchisq(((m$null.deviance-m$deviance))/(m$deviance/df2),df)
   }
   else {
    ans<-1-pchisq(m$null.deviance-m$deviance,m$df.null-m$df.residual)
   }
   ans 
  }

  G<-function(x,y)
  {
   tt<-table(y,x)
   df<-(dim(tt)[1]-1)*(dim(tt)[2]-1)
   tt.r<-apply(tt,1,sum)
   R<-tt.r[1]
   S<-tt.r[2]
   N<-sum(tt)
   n<-apply(tt,2,sum)
   a1<-(tt[1,]*N)/(R*n)
   a2<-(tt[2,]*N)/(S*n)
   ans <- 2*sum(tt*log(rbind(a1,a2)))
   pval <- 1-pchisq(ans,df=df)
   pval
  }


  if (length(levels(x))==1) {
    pval<-NA
  }
  else {
    if (test=="lrt") {
     if (quantitative)  
      pval<-lrt(glm(y~x,family="gaussian"))
     else  
#      pval<-lrt(glm(y~x,family="binomial"))
       pval<-G(x,y)  
    }
  }
  pval
 }

