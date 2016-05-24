`mma.test` <-
function(y,group,tail=c("right","left","both"))
{

### assign tail and name###

 tail<-match.arg(tail)
 METHOD<-"Mudholkar et al. (1993) test"
 DNAME<-deparse(substitute(y))
 y<-y[!is.na(y)]
 group<-group[!is.na(y)]

### stop the code if the length of y does not match the length of group ###
 
 if (length(y)!=length(group))
 {
  stop("the length of the data (y) does not match the length of the group")
 }

### sort the order just in case the input is not sorted by group ###

 reorder<-order(group)
 group<-group[reorder]
 y<-y[reorder]

### calculate component statistics f (Mudholkar et al., 1993) ###

 s<-tapply(y,group,var)
 n<-tapply(y,group,length)
 v<-n-1
 k<-length(n)
 m<-k-1
 s2<-double(m)
 v2<-double(m)

 for(i in 1:m)
 {
  s1<-s[1:i]
  v1<-v[1:i]
  s2[i]<-sum(v1*s1)/sum(v1)
  v2[i]<-sum(v1)
 }
 
 s3<-s[2:k]
 v3<-v[2:k]
 f<-s3/s2
 
### create vectors that store log of the component probabilities ###

 p<-double(m)
 q<-double(m)
 
### calculate log of the component probabilities ###

 if(tail=="right")
 {
  METHOD<-paste(METHOD,"(right-tailed)")

  for(j in 1:m)
  {
   p[j]<-pf(f[j],v3[j],v2[j],lower.tail=FALSE, log.p=TRUE)
   q[j]<-pf(f[j],v3[j],v2[j],lower.tail=TRUE, log.p=TRUE)
  }
 }
 
 else if(tail=="left")
 {
  METHOD<-paste(METHOD,"(left-tailed)")

  for(j in 1:m)
  {
   p[j]<-pf(f[j],v3[j],v2[j],lower.tail=TRUE, log.p=TRUE)
   q[j]<-pf(f[j],v3[j],v2[j],lower.tail=FALSE, log.p=TRUE)
  }
 }
 
 else
 {tail ="both";
  METHOD<-paste(METHOD,"(two-tailed)")

  for(j in 1:m)
  {
   r<-pf(f[j],v3[j],v2[j],lower.tail=TRUE, log.p=TRUE)
   s<-pf(f[j],v3[j],v2[j],lower.tail=FALSE, log.p=TRUE)
   
   if(r<log(1/2))
   {
    p[j]<-r+log(2)
    q[j]<-log(1-exp(p[j]))
   }
   
   else
   {
    p[j]<-s+log(2)
    q[j]<-log(1-exp(p[j]))
   }
  }
 }
 
### combine p-values using four p-value combining methods (T,F,N,L) ###

 psiT<-exp(min(p))
 psiF<--2*sum(p)
 psiN<--sum(qnorm(p,lower.tail=TRUE,log.p=TRUE))
 A<-pi^2*m*(5*m+2)/(15*m+12)
 psiL<--A^(-1/2)*sum(p-q)
 
 psiT.pvalue<-1-(1-psiT)^m
 psiF.pvalue<-pchisq(psiF,df=2*m,lower.tail=FALSE)
 psiN.pvalue<-pnorm(psiN,mean=0,sd=sqrt(m),lower.tail=FALSE)
 psiL.pvalue<-pt(psiL,df=5*m+4,lower.tail=FALSE)
 
### display output ###

 PSIT=psiT
 PSIF=psiF
 PSIN=psiN
 PSIL=psiL

 names(PSIT)="Test Statistic (T)"
 names(PSIF)="Test Statistic (F)"
 names(PSIN)="Test Statistic (N)"
 names(PSIL)="Test Statistic (L)"

 list(
 T=structure(list(statistic = PSIT, p.value=psiT.pvalue, data.name = DNAME, method = METHOD), class = "htest"), 
 F=structure(list(statistic = PSIF, p.value=psiF.pvalue, data.name = DNAME, method = METHOD), class = "htest"),
 N=structure(list(statistic = PSIN, p.value=psiN.pvalue, data.name = DNAME, method = METHOD), class = "htest"), 
 L=structure(list(statistic = PSIL, p.value=psiL.pvalue, data.name = DNAME, method = METHOD), class = "htest"))
}

