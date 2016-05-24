`robust.mmm.test` <-
function(y,group,tail=c("right","left","both"))
{

### stop the code if the length of y does not match the length of group ###

 if (length(y)!=length(group))
 {
  stop("the length of the data (y) does not match the length of the group")
 }

### assign tail and name ###

 tail<-match.arg(tail)
 DNAME<-deparse(substitute(y))
 y<-y[!is.na(y)]
 group<-group[!is.na(y)]
 METHOD<-"Mudholkar et al. (1995) test"

### sort the order just in case the input is not sorted by group ###
 
 reorder<-order(group)
 group<-group[reorder]
 y<-y[reorder]

### initialize variables ###

 N<-length(y)
 n<-tapply(y,group,length)
 k<-length(n)
 m<-k-1
 s<-tapply(y,group,var)
 v<-double()
 start.pos<-c(1)

### calculate statistic ###

 for(i in 1:m)
 { 
  j<-sum(n[1:i])+1
  start.pos<-c(start.pos,j)
 }

 for(i in 1:k)
 {
  start<-start.pos[i]
  end<-start+n[i]-1
  suby<-y[start:end]
  suby.length<-length(suby)

  for(j in 1:suby.length)
  {
   subsuby<-suby[-j]
   jackvar<-var(subsuby)
   v<-c(v,jackvar)
  }
 }
 n.group<-n[group]
 s.group<-s[group]

 w<-n.group*log(s.group)-(n.group-1)*log(v)
 w.mean<-tapply(w,group,mean)
 w.mean.group<-w.mean[group]
 q.sq<-sum((w-w.mean.group)^2)/(N-k)
 yi<-double(m)

 for(i in 1:m)
 {
  ni<-n[1:i]
  j<-i+1
  ni2<-n[1:j]
  wi.mean<-w.mean[1:i]
  yi[i]<-(sum(ni)*w.mean[i+1]-sum(ni*wi.mean))/sqrt((n[i+1])^(-1)*sum(ni)*sum(ni2))
 }

 t<-double(m)
 p<-double(m)
 q<-double(m)

 for(i in 1:m)
 {
  subyi<-0
  j<-i-1

  if(i>1) subyi<-yi[1:j]

  subyi.sq<-(subyi)^2
  t[i]<-yi[i]/sqrt(((N-k)*q.sq+sum(subyi.sq))/(N-k+i-1))
  statistic<-t[i]

### calculate log of the component probabilities ###

  df<-N-k+i-1

  if (tail == "left") 
  {
   if(i==1) METHOD<-paste(METHOD,"(left-tailed)")
   p.value = pt(statistic,df,lower.tail=TRUE)
   log.p.value = pt(statistic,df,lower.tail=TRUE,log.p=TRUE)
   log.q.value = pt(statistic,df,lower.tail=FALSE,log.p=TRUE)
  }

  else if (tail == "right") 
  {
   if(i==1) METHOD<-paste(METHOD,"(right-tailed)")
   p.value = pt(statistic,df,lower.tail=FALSE)
   log.p.value = pt(statistic,df,lower.tail=FALSE,log.p=TRUE)
   log.q.value = pt(statistic,df,lower.tail=TRUE,log.p=TRUE)
  }

  else 
  {tail = "both";
   if(i==1) METHOD<-paste(METHOD,"(two-tailed)")
   p.value = pt(statistic,df,lower.tail=TRUE)
   log.p.value = pt(statistic,df,lower.tail=TRUE,log.p=TRUE)
   log.q.value = pt(statistic,df,lower.tail=FALSE,log.p=TRUE)

   if(p.value >= 0.5) 
   {
    p.value = pt(statistic,df,lower.tail=FALSE)
    log.p.value = pt(statistic,df,lower.tail=FALSE,log.p=TRUE)
    log.q.value = pt(statistic,df,lower.tail=TRUE,log.p=TRUE)
   }

   p.value = p.value*2
   log.p.value = log.p.value+log(2)
   log.q.value = log(1-p.value)
  }

  p[i]<-log.p.value
  q[i]<-log.q.value
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

