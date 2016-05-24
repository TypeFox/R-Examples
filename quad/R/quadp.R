quadp <-
function(y,A,mycoef){
 n=length(y)
 m=dim(A)[1]

 Dprime=matrix(0,n,n)
 Dprime=matrix(0,n,n)
 for (i in (1:n)){
   for (j in (1:n)){
   Dprime[i,j]=y[i]*y[j]-.5*y[i]^2-.5*y[j]^2
   }
 }
 Sy<-sumfun(Dprime)
 Py<-lincombfun(Sy,mycoef)

  ###############

 #C=A/m
 C=A
 C=C-diag(diag(C))
 Sx<-sumfun(C)
 Px<-lincombfun(Sx,mycoef)
 temp=momentfun(Px,Py,n,mycoef)
 m1=temp$first
 m2=temp$second
 m3=temp$third
 m4=temp$fourth
 mean.q=m1
 variance.q=m2-m1^2
 skewness.q=(m3-3*m1*variance.q-m1^3)/variance.q^(3/2)
 kurtosis.q=(m4-4*m1*m3-4*m1^4+6*m1^2*m2+m1^4)/variance.q^2-3
 moments <- c(mean=mean.q,variance=variance.q,skewness=skewness.q,kurtosis=kurtosis.q+3)
 yvar=var(y)
 # V = colSums((x%*%(as.matrix(y)))^2)/(yvar*(n-1)/n)
 V = sum(C*Dprime)
 p.quad=ppearson(V,moments=moments,lower.tail=F)  # this is the right-tailed approximate p-value
 return(list(stat=V,p=p.quad))
 }
