.lrtB<-function(n.fp,n,x,pM,pMd,no.con){

n.poly<-length(pM)

no<-1:n.fp
no.t<-no
info<-array(NA,c(choose(n.poly,n.fp),n.fp))
likelihood<-array(NA,choose(n.poly,n.fp))
mle<-array(NA,choose(n.poly,n.fp))
lr2<-array(NA,choose(n.poly,n.fp))
pMc.v<-array(NA,choose(n.poly,n.fp))


for (i in 1:choose(n.poly,n.fp)){

 check<-0

 if (n.fp>1){

 info[i,]<-no.t
 diff.o<-array()
 diff.ode<-array(NA,c(n.fp,n.fp))

 for (k in 1:n.fp){
  diff.o[k]<-x[no.t[k]]/n[no.t[k]]-pM[no.t[k]]
  for (l in 1:n.fp){
   if (k==l) diff.ode[k,l]<-1
   else diff.ode[k,l]<-(pMd[no.t[k],no.t[l]]-pM[no.t[k]]*pM[no.t[l]])/(pM[no.t[l]]*(1-pM[no.t[l]]))
  } #l
 } #k

if (!any(is.na(diff.ode))){
 if (abs(det(diff.ode))>1e-50){

  diff<- try(solve(diff.ode,tol=1e-100) %*% diff.o)

  temp<-array(NA,c(n.fp,1))
  for (k in 1:n.fp) temp[k]<-pM[no.t[k]]
  if(all(diff<=1-temp & diff>=-temp)){
   check<-1
   pMc<-pM
   for (k in 1:n.poly){
    if (all(k!=no.t)){
     for (l in 1:n.fp) pMc[k]<-pMc[k]+diff[l]*(pMd[k,no.t[l]]-pM[k]*pM[no.t[l]])/(pM[no.t[l]]*(1-pM[no.t[l]]))
    }
    else pMc[k]<-x[k]/n[k]
   }
   if (any(is.na(pMc)) | any(pMc<0 | pMc>1)) check<-2
  }


 } #if det
} #if !any

 z<-sum(no.t==n.poly-n.fp+no)
 if (z!=n.fp) {
  no.t[n.fp-z]<-no.t[n.fp-z]+1
  if (z>0) for (j in (n.fp-z+1):n.fp) no.t[j]<-no.t[j-1]+1
 }


 } # if n.fp>1


 if (n.fp==1){

  check<-1
  pMc<-pM
  pMc[i]<-x[i]/n[i]
  pd<-pM[i]
  pdc<-pMc[i]

  for (j in 1:n.poly){
   if (j!=i){
    for (k in 1:n.fp){
     pMc[j]<-pMc[j]+(pdc[k]-pd[k])*(pMd[j,i]-pM[j]*pd[k])/(pd[k]*(1-pd[k]))
    }
   } #if
  } #j

 info[i,]<-i

 } # if n.fp==1

 if (check==1){
  likelihood[i]<-sum(x*log(pMc)+(n-x)*log(1-pMc))
  mle[i]<-sum(x*log(x/n)+(n-x)*log(1-x/n))

  pMc.v<-.varEst(i,n.fp,n.poly,x/n,pM,pMd,no.con,n[1])  
  lr2[i]<-sum(-2*(x*log(pMc)+(n-x)*log(1-pMc)-(x*log(x/n)+(n-x)*log(1-x/n))) * x*(1-x/n) /(x*(1-x/n)+pMc.v))
 }

} #i

cbind(info,pchisq((lr2),(n.poly-n.fp)),lr2,mle,likelihood)

}
