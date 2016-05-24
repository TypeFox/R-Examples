error.rates<-function(H0,Z,pMc,geno,no.ca,no.con=nrow(geno),sim.no=1000){

pM<-allele.freq(geno)
pMd<-hap.freq(geno)
n.poly<-ncol(geno)/2
J<-n.poly

result.t1<-array(NA,sim.no)
result.t2<-array()
resultLL<-array(NA,sim.no)

for (k in 1:sim.no){

n<-array()
x<-array()
for (i in 1:J){
n[i]<-no.ca
x[i]<-rbinom(1,n[i],pMc[i])
}

pMdt<-array(NA,c(n.poly,n.poly))

while(any(is.na(pMdt))){

n.con<-array()
x.con<-array()
for (i in 1:J){
n.con[i]<-no.con
x.con[i]<-rbinom(1,n.con[i],pM[i])
}

pMt<-x.con/n.con
pMdt<-array(NA,c(n.poly,n.poly))
pMd.lbt<-array(NA,dim=c(n.poly,n.poly))
pMd.ubt<-array(NA,dim=c(n.poly,n.poly))
for (i in 1:n.poly) pMdt[i,i]<--1

for (i in 1:(n.poly-1)){
 for (j in (i+1):n.poly){

 if (pMt[j]<=0.5 && pMt[i]<=0.5){
  pMd.lbt[i,j]<-0
  pMd.ubt[i,j]<-min(pMt[j],pMt[i])
 }
 if (pMt[j]<=0.5 && pMt[i]>0.5){
  pMd.lbt[i,j]<-pMt[i]-1+max(pMt[j],1-pMt[i])
  pMd.ubt[i,j]<-pMt[j]
 }
 if (pMt[j]>0.5 && pMt[i]<=0.5){
  pMd.lbt[i,j]<-pMt[j]-1+max(1-pMt[j],pMt[i])
  pMd.ubt[i,j]<-pMt[i]
 }
 if (pMt[j]>0.5 && pMt[i]>0.5){
  pMd.lbt[i,j]<-pMt[i]+pMt[j]-1
  pMd.ubt[i,j]<-1-max(1-pMt[j],1-pMt[i])
 }

 pMdt[i,j]<-rbinom(1,n.con[i],pMd[i,j])/n.con[i]
 if (pMd.ubt[i,j]-pMd.lbt[i,j]>=1/n.con[i]){

  if(pMdt[i,j]>=pMd.lbt[i,j] & pMdt[i,j]<=pMd.ubt[i,j]){ 
   pMdt[j,i]<-pMdt[i,j]
  }
  else{       ## pMdt<lbt | pMdt>ubt
   return.no<-0
   while((pMdt[i,j]<pMd.lbt[i,j] & return.no<20) | (pMdt[i,j]>pMd.ubt[i,j] & return.no<20)){ 
    pMdt[i,j]<-rbinom(1,n.con[i],pMd[i,j])/n.con[i]
    return.no<-return.no+1
   }
   if (pMdt[i,j]>=pMd.lbt[i,j] & pMdt[i,j]<=pMd.ubt[i,j]){
    pMdt[j,i]<-pMdt[i,j]
   }
  }

 }
 else {            #ubt-lbt<1/n.con
  temp.array<-1:n.con[i]/n.con[i]
  temp.value<-temp.array[temp.array<=pMd.ubt[i,j] & temp.array>=pMd.lbt[i,j]]
  if(any(temp.value)){
   pMdt[i,j]<-temp.value
   pMdt[j,i]<-pMdt[i,j]
  }
 }

 } #j
} #i

} # while check

t2<-array()
for (l in 1:Z){

  z<-.lrtB(l,n,x,pMt,pMdt,no.con)

  if (l<Z){
   t2<-c(t2,z[l+1])
  }
  
  if (l==Z){
   result.t1[k]<-z[H0,l+1]
   t2<-c(t2,z[-H0,l+1])  # need to change for each model
  }

} #for

t2<-t2[-1]
result.t2<-c(result.t2,t2)

if(!is.na(result.t1[k])){
 comp.min<-t2[!is.na(t2)]
 if(is.logical(comp.min)) comp.min<-0
 comp.min<-min(comp.min)  
 if(result.t1[k]<=comp.min) resultLL[k]<-1
 else resultLL[k]<-0
}

} #sim

result.t2<-result.t2[-1]

T1.05<-sum(result.t1>0.95,na.rm=TRUE)/sim.no
T1.01<-sum(result.t1>0.99,na.rm=TRUE)/sim.no
T2.05<-sum(result.t2<0.95,na.rm=TRUE)/nrow(as.matrix(result.t2))
T2.01<-sum(result.t2<0.99,na.rm=TRUE)/nrow(as.matrix(result.t2))


result<-c(T1.05,T1.01,T2.05,T2.01,
sum(resultLL,na.rm=TRUE)/sim.no,
sum(is.na(resultLL))/sim.no
)

result

}
