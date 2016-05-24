hap.freq<-function(geno){

n.poly<-ncol(geno)/2
pM<-allele.freq(geno)

pMd<-array(NA,c(n.poly,n.poly))
pMd.lb<-array(NA,dim=c(n.poly,n.poly))
pMd.ub<-array(NA,dim=c(n.poly,n.poly))

for (i in 1:(n.poly-1)){
 for (j in (i+1):n.poly){
 
  geno.t<-cbind(geno[,(2*i-1):(2*i)],geno[,(2*j-1):(2*j)])
  temp<-haplo.em(geno.t)

  if (nrow(as.matrix(temp$hap.prob))==4){
  pMd[i,j]<-temp$hap.prob[1]
  pMd[j,i]<-temp$hap.prob[1]
  }

  if (nrow(as.matrix(temp$hap.prob))<4){
   a<-sort(unique.default(c(geno[,(2*i-1)],geno[,(2*i)])))
   b<-sort(unique.default(c(geno[,(2*j-1)],geno[,(2*j)])))
   if (temp$haplotype[1,1]==a[1] && temp$haplotype[1,2]==b[1]){
    pMd[i,j]<-temp$hap.prob[1]
    pMd[j,i]<-temp$hap.prob[1]
   }
   else {
    pMd[i,j]<-0
    pMd[j,i]<-0
   }
  }

 }
}

pMd

}
