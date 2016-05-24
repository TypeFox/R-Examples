allele.freq<-function(geno){

n.poly<-ncol(geno)/2
pM<-array()

for (i in 1:n.poly){
 geno.t<-geno[,(2*i-1):(2*i)]
 if (sum(geno.t[geno.t==4])>0) pM[i]<-sum(geno.t[geno.t==4])/(4*nrow(geno)*2)
 if (sum(geno.t[geno.t==3])>0) pM[i]<-sum(geno.t[geno.t==3])/(3*nrow(geno)*2)
 if (sum(geno.t[geno.t==2])>0) pM[i]<-sum(geno.t[geno.t==2])/(2*nrow(geno)*2)
 if (sum(geno.t[geno.t==1])>0) pM[i]<-sum(geno.t[geno.t==1])/(nrow(geno)*2)
}

pM

}
