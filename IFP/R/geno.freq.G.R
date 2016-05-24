
geno.freq<-function(genoG){

nAA<-0
nAa<-0
naa<-0
npop<-ncol(genoG)/2
for(i in 1:npop){
 nAA<-nAA+(rowSums(genoG[,c(2*i-1,2*i)])==0)
 nAa<-nAa+(rowSums(genoG[,c(2*i-1,2*i)])==1) 
 naa<-naa+(rowSums(genoG[,c(2*i-1,2*i)])==2) 
}
cbind(nAA/npop,nAa/npop,naa/npop)

}
