genotype<-function(genoG){

npop<-ncol(genoG)/2
gent<-array(NA,c(nrow(genoG),npop))
for(i in 1:npop){
 gent[,i]<-rowSums(genoG[,c(2*i-1,2*i)])
}
gent

}


