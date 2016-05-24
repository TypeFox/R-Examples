parents.g <-
function(gene, g, allelefreq){
## input gene as a vector of (1,2,3) : AA==1, Aa==2, aa==3
if(length(allelefreq)==1) allelefreq <- c(allelefreq,0)
qgene <- allelefreq
AAq <- qgene^2
Aaq <- 2*qgene*(1-qgene)
aaq <- (1-qgene)^2

tmp.g <- matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3, 3,1, 3,2, 3,3), ncol=2, byrow=T)
 
tmp.prob<-c(AAq[g]^2, 0.5*AAq[g]*Aaq[g], AAq[g]*aaq[g], 
		      0.5*Aaq[g]*AAq[g], 0.25*Aaq[g]^2, 0.5*Aaq[g]*aaq[g], 
		      aaq[g]*AAq[g], 0.5*aaq[g]*Aaq[g], aaq[g]^2)

tmp <- sample(1:9, 1, prob=Pgene(gene, 1:9)*tmp.prob)

return(tmp.g[tmp,])
}
