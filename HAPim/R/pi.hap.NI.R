`pi.hap.NI` <-
function(res.structure,cor.hap){
Bin.type=res.structure$Bin.type
Index=res.structure$Index
Hap=res.structure$Hap
nbre.type=length(Bin.type[,1])
nbre.marq=length(Bin.type[1,])
nbre.hap=length(Hap[,1])
corresp=cor.hap$corresp
haplo.obs=cor.hap$assoc
pi.hap=list()

pi.hap[[1]]=rep(0,nbre.hap)
nbe.obs.sansNA=sum(corresp[,2]==1)
for(ih in haplo.obs){
pi.hap[[1]][ih]=sum(corresp[,2]==1 & corresp[,3]==ih)/nbe.obs.sansNA
}




pi.hap[[nbre.type]]=rep(1,nbre.hap)
for(it in 2:(nbre.type-1)){
#Bin=Bin.type[it,]
#Zero=(1:nbre.marq)[Bin==0]
nc=ncol(Index[[it]])
nr=nrow(Index[[it]])
for(j in 1:nc){
pi.hap[[it]][Index[[it]][,j]]=rep(sum(pi.hap[[1]][Index[[it]][,j]]),nr)
}
}

pi.hap
}

