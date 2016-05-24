rarc <-
function(matrix, samplesize=NULL, nrandom=99)
 {
  sSRobs<- function(D,d,ss)
  {E<-D[d,]
  F<-E[1:ss,]
  return(c(specpool(F)$Species,sum(F)))
  }

if(is.null(nrandom)==TRUE | nrandom<10) nrandom<-99
a<-as.matrix(matrix)
if(is.null(samplesize)==TRUE || is.vector(samplesize)==FALSE) samplesize<-seq(1:dim(a)[1])
sortie<-matrix(NA, ncol=2,nrow=length(samplesize))

for (i in 1:(length(samplesize))) {
  sb<-boot(a, sSRobs, R=nrandom, ss=samplesize[i])
  sortie[i,1]<-mean(sb$t[,1])
  sortie[i,2]<-mean(sb$t[,2])
}
sortie<-as.data.frame(sortie);sortie$sample<-samplesize
names(sortie)[1:3]<-c("richness","individuals","samples")
cat("done", "\n")
return(sortie)
}