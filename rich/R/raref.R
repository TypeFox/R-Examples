raref <-
function(matrix, dens, nrandom=99)
 {
  sSRobs<- function(D,d,ss)
    {E<-D[d,]
    F<-E[1:ss,]
    return(c(specpool(F)$Species,sum(F)))
    }
cat("computing...", "\n")
a<-as.matrix(matrix)
if(is.null(nrandom)==TRUE | nrandom<10) nrandom<-99
if(is.null(dens)==TRUE | dens<=0) stop("invalid dens value")
if (any(is.na(a))) stop("na entries in table")

samplesize<-seq(from=1, to=dim(a)[1],by=1)
sortie<-matrix(NA, ncol=2,nrow=dim(a)[1])

for (i in 1:(length(samplesize))) {
  sb<-boot(a, sSRobs, R=nrandom, ss=samplesize[i])
  sortie[i,1]<-mean(sb$t[,1])
  sortie[i,2]<-mean(sb$t[,2])
}
sortie<-as.data.frame(sortie);sortie$sample<-samplesize
names(sortie)[1:2]<-c("nbsp","ind")

if(dens<=min(sortie$ind) | dens>=max(sortie$ind)) stop("invalid dens value")
indo<-max(which(sortie$ind<=dens))

dmin<-sortie[indo,]
adj<-sortie[indo+1,2]-sortie[indo,2]
opp<-sortie[indo+1,1]-sortie[indo,1]
hyp<-sqrt(adj*adj+opp*opp)
alpha<-asin(opp/hyp)
ad_prime<-dens-sortie[indo,2]
opp_prime<-tan(alpha) * ad_prime
Sinterp<-sortie[indo,1]+opp_prime
cat("done", "\n")
return(list(rar=sortie, Sinterp=c(dens,Sinterp)))}