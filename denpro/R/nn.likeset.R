nn.likeset<-function(dendat,radmat,k,p=0.1,lambda=NULL)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

volunitball<-volball(1,d)

radit<-radmat[,k]
evat<-k/(n*radit^d*volunitball)
if (is.null(lambda)){
  maksi<-max(evat,na.rm=TRUE)
  lambda<-p*maksi
}
grt<-(evat>=lambda)

#dendatsub<-dendat[grt,]

return(grt)
}

