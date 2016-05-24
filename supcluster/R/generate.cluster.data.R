generate.cluster.data=function(ratio,npats=80,clusts=c(12,8,12,12,6),
                sig=1,gamma=1,beta=c(-5,-2.5,0,2.5,5)){
  nclusters=length(clusts)
 ngenes=sum(clusts)
#create jmatrix
 jmatrix=matrix(0,ngenes,nclusters)
 start=1
 for (i in 1:nclusters){
   jmatrix[start:(start+clusts[i]-1),i]=1
   start=start+clusts[i]
 }
parms=list(jmatrix=jmatrix,sig=sig,gamma=gamma,beta=beta,tau=ratio*sig)
 sigma<-function(xx){
    J=xx$jmatrix
    mat2=xx$tau*J%*%xx$beta
    Sigma=rbind(cbind(diag(xx$sig,ngenes,ngenes)+
                        tcrossprod(xx$tau*J,J),
                      mat2),cbind(t(mat2),
                                  xx$tau*sum(xx$beta^2)+xx$gamma))
    return(Sigma)
  }

dat=rmvnorm(npats,sigma=sigma(parms))
colnames(dat)<-c(paste('x',1:ngenes,sep=""),"outcome")
dat[1:npats,1:ngenes]=exp(dat[1:npats,1:ngenes])
return(dat)
} 
