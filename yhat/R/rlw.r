rlw <- function (dataMatrix, dv, ivlist){

ivlist <- unlist(ivlist)
ilist<-match(ivlist,colnames(dataMatrix))
ilist<-na.omit(ilist)
ilist<-colnames(dataMatrix)[ilist]
dataMatrix<-na.omit(dataMatrix[,c(dv,ilist)])

k<-length(ivlist)
ds<-matrix(ncol=k+1,nrow=nrow(dataMatrix))
ds[,1]<-dataMatrix[,dv]
for (i in 1:k){
  ds[,i+1]<-dataMatrix[,ivlist[i]]
}
colnames(ds)<-c(dv,ivlist)

rxx<-cor(ds[,2:(k+1)])
rxy<-cor(ds[,2:(k+1)],ds[,1])
evm<-eigen(rxx)
ev<-evm$values
evec<-evm$vectors
d<-diag(ev)
delta<-sqrt(d)
lambda<-evec %*% delta %*% t(evec)
lambdasq<-lambda^2
beta<-solve(lambda) %*% rxy
rawrl<-lambdasq %*% beta^2
 
return(rawrl)
}

