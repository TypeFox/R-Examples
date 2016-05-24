`similarity` <-
function(A) {
nc<-ncol(A)
nf <- nrow(A)
matriz<-rep(0,nf*nf)
dim(matriz)<-c(nf,nf)
dimnames(matriz)<-list(row.names(A),row.names(A))
for (k1 in 1:(nf-1)) {
for (k2 in (k1+1):nf) {
mm<- na.omit(c(A[k1,] == A[k2,]))
npar<-length(mm)
sii <- sum(mm)
# npares<- min( length(na.omit(A[k1,])),length(na.omit(A[k2,])) )
matriz[k2,k1]<- sii/npar # chequear el numero de posibilidades
}
}
distance<-as.dist(matriz)
return(distance)
}

