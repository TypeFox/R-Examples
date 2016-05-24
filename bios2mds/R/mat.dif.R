mat.dif <- function (align1, align2, gap = FALSE, aa.strict = FALSE, sqrt = FALSE) {

  #prepare align
N<-matrix(as.vector(unlist(align1)),ncol=length(align1[[1]]),byrow=T)
rownames(N)<-names(align1)
M<-matrix(as.vector(unlist(align2)),ncol=length(align2[[1]]),byrow=T)
rownames(M)<-names(align2)   
D<-matrix(0,nrow=length(N[,1]),ncol=length(M[,1]))
# prepare substitution matrix
aaName<-union(names(table(N)),names(table(M)))
if(sum(aaName=="-")>0){
aaName<-aaName[-which(aaName=="-")]
}
if(aa.strict){
	aaName<-c("A","R","N","C","D","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
}
aaName<-c(aaName,"-")
s<-matrix(0,ncol=length(aaName),nrow=length(aaName))
diag(s)<-1
G<-c(rep(0,(length(aaName)-1)),1)
names(G)<-aaName
if(!gap){
	s[length(aaName),length(aaName)]<-0
	D<-D+ length(M[1,]) +(matrix(G[N],nrow=length(N[,1]))%*%t(matrix(G[M],nrow=length(M[,1]))))
}else{
D<-D + length(M[1,])
}
colnames(s)<-c(aaName)
rownames(s)<-c(aaName)



# Length matrix (D) preparation with gap correction


# S , T and Sr matrix preparation
S<-matrix(0,ncol=length(M[,1]),nrow=length(N[,1]))

#calcul
for(i in 1:length(N[1,])){
if(!gap){
      D[which(N[,i]=="-"),]<-D[which(N[,i]=="-"),]-1
      D[,which(M[,i]=="-")]<-D[,which(M[,i]=="-")]-1
}
      # Score between N and M for i position
      S<-s[N[,i],M[,i]]+S
}

# Correction for sequence distance > 1, V > 1 where Sr is too big
V<-round(1-(S/D),3)
rownames(V)<-rownames(N)
colnames(V)<-rownames(M)
    if(sqrt==FALSE)
        return (V)
    else
        return (sqrt(V))

}
