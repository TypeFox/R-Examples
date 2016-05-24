mat.dis <- function (align1, align2, sub.mat.id = "PAM250",  sqrt = FALSE) {

  if (!exists("sub.mat"))
    data("sub.mat", package = "bios2mds", verbose= FALSE)

  if (!is.element(sub.mat.id, names(sub.mat)))
    stop("sub.mat does not contain sub.mat.id")
       

# prepare substitution matrix
s<-matrix(0,ncol=21,nrow=21)
s[1:20,1:20]<-sub.mat[[sub.mat.id]]
aaName<-colnames(sub.mat[[sub.mat.id]])
colnames(s)<-c(aaName,"-")
rownames(s)<-c(aaName,"-")
  #prepare align
N<-matrix(as.vector(unlist(align1)),ncol=length(align1[[1]]),byrow=T)
rownames(N)<-names(align1)
M<-matrix(as.vector(unlist(align2)),ncol=length(align2[[1]]),byrow=T)
rownames(M)<-names(align2)    


# matrix to gap filter
G<-matrix(1,ncol=21,nrow=21)
colnames(G)<-c(aaName,"-")
rownames(G)<-c(aaName,"-")
G[21,]<-c(rep(0,length(aaName)),1)
for(i in 1:20){
      G[i,][i]<-0
}


D<-matrix(0,nrow=length(N[,1]),ncol=length(M[,1]))
Fn<-lapply(1:length(aaName),function(i){D})
Fm<-Fn
# Length matrix (D) preparation with gap correction
D<-D+ length(M[1,]) +(matrix(G[21,N],nrow=length(N[,1]))%*%t(matrix(G[21,M],nrow=length(M[,1]))))

# S , T and Sr matrix preparation
S<-matrix(0,ncol=length(M[,1]),nrow=length(N[,1]))
T<-matrix(0,ncol=length(M[,1]),nrow=length(N[,1]))
Sr<-matrix(0,ncol=length(M[,1]),nrow=length(N[,1]))
#calcul
for(i in 1:length(N[1,])){
      D[which(N[,i]=="-"),]<-D[which(N[,i]=="-"),]-1
      D[,which(M[,i]=="-")]<-D[,which(M[,i]=="-")]-1
      # Score between N and M for i position
      S<-s[N[,i],M[,i]]+S
# Max Score between N and M for i position with suppression of max if presence of gaps
      lambda<-(matrix(diag(s[N[,i],N[,i]]),ncol=length(M[,i]),nrow=length(N[,i]))+
t(matrix(diag(s[M[,i],M[,i]]),ncol=length(N[,i]),nrow=length(M[,i]))))
      lambda[which(N[,i]=="-"),]<-0
      lambda[,which(M[,i]=="-")]<-0
      T<-T+lambda
# Count of aa use pair to pair between N and M sequences
      for(j in 1:20){
	    Fn[[j]]<-Fn[[j]]+(0^(matrix(G[j,N[,i]],ncol=length(M[,i]),nrow=length(N[,i]))+t(matrix(G[21,M[,i]],ncol=length(N[,i]),nrow=length(M[,i])))))
	    Fm[[j]]<-Fm[[j]]+(0^(t(matrix(G[j,M[,i]],ncol=length(N[,i]),nrow=length(M[,i])))+matrix(G[21,N[,i]],ncol=length(M[,i]),nrow=length(N[,i]))))
      }
}
# Score for random sequences
for(i in 1:20){
for(j in 1:20){
  Sr<-Sr+(Fn[[i]]*Fm[[j]]*s[i,j])
}
}
S<-S
T<-T/2
Sr<-Sr/D
# Correction for sequence distance > 1, V > 1 where Sr is too big
V<-round(1-((S-Sr)/(T-Sr)),3)
V[which(V>=1)]<-0.999
rownames(V)<-rownames(N)
colnames(V)<-rownames(M)
    if(sqrt==FALSE)
        return (V)
    else
        return (sqrt(V))

}
