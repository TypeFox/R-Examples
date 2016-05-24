pcasup3 <-
function(X,n,m,p){			

X=as.matrix(X)
ea=eigen(X%*%t(X))
cat("PCASUP: eigenvalues mode A",fill=TRUE)
la=ea$values

mat=matrix(,length(la),2)
mat[,1]=la
mat[,2]=cumsum(la)/sum(la)*100
labComp=paste("Comp.",1:length(la),sep="")
rownames(mat)=labComp
colnames(mat)=c("Eigenvalue","Fit(%)")
print(round(mat,digits=2))
names(la)=labComp
rownames(ea$vectors)=labComp
colnames(ea$vectors)=labComp


X=permnew(X,n,m,p)
eb=eigen(X%*%t(X))
cat("PCASUP: eigenvalues mode B",fill=TRUE)
lb=eb$values

mat=matrix(,length(lb),2)
mat[,1]=lb
mat[,2]=cumsum(lb)/sum(lb)*100
labComp=paste("Comp.",1:length(lb),sep="")
rownames(mat)=labComp
colnames(mat)=c("Eigenvalue","Fit(%)")
print(round(mat,digits=2))
names(lb)=labComp
rownames(eb$vectors)=labComp
colnames(eb$vectors)=labComp

X=permnew(X,m,p,n)
ec=eigen(X%*%t(X))
cat("PCASUP: eigenvalues mode C",fill=TRUE)
lc=ec$values

mat=matrix(,length(lc),2)
mat[,1]=lc
mat[,2]=cumsum(lc)/sum(lc)*100
labComp=paste("Comp.",1:length(lc),sep="")
rownames(mat)=labComp
colnames(mat)=c("Eigenvalue","Fit(%)")
print(round(mat,digits=2))
names(lc)=labComp
rownames(ec$vectors)=labComp
colnames(ec$vectors)=labComp

out=list()
out$A=ea$vectors
out$B=eb$vectors
out$C=ec$vectors
out$la=la
out$lb=lb
out$lc=lc
return(out)
}
