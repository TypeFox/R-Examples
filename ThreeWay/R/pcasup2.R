pcasup2 <-
function(X,n,m,p,model){			

X=as.matrix(X)
out=list()
if (model!=3){
	ea=eigen(X%*%t(X))
	cat("PCASUP: eigenvalues mode A",fill=TRUE)
	la=ea$values
	labComp=paste("Comp.",1:length(la),sep="")
	mat=matrix(,length(la),2)
	mat[,1]=la
	mat[,2]=cumsum(la)/sum(la)*100
	rownames(mat)=labComp
	colnames(mat)=c("Eigenvalue","Fit(%)")
	print(round(mat,digits=2))
	names(la)=labComp
	rownames(ea$vectors)=labComp
	colnames(ea$vectors)=labComp
	out$A=ea$vectors
	out$la=la
}
X=permnew(X,n,m,p)
if (model!=2){
	eb=eigen(X%*%t(X))
	cat("PCASUP: eigenvalues mode B",fill=TRUE)
	lb=eb$values
	labComp=paste("Comp.",1:length(lb),sep="")
	mat=matrix(,length(lb),2)
	mat[,1]=lb
	mat[,2]=cumsum(lb)/sum(lb)*100
	rownames(mat)=labComp
	colnames(mat)=c("Eigenvalue","Fit(%)")
	print(round(mat,digits=2))
	names(lb)=labComp
	rownames(eb$vectors)=labComp
	colnames(eb$vectors)=labComp
	out$B=eb$vectors
	out$lb=lb
}
X=permnew(X,m,p,n)
if (model!=1){
	ec=eigen(X%*%t(X))
	cat("PCASUP: eigenvalues mode C",fill=TRUE)
	lc=ec$values
	labComp=paste("Comp.",1:length(lc),sep="")
	mat=matrix(,length(lc),2)
	mat[,1]=lc
	mat[,2]=cumsum(lc)/sum(lc)*100
	rownames(mat)=labComp
	colnames(mat)=c("Eigenvalue","Fit(%)")
	print(round(mat,digits=2))
	names(lc)=labComp
	rownames(ec$vectors)=labComp
	colnames(ec$vectors)=labComp
	out$C=ec$vectors
	out$lc=lc
}
return(out)
}
