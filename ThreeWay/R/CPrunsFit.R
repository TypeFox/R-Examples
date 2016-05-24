CPrunsFit <-
function(X,n,m,p,maxC){

X=as.matrix(X)
if (maxC>max(n,m,p)){
	maxC=max(n,m,p)
}
o=0
for (r in 1:maxC){
	o=o+1
	A=matrix(rnorm(n*r),ncol=r)
	B=matrix(rnorm(m*r),ncol=r)
	C=matrix(rnorm(p*r),ncol=r)
	PAR=CPfuncrep(X,n,m,p,r,1,1,1,0,1e-6,5000,A,B,C)	
	fp=PAR$fp
	if (o==1){
		out=matrix(c(r,r,r,fp),1,4)
	} else{
		out=rbind(out,c(r,r,r,fp))
	}
}
if (nrow(out)>1){
	s=colSums(t(out[,1:3]))
} else{
	s=sum(out[,1:3])
}
out=cbind(out,s)
rownames(out)=rep("CP",length=dim(out)[1])
colnames(out)=c("S","S","S","Fit (%) ","S+S+S")
return(out)
}
