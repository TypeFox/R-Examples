jointplotgen <-
function(K,A,B,C,fixmode,fixunit,laba,labb,labc){		

r1=ncol(A)	
n=nrow(A)	
r2=ncol(B)		
m=nrow(B)
r3=ncol(C)
p=nrow(C)
	
# adjustment of aspectratio for printer 
if (fixmode==1){
	Gmat=matrix(t(K[fixunit,]),r2,r3)
    Smat=B%*%Gmat%*%t(C)
	SVD=svd(Smat)
	fit=SVD$d[1]^2+SVD$d[2]^2				
	fit=fit/SUM(Smat)$ssq*100
	Bmat=(m/p)^.25*SVD$u[,1:2]%*%diag(SVD$d[1:2])^.5
	Cmat=(p/m)^.25*SVD$v[,1:2]%*%diag(SVD$d[1:2])^.5
	xmax=max(Bmat,Cmat)+.1
	xmin=min(Bmat,Cmat)-.1
	X=Bmat
	id=labb
	plot(X[,1],X[,2],xlim=c(xmin,xmax),ylim=c(xmin,xmax),type='n',xlab="",ylab="")
	text(X[,1],X[,2],id,pos=1,cex=0.8)			
	X=Cmat
	id=labc
	text(X[,1],X[,2],id,pos=1,cex=0.8)	
}
if (fixmode==2){
	K=permnew(K,r1,r2,r3)
	Gmat=matrix(t(K[fixunit,]),r3,r1)
	Smat=C%*%Gmat%*%t(A)
	SVD=svd(Smat)
	fit=SVD$d[1]^2+SVD$d[2]^2
	fit=fit/SUM(Smat)$ssq*100
	Cmat=(p/n)^.25*SVD$u[,1:2]%*%diag(SVD$d[1:2])^.5
	Amat=(n/p)^.25*SVD$v[,1:2]%*%diag(SVD$d[1:2])^.5
	xmax=max(Cmat,Amat)+.1
	xmin=min(Cmat,Amat)-.1
	X=Cmat
	id=labc
	plot(X[,1],X[,2],xlim=c(xmin,xmax),ylim=c(xmin,xmax),type='n',xlab="",ylab="")
	text(X[,1],X[,2],id,pos=1,cex=0.8)			
	X=Amat
	id=laba
	text(X[,1],X[,2],id,pos=1,cex=0.8)	
}
if (fixmode==3){
	K=permnew(K,r1,r2,r3)
	K=permnew(K,r2,r3,r1)
	Gmat=matrix(t(K[fixunit,]),r1,r2)
	Smat=A%*%Gmat%*%t(B)
	SVD=svd(Smat)
	fit=SVD$d[1]^2+SVD$d[2]^2
	fit=fit/SUM(Smat)$ssq*100
	Amat=(p/n)^.25*SVD$u[,1:2]%*%diag(SVD$d[1:2])^.5
	Bmat=(n/p)^.25*SVD$v[,1:2]%*%diag(SVD$d[1:2])^.5
	xmax=max(Amat,Bmat)+.1
	xmin=min(Amat,Bmat)-.1
	X=Amat
	id=laba
	plot(X[,1],X[,2],xlim=c(xmin,xmax),ylim=c(xmin,xmax),type='n',xlab="",ylab="")
	text(X[,1],X[,2],id,pos=1,cex=0.8)			
	X=Bmat
	id=labb
	text(X[,1],X[,2],id,pos=1,cex=0.8)	
}
return(fit)
}
