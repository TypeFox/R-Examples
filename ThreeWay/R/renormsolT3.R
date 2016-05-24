renormsolT3 <-
function(A,B,C,G,mode){

if (mode==1){
	ss=diag(SUM(G)$row^.5,nrow=ncol(A))		
	A=A%*%ss
	G=solve(ss)%*%G				
}
if (mode==2){
	G=permnew(G,ncol(A),ncol(B),ncol(C))
	ss=diag(SUM(G)$row^.5,nrow=ncol(B))
	B=B%*%ss
	G=solve(ss)%*%G
	G=permnew(G,ncol(B),ncol(C),ncol(A))
	G=permnew(G,ncol(C),ncol(A),ncol(B))
}
if (mode==3){
	G=permnew(G,ncol(A),ncol(B),ncol(C))
	G=permnew(G,ncol(B),ncol(C),ncol(A))
	ss=diag(SUM(G)$row^.5,nrow=ncol(C))
	C=C%*%ss
	G=solve(ss)%*%G
	G=permnew(G,ncol(C),ncol(A),ncol(B))
}
out=list()
out$A=A
out$B=B
out$C=C
out$H=G
return(out)
}
