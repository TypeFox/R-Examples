renormsolCP <-
function(A,B,C,mode){

r=ncol(A)
if (mode==1){
	ssb=diag(SUM(B)$col^(-.5),nrow=r)
	ssc=diag(SUM(C)$col^(-.5),nrow=r)
	B=B%*%ssb
	C=C%*%ssc
	A=A%*%solve(ssb)%*%solve(ssc)
}
if (mode==2){
	ssa=diag(SUM(A)$col^(-.5),nrow=r)		
	ssc=diag(SUM(C)$col^(-.5),nrow=r)
	A=A%*%ssa
	C=C%*%ssc
	B=B%*%solve(ssa)%*%solve(ssc)
}
if (mode==3){
	ssa=diag(SUM(A)$col^(-.5),nrow=r)		
	ssb=diag(SUM(B)$col^(-.5),nrow=r)
	A=A%*%ssa
	B=B%*%ssb
	C=C%*%solve(ssa)%*%solve(ssb)
}
out=list()
out$A=A
out$B=B
out$C=C
return(out)
}
