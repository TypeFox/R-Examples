SUM <-
function(A){

A=as.matrix(A)
nc=ncol(A)
nr=nrow(A)
sc=c(array(0,nc))
sr=c(array(0,nr))
mc=c(array(0,nc))
mr=c(array(0,nr))
max=rep(0,nc)
min=rep(0,nc)
Vmaxc=max.col(t(A),ties.method="first")
Vmaxr=max.col(A,ties.method="first")
Vminc=vector(mode="numeric",nc)
Vminr=vector(mode="numeric",nr)
cont=0
for (k in 1:nc){
	cont=cont+1
	sc[k]=sum(A[,k]^2)
	mc[k]=sum(A[,k])/nr
	min[k]=min(A[,k])
	max[k]=max(A[,k])
	Vminc[k]=which.min(A[,k])
	k=k+1
}
for (k in 1:nr){
	sr[k]=sum(A[k,]^2)
	mr[k]=sum(A[k,])/nc
	Vminr[k]=which.min(A[k,])
	k=k+1
}
cumC=A
if (nr>1){
	for (k in 1:(nr-1)){
		cumC[k+1,]=A[k+1,]+cumC[k,]
		k=k+1
	}
}
cumR=(A)
if (nc>1){
	for (k in 1:(nc-1)){
		cumR[,k+1]=A[,k+1]+cumR[,k]
		k=k+1
	}
}
t=sum(A^2)

out=list()
out$row=sr
out$col=sc
out$mr=mr
out$mc=mc
out$minc=min
out$maxc=max
out$valueMinr=Vminr
out$valueMinc=Vminc
out$valueMaxr=Vmaxr
out$valueMaxc=Vmaxc
out$ssq=t
out$cumsumr=cumR
out$cumsumc=cumC
return(out)
}
