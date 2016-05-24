CPfitpartitioning <-
function(Xprep,n,m,p,A,B,C,laba,labb,labc){	

narg=nargs()
if (narg<8){
	laba=paste("a",1:n,sep="")
	labb=paste("b",1:m,sep="")
	labc=paste("c",1:p,sep="")
}

Xprep=as.matrix(Xprep)
r=ncol(A)
H=matrix(0,r,r^2)       # superidentity 3-way array
for (ii in 1:r){
	H[ii,(ii-1)*r+ii]=1
}

Res=Xprep-(A%*%H%*%kronecker(t(C),t(B)))
# A-mode
fitA=100-colSums(t(Res)^2)/colSums(t(Xprep)^2)*100		
# B-mode
Res=permnew(Res,n,m,p)
Xprep=permnew(Xprep,n,m,p)
fitB=100-colSums(t(Res)^2)/colSums(t(Xprep)^2)*100
# C-mode
Res=permnew(Res,m,p,n)
Xprep=permnew(Xprep,m,p,n)
fitC=100-colSums(t(Res)^2)/colSums(t(Xprep)^2)*100
Xprep=permnew(Xprep,p,n,m)
Res=permnew(Res,p,n,m)
fitA=cbind(fitA)
fitB=cbind(fitB)
fitC=cbind(fitC)
rownames(fitA)=laba
rownames(fitB)=labb
rownames(fitC)=labc
out=list()
out$fitA=fitA
out$fitB=fitB
out$fitC=fitC
return(out)
}
