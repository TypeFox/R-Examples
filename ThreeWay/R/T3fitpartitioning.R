T3fitpartitioning <-
function(Xprep,n,m,p,AS,BT,CU,K,renormmode,laba,labb,labc){	

narg=nargs()
if (narg<10){
	laba=paste("a",1:n,sep="")
	labb=paste("b",1:m,sep="")
	labc=paste("c",1:p,sep="")
}

Xprep=as.matrix(Xprep)
r1=ncol(AS)
r2=ncol(BT)
r3=ncol(CU)
Res=Xprep-(AS%*%K%*%kronecker(t(CU),t(BT)))
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

BCcontr=NULL
ACcontr=NULL
ABcontr=NULL
if (renormmode==1){
	# compute fit contribution to total fit% of each B&C-mode component combination
	F=AS%*%K
	BCcontr=colSums(F^2)/SUM(Xprep)$ssq*100
	strCB=noquote(vector(mode="character",length=r2*r3))
	i=1
	for (k in 1:r3){
		for (j in 1:r2){
			strCB[i]=noquote(paste(" B",as.character(j),"xC",as.character(k),sep=""))
			i=i+1
		}
	}
	BCcontr=cbind(BCcontr)
	rownames(BCcontr)=strCB
}
if (renormmode==2){
	# compute fit contribution to total fit% of each A&C-mode component combination
	F=BT%*%permnew(K,r1,r2,r3)
	ACcontr=colSums(F^2)/SUM(Xprep)$ssq*100
	strAC=noquote(vector(mode="character",length=r1*r3))
	i=1
	for (k in 1:r1){
		for (j in 1:r3){
			strAC[i]=noquote(paste(" C",as.character(j),"xA",as.character(k),sep=""))
			i=i+1
		}
	}
	ACcontr=cbind(ACcontr)
	rownames(ACcontr)=strAC
}
if (renormmode==3){
	# compute fit contribution to total fit% of each A&B-mode component combination
	F=CU%*%permnew(permnew(K,r1,r2,r3),r2,r3,r1)
	ABcontr=colSums(F^2)/SUM(Xprep)$ssq*100
	strBA=noquote(vector(mode="character",length=r1*r2))
	i=1
	for (k in 1:r2){
		for (j in 1:r1){
			strBA[i]=noquote(paste(" A",as.character(j),"xB",as.character(k),sep=""))
			i=i+1
		}
	}
	ABcontr=cbind(ABcontr)
	rownames(ABcontr)=strBA
}
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
out$ABcontr=ABcontr
out$ACcontr=ACcontr
out$BCcontr=BCcontr
return(out)
}
