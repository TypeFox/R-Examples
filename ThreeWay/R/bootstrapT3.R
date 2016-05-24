bootstrapT3 <-
function(X,A,B,C,G,n,m,p,r1,r2,r3,conv,centopt,normopt,optimalmatch,laba,labb,labc){

narg=nargs()
if (narg<16){
	laba=paste("a",1:n,sep="")
	labb=paste("b",1:m,sep="")
	labc=paste("c",1:p,sep="")
}
X=as.matrix(X)
ccc=0
while (ccc==0){
    cat("How many bootstrap samples do you want to draw? (1000 default)",fill=TRUE)
	n_boot=scan("",n=1)
	if (length(n_boot)==0) {
		n_boot=1000
	}
    if ((floor(n_boot)-n_boot)!=0){
        cat("Error: Select the number of bootstrap samples,fill=TRUE")
    } else{ccc=1}
}

# setting up supermatrices for collecting bootstrap matrices
BB=matrix(0,m*r2,n_boot)
CC=matrix(0,p*r3,n_boot)
GG=matrix(0,r1*r2*r3,n_boot)
fpfp=matrix(0,1,n)
BBB=matrix(0,m*r2,n)
CCC=matrix(0,p*r3,n)
GGG=matrix(0,r1*r2*r3,n)
fpfpfp=matrix(0,1,n)
for (i in 1:n_boot){
	if (i/50-floor(i/50)==0){
		cat(paste("Bootstrap run ", i),fill=TRUE)
	}
    b=t(ceiling((matrix(runif(n*1,0,1),n,1))*n))        # bootstrap sample
    Xb=X[b,]
    n_bootsample=n
    Astart=A

    # preprocessing
    Xprepb=Xb
    cc=centopt
    if ((cc==1) | (cc==12) | (cc==13)){
        Xprepb=cent3(Xprepb,n_bootsample,m,p,1)
    }
    if ((cc==2) | (cc==12) | (cc==23)){
        Xprepb=cent3(Xprepb,n_bootsample,m,p,2)
    }
    if ((cc==3) | (cc==13) | (cc==23)){
        Xprepb=cent3(Xprepb,n_bootsample,m,p,3)
    }
    cc=normopt
    if (cc==1){
        Xprepb=norm3(Xprepb,n_bootsample,m,p,1)
    }
    if (cc==2 ){
        Xprepb=norm3(Xprepb,n_bootsample,m,p,2)
    }
    if (cc==3){
        Xprepb=norm3(Xprepb,n_bootsample,m,p,3)
    }

    # Tucker3, use two runs: rational start and sample solution start
    start=2                                   # start from sample solution
    T3ABKrep=T3funcrep(Xprepb,n_bootsample,m,p,r1,r2,r3,start,conv,Astart,B,C,G)
    start=0                                   # rational start
    T3ABKrep2=T3funcrep(Xprepb,n_bootsample,m,p,r1,r2,r3,start,conv,Astart,B,C,G)
	# resize the A-mode if n>m*p

    if (T3ABKrep2$f<T3ABKrep$f){
        fb=T3ABKrep2$f
        Bb=T3ABKrep2$B
        Cb=T3ABKrep2$C
        Gb=T3ABKrep2$H
        fpb=T3ABKrep2$fp
    } else{
        fb=T3ABKrep$f
        Bb=T3ABKrep$B
        Cb=T3ABKrep$C
        Gb=T3ABKrep$H
        fpb=T3ABKrep$fp
    }

    if (optimalmatch==0){
        # matching via orthogonal rotation towards full solution
        #minimize || Bb T - B ||, || Cb U - C ||, || S' Gb - G ||.
        # using svd B'Bb = PDQ', then T=QP', etc.
        SVD=svd(t(B)%*%Bb)
        T=SVD$v%*%t(SVD$u)
        SVD=svd(t(C)%*%Cb)
        U=SVD$v%*%t(SVD$u)
        Bb=Bb%*%T
        Cb=Cb%*%U
        Gb=Gb%*%kronecker(U,T)
        SVD=svd(G%*%t(Gb))
        S=SVD$v%*%t(SVD$u)
        Gb=t(S)%*%Gb
    }
    if (optimalmatch==1){
        # matching via optimal transformation towards full solution
        # minimize || Bb T - B ||, || Cb U - C ||, || S Gb - G ||.
        # using linear regressions
        T=solve(t(Bb)%*%Bb)%*%t(Bb)%*%B
        U=solve(t(Cb)%*%Cb)%*%t(Cb)%*%C
        Bb=Bb%*%T
        Cb=Cb%*%U
        Gb=Gb%*%kronecker(solve(t(U)),solve(t(T)))
        S=G%*%t(Gb)%*%solve(Gb%*%t(Gb))
        Gb=S%*%Gb
    }

    # collect bootstrap solutions
	BB[,i]=Bb
	CC[,i]=Cb
	GG[,i]=Gb
	fpfp[i]=fpb
}
BB=t(BB)
CC=t(CC)
GG=t(GG)
se_B=matrix(0,m,r2)
nB=m*r2
for (i in 1:nB){
    se_B[i]=sqrt(var(BB[,i])*(nB-1)/nB)
}
nC=p*r3
se_C=matrix(0,p,r3)
for (i in 1:nC){
    se_C[i]=sqrt(var(CC[,i])*(nC-1)/nC)
}
nG=r1*r2*r3
se_G=matrix(0,r1,r2*r3)
for (i in 1:nG){
    se_G[i]=sqrt(var(GG[,i])*(nG-1)/nG)
}
se_fp=sqrt(var(fpfp)*(length(fpfp)-1)/length(fpfp))

m_B=matrix(SUM(BB)$mc,m,r2)
m_C=matrix(SUM(CC)$mc,p,r3)
m_G=matrix(SUM(GG)$mc,r1,r2*r3)
m_fp=mean(fpfp)

# make labels 
labCompA=paste("A",1:r1,sep="")
strB=noquote(vector(mode="character",length=2*r2))
strC=noquote(vector(mode="character",length=2*r3))
strG=noquote(vector(mode="character",length=2*r2*r3))
i=1
for (j in 1:r2){
	strB[i]=noquote(paste("LB B",as.character(j), sep=""))
	i=i+1
	strB[i]=noquote(paste("UB B",as.character(j), sep=""))
	i=i+1
}
i=1
for (k in 1:r3){
	strC[i]=noquote(paste("LB C",as.character(k),sep=""))
	i=i+1
	strC[i]=noquote(paste("UB C",as.character(k),sep=""))
	i=i+1
}
i=1
for (k in 1:r3){
	for (j in 1:r2){
		strG[i]=noquote(paste("LB B",as.character(j),"xC",as.character(k),sep=""))
		i=i+1
		strG[i]=noquote(paste("UB B",as.character(j),"xC",as.character(k),sep=""))
		i=i+1
	}
}

lo_B=matrix(percentile95(BB)$lo,m,r2)
up_B=matrix(percentile95(BB)$up,m,r2)
Bint=ccmat(lo_B,up_B)
rownames(Bint)=labb
colnames(Bint)=strB
lo_C=matrix(percentile95(CC)$lo,p,r3)
up_C=matrix(percentile95(CC)$up,p,r3)
Cint=ccmat(lo_C,up_C)
rownames(Cint)=labc
colnames(Cint)=strC
lo_G=matrix(percentile95(GG)$lo,r1,r2*r3)
up_G=matrix(percentile95(GG)$up,r1,r2*r3)
Gint=ccmat(lo_G,up_G)
rownames(Gint)=labCompA
colnames(Gint)=strG
lo_fp=quantile(fpfp,0.025) 
up_fp=quantile(fpfp,0.975)
fpint=c(lo_fp,up_fp)
names(fpint)=c("LB Fit (%)", "UB Fit (%)")

out=list()
out$Bint=Bint
out$Cint=Cint
out$Gint=Gint
out$fpint=fpint
return(out)
}
