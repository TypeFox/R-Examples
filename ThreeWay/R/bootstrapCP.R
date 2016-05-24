bootstrapCP <-
function(X,A,B,C,n,m,p,r,ort1,ort2,ort3,conv,centopt,normopt,scaleopt,maxit,laba,labb,labc){

narg=nargs()
if (narg<17){
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
BB=matrix(0,m*r,n_boot)
CC=matrix(0,p*r,n_boot)
fpfp=matrix(0,1,n)
BBB=matrix(0,m*r,n)
CCC=matrix(0,p*r,n)
fpfpfp=matrix(0,1,n)
for (i in 1:n_boot){                # +n
	if (i/50-floor(i/50)==0){
		cat(paste("Bootstrap run ", i),fill=TRUE)
	}
    b=t(ceiling((matrix(runif(n*1,0,1),n,1))*n))        # bootstrap sample
    Xb=X[b,]
    n_bootsample=n

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

    # Parafac, use two runs: rational start and sample solution start
    start=2 # start from sample solution
    CPmatBC=CPfuncrep(Xprepb,n_bootsample,m,p,r,ort1,ort2,ort3,start,conv,maxit,A,B,C)
    start=0 # rational start
    CPmatBC2=CPfuncrep(Xprepb,n_bootsample,m,p,r,ort1,ort2,ort3,start,conv,maxit,A,B,C)

    if (CPmatBC2$f<CPmatBC$f){
        fb=CPmatBC2$f
        Ab=CPmatBC2$A
        Bb=CPmatBC2$B
        Cb=CPmatBC2$C
        fpb=CPmatBC2$fp
    } else{
        fb=CPmatBC$f
        Ab=CPmatBC$A
        Bb=CPmatBC$B
        Cb=CPmatBC$C
        fpb=CPmatBC$fp
    }

	if ((scaleopt==1) | (scaleopt==2) | (scaleopt==3)){
		RNsol=renormsolCP(A,B,C,scaleopt)
		Ab=RNsol$A
		Bb=RNsol$B
		Cb=RNsol$C
	}
	
    # permutation
    if (r>1){
		coeff=rep(0,factorial(r))
		permutation=perms(r)
		Bbperm=matrix(,factorial(r),r)
		Cbperm=matrix(,factorial(r),r)
		for (j in 1:nrow(permutation)){
			Bbperm=Bb[,permutation[j,]]
			Cbperm=Cb[,permutation[j,]]
			for (r1 in 1:r){
				coeff[j]=coeff[j]+abs(phi(B[,r1],Bbperm[,r1]))*abs(phi(C[,r1],Cbperm[,r1]))
			}
		}
		Bopt=Bb[,permutation[which(coeff==max(coeff)),]]
		Copt=Cb[,permutation[which(coeff==max(coeff)),]]
	} else{
		Bopt=Bb
		Copt=Cb
	}
	# reflection
    for (r1 in 1:r){
        if (phi(Bopt[,r1],B[,r1])<0){
		 	Bopt[,r1]=-Bopt[,r1]
			}
        if (phi(Copt[,r1],C[,r1])<0){
			Copt[,r1]=-Copt[,r1]
    	}
	}

    # collect bootstrap solutions
	BB[,i]=Bopt
	CC[,i]=Copt
	fpfp[i]=fpb
}
BB=t(BB)
CC=t(CC)
se_B=matrix(0,m,r)
nB=m*r
for (i in 1:nB){
    se_B[i]=sqrt(var(BB[,i])*(nB-1)/nB)
}
nC=p*r
se_C=matrix(0,p,r)
for (i in 1:nC){
    se_C[i]=sqrt(var(CC[,i])*(nC-1)/nC)
}
se_fp=sqrt(var(fpfp)*(length(fpfp)-1)/length(fpfp))

m_B=matrix(SUM(BB)$mc,m,r)
m_C=matrix(SUM(CC)$mc,p,r)
m_fp=mean(fpfp)

# make labels 
str=noquote(vector(mode="character",length=2*r))
i=1
for (j in 1:r){
	str[i]=noquote(paste("LB Comp.",as.character(j), sep=""))
	i=i+1
	str[i]=noquote(paste("UB Comp.",as.character(j), sep=""))
	i=i+1
}
lo_B=matrix(percentile95(BB)$lo,m,r)
up_B=matrix(percentile95(BB)$up,m,r)
Bint=ccmat(lo_B,up_B)
rownames(Bint)=labb
colnames(Bint)=str
lo_C=matrix(percentile95(CC)$lo,p,r)
up_C=matrix(percentile95(CC)$up,p,r)
Cint=ccmat(lo_C,up_C)
rownames(Cint)=labc
colnames(Cint)=str
lo_fp=quantile(fpfp,0.025) 
up_fp=quantile(fpfp,0.975)
fpint=c(lo_fp,up_fp)
names(fpint)=c("LB Fit (%)", "UB Fit (%)")

out=list()
out$Bint=Bint
out$Cint=Cint
out$fpint=fpint
return(out)
}
