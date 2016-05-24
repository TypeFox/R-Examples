CPfunc <-
function(X,n,m,p,r,ort1,ort2,ort3,start,conv,maxit,A,B,C){

X=as.matrix(X)
ftiter=matrix(0,maxit/10,2)
mintripcos=0
cputime=system.time({
    ssx=sum(X^2)

    if (start==0){
		# rational starts via eigendecompositions, or random if r is too big
        if (n>=r){
            AUT=eigen(X%*%t(X))
            A=AUT$vectors[,1:r]
		} else{
			A=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
	        A=A[1:n,]
        }
        Z=permnew(X,n,m,p)          # yields m x p x n array
        if (m>=r){
            AUT=eigen(Z%*%t(Z))
            B=AUT$vectors[,1:r]
        } else{
            B=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
            B=B[1:m,]
        }
        Z=permnew(Z,m,p,n)          # yields p x n x m
        if (p>=r){
            AUT=eigen(Z%*%t(Z))
            C=AUT$vectors[,1:r]
        } else{
            C=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
            C=C[1:p,]
        }
    }
    if (start==1){
        if (n>=r){
            A=orth(matrix(runif(n*r,0,1),nrow=n)-0.5)
        } else{
            A=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
            A=A[1:n,]
        }
        if (m>=r){
            B=orth(matrix(runif(m*r,0,1),nrow=m)-0.5)
        } else{
            B=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
            B=B[1:m,]
        }
        if (p>=r){
            C=orth(matrix(runif(p*r,0,1),nrow=p)-0.5)
        } else{
            C=orth(matrix(runif(r*r,0,1),nrow=r)-0.5)
            C=C[1:p,]
        }
    }

	H=matrix(0,r,r^2)       # superidentity 3-way array
	for (ii in 1:r){
		H[ii,(ii-1)*r+ii]=1
	}
	H1=permnew(H,r,r,r)
	H1=permnew(B%*%H1,m,r,r)
	H1=permnew(C%*%H1,p,r,m)  # gives r x m*p array  H%*%t(C)xt(B)
	f=sum((X-A%*%H1)^2)
	cat(paste("Candecomp/Parafac function value at Start is ",f, sep=" "),fill=TRUE)
	fold=f+2*conv*f
	iter=0
	BB=t(B)%*%B
	CC=t(C)%*%C
	
	while ((fold-f>conv*f | iter<2) & (f>conv^2) & (iter<maxit)){
		fold=f
		#   Update A
		#   first compute XF = X H C'xB' in FF' (implicitly)
		Z1=permnew(X,n,m,p)
		Z1=permnew(t(B)%*%Z1,r,p,n)
		Z1=permnew(t(C)%*%Z1,r,n,r)   # gives n x r^2 array
		XF=Z1%*%t(H)
		if (ort1==1){
			FF=BB*CC
			A=XF%*%solve(FF)
		}
		if (ort1==2){
			SVD=svd(XF)
			A=SVD$u%*%t(SVD$v)
		}
		if (ort1==3){
			FF=BB*CC
			SVD=svd(XF-matrix(1,n,1)%*%apply(XF,2,mean))
			A=SVD$u%*%t(SVD$v)+matrix(1,n,1)%*%apply(XF,2,mean)%*%solve(FF)
		}
		AA=t(A)%*%A

		# Update B
		# first compute XF = X H A'xC' in FF' (implicitly)
		Z=permnew(X,n,m,p)
		Z1=permnew(Z,m,p,n)
		Z1=permnew(t(C)%*%Z1,r,n,m)
		Z1=permnew(t(A)%*%Z1,r,m,r)   # gives m x r^2 array
		XF=Z1%*%t(H)
		if (ort2==1){
			FF=AA*CC
			B=XF%*%solve(FF)
		}
		if (ort2==2){
			SVD=svd(XF)
			B=SVD$u%*%t(SVD$v)
		}
		if (ort2==3){
			FF=AA*CC
			SVD=svd(XF-matrix(1,m,1)%*%apply(XF,2,mean))
			B=SVD$u%*%t(SVD$v)+matrix(1,m,1)%*%apply(XF,2,mean)%*%solve(FF)
		}
		BB=t(B)%*%B

		# Update C
		# first compute XF = X H B'xA' in FF' (implicitly)
		Z=permnew(Z,m,p,n)
		Z1=permnew(Z,p,n,m)
		Z1=permnew(t(A)%*%Z1,r,m,p)
		Z1=permnew(t(B)%*%Z1,r,p,r)    #  gives p x r^2 array
		XF=Z1%*%t(H)
		if (ort3==1){
			FF=AA*BB
			C=XF%*%solve(FF)
		}
		if (ort3==2){
			SVD=svd(XF)
			C=SVD$u%*%t(SVD$v)
		}
		if (ort3==3){
			FF=AA*BB
			SVD=svd(XF-matrix(1,p,1)%*%apply(XF,2,mean))
			C=SVD$u%*%t(SVD$v)+matrix(1,p,1)%*%apply(XF,2,mean)%*%solve(FF)
		}
		CC=t(C)%*%C

		# Evaluate
		if (ort3==1){
			f=ssx-tr(CC%*%FF)
		} else{
			H1=permnew(H,r,r,r)
			H1=permnew(B%*%H1,m,r,r)
			H1=permnew(C%*%H1,p,r,m)  # gives r x m*p array  H*t(C)xt(B)
			f=sum((X-A%*%H1)^2)

		}

		iter=iter+1
		if ((iter%%10)==0){
			tripcos=min(phi(A,A)*phi(B,B)*phi(C,C))
			if (iter==10){
				mintripcos=tripcos
			}
			if (tripcos<mintripcos){
				mintripcos=tripcos
			}
			if ((iter%%1000)==0){
				cat(paste("Minimal Triple cosine =",tripcos,sep=" "),fill=TRUE)
			}
			ftiter[iter/10,]=c(f, tripcos)
		}
		if ((iter%%50)==0){
			cat(paste("f=",f,"after",iter,"iters; diff.=",(fold-f), sep=" "),fill=TRUE)
		}
	}
})
fp=100-100*f/ssx
tripcos=min(phi(A,A)*phi(B,B)*phi(C,C))
names(tripcos)=c("Minimal triple cosine")
if (iter<10){
    mintripcos=tripcos
}
cat(paste("Candecomp/Parafac function value is",f,"after",iter,"iterations", sep=" "),fill=TRUE)
cat(paste("Fit percentage is",fp,"%",sep=" "),fill=TRUE)
cat(paste("Procedure used",(round(cputime[1],2)),"seconds", sep=" "),fill=TRUE)
out=list()
out$A=A
out$B=B
out$C=C
out$f=f
out$fp=fp
out$iter=iter
out$tripcos=tripcos
out$mintripcos=mintripcos
out$ftiter=ftiter
out$cputime=cputime[1]
return(out)
}
