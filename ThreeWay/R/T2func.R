T2func <-
function(X,n,m,p,r1,r2,r3,start,conv,model,A,B,C,H){

X=as.matrix(X)
if (model==1){
	C=diag(r3)
}
if (model==2){
	B=diag(r2)
}
if (model==3){
	A=diag(r1)
}
cputime=system.time({
	# initialize A, B and C
	ss=sum(X^2)
	dys=0

	if (start==0){
		cat("Rational ORTHONORMALIZED start",fill=TRUE)
		# rational starts via eigendecompositions
		if (model!=3){
			EIG=eigen(X%*%t(X))	
			A=EIG$vectors[,1:r1]
		}
		Z=permnew(X,n,m,p)		# yields m x p x n array
		if (model!=2){
			EIG=eigen(Z%*%t(Z))
			B=EIG$vectors[,1:r2]
		}
		Z=permnew(Z,m,p,n)		# yields p x n x m array
		if (model!=1){
			EIG=eigen(Z%*%t(Z))
			C=EIG$vectors[,1:r3]
		}
	}

	if (start==1){
		cat("Random ORTHONORMALIZED starts",fill=TRUE)
		if (model!=3){
			if (n>=r1){
				A=orth(matrix(runif(n*r1,0,1),n,r1)-.5)
			} else{
				A=orth(matrix(runif(r1*r1,0,1),r1,r1)-.5)
				A=A[1:n,]
			}
		}
		if (model!=2){
			if (m>=r2){
				B=orth(matrix(runif(m*r2,0,1),m,r2)-.5)
			} else{
				B=orth(matrix(runif(r2*r2,0,1),r2,r2)-.5)
				B=B[1:m,]
			}
		}
		if (model!=1){
			if (p>=r3){
				C=orth(matrix(runif(p*r3,0,1),p,r3)-.5)
			} else{
				C=orth(matrix(runif(r3*r3,0,1),r3,r3)-.5)
				C=C[1:p,]
			}
		}
	}

	# Update Core
	if (start!=2){
		Z=permnew(t(A)%*%X,r1,m,p)
		Z=permnew(t(B)%*%Z,r2,p,r1)
		H=permnew(t(C)%*%Z,r3,r1,r2)    
	}

	# Evaluate f
	if (start==2){
		Z=B%*%permnew(A%*%H,n,r2,r3)
		Z=C%*%permnew(Z,m,r3,n)
		Z=permnew(Z,p,n,m)			# Z = Xhat, nxmxp
		f=sum((X-Z)^2)				# use full formula, taking into account possibility of nonoptimal core in start
	} else{
		f=ss-sum(H^2)
	}

	cat(paste("Tucker2 function value at start is ",f),fill=TRUE)
	iter=0
	fold=f+2*conv*f
	while (fold-f>f*conv){
		iter=iter+1
		fold=f

		if (model!=3){
			# update A   (Z=X*C'x B' - GS Z*Z'*A)
			Z=permnew(X,n,m,p)
			Z=permnew(t(B)%*%Z,r2,p,n)
			Z=permnew(t(C)%*%Z,r3,n,r2)			 # yields n x r2 x r3 array
			A=qr.Q(qr(Z%*%(t(Z)%*%A)),complete=FALSE)
		}

		if (model!=2){
			# update B
			Z=permnew(X,n,m,p)
			Z=permnew(Z,m,p,n)
			Z=permnew(t(C)%*%Z,r3,n,m)
			Z=permnew(t(A)%*%Z,r1,m,r3)			 # yields m x r3 x r1 array
			B=qr.Q(qr(Z%*%(t(Z)%*%B)),complete=FALSE)
		}
		
		if (model!=1){
			# update C
			Z=permnew(t(A)%*%X,r1,m,p)
			Z=permnew(t(B)%*%Z,r2,p,r1)			 # yields p x r1 x r2 array
			C=qr.Q(qr(Z%*%(t(Z)%*%C)),complete=FALSE)
		}
		
		# Update Core
		Z=permnew(t(A)%*%X,r1,m,p)
		Z=permnew(t(B)%*%Z,r2,p,r1)
		H=permnew(t(C)%*%Z,r3,r1,r2)

		# Evaluate f
		f=ss-sum(H^2)
		if ((iter%%10)==0){					
			cat(paste("Tucker2 function value after iteration ",iter," is ",f),fill=TRUE)
		}
	}
})
ss=sum(X^2)
fp=100*(ss-f)/ss

# compute "intrinsic eigenvalues"
# eigenvalues for A-mode:
La=H%*%t(H)
Y=permnew(H,r1,r2,r3)
Lb=Y%*%t(Y)
Y=permnew(Y,r2,r3,r1)
Lc=Y%*%t(Y)
cat(paste("Tucker2 function value is",f,"after",iter,"iterations", sep=" "),fill=TRUE)
cat(paste("Fit percentage is",fp,"%",sep=" "),fill=TRUE)
cat(paste("Procedure used",(round(cputime[1],2)),"seconds", sep=" "),fill=TRUE)

out=list()
out$A=A
out$B=B
out$C=C
out$H=H
out$f=f
out$fp=fp
out$iter=iter
out$cputime=cputime[1]
out$La=La
out$Lb=Lb
out$Lc=Lc
return(out)
}
