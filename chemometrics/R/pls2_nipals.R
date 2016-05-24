pls2_nipals <-
function (X,Y,a,it=50,tol=1e-8,scale=FALSE) 
{
Xh <- scale(X,center=TRUE,scale=scale)		#mean-centering of data matrix X
Yh <- scale(Y,center=TRUE,scale=scale)		#mean-centering of data matrix Y

T <- NULL
W <- NULL
Q <- NULL
U <- NULL
P <- NULL
D <- NULL
C <- NULL
W <- NULL
for (h in 1:a){
	nr <- 0
	uh <- Yh[,1]		#starting value for uh is 1st column of Yh
	ende <- FALSE
	#inner steps of PLS2 algorithm
	while (!ende){
		nr <- nr+1
		wh <- t(Xh)%*%uh				#LS regression for wh
		wh <- wh/as.vector(sqrt(t(wh)%*%wh))		#normalization of wh
		th <- Xh%*%wh				 	#LS regression for th
		ch <- t(Yh)%*%th				#LS regression for ch
		ch <- ch/as.vector(sqrt(t(ch)%*%ch))		#normalization of ch
		uhnew <- Yh%*%ch 				#LS regression for uh
		deltau <- uhnew-uh
		unorm <- as.numeric(sqrt(t(deltau)%*%deltau))

		if (unorm<tol) {
			ende <- TRUE
		}
		else if (nr > it) {	#too many iterations
			ende <- TRUE
			cat("\nWARNING! Iteration stop in h=",h," without convergence!\n\n")
		}
		uh <- uhnew
	}
	ph <- t(Xh)%*%th/as.vector(t(th)%*%th)			#LS regression for ph
	qh <- t(Yh)%*%uh/as.vector(t(uh)%*%uh)			#LS regression for qh
	dh <- t(uh)%*%th/as.vector(t(th)%*%th)			#LS regression for dh
	Xh <- Xh - th%*%t(ph)					#calculate new Xh
	Yh <- Yh - (th%*%t(ch))*as.vector(dh)			#calculate new Yh

	T <- cbind(T,th)	#build matrix T
	Q <- cbind(Q,qh)	#build matrix Q
	U <- cbind(U,uh)	#build matrix U
	P <- cbind(P,ph)	#build matrix P
	D <- c(D,dh)		#build vector of diagonal elements of matrix D
	C <- cbind(C,ch)	#build matrix C
	W <- cbind(W,wh)	#build matrix W

	# final regression coefficients:
	B <- W%*%solve(t(P)%*%W)%*%t(C)
}
list(P=P,T=T,Q=Q,U=U,D=D,W=W,C=C,B=B)
}

