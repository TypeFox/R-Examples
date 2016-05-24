pls1_nipals <-
function (X,y,a,it=50,tol=1e-8,scale=FALSE) 
{
Xh <- scale(X,center=TRUE,scale=scale)		#mean-centering of data matrix X
yh <- scale(y,center=TRUE,scale=scale)		#mean-centering of data y

T <- NULL
P <- NULL
C <- NULL
W <- NULL
for (h in 1:a){
	wh <- t(Xh)%*%yh				#LS regression for wh
	wh <- wh/as.vector(sqrt(t(wh)%*%wh))		#normalization of wh
	th <- Xh%*%wh				 	#LS regression for th
	ch <- as.numeric(t(yh)%*%th)			#LS regression for ch
	ch <- ch/as.vector(sqrt(t(th)%*%th))		#normalization of ch
	#ch <- ch/as.vector(sqrt(t(ch)%*%ch))		#normalization of ch
	ph <- t(Xh)%*%th/as.vector(t(th)%*%th)		#LS regression for ph
	Xh <- Xh - th%*%t(ph)				#calculate new Xh
	yh <- yh - th*ch				#calculate new yh
	T <- cbind(T,th)	#build matrix T
	P <- cbind(P,ph)	#build matrix P
	C <- c(C,ch)		#build vector C
	W <- cbind(W,wh)	#build matrix W
}
	# final regression coefficients:
	b <- W%*%solve(t(P)%*%W)%*%C
list(P=P,T=T,W=W,C=C,b=b)
}

