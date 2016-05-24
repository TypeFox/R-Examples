bothsidesmodel.hotelling <-
function(x,y,z,rows,cols) {
	bsm <- bothsidesmodel(x,y,z)
	lstar <- length(cols)
	pstar <- length(rows)
	nu <- bsm$df[1]
	bstar <- bsm$Beta[rows,cols]
	if(lstar==1) bstar <- matrix(bstar,ncol=1)
	if(pstar==1) bstar <- matrix(bstar,nrow=1)
	W.nu <- bsm$Sigmaz[cols,cols]
	cx <- solve(t(x)%*%x)
	B <- t(bstar)%*%solve(cx[rows,rows])%*%bstar
	t2 <- tr(solve(W.nu)%*%B)
	f <- (nu-lstar+1)*t2/(lstar*pstar*nu)
	df <- c(lstar*pstar,nu-lstar+1)
	W <- W.nu*nu
	lambda <- ifelse(lstar==1,W/(W+B),det(W)/det(W+B))
	chis <- -(nu-(lstar-pstar+1)/2)*log(lambda)
	Hotelling <- list(T2 = t2, F = f, df = df,pvalue = 1-pf(f,df[1],df[2]))
	Wilks <- list(Lambda=lambda,Chisq=chis,df=df[1],pvalue=1-pchisq(chis,df[1]))
	list(Hotelling = Hotelling,Wilks = Wilks)
}
