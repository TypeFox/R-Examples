StratSel <-
function(formula,corr=TRUE, Startval, optim.method="BFGS", data, ...){
	user.supplied.startval <- as.numeric(missing(Startval)!=TRUE)

	mf <- match.call(expand.dots = FALSE)
	
	#print(mf)
	m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
	#print(m)
	mf <- mf[c(1, m)]

	f <- Formula(formula)
	mf[[1]] <- as.name("model.frame")
	mf$formula <- f
	mf <- eval(mf, parent.frame())
	#print(names(mf))
		
	Y <- model.response(mf)
	X11 <- model.matrix(f, data = mf, rhs = 1)
	X14 <- model.matrix(f, data = mf, rhs = 2)
	X24 <- model.matrix(f, data = mf, rhs = 3)

	#print(colnames(X24))

	ifelse(dim(as.matrix(Y))[2]==2, ymatrix <- TRUE, ymatrix <- FALSE)
	if(ymatrix==TRUE){
		yy <- rep(NA,dim(Y)[1])
		yy[Y[,1]==0] <- 1	
		yy[Y[,2]==0] <- 3
		yy[Y[,2]==1] <- 4
		}
	if(ymatrix==FALSE) yy <- Y
	
	dim.x11 <- dim(as.matrix(X11))
	dim.x14 <- dim(as.matrix(X14))
	dim.x24 <- dim(as.matrix(X24))

	
	Startval <- gen.Startval(Startval, user.supplied.startval, corr, ys=yy, xs11=X11, xs14=X14, xs24=X24, dim.x11=dim.x11, dim.x14=dim.x14, dim.x24=dim.x24)
	
	#print(Startval)
	
	if (corr==TRUE) out.struc <- optim(Startval, logLikStratSel, x11=X11, x14=X14, x24=X24, y=Y, hessian=TRUE, method=optim.method, ...)
	if (corr==FALSE) out.struc <- optim(Startval, logLikStrat, x11=X11, x14=X14, x24=X24, y=Y, hessian=TRUE, method=optim.method, ...)

	#Names

		if (corr==TRUE)	names.out <- c("u11 constant", paste("u11",colnames(X11)[-1]), "u14 constant" ,paste("u14",colnames(X14)[-1]),"u24 constant", paste("u24",colnames(X24)[-1]), "rho (errors)")

		if (corr==FALSE)	names.out <- c("u11 constant", paste("u11",colnames(X11)[-1]), "u14 constant" ,paste("u14",colnames(X14)[-1]),"u24 constant", paste("u24",colnames(X24)[-1]))


		
	names(out.struc$par) <-	names.out 
	colnames(out.struc$hessian) <- names.out
	
	vcov <- solve(out.struc$hessian)
	vcov[length(out.struc$par),length(out.struc$par)]<- fetch.rho.v(v=vcov, b=out.struc$par) 
	colnames(vcov) <- names.out  
	
	beta <- as.matrix(out.struc$par)
	beta <- fetch.rho.b(b=beta)

	ll <- out.struc$value

	nits <- out.struc$counts[2]
	conv <- out.struc$convergence
	
	#names(nits) <- "Number of Iterations"
	
	#print(nits)
		
	if (corr==TRUE)	est <- list(coefficients=beta,
			vcov=vcov,
			corr= out.struc$par[length(out.struc$par)],
			df = nrow(X11),
			logLik=-ll,
			model = mf,
			call=match.call(),
			formula= f,
			y=yy,
			DIM=c(dim.x11,dim.x14,dim.x24),
			nits=as.numeric(nits),
			conv=conv)

	if (corr==FALSE) {	est <- list(coefficients=beta,
			vcov=vcov,
			corr= FALSE,
			df = nrow(X11),
			logLik=-ll,
			model = mf,
			call=match.call(),
			formula= f,
			y=yy,
			DIM=c(dim.x11,dim.x14,dim.x24),
			nits=as.numeric(nits),
			conv=conv)
			}
	
	class(est) <- "StratSel"
	
	return(est)

}
