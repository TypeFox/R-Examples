penlogreg <- function(X,Y,tau.o=1,
	a0=NULL,b0=NULL,
	pen.method=c("vectors","simple","none"),
	progress=FALSE)
{
	pen.method <- match.arg(pen.method)
	# Preliminaries:
	n <- dim(X)[1]
	d <- dim(X)[2]
	L <- max(Y)
	# Estimation of standard deviations
	# (within groups) etc.:
	IL <- diag(rep(1,L))
	YY <- IL[Y,]
	nvec <- colSums(YY)
	mu.hat <- diag(1/nvec) %*% crossprod(YY, X)
	Xc <- X - mu.hat[Y,]
	if (is.null(a0) || is.null(b0))
	{
		Sigma.hat <- crossprod(Xc, Xc) / (n - L)
		std.X <- sqrt(diag(Sigma.hat))
		# Determine starting values for
		# logistic regression:
		nvec <- table(Y)
		b0 <- qr.solve(Sigma.hat,t(mu.hat))
		a0 <- log(nvec/n) - diag(mu.hat %*% b0)
		a0 <- a0 - mean(a0)
		b0 <- b0 - tcrossprod(rowMeans(b0), rep(1,L))
		Theta <- rbind(a0,b0)
		if (pen.method!="simple")
		{
			Theta <- Theta -
				tcrossprod(rowMeans(Theta), rep(1,L))
		}
		else
		{
			Theta <- Theta -
				tcrossprod(apply(Theta,1,FUN=median), rep(1,L))
			Theta[1,] <- Theta[1,] - mean(Theta[1,])
		}
	}
	else
	{
		std.X <- sqrt(colSums(Xc^2)/(n-L))
		Theta <- rbind(a0,b0)
	}
	
	# Augmented data matrix (feature vectors):
	XX <- cbind(rep(1,n),X)
	dd <- d+1
	tau <- c(0,tau.o*std.X)
	tmp <- switch(pen.method,
		vectors = LR1.full(XX,Y,Theta,tau),
		simple  = LR2.full(XX,Y,Theta,tau),
		none    = LR0.full(XX,Y,Theta))
	LR <- tmp$LR
	# regularize Hessian matrix:
	nu <- mean(diag(tmp$HLR))*10^(-5)
	delta <- qr.solve(tmp$HLR+diag(nu,L*dd,L*dd),tmp$GLR)
	dirderiv <- sum(tmp$GLR*delta)

	iter1 <- 0
	if (progress)
	{
		txt <- paste('Iteration ',
			as.character(iter1),
			': LR = ',
			as.character(round(LR,digits=5)))
		print(txt,quote=FALSE)
	}
	while (iter1 < 200 && dirderiv > 10^(-7))
	{
		iter1 <- iter1+1
		Theta.new <- Theta - matrix(delta,nrow=dd,ncol=L)
		
		LR.new <- switch(pen.method,
			vectors = LR1.only(XX,Y,Theta.new,tau),
			simple  = LR2.only(XX,Y,Theta.new,tau),
			none    = LR0.only(XX,Y,Theta.new))
		Theta.mid <- (Theta + Theta.new)/2
		LR.mid <- switch(pen.method,
			vectors = LR1.only(XX,Y,Theta.mid,tau),
			simple  = LR2.only(XX,Y,Theta.mid,tau),
			none    = LR0.only(XX,Y,Theta.mid))
		iter2 <- 0
		while ((LR.new > LR || LR.mid < LR.new) && iter2 < 20)
		{
			iter2 <- iter2+1
			Theta.new <- Theta.mid
			LR.new <- LR.mid
			Theta.mid <- (Theta + Theta.new)/2
			LR.mid <- switch(pen.method,
				vectors = LR1.only(XX,Y,Theta.mid,tau),
				simple  = LR2.only(XX,Y,Theta.mid,tau),
				none    = LR0.only(XX,Y,Theta.mid))
		}
		if (LR.new < LR)
		{
			Theta <- Theta.new
			tmp <- switch(pen.method,
				vectors = LR1.full(XX,Y,Theta,tau),
				simple  = LR2.full(XX,Y,Theta,tau),
				none    = LR0.full(XX,Y,Theta))
			LR <- tmp$LR
			# regularize Hessian matrix:
			nu <- mean(diag(tmp$HLR))*10^(-5)
			delta <- qr.solve(tmp$HLR+diag(nu,L*dd,L*dd),tmp$GLR)
			dirderiv <- sum(tmp$GLR*delta)
		}
		else
		{
			dirderiv <- 0
		}
		if (progress)
		{
			txt <- paste('Iteration ',
				as.character(iter1),
				': LR = ',
				as.character(round(LR,digits=5)))
			print(txt,quote=FALSE)
		}
	}
	return(list(a0=a0,b0=b0,std.X=std.X,
		a=Theta[1,],b=Theta[2:dd,],
		PM=tmp$PM))
}
