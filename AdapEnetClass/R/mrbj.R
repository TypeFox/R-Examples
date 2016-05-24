mrbj <-
function(formula, data, subset, trace=FALSE, gehanonly=FALSE, cov=FALSE, 
na.action=na.exclude, residue=FALSE, mcsize=100)
    {
	lss.betag<-function(x,y,delta,z)
	{
	row=ncol(x)
	col=ncol(z)
	betagm<-matrix(0,ncol=col,nrow=row)
	dimnum<-dim(x)
	n1<-dimnum[1]
	n2<-dimnum[2]	
	yy0<-rep(y,rep(n1,n1))
	delta1<-rep(delta,rep(n1,n1))
	yy1<-rep(y,n1)
	yy2<-delta1*(yy0-yy1)	
	xx0<-matrix(rep(as.vector(x),rep(n1,n1*n2)),nrow=n1*n1)
	xx1<-t(matrix(rep(as.vector(t(x)),n1),nrow=n2))
	xx2<-xx0-xx1
	
	for(i in 1:col)
	{
		zz=rep(z[,i],rep(n1,n1))*rep(z[,i],n1)
		xxdif<-xx2*zz*delta1
		xnew<-apply(xxdif,2,sum)
		xnew<-rbind(xxdif)
		yynew<-c(yy2*zz)
            fit <- Enet.wls(xnew, yynew, delta1)$beta
    		betagm[,i] <- fit
	}
	betagm
	}

	eps <- .Machine$double.eps^(2/3)
	call <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
	mf <- mf[c(1, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, sys.parent())
	Terms <- attr(mf, "terms")
	xvars <- as.character(attr(Terms, "variables"))
	yvar <- attr(Terms, "response")
	if((yvar <- attr(Terms, "response")) > 0)
		xvars <- xvars[ - yvar]
	else xlevels <- NULL
	y <- model.extract(mf, "response")
	x <- model.matrix(Terms, mf)

	if(all(x[, 1] == 1))
		x <- x[, -1]	
		
	if(ncol(as.matrix(y)) != 2)
		stop("Response must be a right-censored survival object!")

	nobs <- nrow(y)
	#nvar <- ncol(x)
      nvar1 <- ncol(x)

	fit <- list(converged = FALSE, gehanonly=gehanonly, cov=cov, mcsize=mcsize)

	fit$call <- call
	fit$nobs <- nobs
	fit$censored <- nobs - sum(y[,2])
	fit$niter <- 0
	fit$printkm <- residue
	if(gehanonly)
		{
	z <- matrix(rexp(nobs*mcsize), ncol=mcsize)
	zdummy <- matrix(rep(1,nobs), ncol=1)

	beta <- lss.betag(x, y[,1], y[,2], zdummy)
	betastar <- lss.betag(x, y[,1], y[,2], z)

	fit$betag <- beta
      beta <- lss.betag(x, y[,1], y[,2], zdummy)
	betastar <- lss.betag(x, y[,1], y[,2], z)
	fit$cnames <- dimnames(x)[[2]]
      nvar <- ncol(x)	
      bbar <- apply(betastar, 1, mean)
	tmp <- betastar - bbar
	fit$gehancov <- tmp %*% t(tmp)/(mcsize - 1)
	fit$gehansd <- sqrt(diag(fit$gehancov))
	fit$gehanzvalue <- beta/fit$gehansd
	fit$gehanpvalue <- (1 - pnorm(abs(fit$gehanzvalue))) * 2
	dimnames(fit$gehancov) <- list(fit$cnames,fit$cnames)
	if(trace)
		cat("\nbetag: ", format(beta), "\n\n")
       	}
	gnet<-Enet.wls(x, y[,1], y[,2])
	fit$enet <- gnet$beta
	fit$fit <- gnet$fit
	fit
     }
