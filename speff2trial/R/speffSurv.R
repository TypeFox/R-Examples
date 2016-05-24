speffSurv <- function(formula, data, force.in=NULL, nvmax=9, 
	method=c("exhaustive", "forward", "backward"), 
	optimal=c("cp", "bic", "rsq"), trt.id, conf.level=0.95, fixed=FALSE){
	if (missing(trt.id)) stop("Treatment indicator in 'trt.id' is missing.")
	method <- match.arg(method)
	optimal <- match.arg(optimal)
	mf <- match.call()
	mf$trt.id <- mf$method <- mf$conf.level <- mf$optimal <- mf$fixed <- NULL
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())

	X <- as.matrix(model.matrix(terms(formula), mf)[,-1])
	survtime <- model.response(mf)
	U <- survtime[,1]
	d <- survtime[,2]
	N <- sum(n <- table(ind <- as.factor(data[,trt.id])))
	levels(ind) <- 0:1
	ind <- as.numeric(as.vector(ind))
	pi <- n[2]/N

	Y <- function(u) ifelse(U>=u, 1, 0)
	Ymat <- function(u) drop(U>=matrix(u,N,length(u),byrow=TRUE))
	Yi <- function(idx,u) U[idx]>=u
	Ytx <- function(u,tx) (U>=u & ind==tx)
	aBar <- function(time,beta){
		drop(crossprod(ind*exp(beta*ind),Ymat(time)) / crossprod(exp(beta*ind),Ymat(time)))
	}
	Xi <- function(idx) matrix(X[idx,],N,NCOL(X),byrow=TRUE)
	wBar <- function(time,tx){
		t(t(crossprod(X,Ymat(time)*(ind==tx))) / pmax(drop(crossprod(Ymat(time),(ind==tx))),1))
	}
	Ktx <- function(time,tx){
		ntimes <- length(summary(km[tx+1])$time)
		idx <- apply(matrix(time,ntimes,length(time),byrow=TRUE)>=summary(km[tx+1])$time,2,sum)
		stime <- ifelse(idx==0,1,summary(km[tx+1])$surv[pmax(idx,1)])
		ifelse(stime>0,stime,summary(km[tx+1])$surv[pmax(idx-1,1)]/2)
	}
	K <- function(time,tx){
		K0 <- Ktx(time,0)
		K1 <- Ktx(time,1)
		ifelse(tx==0,K0,K1)
	}
	KY <- function(beta){pmax(drop(crossprod(exp(beta*ind),Ymat(U))),1)}
	mHat <- function(beta){ d*(ind - aBar(U,beta)) - sapply(1:N, function(i){ sum(d*(ind[i]- aBar(U,beta))*
	exp(beta*ind[i])*Yi(i,U)/KY(beta)) }) }
	estfnc <- function(beta){ sum(d*(ind-aBar(U,beta)) - (ind-pi)*f - g) }
	
	fitPH <- coxph(as.formula(paste(c(formula[[2]]),"~",trt.id)), data=data)
	betaPH <- fitPH$coef
	varbetaPH <- fitPH$var
	km <- survfit(as.formula(paste("Surv(U,1-d)~",trt.id)), data=data)
	KY1 <- sapply(1:N, function(k) max(sum(Ytx(U[k],1)),1))
	KY0 <- sapply(1:N, function(k) max(sum(Ytx(U[k],0)),1))
	m <- mHat(betaPH)                    
	KM1 <- Ktx(U,1)
	KM0 <- Ktx(U,0)
	K <- K(U,ind)
	wbar1 <- wBar(U,1)
	wbar0 <- wBar(U,0)
	H <- lapply(1:N, function(i){ (1-d[i])*(X[i,]-wBar(U[i],ind[i]))/Ktx(U[i],ind[i]) -
	ind[i]*drop((t(Xi(i))-wbar1) %*% ( (1-d)*ind*Yi(i,U)/(KM1*KY1) )) -
	(1-ind[i])*drop((t(Xi(i))-wbar0) %*% ( (1-d)*(1-ind)*Yi(i,U)/(KM0*KY0) )) })
	H <- t(do.call("cbind",H))
	colnames(H) <- colnames(X)
	if (fixed){
		W <- X
		Q <- H
	} else {
		respRnd <- (ind-pi)*m
		modRnd <- modSearch(as.formula(paste("respRnd~",paste(c(trt.id,colnames(X)),collapse="+"))),
		X, respRnd, "quantitative", method, optimal, force.in, nvmax)
		modCens <- modSearch(as.formula(paste("m~",c(formula[[3]]))), H, m, "quantitative", 
		method, optimal, force.in, nvmax)
		W <- X[,modRnd$names]
		Q <- H[,modCens$names]
	}

	f <- W %*% solve(pi*(1-pi)*crossprod(W)) %*% crossprod(W,(ind-pi)*m)
	g <- Q %*% solve(crossprod(Q)) %*% crossprod(Q,m)
	betaAUG <- uniroot(estfnc, c(-10,10))$root
	abarAUG <- aBar(U,betaAUG)
	varbetaAUG <- sum((mHat(betaAUG)-(ind-pi)*f-g)^2)/sum(d*(1-abarAUG)*abarAUG)^2
	
	fits <- list(beta=c(betaPH,betaAUG))
	class(fits) <- "speffSurv"
	fits$varbeta <- c(varbetaPH, varbetaAUG)
	names(fits$beta) <- names(fits$varbeta) <- c("Prop Haz", "Speff")
	if (!fixed) fits$formula <- list(rndSpace=formula(modRnd$mod), censSpace=formula(modCens$mod))
	fits$fixed <- fixed
	fits$conf.level <- conf.level
	fits$method <- method
	fits$n <- n
	fits
}

