TML.noncensored<-function(formula, data, errors = "Gaussian", cu = NULL, initial = "S",
  otp = "fixed", cov = "parametric", input = NULL, control = list(), ...)
{
	if(!(errors %in% c("Gaussian", "logWeibull")))
		stop(gettextf("Errors distribution should be Gaussian or log-Weibull"))
	if(!(initial %in% c("S", "input")))
    stop(gettextf("initial should be S or input"))
  if(!(otp %in% c("adaptive", "fixed")))
    stop(gettextf("otp should be fixed or adaptive"))
	if(!(cov %in% c("no", "parametric", "nonparametric")))
    stop(gettextf("otp should be one of no, parametric or nonparametric"))
  call <- match.call()
	if(missing(data))
		data <- environment(formula)
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data"), names(mf), 0)
	mf <- mf[c(1, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	y <- model.response(mf, "any")
	if(length(dim(y)) == 1){
		nm <- rownames(y)
		if(!is.null(nm))
			names(y) <- nm
	}
	X <- if(!is.empty.model(mt))
		model.matrix(mt, mf, contrasts)
	else matrix(, NROW(y), 0)
	n <- length(y)
	np <- ncol(X)
	
	control <- do.call("TML.noncensored.control", control)
	iv <- control$iv
	nrep <- control$nrep
	gam <- control$gam
	nitmon <- control$nitmon
        maxit <- control$maxit
        tol  <- control$tol
	fastS <- control$fastS
	seed <- control$seed

	if(errors == "Gaussian"){
		if(is.null(cu))
			cu <- 2.5
		cl <- -cu
		ips <- 2
		xk <- 1.54764
		beta <- 0.5
		if (initial == "S") {
			if (fastS) {
				set.seed(seed)
				zctrl <- lmrob.control(...)
				z <- lmrob.S(X, y, zctrl)
				th0 <- z$coef
				v0 <- z$scale
                                nit0 <- z$k.iter
			}
			else {
				z <- dfcomn2(ipsi = 4, xk = xk, beta = beta)
				if (np <= 2 & n <= 500)
					iopt <- 3
				else iopt <- 1
				if (nrep != 0)
					iopt <- 2
				z <- hysest(X, y, np + 1, iopt = iopt, intch = 1,
				nrep = nrep, tols = tol, tolr = tol, iseed = seed)
				th0 <- z$theta[1:np]
				v0 <- z$smin
                                nit0 <- if (z$ierr<=1) 0 else nrep
			}
		}
		if (initial == "input") {
			z <- input
			th0 <- z$theta
			v0 <- z$sigma
		}
		namat <- matrix(NA, nrow = np, ncol = np)
		nares <- list(th0 = th0, v0 = v0, th1 = rep(NA, np), v1 = NA,
			tl = NA, tu = NA, V0 = NA, COV = namat,
			V1 = NA, alpha = NA, tn = NA, beta = NA, wi = rep(NA,n))
		re <- y - as.vector(X %*% as.matrix(th0))
		rs <- re/v0
		tp <- adaptn(sort(rs), cl, cu, option = otp)
		if (is.na(tp$tu))
			return(nares)
		wi <- tPsin(rs, tp$tl, tp$tu)
		yr <- y[wi != 0 | rs == 0]
		Xr <- X[wi != 0 | rs == 0, , drop = F]
		tp$tn <- length(yr)
		z <- MLnp(Xr, yr, iv, tp)
		fitted.values <- X %*% z$th1
		residuals <- y - X %*% z$th1
		dn <- colnames(X)
		if(is.null(dn))
			dn <- paste("x", 1:np, sep = "")
		names(th0) <- dn
		names(z$th1) <- dn
		res <- list(th0 = th0, v0 = v0, nit0=nit0, th1 = z$th1, v1 = z$v1, nit1=0, tl = tp$tl, 
                       tu = tp$tu, alpha = tp$alpha, tn = tp$tn, beta = tp$beta, weights = (wi != 0) * 1, 
      residuals = residuals, fitted.values = fitted.values, call = call, 
      formula = formula, terms = mt, data = data, errors = errors)
		if (cov != "no") {
			l <- cl
			u <- cu
			if (otp == "adaptive") {
				l <- tp$tl
				u <- tp$tu
			}
			if (cov == "nonparametric")
				K <- Cov2.n(X, y, u, z$th1, z$v1, opt = "averages",xk = xk)
			if (cov == "parametric")
				K <- CovE.n(X, y, u, z$th1, z$v1)
			res <- c(res, list(COV = K$CV1))
		}
	}
	if(errors == "logWeibull"){
		if(is.null(cu))
			cu <- 1.855356
		cl <- Izero(cu)
		ips <- 2
		xk <- 1.717817
		beta <- 0.5
    if (initial == "S"){
			if (fastS) {
				set.seed(seed)
				zctrl <- lmrob.control(...)
				z <- lmrob.S(X, y, zctrl)
				th0 <- z$coef
				v0 <- z$scale
                                nit0 <- z$k.iter
			}
			else {
				if (np <= 2 & n <= 500)
					iopt <- 3
				else iopt <- 1
				if (nrep != 0)
					iopt <- 2
				z <- dfcomn2(ipsi = 4, xk = xk, beta = beta)
				z <- hysest(X, y, nq = np + 1, iopt = iopt, intch = 1,
					nrep = nrep, tols = tol, tolr = tol, iseed = seed)
				th0 <- z$theta[1:np]
				v0 <- z$smin
				b0 <- -0.1352
				th0[1] <- th0[1] - b0 * v0
                                nit0 <- if (z$ierr<=1) 0 else nrep
			}
		}
		if (initial == "input") {
			z <- input
			v0 <- z$v
			th0 <- z$tau
		}
		namat <- matrix(NA, nrow = np, ncol = np)
		nares <- list(th0 = th0, v0 = v0, th1 = rep(NA, np), v1 = NA,
			tl = NA, tu = NA, V0 = NA, COV = namat,
			V1 = NA, alpha = NA, tn = NA, beta = NA, wi = rep(NA,n))
		re <- y - as.vector(X %*% as.matrix(th0))
		rs <- re/v0
		tp <- adaptw(sort(rs), cl, cu, otp)
		if (is.na(tp$tu))
			return(nares)
		wi <- tPsiw(rs, tp$tl, tp$tu)
		yr <- y[wi != 0 | rs == 0]
		Xr <- X[wi != 0 | rs == 0, , drop = FALSE]
		tp$tn <- length(yr)
		z <- MLwp(Xr, yr, th0, v0, iv, n, tp, gamm = gam, maxit,tol,
			nitmon)
		fitted.values <- X %*% z$th1
		residuals <- y - X %*% z$th1
		dn <- colnames(X)
		if(is.null(dn))
			dn <- paste("x", 1:np, sep = "")
		names(th0) <- dn
		names(z$th1) <- dn
		res <- list(th0 = th0, v0 = v0, nit0=nit0, th1 = z$th1, v1 = z$v1, nit1 = z$nit, tl = tp$tl, 
			tu = tp$tu, alpha = tp$alpha, tn = tp$tn, beta = tp$beta, weights = (wi != 0) * 1, 
      residuals = residuals, fitted.values = fitted.values, call = call, 
      formula = formula, terms = mt, data = data, errors = errors)
		if (cov != "no") {
			l <- cl
			u <- cu
			if (otp == "adaptive") {
				l <- tp$tl
				u <- tp$tu
			}
			if (cov == "nonparametric")
				K <- Cov2.w(X, y, l, u, z$th1, z$v1, opt = "averages")
			if (cov == "parametric")
				K <- CovE.w(X, y, l, u, z$th1, z$v1)
			res <- c(res, list(COV = K$CV1))
		}
	}  
	class(res) <- "TML"
	res
}

