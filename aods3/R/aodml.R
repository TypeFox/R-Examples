aodml <- function(formula, data, family = c("bb", "nb"),
                  link = c("logit", "cloglog", "probit"), 
                  phi.formula = ~ 1, 
                  phi.scale = c("identity", "log", "inverse"),
                  phi.start = NULL, fixpar = list(), hessian = TRUE,
                  method = c("BFGS", "Nelder-Mead", "CG", "SANN"),
                  control = list(maxit = 3000, trace = 0), ...) {
	
  call <- match.call(expand.dots = FALSE)
	fam <- match.arg(family)
	phi.scale <- match.arg(phi.scale)
  method <- match.arg(method)

	mu.f <- formula
	phi.f <- phi.formula
  if(phi.f != ~ 1)
    phi.f <- update.formula(phi.f, ~ . - 1)

	dat <- data
	
	mf <- model.frame(formula = mu.f, data = dat)
	resp <- model.response(mf)
	
  if(fam == "bb") {
  	link <- match.arg(link)
		fam.glm <- eval(parse(text = paste("binomial(link =", link,")")))
		m <- resp[, 1]
		n <- rowSums(resp)
		y <- m / n
	  }
  
	if(fam == "nb") {
		fam.glm <- "poisson"
		link <- "log"
		y <- as.vector(resp)
	  }
  
  offset <- model.offset(mf)
	
	# model matrices
  X.b <- model.matrix(object = mu.f, data = dat)
  X.phi <- model.matrix(object = phi.f, data = dat)
  nbb <- ncol(X.b)
	nbphi <- ncol(X.phi)
	# total number of parameters (including the eventual ones that are set as constant)
	nbpar <- nbb + nbphi
	
	# Parameters set as constant (fixpar)
	# Argument fixpar must be a list with two components:
	# the rank(s) and the value(s) of the parameters that are set as constant 
	id.par <- id.est <- 1:nbpar
	fixp <- FALSE
	id.fix <- theta.fix <- NULL
	if(!is.null(unlist(fixpar))) {
		fixp <- TRUE
		id.fix <- fixpar[[1]]
		theta.fix <- fixpar[[2]]
	  id.est <- id.par[-id.fix]
	  }

	# starting values of the parameters
  fm <- glm(formula = mu.f, family = fam.glm, data = dat)
	
  if(is.null(phi.start)) {
		z <- 0.1
		val <- switch(phi.scale, "identity" = z, "log" = log(z), "inverse" = 1 / z)
	  }
  else
		val <- phi.start
  
  phi.start <- rep(val, nbphi)
	param.start <- c(coef(fm), phi.start)

	# logL  
  dbetabin <- function(m, n, mu, k, log = FALSE) {
    a <- (k - 1) * mu
    b <- (k - 1) * (1 - mu)
    if(log)
      lbeta(a + m, b + n - m) - lbeta(a, b) + lchoose(n, m)
    }
  
  m.logL <- function(theta){
  	
  	param <- vector(length = nbpar)
  	param[id.fix] <- theta.fix
  	param[id.est] <- theta
  	
  	b <- param[1:nbb]
    nu <- as.vector(X.b %*% b)
    mu <- if(is.null(offset)) invlink(nu, type = link) else invlink(nu + offset, type = link)
    vecphi <- as.vector(X.phi %*% param[(nbb + 1):(nbb + nbphi)])
    
    k <- switch(phi.scale, identity = 1 / vecphi, log = 1 / exp(vecphi), inverse = vecphi)

    l <- vector(length = length(mu))
    
    if(fam == "bb")
      l <- dbetabin(m = m, n = n, mu = mu, k = k, log = TRUE)
    
    if(fam == "nb")
      l <- dnbinom(x = y, mu = mu, size = k, log = TRUE)
      
    -sum(l)
    
    }

	# fit
  res <- optim(par = param.start[id.est], fn = m.logL, hessian = hessian, 
               method = method, control = control, ...)
  
	## Results
  
	#param
	theta <- res$par
  param <- vector(length = nbpar)
  param[id.fix] <- theta.fix
  param[id.est] <- theta
  namb <- colnames(X.b)
  namphi <- paste("phi", colnames(X.phi), sep = ".")
  names(param) <- c(namb, namphi)
  
  b <- param[seq(along = namb)]
  phi <- param[-seq(along = namb)]

  # varparam
  is.singular <- function(X) qr(X)$rank < nrow(as.matrix(X))
	varparam <- matrix(NA, ncol = nbpar, nrow = nbpar)
  singular.H0 <- NA
	if(hessian) {
    H0 <- res$hessian
    singular.H0 <- is.singular(H0)
    if(!singular.H0)
			varparam[id.est, id.est] <- qr.solve(H0)
		else
			warning("The hessian matrix was singular.\n")
	  }
  dimnames(varparam)[[1]] <- dimnames(varparam)[[2]] <- names(param)
  
  # log-likelihood contributions
  nu <- as.vector(X.b %*% b)
  mu <- if(is.null(offset)) invlink(nu, type = link) else invlink(nu + offset, type = link)
  vecphi <- as.vector(X.phi %*% phi)
  k <- switch(phi.scale, identity = 1 / vecphi, log = 1 / exp(vecphi), inverse = vecphi)
  l <- lmax <- vector(length = length(mu))

  if(fam == "bb") {
    l <- dbetabin(m = m, n = n, mu = mu, k = k, log = TRUE)
    lmax <- dbetabin(m = m, n = n, mu = y, k = k, log = TRUE)
    }
  
  if(fam == "nb") {
    l <- dnbinom(x = y, mu = mu, size = k, log = TRUE)
    lmax <- dnbinom(x = y, mu = y, size = k, log = TRUE)
    }  
  
  l[is.na(l)] <- 0
  lmax[is.na(lmax)] <- 0  
	
  # other results
	# if fixpar is not null, df.model is lower than nbpar
	df.model <- length(theta)
  logL <- -res$value
  df.residual <- length(y) - df.model
  iterations <- res$counts
  code <- res$convergence
  msg <- if(!is.null(res$message)) res$message else character(0)
		
  structure(
  	list(
  		call = call, family = fam, link = link, dat = dat,
  		formula = mu.f, phi.formula = phi.f, phi.scale = phi.scale,
  		X.b = X.b, X.phi = X.phi,
  		resp = resp, offset = offset,
  		param = param, b = b, phi = phi,
    	varparam = varparam,
  		nbpar = nbpar, df.model = df.model, df.residual = df.residual,
  		logL = logL, l = l, lmax = lmax,
  		iterations = iterations, code = code, msg = msg,
      singular.hessian = singular.H0),
  	class = "aodml")
  }
