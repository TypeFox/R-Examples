#' @importFrom stats sd quantile prcomp splinefun model.matrix model.frame model.response qlogis
#' @importFrom stats model.weights .getXlevels runif logLik nobs lm terms delete.response
#' @importFrom splines spline.des
#' @importFrom survival Surv is.Surv

BS <- function (x, df, knots, degree, Boundary.knots = range(x), deriv = 0){
  
  n <- length(x)
  xl <- (x < (Bl <- Boundary.knots[1]))
  xr <- (x > (Br <- Boundary.knots[2]))
  x.l <- x[xl]; x.r <- x[xr]
  x <- x[(xin <- !(xout <- (xl | xr)))]
  flag <- FALSE
  
  ord <- 1 + (degree <- as.integer(degree))
  if(ord <= 1){stop("'degree' must be integer >= 1")}
  if(missing(knots)){
    nIknots <- df - ord + 1
    if(nIknots < 0){
      nIknots <- 0
      stop("'df' cannot be smaller than 'degree'")
    }
    knots <- if(nIknots > 0){
      knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
                         2)[-c(1, nIknots + 2)]
      stats::quantile(x, knots)
    }
  }
  Aknots <- sort(c(rep(Boundary.knots, ord), knots))
  if(length(x) == 0){x <- mean(Boundary.knots); flag <- TRUE}
  basis <- splines::spline.des(Aknots, x, ord, derivs = rep(deriv, length(x)))$design
  df <- ncol(basis <- basis[, -1L, drop = FALSE])
  if(flag){basis <- NULL}
  
  if(any(xout)){
    k <- splines::spline.des(Aknots, Bl, ord, derivs = 1)$design
    k <- k[, -1L, drop = FALSE][1,1]
    temp <- matrix(0, n, df)
    temp[xin,] <- basis; basis <- temp
    if(any(xl)){basis[xl,1] <- (if(deriv == 0) k*(x.l - Bl) else k)}
    if(any(xr)){
      basis[xr,df] <- (if(deriv == 0) 1 + k*(x.r - Br) else k)
      if(df > 1){basis[xr,df - 1] <- (if(deriv == 0) -k*(x.r - Br) else -k)}
    }
  }
  
  a <- list(degree = degree, knots = if (is.null(knots)) numeric(0L) else knots, 
            Boundary.knots = Boundary.knots)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("BS", "basis", "matrix")
  basis
}


#' @export
flexPM <- function(formula, data, weights, df = 3, degree = 3, knots, maxit, tol = 1e-6){

	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "weights", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	
	if(!survival::is.Surv(zyd <- model.response(mf))){stop("the model response must be created with Surv()")}
	if((n <- nrow(zyd)) == 0){stop("zero non-NA cases")}
	type <- attributes(zyd)$type
	if(type == "right"){
		y <- zyd[,1]; z <- rep.int(-Inf,n); d <- zyd[,2]
		type <- (if(any(d == 0)) "cens" else "u")
	}
	else if(type == "counting"){
		z <- zyd[,1]; y <- zyd[,2]; d <- zyd[,3]; type <- "trunc"
		if(all(z == -Inf)){type <- (if(any(d == 0)) "cens" else "u")}
	}
	else{stop("only type = 'right' and type = 'counting' are supported")}
	if(!(any(d == 1))){stop("all observation are censored")}
	
	x <- model.matrix(mt,mf)
	weights <- model.weights(mf)
	if(!is.null(weights)){weights <- weights/sum(weights)*length(y)}
	if(!missing(knots) && !is.null(knots)){
		if(any(knots <= min(y) | knots >= max(y))){stop("'knots' must be in the range of the response variable")}
	}
	
	model <- suppressWarnings(L.estimator(z,y,d,x,weights, df, degree, knots, 
	                     maxit = maxit, tol = tol, type = type))
	fit <- list(converged = model$converged, n.it = model$n.it, 
	            n = model$n, n.free.par = model$n.free.par, 
	            logLik = model$loglik, AIC = model$AIC, BIC = model$BIC
	)
	fit$mf <- mf
	attr(fit$mf, "m") <- model$m
	attr(fit$mf, "s") <- model$s
	attr(fit$mf, "u0") <- model$u0
	attr(fit$mf, "u0inv") <- model$u0inv
	attr(fit$mf, "u") <- model$u
	attr(fit$mf, "uinv") <- model$uinv
	attr(fit$mf, "xystats") <- model$xystats
	attr(fit$logLik, "df") <- model$n.free.par
	fit$call <- cl
	class(fit) <- "flexPM"
	if(!fit$converged){warning("optimization algorithm did not converge")}
  
	fit
}



L.estimator <- function(z,y,d,x,weights, df,degree,knots, theta, maxit = 100, tol = 1e-6, type){

	if(type == "trunc"){loglik <- (if(is.null(weights)) L.loglik.cens.trunc else L.wloglik.cens.trunc)}
	else if(type == "cens"){loglik <- (if(is.null(weights)) L.loglik.cens else L.wloglik.cens)}
	else{loglik <- (if(is.null(weights)) L.loglik else L.wloglik)}

	# x

	v <- qr((x <- cbind(1,x)))
	sel <- v$pivot[1:v$rank]
	x <- x[,sel, drop = FALSE]
	if((q <- ncol(x)) > 2){
		pcx <- stats::prcomp(x[,-1], center = TRUE, scale = TRUE)
		x[,-1] <- pcx$x
	}
	SdX <- apply(x,2,sd); Mx <- colMeans(x)
	Mx[1] <- 0; SdX[1] <- 1
	x <- scale(x, center = Mx, scale = SdX)
	

	# y, z # pretransformation

	n <- length(y)
	ry <- range(y)
	z[notrunc <- which(z < ry[1])] <- -Inf
	P <- c(0.1,0.25,0.5,0.75,0.9)
	qy <- c(2*ry[1] - ry[2], ry[1], quantile(y,P), ry[2], 2*ry[2] - ry[1])
	qq <- qlogis(c(1/10/n, 0.025, P, 0.975, 1 - 1/10/n))
 	u0 <- stats::splinefun(qy,qq, method = "hyman")
 	u0inv <- stats::splinefun(qq,qy, method = "hyman")
 	y <- u0(y)
	if(type == "trunc"){z <- u0(z); z[notrunc] <- -Inf}
	if(!missing(knots) && !is.null(knots)){knots <- u0(knots)}
  
	# y, z # scaling
  
	ry <- range(y)
	min.y <- ry[1]
	max.y <- ry[2]
	y <- (y - min.y)/(max.y - min.y)*10
	z <- (z - min.y)/(max.y - min.y)*10
	if(!missing(knots) && !is.null(knots)){knots <- (sort(knots) - min.y)/(max.y - min.y)*10}
	d1 <- which(d == 1)
	dd <- 1 + d

	# u(.)

	Bk <- c(-1e-5, 10 + 1e-5)
	By <- BS(y, df = df, knots = knots, degree = degree, Boundary.knots = Bk); knots <- attr(By, "knots")
	B1y <- BS(y, df = df, degree = degree, Boundary.knots = Bk, deriv = 1, knots = knots)
	Bz <- (if(type == "trunc") BS(z, df = df, degree = degree, Boundary.knots = Bk, knots = knots) else NULL)

	if((df <- ncol(By)) > 1){
		for(j in (df - 1):1){
			By[,j] <- By[,j] + By[,j + 1]
			B1y[,j] <- B1y[,j] + B1y[,j + 1]
			Bz[,j] <- Bz[,j] + Bz[,j + 1]
		}
	}

	# Starting points ##############################################################################

	if(missing(theta)){
		m0 <- lm(rowSums(By) ~ -1 + x, weights = weights)
		mu <- m0$coef
		sigma <- sd(m0$residuals)*sqrt(3)/pi
		theta <- c(mu, log(sigma), rep(0, q - 1), rep(0, df - 1))
	}

	### global optimization ######################################################################

	if(missing(maxit)){maxit <- 10*length(theta)}
	opt <- suppressWarnings(
		flexPM.newton(c(theta), loglik, maxit = maxit, tol = tol,
		Bz = Bz, By = By, B1y = B1y, x = x, d = d, d1 = d1, dd = dd, 
		weights = weights, notrunc = notrunc)
	)
	phi <- opt$estimate
	theta <- (if(df > 1) c(exp(phi[(2*q + 1):length(phi)]), 1) else 1)
	mu <- phi[1:q]; sigma <- phi[(q + 1):(2*q)]

	# xy stats

	xystats <- list(type = type,
		mu = mu, sigma = sigma, theta = (if(df > 1) phi[(2*q + 1):length(phi)] else NULL),
		sel = sel, rot = (if(q > 2) pcx$rotation else NULL), 
		center0 = (if(q > 2) pcx$center else NULL), scale0 = (if(q > 2) pcx$scale else NULL),
		center = attr(x, "scaled:center"), scale = attr(x, "scaled:scale"),
		knots = knots, df = df, degree = attr(By, "degree"), 
		min.y = min.y, max.y = max.y
	)

	### output

	ll <- -opt$minimum
	r <- length(phi)

	g <- sort(c(-(1:10), seq(0,10, 0.05), (11:20)))
	Bg <- BS(g, df = df, knots = knots, degree = degree, Boundary.knots = Bk)
	if(df > 1){for(j in (df - 1):1){Bg[,j] <- Bg[,j] + Bg[,j + 1]}}
	u <- stats::splinefun(g, c(Bg%*%cbind(theta)), method = "hyman")
	uinv <- stats::splinefun(u(g),g, method = "hyman")
  
	out <- list(
		converged = (opt$n.it < maxit), n.it = opt$n.it, 
		m = c(x%*%cbind(mu)), s = c(exp(x%*%cbind(sigma))) + 0.001,
		u0 = u0, u0inv = u0inv, u = u, uinv = uinv, 
		loglik = ll, n = n, n.free.par = r,
		AIC = 2*(r - ll), BIC = -2*ll + r*log(n),
		xystats = xystats
	)
	out
}

# predict CDF, PDF, and survival (type = "CDF"), p-th quantile(s) (type = "QF"), 
# or simulate data with p ~ U(0,1) (type = "sim").
# 'tol' is used to compute quantiles.

#' @export
predict.flexPM <- function(object, type = c("CDF", "QF", "sim"), newdata, p, ...){
  
  if(is.na(match(type <- type[1], c("CDF", "QF", "sim")))){stop("invalid 'type'")}
  
  mf <- object$mf
  mu <- attr(mf, "m")
  sigma <- attr(mf, "s")
  u0 <- attr(mf, "u0")
  u0inv <- attr(mf, "u0inv")
  u <- attr(mf, "u")
  uinv <- attr(mf, "uinv")
  s <- attr(mf, "xystats")
  
  mt <- terms(mf)
  miss <- attr(mf, "na.action")
  nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:(nrow(mf) + length(miss)))[-miss])
  xlev <- .getXlevels(mt, mf)
  
  if(!missing(newdata)){
    
    if(type == "CDF"){
      yn <- as.character(if(s$type == "trunc") mt[[2]][[3]] else mt[[2]][[2]])
      if(is.na(ind <- match(yn, colnames(newdata))))
      {stop("for 'type = CDF', 'newdata' must contain the y-variable")}
      if(s$type == "cens"){newdata[,as.character(mt[[2]][[3]])] <- 1}
      if(s$type == "trunc"){newdata[,as.character(mt[[2]][[4]])] <- 1; newdata[,as.character(mt[[2]][[2]])] <- -Inf}
    }
    else{mt <- delete.response(mt)}
    if(any(is.na(match(all.vars(mt), colnames(newdata)))))
    {stop("'newdata' must contain all x-variables")}
    mf <- model.frame(mt, data = newdata, xlev = xlev)
    
    if(nrow(mf) == 0){
      nr <- nrow(newdata)
      if(type == "CDF"){
        out <- data.frame(matrix(NA,nr,3))
        colnames(out) <- c("log.f", "log.F", "log.S")
        rownames(out) <- rownames(newdata)
      }
      else if(type == "QF"){
        out <- data.frame(matrix(NA,nr,length(p)))
        colnames(out) <- paste("p",p, sep = "")
        rownames(out) <- rownames(newdata)
      }
      else{out <- rep.int(NA, nr)}
      return(out)
    }
    miss <- attr(mf, "na.action")
    nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:nrow(newdata))[-miss])
    
    x <- model.matrix(mt, mf)
    x <- cbind(1,x)
    x <- x[,s$sel, drop = FALSE]
    if((q <- ncol(x)) > 2)
    {x[,-1] <- scale(x[,-1, drop = FALSE], scale = s$scale0, center = s$center0)%*%s$rot}
    x <- scale(x, center = s$center, scale = s$scale)
    mu <- c(x%*%cbind(s$mu))
    sigma <- c(exp(x%*%cbind(s$sigma))) + 0.001
  }

  n <- length(miss) + length(nomiss)
  if(type == "CDF"){
    
    y0 <- cbind(model.response(mf))[,1 + (s$type == "trunc")]
    y <- u0(y0)
    y <- (y - s$min.y)/(s$max.y - s$min.y)*10 
    z1 <- (u(y) - mu)/sigma
    z2 <- log(1 + exp(-z1))
    z3 <- log(1 + exp(z1))
    if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
    if(any(out3 <- (z3 == Inf))){z3[out3] <- z1[out3]}
    log.F <- log.S <- log.f <- NULL
    log.F[nomiss] <- -z2
    log.S[nomiss] <- -z3
    log.f[nomiss] <- (-z1 - log(sigma) - 2*z2 + 
		  log(pmax(1e-6, u(y, deriv = 1))) + log(pmax(1e-6, u0(y0, deriv = 1))) +  
		  log(10) - log(s$max.y - s$min.y))
    log.F[miss] <- log.S[miss] <- log.f[miss] <- NA
    out <- data.frame(log.f = log.f, log.F = log.F, log.S = log.S)
    
    rownames(out) <- 1:n
    return(out)
  }
  else{
    if(type == "QF"){
      if(missing(p)){stop("please indicate 'p'")}
      if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}
      P <- log(p) - log(1 - p)
      out <- NULL
      for(j in 1:length(p)){
        Qp <- NULL
        Qp[nomiss] <- uinv(mu + sigma*P[j])
        Qp[miss] <- NA  
        out <- cbind(out, u0inv(Qp*(s$max.y - s$min.y)/10 + s$min.y))
      }
      colnames(out) <- paste("p", p, sep = "")
      out <- as.data.frame(out)
      rownames(out) <- 1:n
      return(out)
    }
    else{
      p <- runif(length(mu))
      P <- log(p) - log(1 - p)
      Qp <- NULL
      Qp[nomiss] <- uinv(mu + sigma*P)
      Qp[miss] <- NA   
      Qp <- u0inv(Qp*(s$max.y - s$min.y)/10 + s$min.y)
      return(Qp)
    }
  }
}

# print and summary method
#' @export
print.flexPM <- function(x, ...){
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Use predict() to obtain predictions from the fitted model")
  cat("\n")
}
#' @export
summary.flexPM <- function(object, ...){
  out <- list(call = object$call, n = object$n, n.free.par = object$n.free.par, converged = object$converged,
       n.it = object$n.it, logLik = object$logLik, AIC = object$AIC, BIC = object$BIC)
  class(out) <- "summary.flexPM"
  out
}

#' @export
print.summary.flexPM <- function(x, ...){
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("n. of obs: ", paste(deparse(round(x$n)), sep = " ", collapse = " "), "\n", sep = "")
  cat("n. of free parameters: ", paste(deparse(round(x$n.free.par)), sep = " ", collapse = " "), "\n", sep = "")
  cat("Converged: ", paste(deparse(x$converged), sep = " ", collapse = " "), "\n", sep = "")
  cat("N. of iterations: ", paste(deparse(round(x$n.it)), sep = " ", collapse = " "), "\n", sep = "")
  cat("Log-likelihood: ", paste(deparse(round(as.numeric(x$logLik),1)), sep = " ", collapse = " "), "\n", sep = "")
  cat("AIC: ", paste(deparse(round(x$AIC,1)), sep = " ", collapse = " "), "\n", sep = "")
  cat("BIC: ", paste(deparse(round(x$BIC,1)), sep = " ", collapse = " "), "\n", sep = "")
  cat("Use predict() to obtain predictions from the fitted model")
  cat("\n")
}
#' @export
logLik.flexPM <- function(object, ...){object$logLik}
#' @export
nobs.flexPM <- function(object, ...){object$n}




# Some u(y) is assumed to be logistic
# This is the log-likelihood function, weigthed/unweighted and truncated/non truncated.
# If deriv = 0, only returns the log-likelihood
# If deriv = 1, returns the N*q matrix of first derivatives
# If deriv = 2, returns the gradient and the hessian


L.wloglik <- function(theta, Bz, By, B1y, x,d,d1,dd,weights, notrunc, deriv = 2){

	ColSums <- (if(deriv == 1) I else colSums)
	h <- 1:((df <- ncol(By)) - 1)
	q <- ncol(x)
	mu <- theta[1:q]
	sigma <- theta[(q + 1):(2*q)]
	theta <- (if(df > 1) c(exp(theta[(2*q + 1):length(theta)]), 1) else 1)
	if(any(abs(theta) == Inf)){return(Inf)}

	mu <- c(x%*%cbind(mu))
	sigma <- exp(c(x%*%cbind(sigma))) + 0.001
	if(any(sigma == Inf)){return(Inf)}
	log.sigma <- log(sigma)
	uy <- c(By%*%cbind(theta))
	u1y <- c(B1y%*%cbind(theta))

	# log-likelihood #

	z1 <- (uy - mu)/sigma
	z2 <- log(1 + exp(-z1))
	z3 <- log(1 + exp(z1))
	if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
	if(any(out3 <- (z3 == Inf))){z3[out3] <- z1[out3]}
	log.F <- -z2
	log.S <- -z3
	log.f <- -z1 - log.sigma - 2*z2 + log(u1y)
		
	ll <- log.f
	loglik <- sum(ll*weights)

	if(deriv == 0){return(-loglik)}

	# Gradient

	a <- exp(log.F)*2 - 1
	v <- (sigma - 0.001)/sigma

	d.mu <- ColSums(x*(a*weights/sigma))
	d.sigma <- ColSums(x*(v*(a*z1 - 1)*weights)) # 
	if(df > 1){
		Theta <- t(matrix(theta[h], length(h), length(d)))
		By.h <- By[,h,drop = FALSE]*Theta
		B1y.h <- B1y[,h,drop = FALSE]*Theta
		d.theta <- B1y.h*(1/u1y) - By.h*(a/sigma) #
		d.theta <- ColSums(d.theta*weights)
	}
	else{d.theta <- NULL}
	if(deriv == 1){return(cbind(d.mu,d.sigma,d.theta))}
	gradient <- c(d.mu,d.sigma,d.theta)

	# Hessian
		
	dF_dg <- exp(log.S + log.F - log.sigma)
	dF_dmu <- -dF_dg
	dF_dsigma <- dF_dmu*z1
		
	H_mu <- t(x)%*%(x*(dF_dmu*2*weights/sigma))
	H_sigma <- t(x)%*%(x*(weights*((v^2)*z1*(2*sigma*dF_dsigma - a) + v*(0.001/sigma*(a*z1 - 1)))))
	H_mu_sigma <- t(x)%*%(x*(v*(dF_dsigma*2 - a/sigma)*weights))

	if(df > 1){
		H_theta <- (diag(d.theta, df - 1,df - 1) 
			- t(B1y.h)%*%(B1y.h*(weights/u1y^2)) - t(By.h)%*%(By.h*(weights*2/sigma*dF_dg)))
		H_mu_theta <- t(x)%*%(By.h*(dF_dg*2*weights/sigma))
		H_sigma_theta <- t(x)%*%(By.h*(weights*v*(a/sigma + z1*2*dF_dg)))

		H <- cbind(H_mu, H_mu_sigma, H_mu_theta)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma, H_sigma_theta))
		H <- rbind(H, cbind(t(H_mu_theta), t(H_sigma_theta), H_theta))
	}
	else{
		H <- cbind(H_mu, H_mu_sigma)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma))	
	}


	# finish #

	attributes(loglik) <- list(gradient = -gradient, hessian = -H)
	-loglik
}

L.loglik <- function(theta, Bz, By, B1y, x,d,d1,dd,weights, notrunc, deriv = 2){

	ColSums <- (if(deriv == 1) I else colSums)
	h <- 1:((df <- ncol(By)) - 1)
	q <- ncol(x)
	mu <- theta[1:q]
	sigma <- theta[(q + 1):(2*q)]
	theta <- (if(df > 1) c(exp(theta[(2*q + 1):length(theta)]), 1) else 1)
	if(any(abs(theta) == Inf)){return(Inf)}

	mu <- c(x%*%cbind(mu))
	sigma <- exp(c(x%*%cbind(sigma))) + 0.001
	if(any(sigma == Inf)){return(Inf)}
	log.sigma <- log(sigma)
	uy <- c(By%*%cbind(theta))
	u1y <- c(B1y%*%cbind(theta))


	# log-likelihood #

	z1 <- (uy - mu)/sigma
	z2 <- log(1 + exp(-z1))
	z3 <- log(1 + exp(z1))
	if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
	if(any(out3 <- (z3 == Inf))){z3[out3] <- z1[out3]}
	log.F <- -z2
	log.S <- -z3
	log.f <- -z1 - log.sigma - 2*z2 + log(u1y)
		
	ll <- log.f
	loglik <- sum(ll)

	if(deriv == 0){return(-loglik)}

	# Gradient

	a <- exp(log.F)*2 - 1
	v <- (sigma - 0.001)/sigma

	d.mu <- ColSums(x*(a/sigma))
	d.sigma <- ColSums(x*(v*(a*z1 - 1)))
	if(df > 1){
		Theta <- t(matrix(theta[h], length(h), length(d)))
		By.h <- By[,h,drop = FALSE]*Theta
		B1y.h <- B1y[,h,drop = FALSE]*Theta
		d.theta <- B1y.h*(d/u1y) - By.h*(a/sigma)
		d.theta <- ColSums(d.theta)
	}
	else{d.theta <- NULL}
	if(deriv == 1){return(cbind(d.mu,d.sigma,d.theta))}
	gradient <- c(d.mu,d.sigma,d.theta)

	# Hessian
		
	dF_dg <- exp(log.S + log.F - log.sigma)
	dF_dmu <- -dF_dg
	dF_dsigma <- dF_dmu*z1
		
	H_mu <- t(x)%*%(x*(dF_dmu*2/sigma))
	H_sigma <- t(x)%*%(x*(((v^2)*z1*(2*sigma*dF_dsigma - a) + v*(0.001/sigma*(a*z1 - 1)))))
	H_mu_sigma <- t(x)%*%(x*(v*(dF_dsigma*2 - a/sigma)))

	if(df > 1){
		H_theta <- (diag(d.theta, df - 1, df - 1) 
			- t(B1y.h)%*%(B1y.h*(1/u1y^2)) - t(By.h)%*%(By.h*(2/sigma*dF_dg)))
		H_mu_theta <- t(x)%*%(By.h*(dF_dg*2/sigma))
		H_sigma_theta <- t(x)%*%(By.h*(v*(a/sigma + z1*2*dF_dg)))

		H <- cbind(H_mu, H_mu_sigma, H_mu_theta)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma, H_sigma_theta))
		H <- rbind(H, cbind(t(H_mu_theta), t(H_sigma_theta), H_theta))
	}
	else{
		H <- cbind(H_mu, H_mu_sigma)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma))	
	}


	# finish #

	attributes(loglik) <- list(gradient = -gradient, hessian = -H)
	-loglik
}


L.wloglik.cens <- function(theta, Bz, By, B1y, x,d,d1,dd,weights, notrunc, deriv = 2){

	ColSums <- (if(deriv == 1) I else colSums)
	h <- 1:((df <- ncol(By)) - 1)
	q <- ncol(x)
	mu <- theta[1:q]
	sigma <- theta[(q + 1):(2*q)]
	theta <- (if(df > 1) c(exp(theta[(2*q + 1):length(theta)]), 1) else 1)
	if(any(abs(theta) == Inf)){return(Inf)}

	mu <- c(x%*%cbind(mu))
	sigma <- exp(c(x%*%cbind(sigma))) + 0.001
	if(any(sigma == Inf)){return(Inf)}
	log.sigma <- log(sigma)
	uy <- c(By%*%cbind(theta))
	u1y <- c(B1y%*%cbind(theta))

	# log-likelihood #

	z1 <- (uy - mu)/sigma
	z2 <- log(1 + exp(-z1))
	z3 <- log(1 + exp(z1))
	if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
	if(any(out3 <- (z3 == Inf))){z3[out3] <- z1[out3]}
	log.F <- -z2
	log.S <- -z3
	log.f <- -z1 - log.sigma - 2*z2 + log(u1y)
		
	ll <- log.f
	ll[-d1] <- log.S[-d1]
	loglik <- sum(ll*weights)

	if(deriv == 0){return(-loglik)}

	# Gradient

	a <- exp(log.F)*dd - d
	v <- (sigma - 0.001)/sigma

	d.mu <- ColSums(x*(a*weights/sigma))
	d.sigma <- ColSums(x*(v*(a*z1 - d)*weights))
	if(df > 1){
		Theta <- t(matrix(theta[h], length(h), length(d)))
		By.h <- By[,h,drop = FALSE]*Theta
		B1y.h <- B1y[,h,drop = FALSE]*Theta
		d.theta <- B1y.h*(d/u1y) - By.h*(a/sigma)
		d.theta <- ColSums(d.theta*weights)
	}
	else{d.theta <- NULL}
	if(deriv == 1){return(cbind(d.mu,d.sigma,d.theta))}
	gradient <- c(d.mu,d.sigma,d.theta)

	# Hessian
		
	dF_dg <- exp(log.S + log.F - log.sigma)
	dF_dmu <- -dF_dg
	dF_dsigma <- dF_dmu*z1
		
	H_mu <- t(x)%*%(x*(dF_dmu*dd*weights/sigma))
	H_sigma <- t(x)%*%(x*(weights*((v^2)*z1*(dd*sigma*dF_dsigma - a) + v*(0.001/sigma*(a*z1 - d)))))
	H_mu_sigma <- t(x)%*%(x*(v*(dF_dsigma*dd - a/sigma)*weights))

	if(df > 1){
		H_theta <- (diag(d.theta, df - 1, df - 1) 
			- t(B1y.h)%*%(B1y.h*(weights*d/u1y^2)) - t(By.h)%*%(By.h*(weights*dd/sigma*dF_dg)))
		H_mu_theta <- t(x)%*%(By.h*(dF_dg*dd*weights/sigma))
		H_sigma_theta <- t(x)%*%(By.h*(weights*v*(a/sigma + z1*dd*dF_dg)))

		H <- cbind(H_mu, H_mu_sigma, H_mu_theta)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma, H_sigma_theta))
		H <- rbind(H, cbind(t(H_mu_theta), t(H_sigma_theta), H_theta))
	}
	else{
		H <- cbind(H_mu, H_mu_sigma)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma))	
	}


	# finish #

	attributes(loglik) <- list(gradient = -gradient, hessian = -H)
	-loglik
}

L.loglik.cens <- function(theta, Bz, By, B1y, x,d,d1,dd,weights, notrunc, deriv = 2){

	ColSums <- (if(deriv == 1) I else colSums)
	h <- 1:((df <- ncol(By)) - 1)
	q <- ncol(x)
	mu <- theta[1:q]
	sigma <- theta[(q + 1):(2*q)]
	theta <- (if(df > 1) c(exp(theta[(2*q + 1):length(theta)]), 1) else 1)
	if(any(abs(theta) == Inf)){return(Inf)}

	mu <- c(x%*%cbind(mu))
	sigma <- exp(c(x%*%cbind(sigma))) + 0.001
	if(any(sigma == Inf)){return(Inf)}
	log.sigma <- log(sigma)
	uy <- c(By%*%cbind(theta))
	u1y <- c(B1y%*%cbind(theta))


	# log-likelihood #

	z1 <- (uy - mu)/sigma
	z2 <- log(1 + exp(-z1))
	z3 <- log(1 + exp(z1))
	if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
	if(any(out3 <- (z3 == Inf))){z3[out3] <- z1[out3]}
	log.F <- -z2
	log.S <- -z3
	log.f <- -z1 - log.sigma - 2*z2 + log(u1y)
		
	ll <- log.f
	ll[-d1] <- log.S[-d1]
	loglik <- sum(ll)

	if(deriv == 0){return(-loglik)}

	# Gradient

	a <- exp(log.F)*dd - d
	v <- (sigma - 0.001)/sigma

	d.mu <- ColSums(x*(a/sigma))
	d.sigma <- ColSums(x*(v*(a*z1 - d)))
	if(ncol(By) > 1){
		Theta <- t(matrix(theta[h], length(h), length(d)))
		By.h <- By[,h,drop = FALSE]*Theta
		B1y.h <- B1y[,h,drop = FALSE]*Theta
		d.theta <- B1y.h*(d/u1y) - By.h*(a/sigma)
		d.theta <- ColSums(d.theta)
	}
	else{d.theta <- NULL}
	if(deriv == 1){return(cbind(d.mu,d.sigma,d.theta))}
	gradient <- c(d.mu,d.sigma,d.theta)

	# Hessian
		
	dF_dg <- exp(log.S + log.F - log.sigma)
	dF_dmu <- -dF_dg
	dF_dsigma <- dF_dmu*z1
		
	H_mu <- t(x)%*%(x*(dF_dmu*dd/sigma))
	H_sigma <- t(x)%*%(x*(((v^2)*z1*(dd*sigma*dF_dsigma - a) + v*(0.001/sigma*(a*z1 - d)))))
	H_mu_sigma <- t(x)%*%(x*(v*(dF_dsigma*dd - a/sigma)))

	if(df > 1){
		H_theta <- (diag(d.theta, df - 1, df - 1) 
			- t(B1y.h)%*%(B1y.h*(d/u1y^2)) - t(By.h)%*%(By.h*(dd/sigma*dF_dg)))
		H_mu_theta <- t(x)%*%(By.h*(dF_dg*dd/sigma))
		H_sigma_theta <- t(x)%*%(By.h*(v*(a/sigma + z1*dd*dF_dg)))

		H <- cbind(H_mu, H_mu_sigma, H_mu_theta)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma, H_sigma_theta))
		H <- rbind(H, cbind(t(H_mu_theta), t(H_sigma_theta), H_theta))
	}
	else{
		H <- cbind(H_mu, H_mu_sigma)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma))	
	}


	# finish #

	attributes(loglik) <- list(gradient = -gradient, hessian = -H)
	-loglik
}

L.wloglik.cens.trunc <- function(theta, Bz, By, B1y, x,d,d1,dd,weights, notrunc, deriv = 2){

	ColSums <- (if(deriv == 1) I else colSums)

	h <- 1:((df <- ncol(By)) - 1)
	q <- ncol(x)
	mu <- theta[1:q]
	sigma <- theta[(q + 1):(2*q)]
	theta <- (if(df > 1) c(exp(theta[(2*q + 1):length(theta)]), 1) else 1)
	if(any(abs(theta) == Inf)){return(Inf)}

	mu <- c(x%*%cbind(mu))
	sigma <- exp(c(x%*%cbind(sigma))) + 0.001
	if(any(sigma == Inf)){return(Inf)}
	log.sigma <- log(sigma)
	uy <- c(By%*%cbind(theta))
	u1y <- c(B1y%*%cbind(theta))
	uz <- c(Bz%*%cbind(theta))


	# log-likelihood #

	z1 <- (uy - mu)/sigma
	z2 <- log(1 + exp(-z1))
	z3 <- log(1 + exp(z1))	
	if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
	if(any(out3 <- (z3 == Inf))){z3[out3] <- z1[out3]}

	log.F <- -z2
	log.S <- -z3
	log.f <- -z1 - log.sigma - 2*z2 + log(u1y)

	# truncation
	
	zz <- (uz - mu)/sigma
	log.Fz <- -log(1 + exp(-zz))
	log.Sz <- -log(1 + exp(zz))
	if(any(outFz <- (log.Fz == -Inf))){log.Fz[outFz] <- zz[outFz]}
	if(any(outSz <- (log.Sz == -Inf))){log.Sz[outSz] <- -zz[outSz]}
	
	ll <- log.f
	ll[-d1] <- log.S[-d1]
	loglik <- sum(ll*weights) - sum(log.Sz*weights)

	###
	
	if(deriv == 0){return(-loglik)}
	Bz[notrunc,] <- zz[notrunc] <- 1 # an arbitrary finite value

	# Gradient

	a <- exp(log.F)*dd - d
	az <- exp(log.Fz)
	a_az <- a - az
	az_zz <- az*zz
	v <- (sigma - 0.001)/sigma

	d.mu <- ColSums(x*((a_az)*weights/sigma))
	d.sigma <- ColSums(x*(v*(a*z1 - d - az_zz)*weights))
	if(df > 1){
		Theta <- t(matrix(theta[h], length(h), length(d)))
		By.h <- By[,h,drop = FALSE]*Theta
		Bz.h <- Bz[,h,drop = FALSE]*Theta
		B1y.h <- B1y[,h,drop = FALSE]*Theta
		d.theta <- B1y.h*(d/u1y) - By.h*(a/sigma) + Bz.h*(az/sigma)
		d.theta <- ColSums(d.theta*weights)
	}
	else{d.theta <- NULL}
	if(deriv == 1){return(cbind(d.mu,d.sigma,d.theta))}
	gradient <- c(d.mu,d.sigma,d.theta)

	# Hessian
		
	dF_dg <- exp(log.S + log.F - log.sigma)
	dF_dmu <- -dF_dg
	dF_dsigma <- dF_dmu*z1
	dFz_dg <- exp(log.Sz + log.Fz - log.sigma)
	dFz_dmu <- -dFz_dg
	dFz_dsigma <- dFz_dmu*zz


	H_mu <- t(x)%*%(x*((dF_dmu*dd - dFz_dmu)*weights/sigma))
	H_sigma <- t(x)%*%(x*(weights*((v^2)*(z1*(dd*sigma*dF_dsigma - a) - zz*(sigma*dFz_dsigma - az)) + v*(0.001/sigma*(a*z1 - d - az_zz)))))
	H_mu_sigma <- t(x)%*%(x*(v*(dF_dsigma*dd - dFz_dsigma - a_az/sigma)*weights))

	if(df > 1){
		H_theta <- (diag(d.theta, df - 1, df - 1) 
			- t(B1y.h)%*%(B1y.h*(weights*d/u1y^2)) 
			- t(By.h)%*%(By.h*(weights*dd/sigma*dF_dg)) 
			+ t(Bz.h)%*%(Bz.h*(weights/sigma*dFz_dg)))
		H_mu_theta <- t(x)%*%((By.h*(dF_dg*dd) - Bz.h*dFz_dg)*weights/sigma)
		H_sigma_theta <- t(x)%*%(By.h*(weights*v*(a/sigma + z1*dd*dF_dg)) - Bz.h*(weights*v*(az/sigma + zz*dFz_dg)))

		H <- cbind(H_mu, H_mu_sigma, H_mu_theta)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma, H_sigma_theta))
		H <- rbind(H, cbind(t(H_mu_theta), t(H_sigma_theta), H_theta))
	}
	else{
		H <- cbind(H_mu, H_mu_sigma)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma))	
	}


	# finish #

	attributes(loglik) <- list(gradient = -gradient, hessian = -H)

	-loglik
}

L.loglik.cens.trunc <- function(theta, Bz, By, B1y, x,d,d1,dd,weights, notrunc, deriv = 2){

	ColSums <- (if(deriv == 1) I else colSums)

	h <- 1:((df <- ncol(By)) - 1)
	q <- ncol(x)
	mu <- theta[1:q]
	sigma <- theta[(q + 1):(2*q)]
	theta <- (if(df > 1) c(exp(theta[(2*q + 1):length(theta)]), 1) else 1)
	if(any(abs(theta) == Inf)){return(Inf)}

	mu <- c(x%*%cbind(mu))
	sigma <- exp(c(x%*%cbind(sigma))) + 0.001
	if(any(sigma == Inf)){return(Inf)}
	log.sigma <- log(sigma)
	uy <- c(By%*%cbind(theta))
	u1y <- c(B1y%*%cbind(theta))
	uz <- c(Bz%*%cbind(theta))

	# log-likelihood #

	z1 <- (uy - mu)/sigma
	z2 <- log(1 + exp(-z1))
	z3 <- log(1 + exp(z1))	
	if(any(out2 <- (z2 == Inf))){z2[out2] <- -z1[out2]}
	if(any(out3 <- (z3 == Inf))){z3[out3] <- z1[out3]}

	log.F <- -z2
	log.S <- -z3
	log.f <- -z1 - log.sigma - 2*z2 + log(u1y)

	# truncation
	
	zz <- (uz - mu)/sigma
	log.Fz <- -log(1 + exp(-zz))
	log.Sz <- -log(1 + exp(zz))
	if(any(outFz <- (log.Fz == -Inf))){log.Fz[outFz] <- zz[outFz]}
	if(any(outSz <- (log.Sz == -Inf))){log.Sz[outSz] <- -zz[outSz]}
		
	ll <- log.f
	ll[-d1] <- log.S[-d1]
	loglik <- sum(ll) - sum(log.Sz)

	###
	
	if(deriv == 0){return(-loglik)}
	Bz[notrunc,] <- zz[notrunc] <- 1 # an arbitrary finite value

	# Gradient

	a <- exp(log.F)*dd - d
	az <- exp(log.Fz)
	a_az <- a - az
	az_zz <- az*zz
	v <- (sigma - 0.001)/sigma

	d.mu <- ColSums(x*((a_az)/sigma))
	d.sigma <- ColSums(x*(v*(a*z1 - d - az_zz)))

	if(df > 1){
		Theta <- t(matrix(theta[h], length(h), length(d)))
		By.h <- By[,h,drop = FALSE]*Theta
		Bz.h <- Bz[,h,drop = FALSE]*Theta
		B1y.h <- B1y[,h,drop = FALSE]*Theta
		d.theta <- B1y.h*(d/u1y) - By.h*(a/sigma) + Bz.h*(az/sigma)
		d.theta <- ColSums(d.theta)
	}
	else{d.theta <- NULL}
	if(deriv == 1){return(cbind(d.mu,d.sigma,d.theta))}
	gradient <- c(d.mu,d.sigma,d.theta)

	# Hessian
		
	dF_dg <- exp(log.S + log.F - log.sigma)
	dF_dmu <- -dF_dg
	dF_dsigma <- dF_dmu*z1
	dFz_dg <- exp(log.Sz + log.Fz - log.sigma)
	dFz_dmu <- -dFz_dg
	dFz_dsigma <- dFz_dmu*zz

	H_mu <- t(x)%*%(x*((dF_dmu*dd - dFz_dmu)/sigma))
	H_sigma <- t(x)%*%(x*(((v^2)*(z1*(dd*sigma*dF_dsigma - a) - zz*(sigma*dFz_dsigma - az)) + v*(0.001/sigma*(a*z1 - d - az_zz)))))
	H_mu_sigma <- t(x)%*%(x*(v*(dF_dsigma*dd - dFz_dsigma - a_az/sigma)))

	if(df > 1){
		H_theta <- (diag(d.theta, df - 1, df - 1) 
			- t(B1y.h)%*%(B1y.h*(d/u1y^2)) 
			- t(By.h)%*%(By.h*(dd/sigma*dF_dg)) 
			+ t(Bz.h)%*%(Bz.h*(1/sigma*dFz_dg)))
		H_mu_theta <- t(x)%*%((By.h*(dF_dg*dd) - Bz.h*dFz_dg)/sigma)
		H_sigma_theta <- t(x)%*%(By.h*(v*(a/sigma + z1*dd*dF_dg)) - Bz.h*(v*(az/sigma + zz*dFz_dg)))

		H <- cbind(H_mu, H_mu_sigma, H_mu_theta)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma, H_sigma_theta))
		H <- rbind(H, cbind(t(H_mu_theta), t(H_sigma_theta), H_theta))
	}
	else{
		H <- cbind(H_mu, H_mu_sigma)
		H <- rbind(H, cbind(t(H_mu_sigma), H_sigma))	
	}


	# finish #

	attributes(loglik) <- list(gradient = -gradient, hessian = -H)
	-loglik
}


flexPM.newton <- function(start, f, tol = 1e-5, maxit = 200, ...){

	f0 <- f(start, ..., deriv = 2)
	g <- attr(f0, "gradient")
	h <- attr(f0, "hessian")
	conv <- FALSE
	eps <- 1
	alg <- "nr"

	for(i in 1:maxit){

		if(conv | max(abs(g)) < tol){break}

		####
		
		H1 <- try(chol(h), silent = TRUE)
		if(class(H1) != "try-error"){
			if(alg == "gs"){alg <- "nr"; eps <- 1}
			delta <- chol2inv(H1)%*%g
		}
		else{
			if(alg == "nr"){alg <- "gs"; eps <- 1}
			delta <- g
		}

		####

		f1 <- Inf
		while(f1 > f0){
			new.start <- start - delta*eps
			if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
			f1 <- try(f(new.start, ..., deriv = 0), silent = TRUE)
			eps <- eps*0.5
			if(class(f1) == "try-error" || is.na(f1)){f1 <- Inf}
		}

		if(conv | f0 - f1 < tol){break}
		f1 <- f(new.start, ..., deriv = 2)
		g1 <- attr(f1, "gradient")
		h1 <- attr(f1, "hessian")

		start <- new.start; f0 <- f1; g <- g1; h <- h1
		eps <- min(eps*10,1)
	}

	list(estimate = start, n.it = i, minimum = as.numeric(f0), gradient = g, hessian = h)
}




























