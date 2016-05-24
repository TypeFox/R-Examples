`hglm.default` <-
	function(X = NULL, y = NULL, Z = NULL, family = gaussian(link = identity),
             rand.family = gaussian(link = identity), method = "EQL", conv = 1e-6, maxit = 50, 
			 startval = NULL, fixed = NULL, random = NULL, X.disp = NULL, disp = NULL, 
			 link.disp = "log", X.rand.disp = NULL, rand.disp = NULL, link.rand.disp = "log", 
			 data = NULL, weights = NULL, fix.disp = NULL, offset = NULL, RandC = ncol(Z), 
			 sparse = TRUE, vcovmat = FALSE, calc.like = FALSE, bigRR = FALSE, verbose = FALSE, ...) {

Call <- match.call()
if (is.null(X)) stop("X is missing with no default")
if (max(inherits(X, c("Matrix", "matrix"))) != 1) stop("X is not a sparse/dense matrix.")
x <- as.matrix(X)
y <- as.vector(y)
if (nrow(x) != length(y)) stop("Length of X and y differ. Remove all NA before input to the hglm function.") else nobs <- nrow(x)

### Check method ###
if (!(method) %in% c('EQL', 'EQL1', 'HL11')) stop('Invalid method option.')

### Check data consistency ###
if (is.null(Z)) stop("Design matrix for the random effects is missing.")
if (!is.null(Z)) {
	if (max(inherits(Z, c("Matrix", "matrix"))) != 1) stop("Z is not a sparse/dense matrix.")
	z <- Z
	if (nrow(x) != nrow(z)) stop("Length of X and Z differ. Remove all NA before input to the hglm function.")
	nRand <- cumsum(RandC)
	k <- length(nRand) ### gets number of random effects
	if (max(nRand) != ncol(z)) stop("sum(RandC) differs from number of columns in Z")
	if (k == 1) {
		colidx <- 1:nRand
	} else {
		colidx <- c()
		for (i in k:2) colidx[[i]] <- (nRand[i - 1] + 1):(nRand[i])
		colidx[[1]] <- 1:nRand[1]
	}
} else stop("Random effects are missing with no default.")
if (!is.null(X.disp)) x.disp <- as.matrix(X.disp) else x.disp <- NULL

####Check NA in y, X, Z and X.disp####
#### added by lrn 2015-03-24

if ( sum( is.na( cbind( y, x, Z, x.disp ) ) ) > 0) warning( "Remove all NA before input to the hglm function.", immediate.=TRUE)

### Check prior weights ###
if (is.null(weights)) {
    prior.weights <- rep(1, nobs + nRand[k])
} else {
	if (!is.numeric(weights) || any(weights <= 0)) stop("Weights must be a numeric vector of positive values")
	if (length(weights) < nobs) stop("Length of the weights differ from the length of the data")
    prior.weights <- c(weights, rep(1, nRand[k]))
}

### Check offset ###
if (is.null(offset)) {
	off <- rep(0, nobs)
} else {
	if (!is.numeric(offset) || length(offset) != nobs) stop("Offset must be a numeric vector of the same length as the data")
	off <- as.numeric(offset)
}
### Data consistency checked ###


### Check consistency of GLM family and number of random effects terms
if (class(rand.family) != 'family' & length(rand.family) != length(RandC)) stop('Number of given families and number of random effects terms differ.')
# 130918 xia: bug fixed - length(family) should be length(rand.family)

### Get GLM family and link ###
if (is.character(family)) family <- get(family)
if (is.function(family)) family <- eval(family)
if (class(rand.family) == 'family') {
	if (is.character(rand.family)) rand.family <- get(rand.family)
	if (is.function(rand.family)) rand.family <- eval(rand.family)
	if (rand.family$family == "Gamma") {
		GammaLink <- rand.family$link
		if (GammaLink == 'log') rand.family <- eval(GAMMA(link = 'log')) ### Note this defintion of random Gamma effects
		if (GammaLink == 'identity') rand.family <- eval(GAMMA(link = 'identity'))
		if (GammaLink == 'inverse') rand.family <- eval(GAMMA('inverse'))
		if (!(GammaLink %in% c('log', 'identity', 'inverse'))) rand.family <- eval(GAMMA(link = 'log'))
	}
	if (rand.family$family %in% c("CAR", "SAR")) {
		link.rand.disp <- rand.family$link.rand.disp
		#z <- tcrossprod(z, rand.family$Dvec)
		z <- z %*% as.matrix(rand.family$Dvec) # svd to eigen bug fix by Moudud
		X.rand.disp <- list(model.matrix(~ rand.family$Deigen))
	}
} else {
	X.rand.disp <- c()
	X.rand.disp[[length(rand.family) + 1]] <- c(NA, NA)
	for (i in 1:length(rand.family)) {
		if (is.character(rand.family[[i]])) rand.family[[i]] <- get(rand.family)
		if (is.function(rand.family[[i]])) rand.family[[i]] <- eval(rand.family)
		if (rand.family[[i]]$family == "Gamma") {
			GammaLink <- rand.family[[i]]$link
			if (GammaLink == 'log') rand.family[[i]] <- eval(GAMMA(link = 'log')) ### Note this defintion of random Gamma effects
			if (GammaLink == 'identity') rand.family[[i]] <- eval(GAMMA(link = 'identity'))
			if (GammaLink == 'inverse') rand.family[[i]] <- eval(GAMMA('inverse'))
			if (!(GammaLink %in% c('log', 'identity', 'inverse'))) rand.family[[i]] <- eval(GAMMA(link = 'log'))
		}
		if (rand.family[[i]]$family %in% c("CAR", "SAR")) {
			link.rand.disp <- rand.family[[i]]$link.rand.disp
			#z[,colidx[[i]]] <- tcrossprod(z[,colidx[[i]]], rand.family[[i]]$Dvec)
			z[,colidx[[i]]] <- z[,colidx[[i]]] %*% as.matrix(rand.family[[i]]$Dvec) # svd to eigen bug fix by Moudud
			X.rand.disp[[i]] <- model.matrix(~ rand.family[[i]]$Deigen)
		}
	}
}
### GLM family and link are checked ###
### Only the GLM families (Lee et al. 2006) will pass this test ###

### Get augmented response, psi (Lee et al., 2006) ###
if (class(rand.family) == 'family') {
	if (rand.family$family %in% c("gaussian", "CAR", "SAR")) {
    	psi <- rep(0, nRand[k])
	} else if (rand.family$family == "GAMMA" || rand.family$family == "inverse.gamma") {
    	psi <- rep(1, nRand[k])
	} else if (rand.family$family == "Beta") {
    	psi <- rep(1/2, nRand[k])
	} else {
    	stop(paste("random.family = ", rand.family$family, " is not recognized as a member of the GLM family"))
	}
} else {
	psi <- c()
	for (i in 1:length(rand.family)) {
		if (rand.family[[i]]$family %in% c("gaussian", "CAR", "SAR")) {
			psi <- c(psi, rep(0, RandC[i]))
		} else if (rand.family[[i]]$family == "GAMMA" || rand.family[[i]]$family == "inverse.gamma") {
			psi <- c(psi, rep(1, RandC[i]))
		} else if (rand.family[[i]]$family == "Beta") {
			psi <- c(psi, rep(1/2, RandC[i]))
		} else {
			stop(paste("random.family = ", rand.family[[i]]$family, " is not recognized as a member of the GLM family"))
		}
	}
}

if (!is.character(link.disp)) link.disp <- deparse(link.disp)
if (link.disp == "log") {
	DispFamily <- Gamma(link = "log")
} 
else if (link.disp == "identity") {
    DispFamily <- Gamma(link = "identity")
} 
else if (link.disp == "inverse") {
    DispFamily <- Gamma(link = "inverse")
} 
else stop("link.disp must be a valid link for the Gamma family GLM") 

### bigRR? ###
if (!bigRR) {
if (ncol(z) > nrow(z) + 1 & length(RandC) == 1) {
	if (rand.family$family == "gaussian") {
		cat("NOTE: You are fitting a model with one Gaussian random effect term,\n",
			"and the number of effects (p) is greater than the number of\n",
			"observations (n). Consider turning on the argument 'bigRR' that may\n",
			"speed up a lot if p >> n.\n")
	}
}

### Check starting values ###
g1 <- glm(y ~ x - 1, family = family, weights = weights, offset = off) 
if (!is.null(startval)) {
    if (!is.numeric(startval)) stop("Non-numeric starting value is not allowed.")
    if (length(startval) < ncol(x) + k + sum(RandC)) stop("Too few starting values. See hglm documentation.")
    if ((family$family == "gaussian" || family$family=="Gamma") & length(startval) < ncol(x) + nRand[k] + k + 1) stop("Too few starting values. See the documentation of hglm.")
    b.hat <- startval[1:ncol(x)]
    if (length(startval) > (ncol(x) + k + sum(RandC))) {
    	init.sig.e <- as.numeric(startval[length(startval)])
    } else {
    	init.sig.e <- 1
    }
	init.u <- startval[(ncol(x) + 1):(ncol(x) + sum(RandC))]
	init.sig.u <- as.numeric(startval[(ncol(x) + sum(RandC) + 1)])
	if (min(init.sig.e, init.sig.u) < 1e-4) stop("Unacceptable initial value is supplied for the variance parameter.")
} else {
	### Generate default initial values of the fixed effects via a GLM ###
    b.hat <- as.numeric(coef(g1))
    init.sig.u <- (init.sig.e <- as.numeric(.6*deviance(g1)/g1$df.residual))*.66
    if (k > 1) init.sig.u <- rep(init.sig.u/k, k)
    if (!is.null(fix.disp)) {
    	if (!is.numeric(fix.disp) | fix.disp <= 0) stop("\"fix.disp\" must be numeric and greater than 1e-4.")
    	init.sig.e <- as.numeric(fix.disp)
    }
    if (min(init.sig.u) < 1e-4) {
    	init.sig.u < rep(.1, k)
    	message("0.1 is chosen as the initial values for the dispersion parameter of the random effects.")
    }
    if (init.sig.e < 1e-4) {
    	init.sig.e <- .1
    	message("0.1 is chosen as the initial values for the dispersion parameter of the mean model.")
    }
    init.u <- rep(0, nRand[k])
}

### Create Augmented data ###
if (!is.null(colnames(x))) {
    x.names <- colnames(x)
    colnames(x) <- NULL
} else x.names <- paste("X.", 1:ncol(x), sep = "")
if (!is.null(z)) {
    if (!is.null(colnames(z))) {
    	z.names <- colnames(z)
    	colnames(z) <- NULL
    } else z.names <- paste("Z.", 1:ncol(z), sep = "")
    Augy <- c(y, psi)      
    AugXZ <- cBind(x,z)
    XX1 <- Matrix(0, nrow = nRand[k], ncol = ncol(x))
    ZZ2 <- Diagonal(nRand[k])
    ZZ2 <- cBind(XX1, ZZ2)
    AugXZ <- rBind(AugXZ, ZZ2)
    rm(list = c("XX1", "ZZ2"))
    colnames(AugXZ) <- c(x.names, z.names)
} else {
    Augy <- y
    AugXZ <- x
}

##################################
### Begin Parameter estimation ###
iter <- 1
if (!is.null(z)) phi <- rep(init.sig.u, RandC)  ### Random effects variance
tau <- rep(init.sig.e, nobs)
if (class(rand.family) == 'family') {
	v.i <- rand.family$linkinv(init.u)
} else {
	v.i <- c()
	for (i in 1:k) v.i <- c(v.i, rand.family[[i]]$linkinv(init.u[colidx[[i]]]))
}
eta.i <- as.numeric(x%*%b.hat) + off
eta0 <- eta.i
mu.i <- family$linkinv(eta.i)
dmu_deta <- family$mu.eta(eta.i)
zi <- eta.i - off + (y - mu.i)/dmu_deta
zi <- zi # - HL.correction

if (!is.null(z)) {
	zmi <- psi
	Augz <- c(zi, zmi)
	if (class(rand.family) == 'family') {
		du_dv <- rand.family$mu.eta(psi)
		w <- sqrt(as.numeric(c((dmu_deta^2/family$variance(mu.i))*(1/tau), (du_dv^2/rand.family$variance(psi))*(1/phi)))*prior.weights)
	} else {
		du_dv <- c()
		w <- as.numeric((dmu_deta^2/family$variance(mu.i))*(1/tau))
		for (i in 1:k) {
			du_dv <- c(du_dv, rand.family[[i]]$mu.eta(psi[colidx[[i]]]))
			w <- c(w, as.numeric((du_dv[colidx[[i]]]^2/rand.family[[k]]$variance(psi[colidx[[i]]]))*(1/phi[colidx[[i]]])))
		}
		w <- sqrt(w*prior.weights)
	}
} else {
	w <- sqrt(as.numeric((dmu_deta^2/family$variance(mu.i))*(1/tau))*prior.weights)
} 
n <- length(Augy)
p <- ncol(AugXZ)

if (method == 'HL11') {
	method <- 'EQL1'
	cat('HL11 method option is now renamed as EQL1. Consider updating your code.\n')
}
HL.correction <- 0
bad <- NULL

while (iter <= maxit) {
	#ii <<- iter
	#if (iter == 1) {
	#	mmein <- list(Augy = Augy, AugXZ = AugXZ, starting.delta = c(b.hat, v.i), tau = tau, phi = phi, 
    #                 n.fixed = ncol(x), n.random = nRand[k], weights.sqrt = w, prior.weights, family, 
    #                 rand.family, maxit, sparse = sparse, off = off, tol = 1e-7, colidx = colidx,
	#				 HL.correction = HL.correction)
	#	save(mmein, file = 'mmein.RData')
	#}
    g.mme <- GLM.MME(Augy = Augy, AugXZ = AugXZ, starting.delta = c(b.hat, v.i), tau = tau, phi = phi, 
                     n.fixed = ncol(x), n.random = nRand[k], weights.sqrt = w, prior.weights, family, 
                     rand.family, maxit, sparse = sparse, off = off, tol = 1e-7, colidx = colidx,
					 HL.correction = HL.correction)
    b.hat <- g.mme$b.hat
    eta.i <- g.mme$eta.i
    v.i <- g.mme$v.i
    Augz <- g.mme$Augz
    dev <- g.mme$dev
    hv <- g.mme$hv
    resid <- g.mme$resid
	if (any(dev < 1e-8)) {
		warning("Residuals numerically 0 are replaced by 1e-8")
		dev <- ifelse(dev < 1e-8, 1e-8, dev)
	}
    fv <- g.mme$fv
    hv <- g.mme$hv
	if (any(abs(hv) > 1 - 1e-8)) {
		warning("Hat-values numerically 1 are replaced by 1 - 1e-8")
		hv <- ifelse(abs(hv) > 1 - 1e-8, 1 - 1e-8, hv)
	}
    mu.i <- family$linkinv(eta.i)
    dmu_deta <- family$mu.eta(eta.i)
    if (!is.null(z)) {
		if (class(rand.family) == 'family') {
    		ui <- rand.family$linkinv(v.i)
    		du_dv <- rand.family$mu.eta(v.i)
		} else {
			ui <- du_dv <- c()
			for (i in 1:k) {
				ui <- c(ui, rand.family[[i]]$linkinv(v.i[colidx[[i]]]))
				du_dv <- c(du_dv, rand.family[[i]]$mu.eta(v.i[colidx[[i]]]))
			}
		}
    }
    if (is.null(fix.disp)) {
    	if (is.null(x.disp)) {
			#tmp <<- as.numeric(dev[1:nobs]/(1 - hv[1:nobs]))
    		g11 <- glm((as.numeric(dev[1:nobs]/(1 - hv[1:nobs]))) ~ 1, family = DispFamily, weights = as.numeric((1 - hv[1:nobs])/2))
        	if (length(g11$coef) == 1) {
      			sigma2e <- DispFamily$linkinv(as.numeric(g11$coef[1]))
      		} else {
      			sigma2e <- NULL
      		}
			tau <- as.numeric(g11$fitted.values) 
    	} else {
			#tmp <<- as.numeric(dev[1:nobs]/(1 - hv[1:nobs]))
			g11 <- glm((as.numeric(dev[1:nobs]/(1 - hv[1:nobs]))) ~ x.disp - 1, family = DispFamily, weights = as.numeric((1 - hv[1:nobs])/2))
        	if (length(g11$coef) == 1) {
      			sigma2e <- DispFamily$linkinv(as.numeric(g11$coef[1]))
      		} else {
      			sigma2e <- NULL
			}
    		tau <- as.numeric(g11$fitted.values) ### Error variance updated
    	}
    	disp.fv <- g11$fitted.values
    	disp.resid <- residuals(g11)/sqrt(1 - hatvalues(g11))
    } else {
    	sigma2e <- as.numeric(fix.disp)
    	disp.fv <- disp.resid <- NULL
    }
	VC2 <- vector("list", k)
	devused <- nobs
	
	for (K in 1:k) {
		devtouse <- nobs + nRand[K]
    	devu <- as.numeric(dev[(devused + 1):devtouse])
    	hvu <- as.numeric(1 - hv[(devused + 1):devtouse])

		#tmp2 <<- as.numeric(devu/hvu)
		#XK <<- X.rand.disp
		#w <<- hvu/2
		
		if (is.null(X.rand.disp[[K]])) {
			g12 <- glm(devu/hvu ~ 1, family = Gamma(link = log), weights = hvu/2)
			sigma2u <- exp(as.numeric(g12$coef[1]))
		} else {
			if (class(rand.family) == 'family') {
				if (!(rand.family$family %in% c("CAR", "SAR"))) {
					g12 <- glm(devu/hvu ~ X.rand.disp[[K]] - 1, family = Gamma(link = link.rand.disp), weights = hvu/2)
					sigma2u <- NULL
				} else if (rand.family$family == 'CAR') {
					g12 <<- try(glm(devu/hvu ~ X.rand.disp[[K]] - 1, family = Gamma(link = link.rand.disp), weights = hvu/2), silent = TRUE)
					if (inherits(g12, 'try-error')) {
						for (l in 1:100) {
							ndev <- length(devu)
							idx <- sample(1:ndev, round(.75*ndev))
							g120 <- try(glm(devu[idx]/hvu[idx] ~ X.rand.disp[[K]][idx,] - 1, family = Gamma(link = link.rand.disp), weights = hvu[idx]/2), silent = TRUE)
							if (!inherits(g120, 'try-error')) {
								mustart <- rep(NA, ndev)
								mustart[idx] <- g120$fitted.values
								mustart[is.na(mustart)] <- mean(devu/hvu)
								g12 <- try(glm(devu/hvu ~ X.rand.disp[[K]] - 1, family = Gamma(link = link.rand.disp), weights = hvu/2, mustart = mustart), silent = TRUE)
							}
							if (!inherits(g12, 'try-error')) break
							#if (inherits(g12, 'try-error')) l <- l + 1 else l <- length(devu)
							#print(l)
						}
					}
					if (inherits(g12, 'try-error')) warning('Internal Gamma GLM failed for CAR family!')
					#print(g12)
					sigma2u <- NULL
					CAR.tau <- 1/g12$coef[1]
					CAR.rho <- -g12$coef[2]/g12$coef[1]
				} else {
					g12 <- try(glm(devu/hvu ~ X.rand.disp[[K]] - 1, family = Gamma(link = inverse.sqrt()), weights = hvu/2), silent = TRUE)
					if (inherits(g12, 'try-error')) {
						for (l in 1:100) {
							ndev <- length(devu)
							idx <- sample(1:ndev, round(.75*ndev))
							g120 <- try(glm(devu[idx]/hvu[idx] ~ X.rand.disp[[K]][idx,] - 1, family = Gamma(link = inverse.sqrt()), weights = hvu[idx]/2), silent = TRUE)
							if (!inherits(g120, 'try-error')) {
								mustart <- rep(NA, ndev)
								mustart[idx] <- g120$fitted.values
								mustart[is.na(mustart)] <- mean(devu/hvu)
								g12 <- try(glm(devu/hvu ~ X.rand.disp[[K]] - 1, family = Gamma(link = inverse.sqrt()), weights = hvu/2, mustart = mustart), silent = TRUE)
							}
							if (!inherits(g12, 'try-error')) break
							#if (inherits(g12, 'try-error')) l <- l + 1 else l <- length(devu)
							#print(l)
						}
					}
					if (inherits(g12, 'try-error')) warning('Internal Gamma GLM failed for SAR family!')
					#print(g12)
					sigma2u <- NULL
					SAR.tau <- 1/g12$coef[1]**2
					SAR.rho <- -g12$coef[2]/g12$coef[1]
				}
			} else {
				if (!(rand.family[[K]]$family %in% c("CAR", "SAR"))) {
					g12 <- glm(devu/hvu ~ X.rand.disp[[K]] - 1, family = Gamma(link = link.rand.disp), weights = hvu/2)
					sigma2u <- NULL
				} else if (rand.family[[K]]$family == 'CAR') {
					g12 <- try(glm(devu/hvu ~ X.rand.disp[[K]] - 1, family = Gamma(link = link.rand.disp), weights = hvu/2), silent = TRUE)
					if (inherits(g12, 'try-error')) {
						for (l in 1:100) {
							ndev <- length(devu)
							idx <- sample(1:ndev, round(.75*ndev))
							g120 <- try(glm(devu[idx]/hvu[idx] ~ X.rand.disp[[K]][idx,] - 1, family = Gamma(link = link.rand.disp), weights = hvu[idx]/2), silent = TRUE)
							if (!inherits(g120, 'try-error')) {
								mustart <- rep(NA, ndev)
								mustart[idx] <- g120$fitted.values
								mustart[is.na(mustart)] <- mean(devu/hvu)
								g12 <- try(glm(devu/hvu ~ X.rand.disp[[K]] - 1, family = Gamma(link = link.rand.disp), weights = hvu/2, mustart = mustart), silent = TRUE)
							}
							if (!inherits(g12, 'try-error')) break
							#if (inherits(g12, 'try-error')) l <- l + 1 else l <- length(devu)
							#print(l)
						}
					}
					if (inherits(g12, 'try-error')) warning('Internal Gamma GLM failed for CAR/SAR family!')
					#print(g12)
					sigma2u <- NULL
					CAR.tau <- 1/g12$coef[1]
					CAR.rho <- -g12$coef[2]/g12$coef[1]
				} else {
					g12 <- try(glm(devu/hvu ~ X.rand.disp[[K]] - 1, family = Gamma(link = inverse.sqrt()), weights = hvu/2), silent = TRUE)
					if (inherits(g12, 'try-error')) {
						for (l in 1:100) {
							ndev <- length(devu)
							idx <- sample(1:ndev, round(.75*ndev))
							g120 <- try(glm(devu[idx]/hvu[idx] ~ X.rand.disp[[K]][idx,] - 1, family = Gamma(link = inverse.sqrt()), weights = hvu[idx]/2), silent = TRUE)
							if (!inherits(g120, 'try-error')) {
								mustart <- rep(NA, ndev)
								mustart[idx] <- g120$fitted.values
								mustart[is.na(mustart)] <- mean(devu/hvu)
								g12 <- try(glm(devu/hvu ~ X.rand.disp[[K]] - 1, family = Gamma(link = inverse.sqrt()), weights = hvu/2, mustart = mustart), silent = TRUE)
							}
							if (!inherits(g12, 'try-error')) break
							#if (inherits(g12, 'try-error')) l <- l + 1 else l <- length(devu)
							#print(l)
						}
					}
					if (inherits(g12, 'try-error')) warning('Internal Gamma GLM failed for CAR/SAR family!')
					#print(g12)
					sigma2u <- NULL
					SAR.tau <- 1/g12$coef[1]**2
					SAR.rho <- -g12$coef[2]/g12$coef[1]
				}
			}
		}
    	if (K == 1) {
        	phi[1:nRand[K]] <- g12$fitted.value
    	} else {
        	phi[(nRand[(K - 1)] + 1):nRand[K]] <- g12$fitted.value #rep(sigma2u, RandC[K])
    	}  ### Random effects variance updated
    	SummarySt <- summary(g12, dispersion = 1)
    	VC2[[K]] <- SummarySt$coefficients[,1:2]
		rm(list = c("g12", "SummarySt", "devu", "hvu"))
		devused <- devtouse
	}
	
	#if (is.null(colnames(Z))) { # names simplified, xia 130221
    	names(VC2) <- paste(".|Random", 1:k, sep = "")
    #} else {
    #	CVnames <- unlist(strsplit(z.names[nRand], "\\:"))
    #	names(VC2) <- CVnames[seq(1, length(CVnames), by = 2)]
	#}
	if (verbose) {
		cat("\n-------------------")
		cat("\nIteration", iter, "\n")
		cat("-------------------\n")
		#cat("Fitted residual variance:", sigma2e, "\n")
		#cat("Fitted random effects variance:", exp(as.numeric(unlist(VC2)[seq(1, 2*k, 2)])), "\n")
		cat("Sum of squared linear predictor:", sum(eta.i^2), "\n")
		cat("Convergence:", sum((eta0 - eta.i)^2)/sum(eta.i^2), "\n")
	}
	if (sum((eta0 - eta.i)^2) < conv*sum(eta.i^2) & iter > 1 ) break
    if (!is.null(z)) {
		if (class(rand.family) == 'family') {
			w <- sqrt(as.numeric(c((dmu_deta^2/family$variance(mu.i))*(1/tau), (du_dv^2/rand.family$variance(ui))*(1/phi)))*prior.weights)
		} else {
			w <- as.numeric((dmu_deta^2/family$variance(mu.i))*(1/tau))
			for (i in 1:k) {
				w <- c(w, as.numeric((du_dv[colidx[[i]]]^2/rand.family[[k]]$variance(ui[colidx[[i]]]))*(1/phi[colidx[[i]]])))
			}
			w <- sqrt(w*prior.weights)
		}
	} else {
		w <- sqrt(as.numeric((dmu_deta^2/family$variance(mu.i))*(1/tau))*prior.weights)
	}
	if (method == 'EQL1' & class(rand.family) == 'family' & iter >= 2*(is.null(fix.disp))) if (rand.family$family %in% c('gaussian', 'CAR', 'SAR')) HL.correction <- HL11(fv = fv, w = w, Z = Z, family = family, tau = tau)
	eta0 <- eta.i
    iter <- iter + 1
}
if (method == 'EQL1' & class(rand.family) == 'family') if (!(rand.family$family %in% c('gaussian', 'CAR', 'SAR'))) {
	warning('EQL1 correction is not implemented yet for non-Gaussian random effects. EQL estimates are provided!')
	method <- 'EQL'
}
names(b.hat) <- x.names
if (!is.null(z)) names(ui) <- z.names
fixef <- b.hat                        
if (!is.null(z)) ranef <- ui else ranef <- phi <- NULL
if (class(rand.family) == 'family') {
	if (rand.family$family %in% c('CAR', 'SAR')) ranef <- rand.family$Dvec %*% ranef ## ranef for CAR calculation bug fixed by Lars 2014-01-20
} else {
	for (i in 1:k) {
		if (rand.family[[i]]$family %in% c('CAR', 'SAR')) ranef[colidx[[i]]] <- rand.family[[i]]$Dvec %*% ranef[colidx[[i]]] ## ranef for CAR calculation bug fixed by Lars 2014-01-20
	}
}

sigma6 <- mean(hv[1:nobs]) + 6*sd(hv[1:nobs])
if (max(hv[1:nobs]) > sigma6) {
	bad <- which.max(hv[1:nobs])
}

sigma2u <- numeric(k)
for (K in 1:k) sigma2u[K] <- exp(unlist(VC2[[K]])[1])

val <- list(call = Call, fixef = fixef, ranef = ranef, RandC = RandC, phi = phi, varFix = sigma2e, 
            varRanef = sigma2u, CAR.tau = NULL, CAR.rho = NULL, SAR.tau = NULL, SAR.rho = NULL, iter = iter, 
			Converge = "did not converge", SeFe = NULL, SeRe = NULL,
            dfReFe = NULL, SummVC1 = NULL, SummVC2 = NULL, method = method, dev = dev, hv = hv, 
            resid = resid, fv = fv, disp.fv = disp.fv, disp.resid = disp.resid, link.disp = link.disp, 
			link.rand.disp = link.rand.disp, vcov = NULL, likelihood = NULL, call.rand.family = rand.family,
			null.model = g1, bad = bad)

if (iter < maxit & all(sigma2u/(sum(sigma2u)+sigma2e) < (1-1e-4))) {
	val$Converge <- "converged"
	### Calculate the standard errors of the fixed and random effects ###
	p1 <- 1:p
	QR <- g.mme$qr
	if (sparse == TRUE) {
		qqq <- 1+QR@q
		OPT <- options()
		options(warn = -1)
		QR <- qr.R(QR)
		covmat <- as.matrix(chol2inv(QR))
		covmat[qqq,qqq] <- covmat
		options(OPT)
		SeFeRe <- sqrt(diag(covmat))
	} else {
		QR <- qr.R(QR)
		covmat <- chol2inv(QR)
		SeFeRe <- sqrt(diag(covmat))
	}
	
	if (class(rand.family) == 'family') {
		if (rand.family$family == "CAR") {
			val$CAR.tau = as.numeric(CAR.tau)
			val$CAR.rho = as.numeric(CAR.rho)
		}
		if (rand.family$family == "SAR") {
			val$SAR.tau = as.numeric(SAR.tau)
			val$SAR.rho = as.numeric(SAR.rho)
		}
	} else {
		for (i in 1:k) {
			if (rand.family[[i]]$family == "CAR") {
				val$CAR.tau = as.numeric(CAR.tau)
				val$CAR.rho = as.numeric(CAR.rho)
			}
			if (rand.family[[i]]$family == "SAR") {
				val$SAR.tau = as.numeric(SAR.tau)
				val$SAR.rho = as.numeric(SAR.rho)
			}
		}
	}
	
	val$SeFe <- SeFeRe[1:NCOL(x)]
	val$SeRe <- SeFeRe[(NCOL(x) + 1):length(SeFeRe)]
	### Calculate the deviance degrees of freedom (Lee et al. 2006, pp 193) ###
	Sigma0 <- as.numeric(g.mme$wt)[1:nobs]
	Pd <- sum(diag(covmat%*%crossprod(as.matrix(AugXZ[1:nobs,])*Sigma0)))
	val$dfReFe <- round(nobs - Pd)
    if (is.null(fix.disp)) {
		SummVC1 <- summary(g11, dispersion = 1)
		SummVC1 <- SummVC1$coefficients[,1:2]
		if (!is.null(row.names(SummVC1))) {
			dnames <- row.names(SummVC1)
			row.names(SummVC1) <- sub("x.disp", '', dnames)
		}
		val$SummVC1 <- SummVC1
	} else {
		val$SummVC1 <- fix.disp
	}
	if (!is.null(z)) {
		val$SummVC2 <- VC2
		val$varRanef <- exp(as.numeric(unlist(VC2)[seq(1, 2*k, 2)]))
    }
    
    if (vcovmat) {
    	val$vcov <- Matrix(diag(tau))
    	if (class(rand.family) == 'family') {
    		val$vcov <- val$vcov + tcrossprod(Matrix(z))*val$varRanef
    	} else {
    		val$vcov <- val$vcov + tcrossprod(Matrix(z[,1:RandC[1]]))*val$varRanef[1]
    		for (J in 2:k) {
    			val$vcov <- val$vcov + tcrossprod(Matrix(z[,(nRand[J - 1] + 1):nRand[J]]))*val$varRanef[J]
    		}
    	}
	}
	
	if (calc.like) {
		z <- as.matrix(z)
		val$likelihood <- likelihood(val, y, X, z, family = family, weights = weights)
	}
	
##### Calculate Profile Likelihood ##########
#  Sigma <- Matrix(diag(tau))
#    if (!is.null(z)) {
#     D <- Matrix(diag(phi))
#V <- z%*%D%*%t(z)+Sigma
#V<-tcrossprod((z*sqrt(phi)))+Sigma
#   } else {
#     V <- Sigma
#   }
#   V.inv <- solve(V)
#   temp <- as.numeric(y-x%*%fixef)
#   logdet.V <- as.numeric(determinant(V,log=TRUE))
#profile <- as.numeric(-nobs/2*log(2*pi)-1/2*logdet.V-1/2*t(temp)%*%V.inv%*%temp-
#    1/2*(detterminant((t(x)%*%V.inv%*%x),log=TRUE))+k/2*log(2*pi))
##   #loglihood <- as.numeric(-nobs/2*log(2*pi)-1/2*logdet.V-1/2*t(temp)%*%V.inv%*%temp)
#   if(method=="REML"){
#       val$ProfLogLik<-NULL #as.numeric(profile)
} else {
	val$iter <- iter - 1
}

} else {
	model <- bigRR(formula = NULL, y = y, X = X, Z = Z, data = data, shrink = NULL, weight = NULL,
   				   family = family, lambda = NULL, impute = FALSE, tol.err = 1e-6, 
				   tol.conv = conv, only.estimates = FALSE, GPU = FALSE)
		   
	val <- list(call = Call, fixef = model$beta, ranef = model$u, RandC = RandC, phi = rep(model$lambda, RandC), varFix = model$phi, 
			   varRanef = model$lambda, iter = model$hglm$iter, Converge = model$hglm$Converge, SeFe = model$hglm$SeFe, SeRe = model$hglm$SeRe,
			   dfReFe = model$hglm$dfReFe, SummVC1 = model$hglm$SummVC1, SummVC2 = model$hglm$SummVC2, method = method, dev = model$hglm$dev, hv = NULL, 
			   resid = model$hglm$resid, fv = model$hglm$fv, disp.fv = NULL, disp.resid = NULL, link.disp = NULL, 
			   vcov = model$hglm$vcov, likelihood = model$hglm$likelihood, null.model = model$hglm$null.model)
}

class(val) <- "hglm"
	
return(val)
	
}
