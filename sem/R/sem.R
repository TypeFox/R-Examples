# last modified 2012-11-03 by J. Fox

sem <- function(model, ...){
	if (is.character(model)) class(model) <- "semmod"
	UseMethod("sem", model)
}

sem.semmod <- function(model, S, N, data, raw=identical(na.action, na.pass), obs.variables=rownames(S), 
		fixed.x=NULL, formula= ~ ., na.action=na.omit, robust=!missing(data), debug=FALSE, 
		optimizer=optimizerSem, objective=objectiveML, ...){
	parse.path <- function(path) {                                           
		path.1 <- gsub("-", "", gsub(" ","", path))
		direction <- if (regexpr("<>", path.1) > 0) 2 
				else if (regexpr("<", path.1) > 0) -1
				else if (regexpr(">", path.1) > 0) 1
				else stop(paste("ill-formed path:", path))
		path.1 <- strsplit(path.1, "[<>]")[[1]]
		list(first=path.1[1], second=path.1[length(path.1)], direction=direction)
	}
	any.NA <- FALSE
	unique.patterns <- valid.pattern <- valid <- pattern.number <- valid.data.patterns <- NULL
	if (missing(S)){
		if (missing(data)) stop("S and data cannot both be missing")
		N.all <- nrow(data)
		data <- model.frame(formula, data=data, na.action=na.action)
		data <- model.matrix(formula, data=data)
		colnames(data)[colnames(data) == "(Intercept)"] <- "Intercept"
		S <- if (raw) rawMoments(data, na.rm=TRUE) else {
					data <-  data[, colnames(data) != "Intercept"]
					cov(data, use="complete.obs")
				}
		N <- nrow(data)
		if (N < N.all) warning(N - N.all, " observations removed due to missingness")        
		if (identical(na.action, na.pass) && any(is.na(data))){
				any.NA <- TRUE
				valid <- !is.na(data)
				colnames(valid) <- colnames(data)
			}
			else {
					valid <- matrix(TRUE, nrow(data), ncol(data))
					colnames(valid) <- colnames(data)
			}
		}
	if ((!is.matrix(model)) | ncol(model) != 3) stop ("model argument must be a 3-column matrix")
	startvalues <- as.numeric(model[,3])
	par.names <- model[,2]
	n.paths <- length(par.names)
	heads <- from <- to <- rep(0, n.paths)
	for (p in 1:n.paths){
		path <- parse.path(model[p,1])
		heads[p] <- abs(path$direction)
		to[p] <- path$second
		from[p] <- path$first
		if (path$direction == -1) {
			to[p] <- path$first
			from[p] <- path$second
		}
	}
	ram <- matrix(0, p, 5)
	all.vars <- unique(c(to, from))
	latent.vars <- setdiff(all.vars, obs.variables)
	not.used <- setdiff(obs.variables, all.vars)
	if (length(not.used) > 0){
		rownames(S) <- colnames(S) <- obs.variables
		obs.variables <- setdiff(obs.variables, not.used)
		S <- S[obs.variables, obs.variables]
		warning("The following observed variables are in the input covariance or raw-moment matrix ",
				"but do not appear in the model:\n",
				paste(not.used, collapse=", "), "\n")
	}
	vars <- c(obs.variables, latent.vars)
	pars <- na.omit(unique(par.names))
	ram[,1] <- heads
	ram[,2] <- apply(outer(vars, to, "=="), 2, which)
	ram[,3] <- apply(outer(vars, from, "=="), 2, which)   
	par.nos <- apply(outer(pars, par.names, "=="), 2, which)
	if (length(par.nos) > 0)
		ram[,4] <- unlist(lapply(par.nos, function(x) if (length(x) == 0) 0 else x))
	ram[,5]<- startvalues
	colnames(ram) <- c("heads", "to", "from", "parameter", "start")
	if (!is.null(fixed.x)) fixed.x <- apply(outer(vars, fixed.x, "=="), 2, which)
	n <- length(obs.variables)
	m <- length(all.vars)
	t <- length(pars)
	if (debug) {
		cat("\n observed variables:\n") 
		print(paste(paste(1:n,":", sep=""), obs.variables, sep=""))
		cat("\n")
		if (m > n){ 
			cat("\n latent variables:\n")
			print(paste(paste((n+1):m,":", sep=""), latent.vars, sep=""))
			cat("\n")
		}
		cat("\n parameters:\n") 
		print(paste(paste(1:t,":", sep=""), pars, sep=""))
		cat("\n\n RAM:\n")
		print(ram)
	}
	if (missing(data)) data <- NULL 
	else {
		data <- data[, obs.variables]
		if (!is.null(valid)) {
			valid <- valid[, obs.variables]
			valid.pattern <- apply(valid, 1, function(row) paste(row, collapse="."))
			unique.patterns <- unique(valid.pattern)
			pattern.number <- apply(outer(valid.pattern, unique.patterns, `==`), 1, which)
			valid.data.patterns <- t(sapply(strsplit(unique.patterns, "\\."), as.logical))
		}
	}
	if (identical(objective, objectiveFIML2) || identical(objective, objectiveFIML)){
		message("NOTE: start values computed from preliminary ML fit")
        opt <- options(warn=-1)
        on.exit(options(opt)) # in case of error
		prelim.fit <- sem(ram, S=S, N=N, raw=raw, data=na.omit(data), valid=valid, param.names=pars, var.names=vars, fixed.x=fixed.x,
				semmod=model, robust=robust, debug=debug, ...)
        if (!prelim.fit$convergence) message("NOTE: preliminary ML fit may not have converged")
        options(opt)
		message("NOTE: preliminary iterations, ", prelim.fit$iterations)
		message("NOTE: iterations reported for final fit are post preliminary fit")
		coeff <- coef(prelim.fit)
		rownames(ram) <- rownames(prelim.fit$ram)[1:nrow(ram)]
		ram[names(coeff), 5] <- coeff
	}
	cls <- gsub("\\.", "", deparse(substitute(objective)))
	cls <- gsub("2", "", cls)
	result <- sem(ram, S=S, N=N, raw=raw, data=data, 
			pattern.number=pattern.number, valid.data.patterns=valid.data.patterns,
			param.names=pars, var.names=vars, fixed.x=fixed.x,
			semmod=model, robust=robust, debug=debug, optimizer=optimizer, objective=objective, cls=cls, ...)
	class(result) <- c(cls, "sem")
	result
}

sem.default <- function(model, S, N, raw=FALSE, data=NULL, start.fn=startvalues,
		pattern.number=NULL, valid.data.patterns=NULL,
		use.means=TRUE, param.names, 
		var.names, fixed.x=NULL, robust=!is.null(data), semmod=NULL, debug=FALSE,
		analytic.gradient=!identical(objective, objectiveFIML), warn=FALSE, maxiter=1000, par.size=c("ones", "startvalues"), 
		start.tol=1E-6, optimizer=optimizerSem, objective=objectiveML, cls, ...){
	ord <- function(x) 1 + apply(outer(unique(x), x, "<"), 2, sum)
	is.triangular <- function(X) {
		is.matrix(X) && (nrow(X) == ncol(X)) && 
				(all(0 == X[upper.tri(X)])) || (all(0 == X[lower.tri(X)]))
	} 
	ram <- model
	S <- unclass(S) # in case S is a rawmoment object
	if (nrow(S) > 1 && is.triangular(S)) S <- S + t(S) - diag(diag(S))
	if (!isSymmetric(S)) stop("S must be a square triangular or symmetric matrix")
	if (qr(S)$rank < ncol(S)) warning("S is numerically singular: expect problems")
	if (any(eigen(S, symmetric=TRUE, only.values=TRUE)$values <= 0)) 
		warning("S is not positive-definite: expect problems")
	if ((!is.matrix(ram)) | ncol(ram) != 5 | (!is.numeric(ram)))
		stop ("ram argument must be a 5-column numeric matrix")
	par.size <- if (missing(par.size)) {
				range <- range(diag(S))
				if (range[2]/range[1] > 100) "startvalues" else "ones"
			}
			else match.arg(par.size)
	n <- nrow(S)
	observed <- 1:n
	n.fix <- length(fixed.x)
	if (!is.null(fixed.x)){
		for (i in 1:n.fix){
			for (j in 1:i){
				ram <- rbind(ram, c(2, fixed.x[i], fixed.x[j], 
								0, S[fixed.x[i], fixed.x[j]]))
			}
		}
	}
	m <- max(ram[,c(2,3)])
	missing.variances <- setdiff(1:m, ram[,2][ram[,2] == ram[,3]])
	if (length(missing.variances) > 0) warning(paste("The following variables have no variance or error-variance parameter (double-headed arrow):\n",
						paste(var.names[missing.variances], collapse=", "), 
						"\nThe model is almost surely misspecified; check also for missing covariances.\n"))
	t <- max(ram[,4])
	df <- n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
	if (df < 0) stop(paste("The model has negative degrees of freedom =", df))
	J <- matrix(0, n, m)
	correct <- matrix(2, m, m)
	diag(correct) <- 1
	J[cbind(1:n, observed)]<-1
	par.posn <-  sapply(1:t, function(i) which(ram[,4] == i)[1])
	colnames(ram)<-c("heads", "to", "from", "parameter", "start value")
	rownames(ram)<-rep("",nrow(ram))
	if (length(param.names) > 0) rownames(ram)[par.posn]<-param.names
	fixed <- ram[,4] == 0
	sel.free <- ram[,4]
	sel.free[fixed] <- 1
	one.head <- ram[,1] == 1
	one.free <- which( (!fixed) & one.head )
	two.free <- which( (!fixed) & (!one.head) )
	arrows.1 <- ram[one.head, c(2,3), drop=FALSE]
	arrows.2 <- ram[!one.head, c(2,3), drop=FALSE]
	arrows.2t <- ram[!one.head, c(3,2), drop=FALSE]
	arrows.1.free <- ram[one.free,c(2,3), drop=FALSE]
	arrows.2.free <- ram[two.free,c(2,3), drop=FALSE]
	sel.free.1 <- sel.free[one.free]
	sel.free.2 <- sel.free[two.free]
	unique.free.1 <- unique(sel.free.1)
	unique.free.2 <- unique(sel.free.2)
	rownames(S) <- colnames(S) <- var.names[observed]
	result <- list(var.names=var.names, ram=ram, S=S, J=J, n.fix=n.fix, n=n, N=N, m=m, t=t, raw=raw,
			data=data, semmod=semmod, optimizer=optimizer, objective=objective,
			# remaining values to be supplied after optimization
			coeff=NULL, vcov=NULL, par.posn=NULL, convergence=NULL, iterations=NULL, criterion=NULL, C=NULL, A=NULL, P=NULL,
			adj.obj=NULL, robust.vcov=NULL)
	if (length(param.names)== 0){
		warning("there are no free parameters in the model")
	}
	else {
		if (!is.null(data) && raw && use.means){
			to <- ram[, 2]
			from <- ram[, 3]
			rows <- (from == which(var.names == "Intercept")) & (ram[, 1] == 1) & (ram[, 4] != 0) & (to <= n) & is.na(ram[, 5])
			ram[rows, 5] <- colMeans(data, na.rm=TRUE)[var.names[to[rows]]]
		}
		start <- if (any(is.na(ram[,5][par.posn])))
					start.fn(S, ram, debug=debug, tol=start.tol)
				else ram[,5][par.posn]
		typsize <- if (par.size == "startvalues") abs(start) else rep(1,t)
		model.description <- list(data=data, 
				pattern.number=pattern.number, valid.data.patterns=valid.data.patterns,
				S=S, logdetS=log(det(S)), invS=solve(S), N=N, m=m, n=n, t=t, 
				fixed=fixed, ram=ram, sel.free=sel.free, arrows.1=arrows.1, arrows.1.free=arrows.1.free,
				one.head=one.head, arrows.2t=arrows.2t, arrows.2=arrows.2, arrows.2.free=arrows.2.free, 
				unique.free.1=unique.free.1, unique.free.2=unique.free.2,
				J=J, correct=correct, param.names=param.names, var.names=var.names, observed=observed, raw=raw)
		res <- optimizer(start=start, 
				objective=objective, gradient=analytic.gradient, maxiter=maxiter, debug=debug, par.size=par.size, 
				model.description=model.description, warn=warn, ...)
		ram[par.posn, 5] <- start
		par.code <- paste(var.names[ram[,2]], c("<---", "<-->")[ram[,1]],
				var.names[ram[,3]])
		result$coeff <- res$par
		result$vcov <- res$vcov
		result$par.posn <- par.posn
		result$convergence <- res$convergence
		result$iterations <- res$iterations
		result$criterion <-  res$criterion
		result$C <- res$C
		result$A <- res$A
		result$P <- res$P
		if (!is.na(result$iterations)) if(result$iterations >= maxiter) warning("maximum iterations exceeded")
	}
	if (missing(cls)){
		cls <- gsub("\\.", "", deparse(substitute(objective)))
		cls <- gsub("2", "", cls)
	}
#	if(cls == "objectiveCompiledGLS") 
#			cls <- c(cls, "objectiveGLS")
#	else if(cls == "objectiveCompiledML") 
#			cls <- c(cls, "objectiveML")
	class(result) <- c(cls, "sem")
	if (robust && !is.null(data) && inherits(result, "objectiveML")){
		result$adj.obj <- sbchisq(result, na.omit(data))
		result$robust.vcov <- robustVcov(result, adj.obj=result$adj.obj)
	}
	result
}

vcov.sem <- function(object, robust=FALSE, analytic=inherits(object, "objectiveML") && object$t <= 500, ...) {
	if (robust) return(object$robust.vcov)
	if (!analytic) return(object$vcov)
	if (!inherits(object, "objectiveML")) stop("analytic coefficient covariance matrix unavailable")
	hessian <- function(model){
#		accumulate <- function(A, B, C, D, d) {
#				res <- matrix(0, d^2, d^2)
#			B[1:d, 1:d] %x% A[1:d, 1:d] + matrix(rep(rep(t(C[1:d, 1:d]), 1, each=d), d), d^2, d^2, byrow=TRUE) * matrix(rep(rep((D[1:d, 1:d]), 1, each=d), d), d^2, d^2)
#		}    
		A <- model$A
		P <- model$P
		S <- model$S
		C <- model$C
		J <- model$J
		m <- model$m
		t <- model$t
		I.Ainv <- solve(diag(m) - A) 
		Cinv <- solve(C)    
		AA <- t(I.Ainv) %*% t(J)
		BB <- J %*% I.Ainv %*% P %*% t(I.Ainv)
		CC <- t(I.Ainv) %*% t(J)
		DD <- J %*% I.Ainv
		dF.dBdB <- accumulate(AA %*% Cinv %*% t(AA), t(BB) %*% Cinv %*% BB,
				AA %*% Cinv %*% BB, t(BB) %*% Cinv %*% t(AA), m)                
		dF.dPdP <- accumulate(CC %*% Cinv %*% t(CC), t(DD) %*% Cinv %*% DD,
				CC %*% Cinv %*% DD, t(DD) %*% Cinv %*% t(CC), m)                
		dF.dBdP <- accumulate(AA %*% Cinv %*% t(CC), t(BB) %*% Cinv %*% DD,
				AA %*% Cinv %*% DD, t(BB) %*% Cinv %*% t(CC), m)               
		ram <- model$ram    
		fixed <- ram[, 4] == 0
		sel.free <- ram[, 4]
		sel.free[fixed] <- 0
		one.head <- ram[, 1] == 1
		one.free <- which( (!fixed) & one.head )
		two.free <- which( (!fixed) & (!one.head) )
		two.free.cov <- which((!fixed) & (!one.head) & (ram[, 2] != ram[, 3]))
		arrows.1 <- ram[one.head, c(2, 3), drop=FALSE]
		arrows.2 <- ram[!one.head, c(2, 3), drop=FALSE]
		arrows.2t <- ram[!one.head, c(3, 2), drop=FALSE]
		arrows.1.free <- ram[one.free, c(2, 3), drop=FALSE]
		arrows.2.free <- ram[two.free, c(2, 3), drop=FALSE]
		sel.free.1 <- sel.free[one.free]
		sel.free.2 <- sel.free[two.free]
		unique.free.1 <- unique(sel.free.1)
		unique.free.2 <- unique(sel.free.2)    
		posn.matrix <- matrix(1:(m^2), m, m)    
		posn.free <- c(posn.matrix[arrows.1.free], 
				(m^2) + posn.matrix[arrows.2.free])        
		DBB <- dF.dBdB[posn.matrix[arrows.1.free], 
				posn.matrix[arrows.1.free], drop=FALSE]
		DPP <- dF.dPdP[posn.matrix[arrows.2.free], 
				posn.matrix[arrows.2.free], drop=FALSE]
		DBP <- dF.dBdP[posn.matrix[arrows.1.free],
				posn.matrix[arrows.2.free], drop=FALSE]
#    browser()
		hessian <- rbind( cbind(DBB,    DBP),
				cbind(t(DBP), DPP))    
		n1 <- length(one.free)
		n2 <- length(two.free)
		nn <- rep(c(sqrt(2), sqrt(2)/2), c(n1, n2))
		nn[c(one.free, two.free) %in% two.free.cov] <- sqrt(2)
#    browser()
		hessian <- hessian * outer(nn, nn)
		pars <- ram[, 4][!fixed]
		Z <- outer(1:t, pars, function(x, y) as.numeric(x == y))
		hessian <- Z %*% hessian %*% t(Z)
		par.names <- c(names(one.free), names(two.free))
		par.names <- par.names[par.names != ""]
		rownames(hessian) <- colnames(hessian) <- par.names
		nms <- names(coef(object))
		hessian[nms, nms]
	}
	h <- hessian(object)
	t <- object$t
	N <- object$N
	raw <- object$raw
	param.names <- rownames(h)
	vcov <- matrix(NA, t, t)
	qr.hess <- try(qr(h), silent=TRUE)
	if (class(qr.hess) == "try-error"){
		warning("Could not compute QR decomposition of Hessian.\nOptimization probably did not converge.\n")
	}
	else if (qr.hess$rank < t){
		warning(" singular Hessian: model is probably underidentified.\n")
		which.aliased <- qr.hess$pivot[-(1:qr.hess$rank)]
		attr(vcov, "aliased") <- param.names[which.aliased]
	}
	else {
		vcov <- (2/(N - (!raw))) * solve(h)
		if (any(diag(vcov) < 0)) {
			attr(vcov, "aliased") <- param.names[diag(vcov) < 0]
			warning("Negative parameter variances.\nModel may be underidentified.\n")
		}
	}
	vcov
}

coef.sem <- function(object, standardized=FALSE, ...){
	if (!standardized) return(object$coeff)
	sc <- stdCoef(object, ...)
	names <- as.character(sc[, 1])
	which <- names != " "
	coef <- sc[which, 2]
	names(coef) <- names[which]
	coef
}

# the following auxiliary function is for computing Hessians

accumulate <- function(A, B, C, D, d) {
	B[1:d, 1:d] %x% A[1:d, 1:d] + matrix(rep(rep(t(C[1:d, 1:d]), 1, each=d), d), d^2, d^2, byrow=TRUE) * matrix(rep(rep((D[1:d, 1:d]), 1, each=d), d), d^2, d^2)
}    
