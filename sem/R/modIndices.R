# last modified 2015-06-09 by J. Fox

mod.indices <- function(...){
	.Deprecated("modIndices", package="sem")
	modIndices(...)
}

modIndices <- function(model, ...){
    UseMethod("modIndices")
    }

# incorporates corrections by Michael Culbertson:

modIndices.objectiveML <- function(model, duplicated, deviance=NULL, ...){
	Duplicated <- function(x){
		X <- outer(x, x, FUN="==")
		diag(X) <- NA
		apply(X, 2, any, na.rm=TRUE)
	}
	if (is.null(deviance)) deviance <- deviance(model)
	A <- model$A
	P <- model$P
	S <- model$S
	C <- model$C
	J <- model$J
	m <- model$m
	t <- model$t
	NM <- model$N - (!model$raw)
	vars <- model$var.names    
	I.Ainv <- solve(diag(m) - A) 
	Cinv <- solve(C)    
	AA <- t(I.Ainv) %*% t(J)
	BB <- J %*% I.Ainv %*% P %*% t(I.Ainv)
	correct <- matrix(2, m, m)
	diag(correct) <- 1
	grad.P <- correct * AA %*% Cinv %*% (C - S) %*% Cinv %*% t(AA)
	grad.A <- correct * AA %*% Cinv %*% (C - S) %*% Cinv %*% BB 
	grad <- c(grad.A, grad.P) * NM
	dF.dBdB <- accumulate(AA %*% Cinv %*% t(AA), t(BB) %*% Cinv %*% BB,
			AA %*% Cinv %*% BB, t(BB) %*% Cinv %*% t(AA), m)                
	dF.dPdP <- accumulate(AA %*% Cinv %*% t(AA), AA %*% Cinv %*% t(AA),
			AA %*% Cinv %*% t(AA), AA %*% Cinv %*% t(AA), m)                
	dF.dBdP <- accumulate(AA %*% Cinv %*% t(AA), t(BB) %*% Cinv %*% t(AA),
			AA %*% Cinv %*% t(AA), t(BB) %*% Cinv %*% t(AA), m)    
	correct.BB <- correct.PP <- correct.BP <- matrix(1, m^2, m^2)
	d0 <- as.vector(diag(m) ==0 )
	d1 <- as.vector(diag(m) == 1)
	correct.BB[d0, d0] <- 2
	correct.PP[d1, d1] <- 0.5
	correct.PP[d0, d0] <- 2
	correct.BP[d0, d0] <- 2
	Hessian <- NM*rbind(cbind(dF.dBdB * correct.BB,    dF.dBdP * correct.BP),
			cbind(t(dF.dBdP * correct.BP), dF.dPdP * correct.PP))
	ram <- model$ram   
	fixed <- ram[,4] == 0
	one.head <- ram[,1] == 1
	one.free <- which( (!fixed) & one.head )
	two.free <- which( (!fixed) & (!one.head) )
	arrows.1.free <- ram[one.free,c(2,3)]
	arrows.2.free <- ram[two.free,c(2,3)]
	duplicated <- if (missing(duplicated)) (Duplicated(ram[, 4]) & (!fixed)) else (duplicated & (!fixed))
	one.free.duplicated <- which( (!fixed) & one.head & duplicated)
	two.free.duplicated <- which( (!fixed) & (!one.head) & duplicated)
	arrows.1.free.duplicated <- ram[one.free.duplicated, c(2,3)]
	arrows.2.free.duplicated <- ram[two.free.duplicated, c(2,3)]
	posn.matrix <- matrix(1:(m^2), m, m)
	posn.free <- c(posn.matrix[arrows.1.free], 
			(m^2) + posn.matrix[arrows.2.free]) 
	posn.free.duplicated <- c(posn.matrix[arrows.1.free.duplicated], 
			(m^2) + posn.matrix[arrows.2.free.duplicated])
	hessian <- Hessian[posn.free, posn.free]
	par.no <- ram[ram[, 4] > 0, 4]
	pars <- ram[, 4][!fixed]
	Z <- outer(1:t, pars, function(x, y) as.numeric(x == y))
	hessian <- Z %*% hessian %*% t(Z)
	E.inv <- solve(hessian)                      
	par.change <- mod.indices <- rep(0, 2*(m^2))       
	n.posn.free.1 <- length(posn.free) + 1
	for (i in 1:(2*(m^2))) {
		if ((! i %in% posn.free.duplicated) && (i %in% posn.free || qr(as.matrix(Hessian[c(posn.free, i), c(posn.free, i)]))$rank < n.posn.free.1)){  
			par.change[i] <- mod.indices[i] <- NA  
		} 
		else {
			k <- Hessian[i, i]
			d <- Hessian[i, posn.free]
			d <- sapply(1:t, function(i) sum(d[which(par.no==i)]))
			par.change[i] <- -grad[i]/ (k - d %*% E.inv %*% d)
			mod.indices[i] <- -0.5 * grad[i] * par.change[i]
		}
	}
	par.change[mod.indices > 1.1*deviance] <- NA
	mod.indices[mod.indices > 1.1*deviance] <- NA
	mod.A <- matrix(mod.indices[1:(m^2)], m, m)
	mod.P <- matrix(mod.indices[-(1:(m^2))], m, m)
	par.A <- matrix(par.change[1:(m^2)], m, m)
	par.P <- matrix(par.change[-(1:(m^2))], m, m)
	diag(mod.A) <- diag(par.A) <- NA
	rownames(mod.A) <- colnames(mod.A) <- vars
	rownames(mod.P) <- colnames(mod.P) <- vars
	rownames(par.A) <- colnames(par.A) <- vars
	rownames(par.P) <- colnames(par.P) <- vars
	result <- list(mod.A=mod.A, mod.P=mod.P, par.A=par.A, par.P=par.P)
	class(result) <- "modIndices"
	result
}

summary.modIndices <- function(object, round=2, 
    print.matrices=c('both', 'par.change', 'mod.indices'), ...) {
    print.matrices <- match.arg(print.matrices)
    if (print.matrices != "mod.indices"){ 
        cat("\n Parameter change: A matrix (regression coefficients, row <- column)\n")
        print(object$par.A)
        }
    if (print.matrices != "par.change"){ 
        cat("\n Modification indices: A matrix (regression coefficients, row <- column)\n")
        print(round(object$mod.A, round))
        }
    if (print.matrices != "mod.indices"){ 
        cat("\n Parameter change: P matrix (variances/covariances)\n")
        print(object$par.P)
        }
    if (print.matrices != "par.change"){ 
        cat("\n Modification indices: P matrix (variances/covariances)\n")
        print(round(object$mod.P, round))
        }
    invisible(NULL)
    }

print.modIndices <- function(x, n.largest=5, ...){
    cat("\n", n.largest, "largest modification indices, A matrix (regression coefficients):\n")
    mod.A <- as.vector(x$mod.A)
    names <- rownames(x$mod.A)
    names(mod.A) <- outer(names, names, paste, sep="<-")
    print(rev(sort(mod.A))[1:n.largest])
    cat("\n ", n.largest, "largest modification indices, P matrix (variances/covariances):\n")
    mod.P <- as.vector(x$mod.P)
    names(mod.P) <- outer(names, names, paste, sep="<->")
    print(rev(sort(mod.P[lower.tri(x$mod.P, diag=TRUE)]))[1:n.largest])
    }

modIndices.msemObjectiveML <- function(model, ...){
	deviance <- deviance(model)
	G <- length(model$groups)
	t <- model$t  
	N <- model$N
	raw <- model$raw
	wts <- (N - !raw)/(sum(N) - !raw*G)
	hessian <- matrix(0, t, t)
	rownames(hessian) <- colnames(hessian) <- model$param.names
	result <- vector(G, mode="list")
	names(result) <- paste(model$group, model$groups, sep=": ")
	all.pars <- unlist(lapply(model$ram, function(ram) ram[, 4]))
	unique.pars <- unique(all.pars)
	pars.count <- table(all.pars)
	duplicated.pars <- as.numeric(names(pars.count)[pars.count > 1])
	for (g in 1:G){
		ram <- model$ram[[g]]
		parameters <- ram[, 4]
		duplicated <- parameters %in% duplicated.pars
		unique.pars <- unique(parameters[parameters != 0])
		par.posn <- sapply(unique.pars, function(x) which(x == parameters)[1])
		unique.posn <- which(parameters %in% unique.pars)
		rownames(ram)[unique.posn] <- unique(model$param.names[ram[, 4]])
		ram[unique.posn, 4] <- unlist(apply(outer(unique.pars, parameters, "=="), 2, which))
		mod.g <- list(var.names=model$var.names[[g]], ram=ram,J=model$J[[g]], n.fix=model$n.fix, n=model$n[[g]], 
				N=model$N[g], m=model$m[[g]], t=length(unique.pars),
				coeff=model$coeff[parameters],  criterion=model$group.criteria[g],  S=model$S[[g]], raw=model$raw,
				C=model$C[[g]], A=model$A[[g]], P=model$P[[g]])
		class(mod.g) <- c("objectiveML", "sem")
		result[[g]] <- modIndices(mod.g, duplicated=duplicated, deviance=deviance, ...)
	}
	class(result) <- "msemModIndices"
	result
}

print.msemModIndices <- function(x, ...){
	G <- length(x)
	groups <- names(x)
	for (g in 1:G) {
		cat("\n\n", groups[g], "\n")
		print(x[[g]], ...)
	}
}

summary.msemModIndices <- function(object, ...){
	G <- length(object)
	groups <- names(object)
	for (g in 1:G) {
		cat("\n\n", groups[g], "\n")
		print(summary(object[[g]], ...))
	}
}