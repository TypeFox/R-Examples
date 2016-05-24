# these functions originally by Jarrett Byrnes

# last modified 2015-04-30 by J. Fox

robustVcov <- function(sem.obj, adj.obj, data.obj, use.fit=FALSE, use.ginv=FALSE){
	if (missing(adj.obj) && missing(data.obj)) 
		stop("Need a data or sbchisq object")
	if (missing(adj.obj)) {
		adj.obj <- sbchisq(sem.obj, data.obj, use.fit = use.fit)
	}
	ncases <- sem.obj$N
	hes <- semHessian(adj.obj$w_mat, adj.obj$p_deriv_mat)
	info_m <- try(solve(hes), silent = TRUE)
	if (class(info_m) == "try-error" && use.ginv == TRUE) {
		info_m <- ginv(hes)
		ginvFlag <- TRUE
	}
	acov <- info_m %*% (t(adj.obj$p_deriv_mat) %*% adj.obj$w_mat %*% 
			adj.obj$w_adf %*% adj.obj$w_mat %*% adj.obj$p_deriv_mat) %*% 
		info_m
	acov <- acov/(ncases - 1)
	rownames(acov) <- colnames(acov) <- colnames(adj.obj$p_deriv_mat)
	acov
}

sbchisq <- function(sem.obj, sem.data, adj=1e-04, use.fit=FALSE, use.ginv=FALSE){
	props <- semProps(sem.obj)
	chisq <- props$chisq
	df <- props$df
	w_adf <- adfWmat(sem.data)
	w_mat <- mlWmat(sem.obj, use.fit = use.fit)
	p_deriv_mat <- deltaMatrix(sem.obj)
	ginvFlag <- FALSE
	invMat <- try(solve(t(p_deriv_mat) %*% w_mat %*% p_deriv_mat), 
		silent = TRUE)
	if (class(invMat) == "try-error" && use.ginv == TRUE) {
		invMat <- ginv(t(p_deriv_mat) %*% w_mat %*% p_deriv_mat)
		ginvFlag <- TRUE
	}
	res_u <- w_mat - (w_mat %*% p_deriv_mat %*% invMat %*% t(p_deriv_mat) %*% 
			w_mat)
	ug <- res_u %*% w_adf
	scale_stat <- sum(diag(ug))/df
	chisq.scaled <- chisq/scale_stat
	p.old <- pchisq(chisq, df, lower.tail=FALSE)
	p <- pchisq(chisq.scaled, df, lower.tail=FALSE)
	ret <- list(chisq = chisq, t = sem.obj$t, df = df, p.old = p.old, 
		c = scale_stat, chisq.scaled = chisq.scaled, p = p, w_mat = w_mat, 
		p_deriv_mat = p_deriv_mat, w_adf = w_adf, res_u = res_u, 
		N = length(sem.data[, 1]), ginvFlag = ginvFlag)
	class(ret) <- "adjchisq"
	return(ret)
}

adfWmat <- function(rawdata) {
	names <- colnames(rawdata)
	n <- nrow(rawdata)
	n.col <- ncol(rawdata)
	nc.star <- n.col * (n.col + 1)/2
	nc2 <- n.col^2
	z <- scale(rawdata, center=TRUE, scale=FALSE)
	sc <- vector(nc.star, mode="list")
	outnames <- vector(nc.star, mode="character")
	ind <- combn(n.col + 1, 2)
	ind[2, ] <- ind[2, ] - 1
	for (q in 1:nc.star) {
		i <- ind[1, q]
		j <- ind[2, q]
		outnames[q] <- paste(names[i], names[j], sep="_")
		sc[[q]] <- z[, i] * z[, j]
	}
	names(sc) <- outnames
	adf_mat <- var(data.frame(sc)) * (n - 1)/n
	return(adf_mat)
}

mlWmat <- function(sem.obj, use.fit=FALSE) {
	p <- nrow(sem.obj$C) 
	if (use.fit) {
		An <- sem.obj$C
	}
	else {
		An <- sem.obj$S
	}
	Dp <- Ktrans(p)
	An.inv <- solve(An)
	w_mat <- 0.5 * t(Dp) %*% kronecker(An.inv, An.inv) %*% Dp
	rownames(w_mat) <- vech(matrixNames(sem.obj$C))
	colnames(w_mat) <- rownames(w_mat)
	return(w_mat)
}

deltaMatrix <- function(sem.object, adj=1e-04) {
	p.star <- sem.object$n * (sem.object$n + 1)/2
	pars <- names(sem.object$coeff)
	nparms <- length(pars)
	delta.mat <- matrix(0, nparms, p.star)
	rownames(delta.mat) <- pars
	vars <- sem.object$var.names
	J <- sem.object$J
	m <- sem.object$m
	for (j in 1:nparms) {
		A <- sem.object$A
		P <- sem.object$P
		i <- which(rownames(sem.object$ram) == pars[j])
		from <- sem.object$ram[i, 2]
		to <- sem.object$ram[i, 3]
		path_type <- sem.object$ram[i, 1][1]
		if (path_type == 1){
			AA <- A[cbind(from, to)][1]
			adjust <- abs(AA) * adj
			A[cbind(from, to)] <- AA + adjust
		}
		else {
			PP <- P[cbind(to, from)][1]
			adjust <- PP * adj
			P[cbind(from, to)] <- P[cbind(to, from)] <- PP + adjust
		}
		I.Ainv <- solve(diag(m) - A)
		C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
		delta.mat[j, ] <- (vech(C) - vech(sem.object$C))/adjust
	}
	t(delta.mat)
}

semProps <- function(object){
	N <- object$N
	n <- object$n
	t <- object$t
	n.fix <- object$n.fix
	list(N=N, n=n, t=t, n.fix=n.fix, df=n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2, 
			chisq=object$criterion*(N - (!object$raw)))
}

Ktrans <- function(num.vars){
	duplication.matrix(num.vars)
}

matrixNames <- function(mat, sep="_"){
	outer(rownames(mat), colnames(mat), function(x, y) paste(x, y, sep=sep))
}

semHessian <- function(w.mat, delta.mat){
	ret <- t(delta.mat) %*% w.mat %*% delta.mat
	return(ret)
}
