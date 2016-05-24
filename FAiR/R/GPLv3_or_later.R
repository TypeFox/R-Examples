#     This file is part of FAiR, a program to conduct Factor Analysis in R
#     Copyright 2008 Benjamin King Goodrich
#
#     Some portions of this code are Copyright 1995-2007 R Core Development Team,
#     Copyright 2005-2006 Coen Bernaards and Robert Jennrich, Copyright 2008 William
#     Revelle, Copyright 2007 Jens Oehlschlögel, Copyright 2007 Gilles Raiche,
#     and Copyright 2007 John Fox. In each case, the code is licensed under the GNU
#     General Public License version 2 or later.
#
#     FAiR is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     FAiR is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with FAiR.  If not, see <http://www.gnu.org/licenses/>.

## avoid clobbering existing loadings function()
setMethod("loadings", "ANY",
function(x) {
	## This function is copied from stats:::loadings, which is
	## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+
	return(x$loadings)
})

## gets preliminary primary pattern matrix (Lambda) in EFA models
FAiR_get_loadings <- 
function(Psi, S, q) {
	## This function is slightly modified from stats:::factanal.fit.mle, which is 
        ## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+

	sc <- diag(1/sqrt(Psi))
	Sstar <- sc %*% S %*% sc
	E <- eigen(Sstar, symmetric = TRUE)
	L <- E$vectors[, 1:q, drop = FALSE]
	load <- L %*% diag(sqrt(pmax(E$values[1:q] - 1, 0)), q)
	Lambda <- diag(sqrt(Psi)) %*% load
	sorter <- order(colSums(Lambda^2), decreasing = TRUE)
	sorter <- 1:ncol(Lambda) # temporary
	Lambda[,sorter]
}

## Begin methods for restrictions.factanal, which is an EFA model with the same objective
## function as in stats:::factanal and is essentially a Lawley & Maxwell (1971) algorithm
## Note: bfgs_fitS4 and bfgs_helpS4 are not needed because this is not lexical
setMethod("fitS4", signature(restrictions = "restrictions.factanal", 
				 manifest = "manifest.basic"), definition = 
function(par, restrictions, manifest, lower, mapping_rule) {
	## This function is slightly modified from stats:::factanal.fit.mle, which is 
        ## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+

	factors <- restrictions@factors[1]
	sc <- 1/sqrt(par)
# 	Sstar <- sweep(sweep(manifest@cor, 1, sc, "*", FALSE), 2, sc, "*", FALSE)
	Sstar <- manifest@cor * tcrossprod(sc)
	E <- eigen(Sstar, symmetric = TRUE, only.values = TRUE)
	e <- E$values[-(1:factors)]
	e <- sum(log(e) - e) - factors + nrow(manifest@cor)
	return(-e)
})

setMethod("gr_fitS4", signature(restrictions = "restrictions.factanal",
				    manifest = "manifest.basic"), definition =
function(par, restrictions, manifest, helper, lower) {
	## This function is slightly modified from stats:::factanal.fit.mle, which is 
        ## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+
	q <- restrictions@factors[1]

	Psi <- par
	sc <- 1/sqrt(Psi)
# 	Sstar <- sweep(sweep(manifest@cor, 1, sc, "*", FALSE), 2, sc, "*", FALSE)
	Sstar <- manifest@cor * tcrossprod(sc)
	E <- eigen(Sstar, symmetric = TRUE)
	L <- E$vectors[, 1:q, drop = FALSE]
# 	load <- sweep(L, 2, sqrt(pmax(E$values[1:q] - 1, 0)), FUN = "*")
	load <- t(t(L) * sqrt(pmax(E$values[1:q] - 1, 0))) / sc
# 	load <- sweep(load, 1, sc, FUN = "/")
# 	load <- t(t(load) / sc)
	g <- tcrossprod(load) + diag(Psi) - manifest@cor
	out <- diag(g)/Psi^2
	return(out)
})

setMethod("create_start", signature(restrictions = "restrictions.factanal",
					manifest = "manifest.basic"), definition =
function(number, start, restrictions, manifest, reject = FALSE) {
	out <- apply(restrictions@Domains, 1, FUN = function(x) {
			return(runif(number - ncol(start), min = x[1], max = x[2]))
		})
	return(rbind(t(start), out))
})

setMethod("create_FAobject", signature(restrictions = "restrictions.factanal",
					   manifest = "manifest.basic"), definition =
function(restrictions, manifest, opt, call, scores, lower, analytic) {
	S <- cormat(manifest)
	factors <- restrictions@factors[1]
	Lambda <- FAiR_get_loadings(opt$par, S, factors)

	## This part is slightly modified from GPArotation:::echelon, which is Copyright 
	## 2005-2006 Coen Bernaards and Robert Jennrich and licensed under the GPL V2+
	A1 <- Lambda[1:ncol(Lambda), , drop = FALSE]
	B1 <- t(chol(A1 %*% t(A1)))
	Tmat <- solve(A1, B1)
	beta <- matrix(0, nrow(Lambda), ncol(Lambda))
	beta[1:ncol(Lambda), ] <- B1
	beta[-c(1:ncol(beta)), ] <- Lambda[-c(1:ncol(beta)), , drop = FALSE] %*% Tmat
	beta[upper.tri(beta)] <- 0 # already zero aside from numerical error
	signs <- ifelse(colSums(beta) >= 0, 1, -1)
	if(any(signs != 1)) beta <- sweep(beta, 2, signs, "*")
	rownames(beta) <- rownames(S)
	## end echelon part, begin hackish conversion to restrictions.orthonormal
	free <- beta != 0
	beta <- new("parameter.coef", x = beta, free = free, num_free = sum(free))
	new_restrictions <- make_restrictions(manifest = manifest, beta = beta,
						discrepancy = restrictions@discrepancy)

	stuff <- FAiR_restrictions2FA(new_restrictions, manifest, scores)
	sorter <- stuff$sorter
	sorter <- 1:factors # temporary

	uncertainty <- try(FAiR_uncertainty(new_restrictions, manifest, factors, 
						sorter, NULL, TRUE), FALSE)
	if(!is.list(uncertainty)) {
		warning("there was a problem calculating the variance-covariance matrix",
			" of the parameters")
		cols <- new_restrictions@nvars
		uncertainty <- list(vcov = matrix(NA, nrow = cols, ncol = cols),
				Jacobian = matrix(NA_real_, nrow = nrow(S), ncol = cols))
	}
	
	seeds <- matrix(NA_integer_, nrow = 1, ncol = 2)
	trans_mats <- array(diag(factors), c(factors, factors, 3), 
			dimnames = list(NULL, NULL, c("primary", "reference", "T")))
	FAobject <- new("FA.EFA", loadings = stuff$loadings, scale = manifest@sds,
			correlations = stuff$correlations, trans_mats = trans_mats,
			uniquenesses = stuff$uniquenesses, seeds = seeds,
			restrictions = new_restrictions, Jacobian = uncertainty$Jacobian,
			vcov = uncertainty$vcov, scores = stuff$scores, call = call,
			manifest = manifest, optimization = list(extraction = opt),
			Lambda = stuff$loadings[,,"PP"], rotated = FALSE)
	return(FAobject)
})

## prints factor analysis matrices with nice formatting
FAiR_print.loadings <-
function(mat) {
	## This part is slightly modified from stats:::print.loadings, which is 
	## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+
	fx <- format(round(mat, 3))
	names(fx) <- NULL
	nc <- nchar(fx[1], type = "c")
	fx[!mat] <- paste(rep(" ", nc), collapse = "")
	print(fx, quote = FALSE)
	invisible(NULL)
}

## tees up a DAG for plotting or exporting
FAiR_DAG <-
function (x, labels = NULL, cut = 0.0, simple = FALSE, size = c(8, 6), 
	node.font = c("Helvetica",14), edge.font = c("Helvetica", 10), 
	rank.direction = "RL", digits = 1, level = 1, ...) {

	## This function is slightly modified from psych:::fa.graph, which is 
        ## Copyright 2008 William Revelle and licensed under GPL V2+

	if(!require(Rgraphviz)) {
		stop("You need to install the Rgraphviz package to use this function")
	}
	if(level == 1) factors <- loadings(x)
	else           factors <- loadings(x, level = 2)
	rank.direction <- match.arg(rank.direction)
	num.var <- dim(factors)[1]
	num.factors <- dim(factors)[2]
	if (simple) k <- 1 # k redefined below, why?
	else k <- num.factors
	if(level == 1) {
		vars <- rownames(factors)
		fact <- paste("F", 1:num.factors, sep = "")
		graph.shape <- c(rep("box", num.var), rep("ellipse", num.factors))
		graph.rank <- c(rep("sink", num.var), rep("", num.factors))
	}
	else {
		vars <- paste("F", 1:num.var, sep = "")
		fact <- paste("G", 1:num.factors, sep = "")
		graph.shape <- rep("ellipse", num.var + num.factors)
		graph.rank <- rep("", num.var + num.factors)
	}
	clust.graph <- new("graphNEL", nodes = c(vars, fact), edgemode = "directed")
	names(graph.shape) <- nodes(clust.graph)
	names(graph.rank) <- nodes(clust.graph)
	edge.label <- rep("", num.var * k)
	edge.name  <- rep("", num.var * k)
	names(edge.label) <- seq(1:num.var * k)
	l <- factors
	if (num.factors == 1) {
		for (i in 1:num.var) {
			clust.graph <- addEdge(fact[1], vars[i], clust.graph, 1)
			edge.label[i] <- round(factors[i], digits)
			edge.name[i] <- paste(fact[1], "~", vars[i], sep = "")
		}
	}
	else {
		if (simple) {
			m1 <- matrix(apply(t(apply(l, 1, abs)), 1, which.max), ncol = 1)
			for (i in 1:num.var) {
				clust.graph <- addEdge(fact[m1[i]], vars[i],
							clust.graph, 1)
				edge.label[i] <- round(factors[i, m1[i]], digits)
				edge.name[i] <- paste(fact[m1[i]], "~", vars[i], sep = "")
			}
        	}
		else {
			k <- 1
			for (i in 1:num.var) for (f in 1:num.factors) {
				if (abs(factors[i, f]) > cut) {
					clust.graph <- addEdge(fact[f], vars[i], 
								clust.graph, 1)
					edge.label[k] <- round(factors[i, f], digits)
					edge.name[k] <- paste(fact[f], "~", vars[i], 
									sep = "")
					k <- k + 1
				}
			}
		}
        }
	nAttrs <- list()
	eAttrs <- list()
	if (!is.null(labels)) {
		var.labels <- c(labels, fact)
		names(var.labels) <- nodes(clust.graph)
		nAttrs$label <- var.labels
		names(edge.label) <- edge.name
	}
	names(edge.label) <- edge.name
	nAttrs$shape <- graph.shape
	nAttrs$rank <- graph.rank
	eAttrs$label <- edge.label
	attrs <- list(node = list(shape = "ellipse", fixedsize = FALSE), 
			graph = list(rankdir = rank.direction, fontsize = 10,
			bgcolor = "white"))
	obs.var <- subGraph(vars, clust.graph)
	cluster.vars <- subGraph(fact, clust.graph)
	observed <- list(list(graph = obs.var, cluster = TRUE, attrs = c(rank = "")))
	out <- list(x = clust.graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs,
			attrs = attrs, subGList = observed)
	return(out)
}

## makes x into the nearest positive definite correlation matrix
FAiR_nearPD <- 
function(x, corr = TRUE, eig.tol = 1e-06, conv.tol = 1e-07, 
	posd.tol = sqrt(.Machine$double.eps), do2eigen = TRUE, maxit = 100) {
	## This function is slightly modified from Matrix:::nearPD, which is
	## Copyright 2007 Jens Oehlschlögel and is licensed under the GPL V2+
	e <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
	n <- ncol(x)
	if(e[n] > posd.tol * e[1]) {
		attributes(x)$ev <- 1.0
		return(x)
	}
	U <- x
	U[] <- 0
	X <- x
	iter <- 0
	converged <- FALSE
	conv <- Inf
	inorm <- function(x) max(rowSums(abs(x)))
	while (iter < maxit) {
		Y <- X
		T <- Y - U
		e <- eigen(Y, symmetric = TRUE)
		Q <- e$vectors
		d <- e$values
		D <- diag(d)
		p <- d > eig.tol * d[1]
		X <-      Q[, p, drop = FALSE] %*% D[p, p, drop = FALSE] %*% 
			t(Q[, p, drop = FALSE])
		U <- X - T
		X <- (X + t(X))/2
		if (corr) diag(X) <- 1
		conv <- inorm(Y - X) / inorm(Y)
		converged <- (conv <= conv.tol)
		if(converged) break
		iter <- iter + 1
	}
	X <- (X + t(X))/2
	if (do2eigen) {
		e <- eigen(X, symmetric = TRUE)
		d <- e$values
		Eps <- posd.tol * abs(d[1])
		if (d[n] < Eps) {
			d[d < Eps] <- Eps
			Q <- e$vectors
			o.diag <- diag(X)
			X <- Q %*% (d * t(Q))
			D <- sqrt(pmax(Eps, o.diag)/diag(X))
			X[] <- D * X * rep(D, each = n)
		}
	}
	if (corr) X <- cov2cor(X)
	attributes(X)$ev <- min(d)
	return(X)
}

## throws warning if condition is not satisfied
FAiR_warnifnot <-
function(...) {
	## This function is slightly modified from base:::stopifnot, which is 
        ## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+

	n <- length(ll <- list(...))
	if (n == 0) return(invisible())
	mc <- match.call()
	for (i in 1:n) if(!(is.logical(r <- eval(ll[[i]])) && !any(is.na(r)) && all(r))) {
		ch <- deparse(mc[[i + 1]], width.cutoff = 60)
		if (length(ch) > 1) ch <- paste(ch[1], "....")
		warning(paste(ch, " is not ", if(length(r) > 1) "all ", "TRUE", sep = ""),
			call. = FALSE)
	}
}

## parallel analysis for the 4th order moment matrix, to be used with library(nFactors)
FAiR_parallel_4th <-
function (subject = 100, var = 10, rep = 100, cent = 0.05, mcd = FALSE) {
	## This function is slightly modified from nFactors:::parellel, which is
	## Copyright 2007 Gilles Raiche and is licensed under the GPL V2+
	r <- subject
	c <- var
	y <- matrix(c(1:r * c), nrow = r, ncol = c)
	ycor <- matrix(c(1:c * c), nrow = c, ncol = c)
	evpea <- matrix(NA_real_, nrow = rep, ncol = ncol(y) * (ncol(y) + 1))
	leg.txt <- "Pearson"
	foo <- function(x) {z <- tcrossprod(x); z[lower.tri(z, TRUE)]}
	for (k in c(1:rep)) {
		y <- rnorm(y)
		y <- matrix(y, nrow = r, ncol = c)
		if(mcd) Gamma <- FAiR_ADF_robust(y)
		else {
			y <- t(t(y) - colMeans(y))
			z <- t(apply(y, 1, foo))
			Gamma <- cov.wt(z, method = "ML")
		}
		evpea[k,] <- eigen(cov2cor(Gamma), TRUE, TRUE)$values
	}
	SEcentile <- function(sd, n = 100, p = 0.95) {
		return(sd/sqrt(n) * sqrt(p * (1 - p))/dnorm(qnorm(p)))
	}
	sprob <- c(cent)
	mevpea <- sapply(as.data.frame(evpea), mean)
	sevpea <- sapply(as.data.frame(evpea), sd)
	quant <- function(x, sprobs = sprobs) {
		return(as.vector(quantile(x, probs = sprob)))
	}
	qevpea <- sapply(as.data.frame(evpea), quant)
	sqevpea <- sevpea
	sqevpea <- sapply(as.data.frame(sqevpea), SEcentile, n = rep, p = cent)
	result <- list(eigen = data.frame(mevpea, sevpea, qevpea, sqevpea),
			subject = r, variables = c, centile = cent)
	class(result) <- "parallel"
	return(result)
}

## make a square covariance matrix from triangular input
read.triangular <- 
function(file = "", diag = TRUE, names = paste("X", 1:n, sep=""), ...) {
	## The following is slightly modified from sem:::read.moments, which is
	## Copyright 2007 John Fox and is licensed under the GPL V2+

	elements <- scan(file = file, ...)
	m <- length(elements)
	d <- if (diag) 1 else -1
	n <- floor((sqrt(1 + 8 * m) - d) / 2)
	if (m != n * (n + d) / 2) stop("wrong number of elements")
	if (length(names) != n)   stop("wrong number of variable names")
	X <- diag(n)
	X[upper.tri(X, diag=diag)] <- elements
	rownames(X) <- colnames(X) <- names
	X <- X + t(X)
	diag(X) <- diag(X) / 2
	return(X)
}

## simulate on the basis of a FA object
FA2draws <- 
function(object, nsim = 1001, seed = NULL, 
	covariances = FALSE, standardized = TRUE, ...) {
	## This function is slightly modified from stats:::simulate.lm, which is 
        ## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+
	if(!FAiR_is.FA(object)) {
		stop("object must be an object of class FA or inherit from class FA")
	}

	if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
	if(is.null(seed)) RNGstate <- get(".Random.seed", envir = .GlobalEnv)
	else {
		R.seed <- get(".Random.seed", envir = .GlobalEnv)
		set.seed(seed)
		RNGstate <- structure(seed, kind = as.list(RNGkind()))
		on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
	}
	restrictions_copy <- FAiR_copy_restrictions(object@restrictions)
	out <- restrictions2draws(restrictions = restrictions_copy,
				manifest = object@manifest, vcov = vcov(object),
				nsim = nsim, standardized = standardized,
				covariances = covariances)
	return(out)
}

setMethod("simulate", "FA",
function(object, nsim = 1, seed = NULL, standardized = TRUE, ...) {
	## This function is slightly modified from stats:::simulate.lm, which is 
        ## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+
	if(!FAiR_is.FA(object)) {
		stop("object must be an object of class FA or inherit from class FA")
	}

	if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
	if(is.null(seed)) RNGstate <- get(".Random.seed", envir = .GlobalEnv)
	else {
		R.seed <- get(".Random.seed", envir = .GlobalEnv)
		set.seed(seed)
		RNGstate <- structure(seed, kind = as.list(RNGkind()))
		on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
	}
	draws <- FA2draws(object, nsim = nsim, covariances = TRUE, seed = seed,
				standardized = standardized)
	return(draws[[1]])
})

FAiR_parse_formula <-
function(x, data, na.action) {
	## This function is slightly modified from stats:::factanal, which is 
        ## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+

	mt <- terms(x, data = data)
	if (attr(mt, "response")) stop("response not allowed in formula")
	attr(mt, "intercept") <- 0
	mf <- match.call(expand.dots = FALSE)
	names(mf)[names(mf) == "x"] <- "formula"
	mf[[1]] <- as.name("model.frame")
	mf[!names(mf) %in% c("", "formula", "data")] <- NULL
	mf$na.action <- na.action
	mf <- eval.parent(mf)
# 	na.act <- attr(mf, "na.action")
	z <- model.matrix(mt, mf)
	return(z)
}

## comment out the check when pushing to CRAN
FAiR_set_slot <-
function (object, name, check = TRUE, value) {
	## This function is slightly modified from methods:::slot<-, which is 
        ## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+

#       if(check) value <- checkSlotAssignment(object, name, value)
      .Call("FAiR_do_slot_assign", object, name, value, PACKAGE="FAiR")

}
