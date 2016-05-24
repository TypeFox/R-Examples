#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
##   This file is part of the R package rugarch.
##
##   The R package rugarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rugarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

#----------------------------------------------------------------------------------
# linear.hypothesis from pakage car
# included as a direct import as the new r-forge car package does not include
# this at the time of writing which then prevents the rugarch build if using an
# import from statement.
#----------------------------------------------------------------------------------

# last modified 27 Dec 2008 by J. Fox
is.aliased<-function(model){
	!is.null(alias(model)$Complete)
}

.vcov.default <- function(object, ...){
	stop(paste("there is no vcov() method for models of class", 
					paste(class(object), collapse=", ")))
}

has.intercept.matrix <- function (model, ...) {
	"(Intercept)" %in% colnames(model)
}

.hccm <-function(model, type=c("hc3", "hc0", "hc1", "hc2", "hc4"), ...) {
	#last modified 6 Feb 2003 by J. Fox
	if (!is.null(weights(model))) stop("requires unweighted lm")
	type <- match.arg(type)
	sumry <- summary(model, corr = FALSE)
	s2 <- sumry$sigma^2
	V <- sumry$cov.unscaled
	if (type == FALSE) return(s2*V)
	e <- na.omit(residuals(model))
	X <- model.matrix(model)
	df.res <- df.residual(model)
	n <- length(e)
	h <- hat(X)
	p <- ncol(X)
	factor<-switch(type,
			hc0=1,
			hc1=df.res/n,
			hc2=1 - h,
			hc3=(1 - h)^2,
			hc4=(1 - h)^pmin(4, n*h/p))
	V %*% t(X) %*% apply(X, 2, "*", (e^2)/factor) %*% V
}


.makeHypothesis <- function(cnames, hypothesis, rhs = NULL){
	parseTerms <- function(terms){
		component <- gsub("^[-\\ 0-9\\.]+", "", terms)
		component <- gsub(" ", "", component, fixed=TRUE)
		component
	}
	stripchars <- function(x) {
		x <- gsub(" ", "", x, fixed = TRUE)
		x <- gsub("*", "", x, fixed = TRUE)
		x <- gsub("-", "+-", x, fixed = TRUE)
		x <- strsplit(x, "+", fixed = TRUE)[[1]]
		x <- x[x!=""]
		x
	}
	char2num <- function(x) {
		x[x == ""] <- "1"
		x[x == "-"] <- "-1"
		as.numeric(x)
	}
	constants <- function(x) {
		with.coef <- unique(unlist(sapply(cnames, 
								function(y) which(y == parseTerms(x)))))
		if (length(with.coef) > 0) x <- x[-with.coef]
		x <- if (is.null(x)) 0 else sum(as.numeric(x))
		if (any(is.na(x)))
			stop('The hypothesis "', hypothesis, 
					'" is not well formed: contains bad coefficient/variable names.')
		x
	}
	coefvector <- function(x, y) {
		rv <- gsub(" ", "", x, fixed=TRUE) ==
				parseTerms(y)
		if (!any(rv)) return(0)
		if (sum(rv) > 1) stop('The hypothesis "', hypothesis, 
					'" is not well formed.')
		rv <- sum(char2num(unlist(strsplit(y[rv], x, fixed=TRUE))))
		if (is.na(rv))
			stop('The hypothesis "', hypothesis, 
					'" is not well formed: contains non-numeric coefficients.')
		rv
	}
	if(!is.null(rhs)) rhs <- rep(rhs, length.out = length(hypothesis))
	if (length(hypothesis) > 1)
		return(rbind(Recall(cnames, hypothesis[1], rhs[1]), 
						Recall(cnames, hypothesis[-1], rhs[-1])))
	lhs <- strsplit(hypothesis, "=", fixed=TRUE)[[1]]
	if (is.null(rhs)) {
		if (length(lhs) < 2) rhs <- "0"
		else if (length(lhs) == 2) {
			rhs <- lhs[2]
			lhs <- lhs[1]
		} 
		else stop('The hypothesis "', hypothesis, 
					'" is not well formed: contains more than one = sign.')
	} 
	else {
		if (length(lhs) < 2) as.character(rhs)
		else stop('The hypothesis "', hypothesis,
					'" is not well formed: contains a = sign although rhs was specified.')
	}
	lhs <- stripchars(lhs)
	rhs <- stripchars(rhs)
	rval <- sapply(cnames, coefvector, y = lhs) - sapply(cnames, coefvector, y = rhs)
	rval <- c(rval, constants(rhs) - constants(lhs))
	names(rval) <- c(cnames, "*rhs*")  
	rval
}

.printHypothesis <- function(L, rhs, cnames){
	hyp <- rep("", nrow(L))
	for (i in 1:nrow(L)){
		sel <- L[i,] != 0
		h <- L[i, sel]
		h <- ifelse(h < 0, as.character(h), paste("+", h, sep="")) 
		nms <- cnames[sel]
		h <- paste(h, nms) 
		h <- gsub("-1", "-", h)
		h <- gsub("+1", "+", h, fixed=TRUE)
		h <- gsub("-", " - ", h)
		h <- gsub("+", "  + ", h, fixed=TRUE)
		h <- paste(h, collapse="")
		h <- gsub("  ", " ", h, fixed=TRUE)
		h <- sub("^\\ \\+", "", h)
		h <- sub("^\\ ", "", h)
		h <- sub("^-\\ ", "-", h)
		hyp[i] <- paste(h, "=", rhs[i])
	}
	hyp
}

.linear.hypothesis.default <- function(model, hypothesis.matrix, rhs=NULL, 
		test=c("Chisq", "F"), vcov.=NULL, verbose=FALSE, ...){
	df <- df.residual(model)
	if (is.null(df)) df <- Inf ## if no residual df available
	V <- if (is.null(vcov.)) vcov(model)  
			else if (is.function(vcov.)) vcov.(model) else vcov.
	b <- coef(model)
	if (is.null(b)) stop(paste("there is no coef() method for models of class", 
						paste(class(model), collapse=", ")))
	if (is.character(hypothesis.matrix)) {    
		L <- .makeHypothesis(names(b), hypothesis.matrix, rhs)
		if (is.null(dim(L))) L <- t(L)
		rhs <- L[, NCOL(L)]
		L <- L[, -NCOL(L), drop = FALSE]
		rownames(L) <- hypothesis.matrix
	} 
	else {
		L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix) 
				else hypothesis.matrix
		if (is.null(rhs)) rhs <- rep(0, nrow(L))
	}  
	q <- NROW(L)
	if (verbose){
		cat("\nHypothesis matrix:\n")    
		print(L)
		cat("\nRight-hand-side vector:\n")
		print(rhs)
		cat("\nEstimated linear function (hypothesis.matrix %*% coef - rhs)\n")
		print(drop(L %*% b - rhs))
		cat("\n")
	}
	SSH <- as.vector(t(L %*% b - rhs) %*% solve(L %*% V %*% t(L)) %*% (L %*% b - rhs))
	test <- match.arg(test)
	if (!(is.finite(df) && df > 0)) test <- "Chisq"
	name <- try(formula(model), silent = TRUE)
	if (inherits(name, "try-error")) name <- substitute(model)  
	title <- "Linear hypothesis test\n\nHypothesis:"
	topnote <- paste("Model 1: ", paste(deparse(name), collapse = "\n"), "\n",
			"Model 2: restricted model", sep = "")
	note <- if (is.null(vcov.)) "" 
			else "\nNote: Coefficient covariance matrix supplied.\n"
	rval <- matrix(rep(NA, 8), ncol = 4)
	colnames(rval) <- c("Res.Df", "Df", test, paste("Pr(>", test, ")", sep = ""))
	rownames(rval) <- 1:2
	rval[,1] <- c(df, df+q) 
	if (test == "F") {
		f <- SSH/q
		p <- pf(f, q, df, lower.tail = FALSE)
		rval[2,2:4] <- c(-q, f, p)
	} 
	else {
		p <- pchisq(SSH, q, lower.tail = FALSE)
		rval[2,2:4] <- c(-q, SSH, p)
	}
	if (!(is.finite(df) && df > 0)) rval <- rval[,-1]
	structure(as.data.frame(rval), 
			heading = c(title, .printHypothesis(L, rhs, names(b)), "", topnote, note), 
			class = c("anova", "data.frame"))
}

.linear.hypothesis <- function(model, hypothesis.matrix, rhs=NULL,
		test=c("F", "Chisq"), vcov.=NULL, white.adjust=FALSE, ...){
	if (is.aliased(model)) stop("One or more terms aliased in model.")
	test <- match.arg(test)
	if (identical(white.adjust, TRUE)) white.adjust <- "hc3"
	if (is.null(vcov.) && is.character(white.adjust))
		vcov. <- .hccm(model, type = white.adjust)
	rval <- .linear.hypothesis.default(model, hypothesis.matrix, rhs = rhs,
			test = test, vcov. = vcov., ...)
	if (is.null(vcov.)) {
		rval2 <- matrix(rep(NA, 4), ncol = 2)
		colnames(rval2) <- c("RSS", "Sum of Sq")
		SSH <- rval[2,test]
		if (test == "F") SSH <- SSH * abs(rval[2, "Df"])
		df <- rval[1, "Res.Df"]
		error.SS <- deviance(model)
		rval2[,1] <- c(error.SS, error.SS + SSH * error.SS/df)
		rval2[2,2] <- -diff(rval2[,1])
		rval2 <- cbind(rval, rval2)[,c(1, 5, 2, 6, 3, 4)]
		class(rval2) <- c("anova", "data.frame")
		attr(rval2, "heading") <- attr(rval, "heading")
		rval <- rval2
	}
	rval
}
#----------------------------------------------------------------------------------

# kde2d from package MASS
.bandwidth.nrd = function (x) 
{
    r <- quantile(x, c(0.25, 0.75))
    h <- (r[2L] - r[1L])/1.34
    4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
}

.kde2d = function (x, y, h, n = 25, lims = c(range(x), range(y))) 
{
    nx <- length(x)
    if (length(y) != nx) 
        stop("data vectors must be the same length")
    if (any(!is.finite(x)) || any(!is.finite(y))) 
        stop("missing or infinite values in the data are not allowed")
    if (any(!is.finite(lims))) 
        stop("only finite values are allowed in 'lims'")
    gx <- seq.int(lims[1L], lims[2L], length.out = n)
    gy <- seq.int(lims[3L], lims[4L], length.out = n)
    if (missing(h)) 
        h <- c(.bandwidth.nrd(x), .bandwidth.nrd(y))
    h <- h/4
    ax <- outer(gx, x, "-")/h[1L]
    ay <- outer(gy, y, "-")/h[2L]
    z <- matrix(dnorm(ax), n, nx) %*% t(matrix(dnorm(ay), n, 
        nx))/(nx * h[1L] * h[2L])
    return(list(x = gx, y = gy, z = z))
}



#----------------------------------------------------------------------------------
# From HyperbolicDist of D.Scott
#----------------------------------------------------------------------------------
### Test for whole number, with tolerance for representation
### From post by Tony Plate <tplate_at_acm.org>
is.wholenumber <- function(x, tolerance = .Machine$double.eps^0.5){
	if (!is.numeric(x)){
		return(FALSE)
	} else {
		return(isTRUE(all(abs(x - round(x)) < tolerance)))
	}
}

### Calculate recursion for any generalized hyperbolic distribution
### Christine Yang Dong and David Scott from code by Diethelm Wuertz
ghypMean = function(lambda, alpha, beta, delta, mu) 
{	
	gamma <- sqrt(alpha^2 - beta^2)
	mu + delta * beta * besselRatio(delta * gamma, lambda, 1)/gamma
}

momRecursion = function(order = 12, printMatrix = FALSE) {
	## Description:
	##   Computes the moment coefficients recursively
	## Setting Start Values:
	a <- matrix(rep(0, times = order*order), ncol = order)
	a[1, 1] <- 1
	if (order > 1) {
		a[2, 1] <- 1
	}
	## Compute all Cofficients by Recursion:
	if (order > 1) {
		for (d in 2:order) {
			for (l in 2:d) {
				a[d,l] <- a[d - 1,l - 1] + a[d - 1, l]*(2*l + 1 - d)
			}
		}
	}
	rownames(a) <- paste("order=", 1:order, sep = "")
	colnames(a) <- paste("l=", 1:order, sep = "")
	## Print the matrix:
	if (printMatrix) {
		cat("\n")
		print(a)
		cat("\n")
	}
	for (k in 1:order) {
		L <- trunc((k + 1)/2):k
		M <- 2*L - k
	}
	return(list(a = a[order, L], L = L, M = M,
					lmin = trunc((order + 1)/2)))
}

besselRatio <- function(x, nu, orderDiff, useExpScaled = 700){
	if (x > useExpScaled){
		besselK(x, nu + orderDiff, expon.scaled = TRUE)/
				besselK(x, nu, expon.scaled = TRUE)
	}else{
		besselK(x, nu + orderDiff)/besselK(x, nu)
	}
}

## transfer moments about different locations for any distributions
momChangeAbout = function(order = "all", oldMom, oldAbout, newAbout) {
	if (!is.vector(oldMom)){
		stop("A vector of moments must be supplied")
	}
	if (order == "all") {
		## Compute moment of up to length(oldMom) about location new
		mom <- rep(NA,length(oldMom))
		oldMoment <- c(1,oldMom)
		for (i in 1:length(oldMom)) {
			oldMom <- oldMoment[1:(i+1)]
			binomCoeff <- choose(i, 0:i)
			diffPower <- (oldAbout - newAbout)^(i:0)
			mom[i] <- sum(binomCoeff*diffPower*oldMom)
		}
	} else {
		## Check order is within in the right range
		if (length(oldMom) < order) {
			stop("The length of of the vector oldMom must not be less than the
							value of order")
		}
		if (!is.wholenumber(order)){
			stop("Order must be a whole number")
		}
		if ((order < 0)) {
			stop("Order must be positive")
		}
		## Compute moment of a specific order about location new
		oldMom <- c(1,oldMom)
		oldMom <- oldMom[1:(order+1)]
		binomCoeff <- choose(order, 0:order)
		diffPower <- (oldAbout - newAbout)^(order:0)
		mom <- sum(binomCoeff*diffPower*oldMom)
	}
	
	## Return moment
	return(mom)
}


ghypMom = function(order, lambda, alpha, beta, delta, mu, momType = "raw", about = 0) {
	
	## check order is whole number 
	if (!is.wholenumber(order)){
		stop("Order must be a whole number")
	}
	if ((order < 0)) {
		stop("Order must be positive")
	} 
	
	## check momType
	momType <- as.character(momType)
	momType <- tolower(momType)
	if (momType != "raw" & momType != "central" & momType != "mu") {
		stop ("Unrecognised moment type")
	} 
	
	## unpack parameters	
	gm <- sqrt(alpha^2 - beta^2)
	zeta <- delta*gm
	
	if (order == 0) {
		mom <- 1
	} else {         
		## calculate mu moments     
		muMom <- rep (NA,order)
		for (i in 1:order) {
			a <- momRecursion(order = i) 
			coeff <- a$a              
			betaPow <- a$M        
			deltaPow <- 2*a$L
			zetaPow <- a$L
			lengthZetaPow <- length(zetaPow)
			
			## calculate terms and sum
			muM <- coeff*(delta^deltaPow)*(beta^betaPow)*
					sapply(zetaPow, besselRatio, x = zeta, nu = lambda)/(zeta^zetaPow)
			muMom[i] <- sum(muM)   
		}
	}  
	
	if (about != 0) {                    
		mom <- momChangeAbout(order = order, oldMom = muMom, 
				oldAbout = mu, newAbout = about)
	} else {
		if (momType == "mu") {
			mom = muMom[order]
		} else if (momType == "raw") {
			about <- 0
			mom <- momChangeAbout(order = order, oldMom = muMom, 
					oldAbout = mu, newAbout = about)
		} else if (momType == "central") {
			about <- ghypMean(lambda, alpha, beta, delta, mu)
			mom <- momChangeAbout(order = order, oldMom = muMom, 
					oldAbout = mu, newAbout = about)
		}
	}  
	return(mom)
}
#----------------------------------------------------------------------------------

# From Package combinat
combn = function(x, m, fun = NULL, simplify = TRUE, ...)
{
#       DATE WRITTEN: 14 April 1994          LAST REVISED:  10 July 1995
#       AUTHOR:  Scott Chasalow
#
#       DESCRIPTION:
#             Generate all combinations of the elements of x taken m at a time. 
#             If x is a positive integer,  returns all combinations
#             of the elements of seq(x) taken m at a time.
#             If argument "fun" is not null,  applies a function given
#             by the argument to each point.  If simplify is FALSE,  returns 
#             a list; else returns a vector or an array.  "..." are passed 
#             unchanged to function given by argument fun,  if any.
#       REFERENCE:
#             Nijenhuis, A. and Wilf, H.S. (1978) Combinatorial Algorithms for 
#             Computers and Calculators.  NY:  Academic Press.
#       EXAMPLES:
#             > combn(letters[1:4], 2)
#             > combn(10, 5, min)  # minimum value in each combination
#             Different way of encoding points:
#             > combn(c(1,1,1,1,2,2,2,3,3,4), 3, tabulate, nbins = 4)
#             Compute support points and (scaled) probabilities for a
#             Multivariate-Hypergeometric(n = 3, N = c(4,3,2,1)) p.f.:
#             > table.mat(t(combn(c(1,1,1,1,2,2,2,3,3,4), 3, tabulate,nbins=4)))
#
	if(length(m) > 1) {
		warning(paste("Argument m has", length(m), 
						"elements: only the first used"))
		m <- m[1]
	}
	if(m < 0)
		stop("m < 0")
	if(m == 0)
		return(if(simplify) vector(mode(x), 0) else list())
	if(is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == x)
		x <- seq(x)
	n <- length(x)
	if(n < m)
		stop("n < m")
	e <- 0
	h <- m
	a <- 1:m
	nofun <- is.null(fun)
	count <- nCm(n, m, 0.10000000000000002)
	out <- vector("list", count)
	out[[1]] <- if(nofun) x[a] else fun(x[a], ...)
	if(simplify) {
		dim.use <- NULL
		if(nofun) {
			if(count > 1)
				dim.use <- c(m, count)
		}
		else {
			out1 <- out[[1]]
			d <- dim(out1)
			if(count > 1) {
				if(length(d) > 1)
					dim.use <- c(d, count)
				else if(length(out1) > 1)
					dim.use <- c(length(out1), count)
			}
			else if(length(d) > 1)
				dim.use <- d
		}
	}
	i <- 2
	nmmp1 <- n - m + 1
	mp1 <- m + 1
	while(a[1] != nmmp1) {
		if(e < n - h) {
			h <- 1
			e <- a[m]
			j <- 1
		}
		else {
			h <- h + 1
			e <- a[mp1 - h]
			j <- 1:h
		}
		a[m - h + j] <- e + j
		out[[i]] <- if(nofun) x[a] else fun(x[a], ...)
		i <- i + 1
	}
	if(simplify) {
		if(is.null(dim.use))
			out <- unlist(out)
		else out <- array(unlist(out), dim.use)
	}
	out
}

nCm = function(n, m, tol = 9.9999999999999984e-009)
{
#  DATE WRITTEN:  7 June 1995               LAST REVISED:  10 July 1995
#  AUTHOR:  Scott Chasalow
#
#  DESCRIPTION: 
#        Compute the binomial coefficient ("n choose m"),  where n is any 
#        real number and m is any integer.  Arguments n and m may be vectors;
#        they will be replicated as necessary to have the same length.
#
#        Argument tol controls rounding of results to integers.  If the
#        difference between a value and its nearest integer is less than tol,  
#        the value returned will be rounded to its nearest integer.  To turn
#        off rounding, use tol = 0.  Values of tol greater than the default
#        should be used only with great caution, unless you are certain only
#        integer values should be returned.
#
#  REFERENCE: 
#        Feller (1968) An Introduction to Probability Theory and Its 
#        Applications, Volume I, 3rd Edition, pp 50, 63.
#
	len <- max(length(n), length(m))
	out <- numeric(len)
	n <- rep(n, length = len)
	m <- rep(m, length = len)
	mint <- (trunc(m) == m)
	out[!mint] <- NA
	out[m == 0] <- 1	# out[mint & (m < 0 | (m > 0 & n == 0))] <-  0
	whichm <- (mint & m > 0)
	whichn <- (n < 0)
	which <- (whichm & whichn)
	if(any(which)) {
		nnow <- n[which]
		mnow <- m[which]
		out[which] <- ((-1)^mnow) * Recall(mnow - nnow - 1, mnow)
	}
	whichn <- (n > 0)
	nint <- (trunc(n) == n)
	which <- (whichm & whichn & !nint & n < m)
	if(any(which)) {
		nnow <- n[which]
		mnow <- m[which]
		foo <- function(j, nn, mm)
		{
			n <- nn[j]
			m <- mm[j]
			iseq <- seq(n - m + 1, n)
			negs <- sum(iseq < 0)
			((-1)^negs) * exp(sum(log(abs(iseq))) - lgamma(m + 1))
		}
		out[which] <- unlist(lapply(seq(along = nnow), foo, nn = nnow, 
						mm = mnow))
	}
	which <- (whichm & whichn & n >= m)
	nnow <- n[which]
	mnow <- m[which]
	out[which] <- exp(lgamma(nnow + 1) - lgamma(mnow + 1) - lgamma(nnow - 
							mnow + 1))
	nna <- !is.na(out)
	outnow <- out[nna]
	rout <- round(outnow)
	smalldif <- abs(rout - outnow) < tol
	outnow[smalldif] <- rout[smalldif]
	out[nna] <- outnow
	out
}

#################################################################
# some functions from the shape package of Soetaert/ 
# Only needed for a couple of plots
# equivalent to calling shape::femmecol(100)

femme100 = function(){
	ans = c("#00008F","#000098","#0000A1","#0000AA","#0000B3","#0000BC","#0000C5",
"#0000CE","#0000D7","#0000E0","#0000EA","#0000F3","#0000FC","#0006FF",
"#0011FF","#001BFF","#0025FF","#0030FF","#003AFF","#0044FF","#004FFF",
"#0059FF","#0063FF","#006DFF","#0078FF","#0082FF","#008CFF","#0097FF",
"#00A1FF","#00ABFF","#00B6FF","#00C0FF","#00CAFF","#00D5FF","#00DFFF",
"#00E9FF","#00F3FF","#00FEFF","#09FFF6","#13FFEC","#1EFFE1","#28FFD7",
"#32FFCD","#3DFFC2","#47FFB8","#51FFAE","#5BFFA4","#66FF99","#70FF8F",
"#7AFF85","#85FF7A","#8FFF70","#99FF66","#A4FF5B","#AEFF51","#B8FF47",
"#C2FF3D","#CDFF32","#D7FF28","#E1FF1E","#ECFF13","#F6FF09","#FFFE00",
"#FFF300","#FFE900","#FFDF00","#FFD400","#FFCA00","#FFC000","#FFB600",
"#FFAB00","#FFA100","#FF9700","#FF8C00","#FF8200","#FF7800","#FF6D00",
"#FF6300","#FF5900","#FF4F00","#FF4400","#FF3A00","#FF3000","#FF2500",
"#FF1B00","#FF1100","#FF0600","#FB0000","#F10000","#E70000","#DC0000",
"#D20000","#C80000","#BE0000","#B30000","#A90000","#9F0000","#950000",
"#8A0000","#800000")
return(ans)
}

#  function:    shape::drapecol
.drapecol = function (A, col = femme100(), NAcol = "white", lim = NULL) 
{
	nr <- nrow(A)
	nc <- ncol(A)
	ncol <- length(col)
	AA <- 0.25 * (A[1:(nr - 1), 1:(nc - 1)] + A[1:(nr - 1), 2:nc] + 
				A[2:nr, 1:(nc - 1)] + A[2:nr, 2:nc])
	if (!is.null(lim)) 
		Ar <- lim
	else Ar <- range(AA, na.rm = TRUE)
	rn <- Ar[2] - Ar[1]
	ifelse(rn != 0, drape <- col[1 + trunc((AA - Ar[1])/rn * 
									(ncol - 1))], drape <- rep(col[1], ncol))
	drape[is.na(drape)] <- NAcol
	return(drape)
}
