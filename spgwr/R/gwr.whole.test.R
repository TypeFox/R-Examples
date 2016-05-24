# Copyright 2001-2004 Roger Bivand and Danlin Yu
# 

LMZ.F1GWR.test <- function(x) {
	if(class(x) != "gwr") stop("not a gwr object")
	if (!x$hatmatrix) stop("Fit GWR model with argument hatmatrix=TRUE")
	n <- ncol(x$lhat)
	R <- t(diag(n) - x$lhat) %*% (diag(n) - x$lhat)
	RSSg <- as.vector(t(x$lm$y) %*% R %*% x$lm$y)
	DFg1 <- sum(diag(R))
	DFg2 <- sum(diag(R %*% R))
	RSSo <- sum((x$lm$residuals * sqrt(x$lm$weights))^2)
	DFo <- x$lm$df.residual
	statistic <- (RSSg/DFg1)/(RSSo/DFo)
	names(statistic) <- "F"
	parameter <- c((DFg1^2)/DFg2, DFo)
	names(parameter) <- c("df1", "df2")
	ests <- c(RSSo, RSSg)
	names(ests) <- c("SS OLS residuals", "SS GWR residuals")
	pv <- pf(statistic, (DFg1^2)/DFg2, DFo, lower.tail=TRUE)
	res <- list(statistic=statistic, parameter=parameter, p.value=pv,
		method="Leung et al. (2000) F(1) test", alternative="less", 
		data.name=deparse(substitute(x)), estimates=ests)
	class(res) <- "htest"
	res
}

LMZ.F2GWR.test <- function(x) {
	if(class(x) != "gwr") stop("not a gwr object")
	if (!x$hatmatrix) stop("Fit GWR model with argument hatmatrix=TRUE")
	n <- ncol(x$lhat)
	R <- t(diag(n) - x$lhat) %*% (diag(n) - x$lhat)
	RSSg <- as.vector(t(x$lm$y) %*% R %*% x$lm$y)
	DFg1 <- sum(diag(R))
	DFg2 <- sum(diag(R %*% R))
	RSSo <- sum((x$lm$residuals * sqrt(x$lm$weights))^2)
	DFo <- x$lm$df.residual
	statistic <- ((RSSo-RSSg)/(DFo-DFg1))/(RSSo/DFo)
	names(statistic) <- "F"
	parameter <- c(((DFo-DFg1)^2)/(DFo-2*DFg1+DFg2), DFo)
	names(parameter) <- c("df1", "df2")
	ests <- c(RSSo, (RSSo - RSSg))
	names(ests) <- c("SS OLS residuals", "SS GWR improvement")
	pv <- pf(statistic, ((DFo-DFg1)^2)/(DFo-2*DFg1+DFg2), DFo, 
		lower.tail=FALSE)
	res <- list(statistic=statistic, parameter=parameter, p.value=pv,
		method="Leung et al. (2000) F(2) test", alternative="greater", 
		data.name=deparse(substitute(x)), estimates=ests)
	class(res) <- "htest"
	res
}

BFC99.gwr.test <- function(x) {
	if(class(x) != "gwr") stop("not a gwr object")
	if (!x$hatmatrix) stop("Fit GWR model with argument hatmatrix=TRUE")
	n <- ncol(x$lhat)
	R <- t(diag(n) - x$lhat) %*% (diag(n) - x$lhat)
        trS <- sum(diag(x$lhat))
	RSSg <- as.vector(t(x$lm$y) %*% R %*% x$lm$y)
	DFg1 <- sum(diag(R))
	DFg2 <- sum(diag(R %*% R))
        e <- x$lm$residuals * sqrt(x$lm$weights)
	RSSo <- sum(e^2)
	DFo <- x$lm$df.residual
        statistic <- ((RSSo - RSSg)/(DFo-DFg1))/(RSSg/DFg1)
	names(statistic) <- "F"
#	parameter <- c(((DFo-DFg1)^2)/(DFo-2*DFg1+DFg2), (DFg1^2)/DFg2)
# df1 modified to (tr(R_0 - R_1)^2) / tr((R_0 - R_1)^2) as 
# (DFo-2*DFg1+DFg2) != sum(((1-hatvalues(x$lm)) - diag(R))^2)
# but (DFo-DFg1)^2 == (sum((1-hatvalues(x$lm)) - diag(R)))^2
# reported by Deny Kurniawan 070804
	k <- as.integer(x$lm$qr$rank)
	do.coef <- FALSE
# 120822 RSB avoid cross-package calls
#	res <- .Fortran("lminfl", x$lm$qr$qr, n, n, k, as.integer(do.coef), 
#         x$lm$qr$qraux, wt.res = e, hat = double(n), coefficients = 
#	  if (do.coef) matrix(0, n, k) else double(0), sigma = double(n), 
#          tol = 10 * .Machine$double.eps, DUP = FALSE, 
#          PACKAGE = "stats")[c("hat", "coefficients", "sigma", "wt.res")]
	
	hatvalues <- hatvalues(lm(x$lm$y ~ x$lm$x, weights=x$lm$weights))
	parameter <- c(((DFo-DFg1)^2)/(sum(((1-hatvalues) - diag(R))^2)),
		(DFg1^2)/DFg2)
	names(parameter) <- c("df1", "df2")
	ests <- c((RSSo - RSSg), RSSg)
	names(ests) <- c("SS GWR improvement", "SS GWR residuals")
	pv <- pf(statistic, parameter[1], parameter[2], lower.tail=FALSE)
	res <- list(statistic=statistic, parameter=parameter, p.value=pv,
		method="Brunsdon, Fotheringham & Charlton (1999) ANOVA", 
		data.name=deparse(substitute(x)), estimates=ests,
		alternative="greater")
	class(res) <- "htest"
	res
}

anova.gwr <- function(object, ..., approx=FALSE) {
    if(class(object) != "gwr") stop("not a gwr object")
    if (!object$hatmatrix) stop("Fit GWR model with argument hatmatrix=TRUE")
    n <- ncol(object$lhat)
    R <- t(diag(n) - object$lhat) %*% (diag(n) - object$lhat)
    trS <- sum(diag(object$lhat))
    RSSg <- as.vector(t(object$lm$y) %*% R %*% object$lm$y)
    DFg1 <- sum(diag(R))
    e <- object$lm$residuals * sqrt(object$lm$weights)
    RSSo <- sum(e^2)
    DFo <- object$lm$df.residual
    if (approx) {
        ss <- c(RSSo, (RSSo - RSSg), RSSg)
        df <- c(n-DFo, trS - (n-DFo), n - trS)
    } else {
        ss <- c(RSSo, (RSSo - RSSg), RSSg)
        df <- c(n-DFo, n - DFg1 - (n - DFo), DFg1)
    }
    ms <- ss/df
    ms[1] <- NA
    f <- numeric(3)
    f <- NA
    f[3] <- ms[2]/ms[3]
    table <- data.frame(df, ss, ms, f)
    tlabels <- c("OLS Residuals", "GWR Improvement", "GWR Residuals")
    dimnames(table) <- list(tlabels, c("Df", 
        "Sum Sq", "Mean Sq", "F value"))
    structure(table, heading = paste("Analysis of Variance Table", 
        ifelse(approx, "\napproximate degrees of freedom (only tr(S))", "")),
        class = c("anova", "data.frame"))
}

# implemented from 2002 book pp. 91-2

BFC02.gwr.test <- function(x, approx=FALSE) {
  if(class(x) != "gwr") stop("not a gwr object")
  if (!x$hatmatrix) stop("Fit GWR model with argument hatmatrix=TRUE")
  n <- ncol(x$lhat)
  R <- t(diag(n) - x$lhat) %*% (diag(n) - x$lhat)
  RSSg <- as.vector(t(x$lm$y) %*% R %*% x$lm$y)
  DFg1 <- sum(diag(x$lhat))
  DFg2 <- sum(diag(t(x$lhat) %*% x$lhat))
  if (approx) {
    DFg <- n - DFg1
  } else {
    DFg <- n - (2*DFg1 - DFg2)
  }
  RSSo <- sum((x$lm$residuals * sqrt(x$lm$weights))^2)
  DFo <- x$lm$df.residual
  statistic <- RSSo/RSSg
  names(statistic) <- "F"
  parameter <- c(DFo, DFg)
  names(parameter) <- c("df1", "df2")
  ests <- c(RSSo, RSSg)
  names(ests) <- c("SS OLS residuals", "SS GWR residuals")
  pv <- pf(statistic, parameter[1], parameter[2], lower.tail=FALSE)
  res <- list(statistic=statistic, parameter=parameter, p.value=pv,
    method=paste("Brunsdon, Fotheringham & Charlton (2002, pp. 91-2) ANOVA",
    ifelse(approx, "(approximate degrees of freedom - only tr(S))", "")), 
    data.name=deparse(substitute(x)), estimates=ests,
    alternative="greater")
  class(res) <- "htest"
  res
}

#This section tests the coefficient non-stationarity using Leung et al
#(2000)'s F3 test.
#The Bks are already available through the GWR routine, namely,
#from gwr.b[,1] to gwr.b[,k], with k is the number of independent
#variables (include the intercept)
#The central idea of testing the stationarity for each independent
#variable's coefficient is to:
#Construct an appropriate statistic which can reflect the 
#spatial variation. The variance of the estimated coefficients (betas)
#of one independent varialbe (over the space) is a natural choice and
#is used in the GWR book, Leung et al. (2000)'s paper etc., and will
#be coded here.

#The test follows: if the coefficients of one independent variable is
#spatial stationary, that is to say, all the betas for that independent
#variable across the space should be equal. Statistically, the Null
#hypothesis for this is that the expected values for all the betas is
#the same. Detailed analytical deduction of the test statistics and
#values can be found in Leung et al (2000), page 21-23.

#Here goes the function:

LMZ.F3GWR.test <- function(go) {

	if (class(go) != "gwr") stop ("Not a GWR object")
	if (!go$hatmatrix) stop("Fit GWR model with argument hatmatrix=TRUE")
	
#	this.call <- match.call()

#	y <- go$lm$y
	x <- go$lm$x
	n <- NROW(x)
	m <- NCOL(x)
	betas <- as.matrix(as(go$SDF, "data.frame")[,(1+(1:m))])
# (3 + (1:m)) bug spotted by Leong Yin-Yee 080614
	delta1 <- go$results$delta1
	delta2 <- go$results$delta2
	sigma2 <- go$results$sigma2
	gweight <- eval(parse(text=go$gweight))
	bw <- go$bandwidth
	adapt <- go$adapt
	coords <- coordinates(go$SDF)
	if (is.null(adapt)) {
		bandwidth <- rep(bw, n)
	} else {
		bandwidth <- gw.adapt(dp=coords, fp=coords, quant=adapt)
	}
	if (any(bandwidth < 0)) stop("Invalid bandwidth")

	#Define a set of eks, as a matrix to hold the eks, the identity matrix and the matrix J
	
	ek <- diag (m)
	iden <- diag (n)
	J <- matrix (1, nrow = n, ncol = n)
	
	#Vk square can be obtained through formula (56) in Leung et al (2000)'s paper, page 22
	
	Vk2 <- numeric (m)
	for (i in 1:m){
		Vk2[i] <- (1/n)*(t(betas[,i])%*%(iden-(1/n)*J)%*%betas[,i])
	}
	
	#A B matrix is needed for evaluating Gammai (the trace of the matrix 
	#in Leung et al (2000)'s paper, page 22). 
	#Define the gamma1, gamma2, numerator's degree of freedom (numdf), Ftest 
	#and p-values to hold these numerics, the denominator's degree of freedom
	#(dendf) is always delta1^2/delta2
	
	gamma1 <- numeric (m)
	gamma2 <- numeric (m)
	numdf <- numeric (m)
	Ftest <- numeric (m)
	pvalue <- numeric (m)
	dendf <- delta1^2/delta2
	
	for (i in 1:m){
		B <- matrix (nrow = n, ncol = n)
		for (j in 1:n){
			wj <- gweight(spDistsN1(coords, coords[j,])^2, 
				bandwidth[j])
			B[j,] <- ek[i,] %*% solve(t(x)%*%diag(wj)%*%x) %*%
				t(x) %*% diag(wj)
		}
		internal <- (1/n)*(t(B)%*%(iden-(1/n)*J)%*%B)
		gamma1[i] <- sum(diag(internal))
		gamma2[i] <- sum(diag(internal)^2)
		numdf[i] <- gamma1[i]^2/gamma2[i]
		Ftest[i] <- (Vk2[i]/gamma1[i])/sigma2
		pvalue[i] <- pf(Ftest[i], numdf[i], dendf, lower.tail=FALSE)
	}
	
	#Report the result
	
	res <- matrix (nrow = m, ncol = 4)
	res[,1] <- Ftest
	res[,2] <- numdf
	res[,3] <- dendf
	res[,4] <- pvalue
	colnames(res) <- c("F statistic", "Numerator d.f.", 
		"Denominator d.f.", "Pr(>)")
	rownames(res) <- colnames(x)
	cat("\nLeung et al. (2000) F(3) test\n\n")
	printCoefmat(res)
	invisible(res)
}


