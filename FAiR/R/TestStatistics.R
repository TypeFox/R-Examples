#     This file is part of FAiR, a program to conduct Factor Analysis in R
#     Copyright 2008 Benjamin King Goodrich
#
#     Some portions of this code are Copyright 2007-2008 Anne Boomsma and Walter Herzog
#     and are licensed under the GPL version 3 or later and some portions of this code
#     are Copyright 2007 John Fox and licensed under the GPL version 2 or later. The
#     modified code is indicated at the top of each function.
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


## This file contains functions that pertain to hypothesis testing, model comparison, etc.
## They are not intended to be accessed directly by rather via model_comparison()

## NOTE: This file is meant to be read with 90 columns

## Satorra-Bentler (2001) scaled-difference statistic (easy way)
FAiR_paired <-
function(M0, M1, how) {
	U0 <- FAiR_make_U(M0, how)
	U1 <- FAiR_make_U(M1, how)
	Ud <- U0 - U1
	Gamma <- FAiR_make_Gamma(M0, how)
	m <- df.residual(M1) - df.residual(M0)
	c <- sum(diag(Ud %*% Gamma)) / m
	Td <- (M0@manifest@n.obs - 1) * ( deviance(M0) - deviance(M1) ) / c
	p.value   <- pchisq(Td, m, lower.tail = FALSE)
	out <- list(statistic = Td, parameter = m, p.value = p.value, how = how)
	class(out) <- "htest"
	return(out)
}

## returns information criteria
FAiR_infocriteria <-
function(FAobject, null_model) {
	if(!FAiR_is.ML(FAobject)) {
		stop("calculating information criteria is possible only for maximum ",
			"likelihood models")
	}
	S <- model.matrix(FAobject, standardized = FALSE)
	n <- nrow(S)
	N <- FAobject@manifest@n.obs
	BIC <- BIC(FAobject)
	llik_saturated <- -0.5 * (N - 1) * c(determinant(S)$modulus + n)
	df_0 <- 0.5 * n * (n + 1)
	BIC_saturated <- -2 * llik_saturated + df_0 * log(N)
	BIC_null <- BIC(null_model)

	llik <- c(logLik(FAobject))
	vcov <- vcov(FAobject)
	infomatrix <- chol2inv(chol(vcov * N))
	SIC <- -llik + 0.5 * c(determinant(N * infomatrix)$modulus)
# 	SIC <- -llik + 0.5 * (ncol(infomatrix) * log(N) + determinant(infomatrix))$modulus
	llik <- c(logLik(null_model))
	vcov <- vcov(null_model)
	infomatrix <- chol2inv(chol(vcov * N))
	SIC_null <- -llik + 0.5 * c(determinant(N * infomatrix)$modulus)
	infomatrix <- FAiR_Browne1974(S)
	SIC_saturated <- -llik_saturated + 0.5 * c(determinant(N * infomatrix)$modulus)

	return(list(BIC = BIC, BIC_saturated = BIC_saturated, BIC_null = BIC_null,
			SIC = SIC, SIC_saturated = SIC_saturated, SIC_null = SIC_null))
}

## one of the correction factors suggested in Swain (1975)
FAiR_Swain <-
function(FAobject) {
	## The following is modified from http://www.ppsw.rug.nl/~boomsma/swain.R, which
	## is Copyright 2007-2008 Anne Boomsma and Walter Herzog and is licensed under the
	## GPL V3+

	# Equation numbers refer to Herzog, Boomsma and Reinecke (2007)
	d <- df.residual(FAobject)              # degrees of freedom
	p <- nrow(model.matrix(FAobject))       # number of manifest variables
	n <- FAobject@manifest@n.obs - 1        # one less the number of observations
	t <- p * (p + 1) / 2 - d                # number of parameters
	q <- ( sqrt(1 + 8 * t) - 1 ) / 2        # (11) 
	num <- p * (2 * p^2 + 3 * p - 1) - q * (2 * q^2 + 3 * q - 1)   # numerator in (10)
	den <- 12 * d * n                       # denominator in (10)
	s <- 1 - num / den                      # Swain's correction factor s (10)
	return(s)
}

## the correction factor suggested in Bartlett (1950)
FAiR_Bartlett <-
function(FAobject) {
	## This function is slightly modified from stats:::factanal, which is 
        ## Copyright 1995-2007 R Core Development Team and licensed under GPL V2+
	n <- FAobject@manifest@n.obs
	p <- nrow(model.matrix(FAobject))
	factors <- FAobject@restrictions@factors[1]
	b <- n - 1 - (2 * p + 5) / 6 - (2 * factors) / 3
	b <- b / (n - 1)
	if(!FAiR_is.EFA(FAobject)) {
		warning("Bartlett correction is intended for EFA models only")
	}
	return(b)
}

## classic test of the null hypothesis that the model holds in the population
FAiR_test_exact_fit <-
function(FAobject, correction) {
	chisq <- deviance(FAobject)
	     if(correction == "swain")    scale_factor <- FAiR_Swain(FAobject)
	else if(correction == "bartlett") scale_factor <- FAiR_Bartlett(FAobject)
	else scale_factor <- 1
	statistic <- scale_factor * chisq
	names(statistic) <- paste("T (", switch(correction, swain = "Swain",
						bartlett = "Bartlett", none = "No"),
					"correction )")
	parameter <- df.residual(FAobject)
	names(parameter) <- "df"
	p.value   <- pchisq(statistic, parameter, lower.tail = FALSE)
	null.value <- 0
	names(null.value) <- "discrepancy"
	out <- list(statistic = statistic, parameter = parameter, p.value = p.value,
			null.value = null.value, method = "Test of Exact Fit",
			alternative = "greater")
	class(out) <- "htest"
	return(out)
}

## Make Bentler's "U" matrix
FAiR_make_U <-
function(FAobject, how) {
	Gamma <- FAiR_make_Gamma(FAobject, how)
	W <- chol2inv(chol(Gamma))
	W_chol <- chol(W)
	Jacobian <- FAobject@Jacobian
	WJacobian <- W %*% Jacobian
	middle <- t(chol(chol2inv(chol(crossprod(W_chol %*% FAobject@Jacobian)))))
	U <- W - tcrossprod(WJacobian %*% middle)
	return(U)
}

## Make Bentler's "Gamma" various ways
FAiR_make_Gamma <-
function(FAobject, how) {
	C <- fitted(FAobject, standardized = FALSE, reduced = FALSE)
	if(how == "normal") Gamma <- FAiR_Browne1974(C)
	else                Gamma <- chol2inv(chol(FAiR_make_W(FAobject)))
	return(Gamma)
}

## Browne (1984) equation 2.20b
FAiR_test_Browne84 <-
function(FAobject, how) {
	n <- FAobject@manifest@n.obs - 1
	S <- model.matrix(FAobject, standardized = FALSE)
	s <- S[lower.tri(S, FAobject@manifest@diag)]
	C <- fitted(FAobject, standardized = FALSE, reduced = FALSE)
	c <- C[lower.tri(C, FAobject@manifest@diag)]
	U <- FAiR_make_U(FAobject, how) # U is singular
	statistic <- n * (s - c) %*% U %*% (s - c)
	names(statistic) <- "T ( Browne )"
	parameter <- df.residual(FAobject)
	names(parameter) <- "df"
	p.value   <- pchisq(statistic, parameter, lower.tail = FALSE)
	null.value <- 0
	names(null.value) <- "discrepancy"
	out <- list(statistic = statistic, parameter = parameter, p.value = p.value,
			null.value = null.value, method = "Test of Exact Fit",
			alternative = "greater")
	class(out) <- "htest"
	return(out)
}

## Sartorra and Bentler (1994) test
FAiR_test_SB94 <-
function(FAobject, correction, how) {
	statistic <- FAiR_test_exact_fit(FAobject, correction)$statistic
	d <- df.residual(FAobject)
	Gamma <- FAiR_make_Gamma(FAobject, how = "mle")
	U <- FAiR_make_U(FAobject, how)
	UGamma <- U %*% Gamma
	trace <- sum(diag(UGamma))
	parameter <- trace^2 / sum(diag(UGamma %*% UGamma))
	statistic <- parameter * statistic / trace
	p.value   <- pchisq(statistic, parameter, lower.tail = FALSE)
	out <- list(statistic = statistic, parameter = parameter, p.value = p.value,
			correction = correction)
	class(out) <- "htest"
	return(out)
}

## Yuan and Bentler (1998) test, in BJMSP
FAiR_test_YB98 <-
function(FAobject, statistic = NULL) {
	n <- FAobject@manifest@n.obs - 1
	if(is.null(statistic)) statistic <- deviance(FAobject)
	T_YB <- statistic / (1 + statistic / n)
	names(T_YB) <- "T_{YB}"
	parameter <- df.residual(FAobject)
	names(parameter) <- "df"
	p.value   <- pchisq(T_YB, parameter, lower.tail = FALSE)
	null.value <- 0
	names(null.value) <- "discrepancy"
	out <- list(statistic = T_YB, parameter = parameter, p.value = p.value,
			null.value = null.value, alternative = "greater",
			method = "Yuan and Bentler (1998) Test of Exact Fit")
	class(out) <- "htest"
	return(out)
}

## Hotelling's T^2 test
FAiR_test_T2 <-
function(FAobject) {
	if(FAobject@restrictions@discrepancy != "ADF") {
		stop("T^2 statistic requires the ADF discrepancy function")
	}
	N <- FAobject@manifest@n.obs
	fits <- FAobject@optimization@extraction$value
	T2 <- N * fits[length(fits)]
	num_df <- df.residual(FAobject)
	denom_df <- N - num_df
	statistic <- ( T2 / (N - 1) ) * ( denom_df / num_df )
	names(statistic) <- "Hotelling"
	p.value <- pf(statistic, num_df, denom_df, lower.tail = FALSE)
	parameter <- c(num_df, denom_df)
	names(parameter) <- c("df_numerator", "df_denominator")
	null.value <- 0
	names(null.value) <- "discrepancy"
	out <- list(statistic = statistic, parameter = parameter, p.value = p.value,
			null.value = null.value, alternative = "greater",
			method = "T^2 (Hotelling) Test of Exact Fit")
	class(out) <- "htest"
	return(out)
}

## Root Mean Square Error of Approximation
FAiR_RMSEA <-
function(FAobject, statistic, conf.level) {
	## The following is modified from sem:::summary.sem, which is
	## Copyright 2007 John Fox and is licensed under the GPL V2+

	N <- FAobject@manifest@n.obs
	df <- df.residual(FAobject)
	RMSEA <- sqrt(max(statistic / ( (N - 1) * df ) - 1 / (N - 1), 0))
	tail <- (1 - conf.level) / 2
	max <- N	
	while (max > 1) {
		res <- optimize(function(lam) (tail - pchisq(statistic, df, ncp = lam))^2,
                 		interval = c(0, max))
		if (sqrt(res$objective) < tail/100) break
		max <- max / 2
	}
	lam.U <- if (max <= 1) NA_real_ else res$minimum
	max <- max(max, 1)
	while (max > 1) {
		res <- optimize(function(lam) (1 - tail - pchisq(statistic, df, 
								ncp = lam))^2,
                 		interval = c(0, max))
		if (sqrt(res$objective) < tail/100) break
		max <- max / 2
	}
	lam.L <- if (max <= 1) NA_real_ else res$minimum
	RMSEA.U <- sqrt(lam.U / ( (N - 1) * df) )
	RMSEA.L <- sqrt(lam.L / ( (N - 1) * df) )
	ncp <- 0.05^2 * (N - 1) * df
	parameter <- df.residual(FAobject)
	names(parameter) <- "df"
	p.value <- pchisq(statistic, parameter, ncp, lower.tail = FALSE)
	null.value <- 0.05
	names(null.value) <- "discrepancy"
	conf.int <- c(RMSEA.L, RMSEA.U)
	attr(conf.int, "conf.level") <- conf.level
	estimate <- RMSEA
	names(estimate) <- "RMSEA"
	out <- list(statistic = statistic, parameter = parameter, p.value = p.value,
			null.value = null.value, method = "Test of Close Fit",
			alternative = "greater", conf.int = conf.int, estimate = estimate)
	class(out) <- "htest"
	return(out)
}

## Root Mean Square Error of Approximation
FAiR_RMSEA <-
function(FAobject, statistic, conf.level) {
	## The following is modified from http://www.ppsw.rug.nl/~boomsma/swain.R, which
	## is Copyright 2007-2008 Anne Boomsma and Walter Herzog and is licensed under the
	## GPL V3+
	n <- FAobject@manifest@n.obs - 1
	d <- df.residual(FAobject)
	RMSEA <- sqrt(max((statistic - d) / (n * d), 0))
	foo <- function(l) {
		ifelse(l > 0, pchisq(statistic, d, l) -     (1 - conf.level) / 2, 1)
	}
	BIGINT <- .Machine$integer.max
	upper  <- uniroot(foo, interval = c(0,BIGINT))$root
	RMSEAu <- sqrt( max(upper / (n * d), 0) )
	foo <- function(l) {
		ifelse(l > 0, pchisq(statistic, d, l) - 1 + (1 - conf.level) / 2, 1)
	}
	lower  <- uniroot(foo, interval = c(0,upper))$root
	RMSEAl <- sqrt( max(lower / (n * d), 0) )

	L <- 0.05^2 * n * d
	names(d) <- "df"
	p.value <- pchisq(statistic, d, L, lower.tail = FALSE)
	null.value <- 0.05
	names(null.value) <- "discrepancy"
	conf.int <- c(RMSEAl, RMSEAu)
	attr(conf.int, "conf.level") <- conf.level
	estimate <- RMSEA
	names(estimate) <- "RMSEA"
	out <- list(statistic = statistic, parameter = d, p.value = p.value,
			null.value = null.value, method = "Test of Close Fit",
			alternative = "greater", conf.int = conf.int, estimate = estimate)
	class(out) <- "htest"
	return(out)
}

## Steiger's gamma
FAiR_Steiger <-
function(FAobject, RMSEA) {
	## The following is modified from http://www.ppsw.rug.nl/~boomsma/swain.R, which
	## is Copyright 2007-2008 Anne Boomsma and Walter Herzog and is licensed under the
	## GPL V3+
	n <- FAobject@manifest@n.obs - 1
	d <- df.residual(FAobject)
	p <- nrow(model.matrix(FAobject))
	statistic <- RMSEA$statistic
	BIGINT <- .Machine$integer.max
	conf.level <- attributes(RMSEA$conf.int)$conf.level
	foo <- function(l) {
		ifelse(l > 0, pchisq(statistic, d, l) -     (1 - conf.level) / 2, 1)
	}
	upper  <- uniroot(foo, interval = c(0,BIGINT))$root
	foo <- function(l) {
		ifelse(l > 0, pchisq(statistic, d, l) - 1 + (1 - conf.level) / 2, 1)
	}
	lower  <- uniroot(foo, interval = c(0,upper))$root
	Gamma1  <- p / (p + 2 * (statistic - d) / n)
	Gamma1l <- p / (p + 2 * upper / n)
	Gamma1u <- p / (p + 2 * lower / n)
	conf.int <- c(Gamma1l, Gamma1u)
	attr(conf.int, "conf.level") <- conf.level
	estimate <- Gamma1
	names(estimate) <- "Gamma_1"
	out <- list(estimate = estimate, conf.int = conf.int, 
			method = "Gamma Fit Index (Steiger)")
	class(out) <- "htest"
	return(out)	
}

## McDonald's Centrality Index
FAiR_McDonald <-
function(FAobject, statistic) {
	n <- FAobject@manifest@n.obs
	d <- df.residual(FAobject)
	estimate <- exp(-0.5 * (statistic - d) / n)
	names(estimate) <- "Index"
	out <- list(estimate = estimate, method = "Centrality Index (McDonald)")
	class(out) <- "htest"
	return(out)
}

## standardized root mean residual 
FAiR_SRMR <-
function(FAobject) {
	R <- rstandard(FAobject)^2
	return(sqrt(mean(R[lower.tri(R)])))
}

## Tucker-Lewis Index
FAiR_TLI <-
function(FAobject, statistic, statistic_0) {
	## The following is modified from http://www.ppsw.rug.nl/~boomsma/swain.R, which
	## is Copyright 2007-2008 Anne Boomsma and Walter Herzog and is licensed under the
	## GPL V3+
	p  <- nrow(model.matrix(FAobject))
	di <- p * (p - 1) / 2  # degrees of freedom of independence model
	d  <- df.residual(FAobject)
	TLI <- (statistic_0 / di - statistic / d) / (statistic_0 / di -1)
	return(TLI)
}

## Bentler's Comparative Fit Index
FAiR_CFI <-
function(FAobject, statistic, statistic_0) {
	## The following is modified from http://www.ppsw.rug.nl/~boomsma/swain.R, which
	## is Copyright 2007-2008 Anne Boomsma and Walter Herzog and is licensed under the
	## GPL V3+
	p  <- nrow(model.matrix(FAobject))
	d  <- df.residual(FAobject)
	di <- p * (p - 1) / 2  # degrees of freedom of independence model
	tau  <- max(statistic   - d , 0)
	taui <- max(statistic_0 - di, 0)
	CFI  <- 1 - tau / taui # Comparative Fit Index (Bentler)
	return(CFI)
}

## Goodness of Fit Index, Joreskog and Sorbom (1984)
FAiR_GFI <-
function(FAobject) {
	## The following is modified from sem:::summary.sem, which is
	## Copyright 2007 John Fox and is licensed under the GPL V2+

	S  <- model.matrix(FAobject, standardized = FALSE)
	n  <- nrow(S)
	C  <- fitted(FAobject, reduced = FALSE, standardized = FALSE)
	invC <- chol2inv(chol(C))
	CSC  <- invC %*% (S - C)
	CSC  <- CSC  %*% CSC
	CS   <- invC %*% S
	CS   <- CS   %*% CS
	GFI  <- 1 - sum(diag(CSC)) / sum(diag(CS))
	return(GFI)
}

## Adjusted GFI, Joreskog and Sorbom (1984)
FAiR_AGFI <-
function(FAobject, GFI = NULL) {
	## The following is modified from sem:::summary.sem, which is
	## Copyright 2007 John Fox and is licensed under the GPL V2+
	if(is.null(GFI)) GFI <- FAiR_GFI(FAobject)
	n  <- nrow(model.matrix(FAobject))
	df <- df.residual(FAobject)
	AGFI <- 1 - (n * (n + 1) / (2 * df)) * (1 - GFI)
	return(AGFI)
}

## Normalized Fit Index, Bentler and Bonett (1980)
FAiR_NFI <-
function(statistic, statistic_0) {
	## The following is modified from sem:::summary.sem, which is
	## Copyright 2007 John Fox and is licensed under the GPL V2+

        NFI  <- (statistic_0 - statistic) / statistic_0
	return(NFI)
}

## Nonnormalized Fit Index, Bentler and Bonett (1980)
FAiR_NNFI <-
function(FAobject, statistic, statistic_0) {
	## The following is modified from sem:::summary.sem, which is
	## Copyright 2007 John Fox and is licensed under the GPL V2+

	p <- nrow(model.matrix(FAobject))
	dfNull <- p * (p - 1) / 2 
	df <- df.residual(FAobject)
        NNFI <- (statistic_0 / dfNull - statistic / df) / (statistic_0 / dfNull - 1)
	return(NNFI)
}

## estimate a model with no correlations
FAiR_independence_model <-
function(FAobject) {
	factors <- FAobject@restrictions@factors[1]
	S <- model.matrix(FAobject, standardized = FALSE)
	     if(FAobject@restrictions@discrepancy == "YWLS") {
			stop("impossible to calculate the fit statistic of the null",
				" model when discrepancy = 'YWLS'")
		}
	else if(FAobject@restrictions@discrepancy == "MLE") {
		if(is(FAobject@restrictions, "restrictions.1storder")) {
			FAobject@restrictions@Phi@x <- diag(factors)
		}
		FAobject@restrictions@beta@x <- matrix(0, nrow(S), factors)
		FAobject@restrictions@Omega@x <- FAobject@manifest@sds
		return(FAobject)
	}
	else {
		new_restrictions <- FAiR_make_independent(FAobject@restrictions)
		cat("\nEstimating independence model ...\n")
		null_model <- Factanal(FAobject@manifest, new_restrictions,
					BFGSburnin = 0, gradient.check = FALSE,
					seeds = FAobject@seeds[1,])
		return(null_model)
	}
}

## RMSEA-based simulated test
FAiR_nonparametric <-
function(draws, FAobject) {
	if(length(dim(draws)) != 3) {
		stop("'draws' must be a three dimensional numeric array")
	}
	if(!is(FAobject, "FA")) {
		stop("'FAobject' must inherit from class FA")
	}

	reproduced <- fitted(FAobject, reduced = FALSE, standardized = FALSE)

	if(FAobject@restrictions@discrepancy == "MLE") {
		constant <- determinant(reproduced)$modulus + nrow(reproduced)
		foo <- function(draw) {
			log_det_C <-   determinant(draw)$modulus
			C_inv <- try(chol2inv(chol(draw)), silent = TRUE)
			if(!is.matrix(C_inv)) return(NA_real_)
			llik  <- log_det_C + as.vector(crossprod(as.vector(reproduced),
								 as.vector(C_inv)))
			return(llik - constant)
		}
	}
	else if(FAobject@restrictions@discrepancy == "YWLS") {
		stop("Simulated test statistic is not valid under the YWLS objective",
			" function")
	}
	else {
		l <- length(FAobject@restrictions@criteria)
		bar <-      FAobject@restrictions@criteria[[l]]
		if(formals(bar)$diag) mark <- lower.tri(reproduced, TRUE)
		else                  mark <- lower.tri(reproduced, FALSE)
		formals(bar)$s <- reproduced[mark]

		foo <- function(draw) {
			C <- draw
			environment(bar) <- environment()
			misfit <- bar()
			return(misfit)
		}
	}

	n <- FAobject@manifest@n.obs - 1
	nvars <- FAobject@restrictions@nvars
	null_dist <- apply(draws, 3, foo)
	null_dist <- null_dist[!is.na(null_dist)]
	null_dist <- sqrt( pmax(null_dist / nvars  - 1 / n, 0) )
						
	statistic <- sqrt(max(deviance(FAobject) / (n * df.residual(FAobject))  - 1/n, 0))
	names(statistic) <- "T_{np}"
	p.value <- mean(statistic < null_dist)
	null.value <- 0
	names(null.value) <- "discrepancy"
	out <- list(statistic = statistic, p.value = p.value, null.value = null.value, 
			method = "Nonparametric Test of Exact Fit",
			alternative = "greater")
	class(out) <- "htest"
	return(out)
}
