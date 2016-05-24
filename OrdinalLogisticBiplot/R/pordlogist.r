# file OrdinalLogisticBiplot/R/pordlogist.R
# copyright (C) 2012-2013 J.L. Vicente-Villardon and J.C. Hernandez
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

pordlogist <- function(y, x, penalization = 0.1, tol = 1e-04, maxiter = 200, show = FALSE) {
	    if (is.matrix(x)) {
        n <- nrow(x)
    }
    else {
        n <- length(x)
    }
    model= pordlogistfit(y,x, penalization = penalization, tol = tol,  maxiter = maxiter, show = show)
    null= pordlogistfit(y,x=NULL, penalization = penalization, tol = tol,  maxiter = maxiter, show = show)
	
	
  model$DevianceNull = null$Deviance
	model$Dif = (model$DevianceNull - model$Deviance)
	model$df = model$nvar
	model$pval = 1 - pchisq(model$Dif, df = model$df)
	model$CoxSnell = 1 - exp(-1 * model$Dif/n)
	model$Nagelkerke = model$CoxSnell/(1 - exp((model$DevianceNull/(-2)))^(2/n))
	model$MacFaden = 1 - (model$Deviance/model$DevianceNull)
	
	class(model) = "pordlogist"

	return(model)

}

summary.pordlogist <- function(object,...) {
  x = object
	npar = x$J + x$nvar - 1
	J = x$J
	cat(paste("Ordinal Logistic regression Model", "with Ridge Penalizations", x$penalization, " and logit link"), "\n\n")
	cat("n: ", x$nobs, "\n")
	cat("logLikelihood: ", x$logLik, "\n")
	cat("Iterations: ", x$iter, "\n\n")
	cat("Coefficients : \n")
	coefs <- matrix(NA, length(x$coefficients), 4, dimnames = list(names(x$coefficients), c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
	std = matrix(sqrt(diag(x$Covariaces)), npar, 1)
	coefs[, 1] <- x$coefficients
	coefs[, 2] <- std[J:npar]
	coefs[, 3] <- x$coefficients/std[J:npar]
	coefs[, 4] <- 2 * (1 - pnorm(abs(coefs[, 3])))
	print(coefs)
	cat("\n\nThreshold Coefficients : \n")
	names = "1|2"
	if(J > 2){
  	for (i in 2:(J - 1)) {
  		name = paste(i, "|", i + 1, sep = "")
  		names = c(names, name)
  	}
	}
	coefs2 <- matrix(NA, length(x$thresholds), 4, dimnames = list(names, c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
	std = matrix(sqrt(diag(x$Covariaces)), npar, 1)
	coefs2[, 1] <- x$thresholds
	coefs2[, 2] <- std[1:(J - 1)]
	coefs2[, 3] <- x$thresholds/std[1:(J-1)]
	coefs2[, 4] <- 2 * (1 - pnorm(abs(coefs2[, 3])))
	print(coefs2)
	cat("\n\nPseudo R-squared : \n")
	cat("Cox & Snell: ", x$CoxSnell, "\n")
	cat("Nagelkerke: ", x$Nagelkerke, "\n")
	cat("MacFaden: ", x$MacFaden, "\n")
}


pordlogistfit <- function(y, x = NULL, penalization = 0.1, tol = 1e-04, maxiter = 200, show = FALSE) {

	n <- length(y)
	if (is.null(x)) {
		p = 0
	} else {
		if (!is.matrix(x)) {
			x = as.matrix(x)
		}
		p <- ncol(x)
	}
	J = max(y)
	npar = J + p - 1
	Y = matrix(0, n, J)
	for (i in 1:n) if (y[i] > 0)
		Y[i, y[i]] = 1
	R = matrix(0, n, J)
	for (i in 1:n) {
		R[i, ] = cumsum(Y[i, ])
	}

  if (J > 2) {
		A = matrix(colSums(R[, 1:J - 1])/n, J - 1, 1)
	} else if (J == 2) {
		A = matrix(sum(R[, 1:J - 1])/n, J - 1, 1)
	} else {
		stop("There is a variable with the same value for all the items. Revise the data set.")
	}

	Beta = matrix(0, p, 1)
	B = rbind(A, Beta)
	err = 0.1
	iter = 0
	while ((err > tol) & (iter < maxiter)) {
		iter = iter + 1

		if (is.null(x))
			ETA = matrix(1, n, 1) %*% t(A)
		else ETA = matrix(1, n, 1) %*% t(A) - x %*% Beta %*% matrix(1, 1, (J - 1))


		PIA = exp(ETA)/(1 + exp(ETA))
		PIA = cbind(PIA, matrix(1, n, 1))
		PI = matrix(0, n, J)
		PI[, 1] = PIA[, 1]
		PI[, 2:J] = PIA[, 2:J] - PIA[, 1:(J - 1)]
		Rho = log(PIA[, 1:(J - 1)]/(PIA[, 2:J] - PIA[, 1:(J - 1)]))
		gRho = log(1 + exp(Rho))
		U = PIA[, 2:J]/(PIA[, 1:(J - 1)] * (PIA[, 2:J] - PIA[, 1:(J - 1)]))
		D = R[, 1:(J - 1)] - R[, 2:J] * (PIA[, 1:(J - 1)]/PIA[, 2:J])
		D2 = PIA[, 1:J] * (1 - PIA[, 1:J])
		P = array(0, c(n, J, npar))
		Q = array(0, c(n, J - 1, npar))
		GR = matrix(0, npar, 1)
		for (k in 1:npar) {
			if (k < J)
				P[, k, k] = matrix(1, n, 1)
			else P[, 1:(J - 1), k] = (-1 * x[, k - (J - 1)]) %*% matrix(1, 1, (J - 1))
			Q[, , k] = P[, 1:(J - 1), k] * D2[, 1:(J - 1)] - P[, 2:J, k] * (PIA[, 1:(J - 1)]/PIA[, 2:J]) * D2[, 2:J]
			GR[k] = sum(D * U * Q[, , k])
		}
		# Expectation of the second derivative, (Fisher's scoring method is used)\n
		HESS = matrix(0, npar, npar)
		for (s in 1:npar) for (k in 1:npar) {
			HESS[s, k] = sum(U * Q[, , s] * Q[, , k])
		}
		HESS = HESS + 2 * diag(J + p - 1) * penalization
		Bnew = B + solve(HESS) %*% (GR - 2 * B * penalization)
		err = sum((Bnew - B)^2)
		B = Bnew
		A = matrix(B[1:(J - 1)], J - 1, 1)
		if (is.null(x))
		Beta = matrix(0, p, 1)
		else
		Beta = matrix(B[J:(J + p - 1)], p, 1)

		if (show)
			print(c(iter, err))
	}

	L = sum(R[, 1:(J - 1)] * Rho - R[, 2:J] * gRho)

	Deviance = -2 * sum(L)

	model <- list()
	model$nobs = n
	model$J = J
	model$nvar = p
	model$fitted.values = PI
	model$pred = matrix(max.col(PI), n, 1)
	model$Covariaces = solve(HESS)
	model$clasif = table(y, model$pred)
	model$PercentClasif = sum(y == model$pred)/n
	model$coefficients = Beta
	model$thresholds = A
	model$logLik = L
	model$penalization = penalization
	model$Deviance = Deviance
	model$iter = iter
	return(model)
}