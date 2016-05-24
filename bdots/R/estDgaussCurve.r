est.dgauss.curve <- function(time, fixations, rho, conc, params = NULL, cor = TRUE) {
	if(is.null(params)) {
		mu    <- find.mu(time, fixations, conc)
		ht    <- find.ht(time, fixations, conc)
		sig1  <- find.sig1(time, fixations, conc)
		sig2  <- find.sig2(time, fixations, conc)
		base1 <- find.base1(time, fixations, conc)
		base2 <- find.base2(time, fixations, conc)
	} else {
		mu    <- params[1]
		ht    <- params[2]
		sig1  <- params[3]
		sig2  <- params[4]
		base1 <- params[5]
		base2 <- params[6]
	}
	if(cor) {
		fit.curve <- tryCatch(gnls(fixations ~ (time < mu) * (exp(-1 * (time - mu) ^ 2
												/ (2 * sig1 ^ 2)) * (ht - base1) + base1) + (mu <= time) * (exp(-1 * (time - mu) ^ 2 /
												(2 * sig2 ^ 2)) * (ht - base2) + base2), 
											start=c(mu = mu, ht = ht, sig1 = sig1, sig2 = sig2, base1 = base1, base2 = base2),
											correlation=corAR1(rho)), error = function(e) NULL)
		if(is.null(fit.curve)) cor <- FALSE
	}
	if(!cor) {
		fit.curve <- tryCatch(gnls(fixations ~ (time < mu) * (exp(-1 * (time - mu) ^ 2
												/ (2 * sig1 ^ 2)) * (ht - base1) + base1) + (mu <= time) * (exp(-1 * (time - mu) ^ 2 /
												(2 * sig2 ^ 2)) * (ht - base2) + base2), 
											start=c(mu = mu, ht = ht, sig1 = sig1, sig2 = sig2, base1 = base1, base2 = base2)), error = function(e) NULL)
	}
	list(fit = fit.curve, cor = cor)
}