# FREGAT (c) 2016 Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

pval.famBT <- function(Z) {

	w <- Z$w
	Z <- P11CholInvCo %*% Z$Z # Cor
	if (ncol(Z) > 1) {
		Z <- Z %*% diag(w) # n x m
		Z <- rowMeans(Z)
	} else {
		Z <- Z * w
	}
	se.beta0 <- 1 / sum(Z * Z)
	se.beta <- se.beta0 * nullmod$total.var
	est.beta <- sum(Z * as.vector(pheno)) * se.beta0
	chi2 <- (est.beta ^ 2) / se.beta
	p <- pchisq(chi2, 1, lower.tail = F)

	#a<-summary(lm(pheno ~ Z - 1))$coef
	return(c(p, est.beta, sqrt(se.beta)))#, m0, m1))

}

pval.MLR <- function(Z) {
#browser()
#	if (missing(geno)) geno <- P11CholInvCo %*% Z$Z  # making centralized and independent genotypes # Cor
	if (H2est) Z$Z <- CholInvCo %*% Z$Z

	fit <- glm(pheno ~ covariate + Z$Z - 1, family = "gaussian")
	if (stat == 'F') {
		p <- anova(fit, test = stat)[3, 6]
	} else {
		p <- anova(fit, test = stat)[3, 5]
	}
	return(p)

}

pval.famFLM <- function(Z) {

#	geno <- CholInvCo %*% Z$Z  # making centralized and independent genotypes # Cor
#	m1 <- dim(geno)[2]
	m1 <- dim(Z$Z)[2]
	model <- model0

	# base condition is m >= kg >= kb
	# previously ensured that
	# if g then kg > kb
	# if kg == kb then !g
	# all fourier bases are odd

	# g:
	# m >= kg => default
	# m < kg => kg <- m
		# if fourier if kg is even => kg <- kg -1
		# ensured that m >= kg 
		# if kg <= kb => MLR
		# if kg > kb => genobasis recalculated

	if (g) {
		if (m1 >= kg0) {
			genobasis <- genobasis0
			J <- J0
		} else {
			if (GVF == 'bspline') {
				order <- min(order0, m1)
				kg <- m1
			} else if (GVF == 'fourier') {
				if (m1 %% 2) kg <- m1 else kg <- (m1 - 1)
			}
			if (kg <= kb0) {  # m1 <= kg = kb, (BB, BF) -> 0B, (FF, FB) -> 0F with kb, betabasis recalculated
				test <- 'MLR'
			} else {  # m1 >= kg > kb, genobasis recalculated
				if (GVF == 'bspline') {
					genobasis <- create.bspline.basis(norder = order, nbasis = kg)
				} else if (GVF == 'fourier') {
					genobasis <- create.fourier.basis(c(0, 1), nbasis = kg)
				}
				J <- inprod(genobasis, betabasis0)
				model <- paste(GVF, kg, '-', BSF, kb0, sep = '')
			}
		}

	# !g:
	# m > kb => default
	# m <= kb => MLR

	} else {
		if (m1 > kb0) {
			betabasis <- betabasis0
		} else {
			test <- 'MLR'
		}
	}

	if (test == 'MLR') {
		test <- 'famFLM'
		return(c(pval.MLR(Z), 'MLR'))
	}

	### Calculation of matrix FI in given positions of a region
#	return(pos)
	if (g) {
		B <- eval.basis(Z$pos, genobasis)
		### Calculation of formula (1) for functional genotypes (depend on ONLY positions)
		FFF <- ginv(t(B) %*% B) %*% t(B)
#		U <- geno %*% t(FFF)
		U <- Z$Z %*% t(FFF)
		UJ <- matrix( U %*% J, ncol = dim(J)[2])
#		browser()
	} else {
		B <- eval.basis(Z$pos, betabasis)
#		UJ <- geno %*% B
		UJ <- Z$Z %*% B
	}
	if (H2est) UJ <- CholInvCo %*% UJ  # making centralized and independent genotypes # Cor
#	fit <- glm(pheno ~ covariate + UJ - 1, family = "gaussian")
#browser()
	fit <- glm(pheno ~ covariate + UJ - 1, family = "gaussian")
	if (stat == 'F') {
		p <- anova(fit, test = stat)[3, 6]
	} else { p <- anova(fit, test = stat)[3, 5] }

	c(p, model)

}
