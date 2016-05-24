#     This file is part of FAiR, a program to conduct Factor Analysis in R
#     Copyright 2008,2013 Benjamin King Goodrich
#
#     FAiR is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     FAiR is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with FAiR.  If not, see <http://www.gnu.org/licenses/>.

## This file contains utility functions that are not meant to be accessed by users

## NOTE: This file is meant to be read with 90 columns

## communalities via brute force partitioning method; see Kano(1990)
FAiR_PACE <- 
function(Sigma, factors) { ## this takes forever in nontrivial problems and isn't used now
	temp.fun <- function(x) {
		the.list <- rep(list(x), factors)
		eg <- expand.grid(the.list)
		eg.sorted <- t(apply(eg, 1, FUN = function(z) {
			if(any(duplicated(z))) return(rep(NA_real_, length(x)))
			y <- x[!(x %in% z)]
			return(c(sort(z), sort(y)))
		}))
		eg <- unique(eg.sorted)
		eg <- eg[complete.cases(eg),1:factors,drop = FALSE]
		return(eg[1:(nrow(eg)/2),])
	}

	combos <- combn(1:nrow(Sigma), 2 * factors, FUN = temp.fun)

	dets <- apply(combos, 3, FUN = function(x) {
		uniques <- unique(c(x))
		apply(x, 1, FUN = function(two) {
			three <- uniques[!(uniques %in% two)]
			det(Sigma[two,three,drop = FALSE])
		})
	})
	communalities <- rep(NA_real_, nrow(Sigma))
	for(i in 1:nrow(Sigma)) {
		included <- apply(combos, 3, FUN = function(x) !(i %in% x))
		included_combos <- combos[,,included,drop = FALSE]
		included_dets <- dets[,included, drop = FALSE]
		best <- which.max(abs(included_dets))
		the.shelf <- ceiling(best / nrow(included_combos))
		the.row <- best %% nrow(included_combos)
		if(the.row == 0) the.row <- nrow(included_combos)
		uniques <- unique(c(included_combos[,,the.shelf]))
		two <- included_combos[the.row,,the.shelf]
		three <- uniques[!(uniques %in% two)]
		
		rho_21 <- Sigma[two, i]
		rho_31 <- Sigma[three,i]
		rho_23 <- Sigma[two,three]

		communalities[i] <- rho_31 %*% solve(rho_23) %*% rho_21
	}
	communalities <- ifelse(communalities >= 1, 1 - .Machine$double.eps, 
			 ifelse(communalities <= 0, 0 + .Machine$double.eps,
				communalities))
	names(communalities) <- colnames(Sigma)
	return(communalities)
}

## this is much faster but not guaranteed to give the "right" answer
## it is used to make starting communalities in Factanal()
FAiR_PACE_by_RGENOUD <- 
function(Sigma, factors, pp = paste(tempfile(), "PACE.txt", sep = ""), seeds, ...) {
	if(missing(seeds)) seeds <- c(formals(genoud)[c("unif.seed", "int.seed")])
	else if(length(seeds) == 1) seeds <- rep(seeds, 2)
	opt.fun <- function(par, first) {
		criterion.1 <-  all(first != par)
		criterion.2 <- length(unique(par)) / length(par)
		if( (criterion.1 == 1) && (criterion.2 == 1) ) {
			criterion.3 <- abs(det(Sigma[par[1:factors], 
						  par[-c(1:factors)], drop = FALSE]))
			criterion.3 <- criterion.3^(1/factors)
			out <- c(1,1,criterion.3)
		}
		else    out <- c(criterion.1, criterion.2, 0.0)
		return( out )
	}
	communalities <- (0.5 * factors / ncol(Sigma)) / diag(solve(Sigma))
	starts <- t(replicate(1000, sort(sample(1:nrow(Sigma), 2 * factors))))
	for(i in 1:nrow(Sigma)) {
		fits <- t(apply(starts, 1, opt.fun, first = i))
		the_range <- range(fits[fits[,1] == 1,3])
		## Do something useful with the_range
		proj.path <- paste(pp, i, sep = "_")
		if(i != 1) seeds <- c(formals(genoud)[c("unif.seed", "int.seed")])
		opt <- genoud(fn = opt.fun, nvars = 2 * factors, max = TRUE,
				Domains = cbind(rep(1.0, 2 * factors), nrow(Sigma)),
				boundary.enforcement = 2, lexical = TRUE, print.level = 0,
				data.type.int = TRUE, project.path = proj.path, 
				first = i, MemoryMatrix = FALSE, P3 = 0, 
				unif.seed = seeds[1], int.seed = seeds[2], 
				starting.values = starts, ...)
		if(any(opt$value < sqrt(.Machine$double.eps))) {
			warning("submatrix almost singular, using SMC")
			next
		}
# 		else starts <- opt$par
		two   <- opt$par[   1:factors,  drop = FALSE]
		three <- opt$par[-c(1:factors), drop = FALSE]

		rho_21 <- Sigma[two,  i]
		rho_31 <- Sigma[three,i]
		rho_23 <- Sigma[two,three, drop = FALSE]

		communalities[i] <- rho_31 %*% solve(rho_23) %*% rho_21	
	}
	communalities <- ifelse(communalities >= 1, 1 - sqrt(.Machine$double.eps), 
			 ifelse(communalities <= 0, 0 + sqrt(.Machine$double.eps),
				communalities))
	 names(communalities) <- colnames(Sigma)
	return(communalities)
}

## could make pop.size starting values for Rotate() but is not used currently
FAiR_make_starts_for_Rotate <- 
function(A, pop.size) {
	n <- nrow(A)
	r <- ncol(A)
	A_swept <- A / sqrt(rowSums(A^2))
	cutoff <- sort(combn(1:n, r, FUN = function(x) 
			det(crossprod(A_swept[x,]))), TRUE)[floor(pop.size / 2)]
	
	good <- combn(1:n, r, simplify = FALSE, 
		FUN = function(x) if(det(crossprod(A_swept[x,])) >= cutoff) 
					return(x) else return(NULL) )
	
	out <- array(NA, dim = c(r, r, 2 * floor(pop.size / 2)))
	count <- 0
	for(i in 1:length(good)) {
		if(is.null(good[[i]])) next
		count <- count + 1
		Tmat <- t(A_swept[good[[i]],])
		Tmat_svd <- svd(Tmat)
		out[,,count] <- Tmat
		out[,,pop.size - count + 1] <- Tmat_svd$u %*% t(Tmat_svd$v)
	}
	signs <- out[r,,]
	out <- matrix(c(out[-r,,]), ncol = pop.size)
	out <- rbind(out, signs)
	return(out)
}

## used to generate an orthogonal starting Tmat of order r in some create_start methods
FAiR_Landahl <- 
function(r) {
	if(r == 2) {
		Tmat <- rbind(rep(1/sqrt(r), 2), c(1/sqrt(r), -1/sqrt(r)))
		return(Tmat)
	}
	k <- (r-1):1
	submatrix <- diag( sqrt(k / (k+1) ), ncol = length(k) )
	Tmat <- rbind( 1/sqrt(r), cbind( submatrix, -1 / sqrt( (k+1) * k) ) )
	for(k in 2:(r-1)) for(j in k:(r-1)) Tmat[k,j] <- Tmat[k,r]
	return(Tmat)
}

## coerces a vector into an oblique transformation matrix, used in FAiR_Rotate()
FAiR_make_Tmat <- 
function(par) {
	factors <- sqrt(length(par))
	end <- factors^2 - factors

	# fill all but the last row of Tmat
	Tmat <- matrix(par[1:end], nrow = factors - 1, ncol = factors)
	colsums <- colSums(Tmat^2)
	Tmat <- rbind(Tmat, NA_real_)   # append another row to Tmat
	if(any(colsums > 1)) {          # decomplexify Tmat
		attributes(Tmat)$too.big <- max(colsums)
		Tmat[factors,] <- ifelse(colsums > 1, 
					par[(end +1):length(par)][colsums > 1],
					1 - colsums^(.5))
		Tmat <- t(t(Tmat) / sqrt(colSums(Tmat^2)))
	}
	else {  # all cells of Tmat real
		attributes(Tmat)$too.big <- -1.0
		Tmat[factors,] <- sqrt(1 - colsums) * # get the signs for the last row
					ifelse(par[(end + 1):length(par)] > 0, 1, -1)
	}
	return(Tmat)
}

## inverts FAiR_make_Tmat
FAiR_Tmat2par <-
function(Tmat) {
	par <- c(Tmat[1:(nrow(Tmat) - 1),], Tmat[nrow(Tmat),])
	return(par)
}

## this is the function that gets lexically minimized by Rotate()
FAiR_Rotate <- 
function(par, A, criteria) {
	out <- rep(0, 2 + length(criteria))
	Tmat <- FAiR_make_Tmat(par)
	out[1] <- attributes(Tmat)$too.big
	if(out[1] != -1) return(out)
	Phi  <- crossprod(Tmat)
	out[2] <- criteria[[1]](Phi)
	if(out[2] != -1) return(out)

	Phi_inv <- chol2inv(chol(Phi))
	rotmat_primary   <- Tmat %*% Phi_inv
	rotmat_reference <- t(t(rotmat_primary) / sqrt(colSums(rotmat_primary^2)))

	Upsilon <- A %*% rotmat_reference
	signs <- sign(colSums(Upsilon))
	out[3] <- -mean(signs)

	# calculate all criteria
	Pi <- A %*% Tmat
	FC <- t(t(Pi) * sqrt(diag(Phi_inv))) * Upsilon
	if( (l <- length(criteria)) > 2) for(i in 2:(l-1)) {
		environment(criteria[[i]]) <- environment()
		out[i + 2L] <- criteria[[i]](EFA = TRUE)
	}
	environment(criteria[[l]]) <- environment()
	last <- criteria[[l]]()
	if(length(last) > 1) out <- c(out[-length(out)], last)
	else out[length(out)] <- last
	return(out)
}

## Bentler's (1969) orthogonal approximation to an oblique Tmat, not currently used
FAiR_Tmat_approx <- 
function(Tmat) {
	Tmat_svd <- svd(Tmat)
	return(Tmat_svd$u %*% t(Tmat_svd$v))
}

## this function does the stronger submatrix rank check for Reiersol (1950)
# FAiR_check_Reiersol <- 
# function(coefs) {
# 	coefs <- !coefs
# 	cS <- colSums(coefs)
# 	if(any(cS < (ncol(coefs) - 1))) return(FAiR_uniquify(0))
# 	if(any(cS > (nrow(coefs) - 2))) return(FAiR_uniquify(0))
# 	if(any(apply(coefs, 1, all)))   return(FAiR_uniquify(0))
# 	for(p in 1:ncol(coefs)) {
# 		subcoefs <- coefs[coefs[,p],, drop = FALSE]
# 		for(k in 1:nrow(subcoefs)) {
# 			subsubcoefs <-  subcoefs[-k,, drop = FALSE]
# 			if(any(apply(subsubcoefs[,-p, drop = FALSE], 2, all))) {
# 				return(FAiR_uniquify(0))
# 			}
# 		}
# 	}
# 	return(-1.0)
# }

## this function does the stronger submatrix rank check for Reiersol (1950)
FAiR_check_Reiersol <- 
function(coefs) {
	coefs <- !coefs
	cS <- colSums(coefs)
	if(any(cS < (ncol(coefs) - 1))) return(FAiR_uniquify(0))
	if(any(cS > (nrow(coefs) - 2))) return(FAiR_uniquify(0))
	if(any(apply(coefs, 1, all)))   return(FAiR_uniquify(0))
	for(p in 1:ncol(coefs)) {
		subcoefs <- coefs[coefs[,p],-p, drop = FALSE]
		cS <- colSums(subcoefs)
		if(any(cS >= (nrow(subcoefs) - 1))) return(FAiR_uniquify(0))
# 		for(k in 1:nrow(subcoefs)) {
# 			subsubcoefs <-  subcoefs[-k,, drop = FALSE]
# 			for(l in 1:ncol(subsubcoefs)) if(all(subcoefs[,l])) return(out)
# 		}
	}
	return(-1.0)
}

## this function does the submatrix rank check for Howe (1955)
FAiR_check_Howe <- 
function(coefs) {
	coefs <- t(!coefs)
	rS <- rowSums(coefs)
	if(any(rS < (nrow(coefs) - 1))) return(FAiR_uniquify(0))
	if(any(rS > (ncol(coefs) - 2))) return(FAiR_uniquify(0))
	out <- -nrow(unique(coefs)) / nrow(coefs)
	if(out > -1) out <- FAiR_uniquify(out)
	return(out)
}

## a constructor for the FA object called at the end of Rotate() or from GPA2FA()
FAiR_opt2FAobject <- 
function(opt, FAobject, criteria, GPA = FALSE) {
	manifest <- FAobject@manifest
	S <- cormat(manifest)
	Lambda <- coef(FAobject)

	if(GPA) {
		FAobject@seeds <- rbind(FAobject@seeds, c(NA_integer_, NA_integer_))
		opt$method <- tolower(opt$method)
		orthogonal <- opt$orthogonal
		criteria <- list(method = opt$method)
		Tmat <- opt$Th
	}
	else {
		Tmat <- FAiR_make_Tmat(opt$par)
		orthogonal <- FALSE
		FAobject@optimization$transformation$opt <- opt
	}

	Phi  <- crossprod(Tmat)
	diag(Phi) <- 1
	rotmat_primary   <- Tmat %*% chol2inv(chol(Phi))
	rotmat_reference <- sweep(rotmat_primary, 2, sqrt(colSums(rotmat_primary^2)), "/")

	Psi <- crossprod(rotmat_reference)
	D   <- diag(diag(crossprod(Tmat, rotmat_reference)))

	PP <- Lambda %*% rotmat_primary
	RS <- Upsilon <- Lambda %*% rotmat_reference
	PS <- Lambda %*% Tmat
	RP <- Lambda %*% rotmat_reference %*% chol2inv(chol(Psi))
	FC <- PP * PS
	sorter <- order(colSums(FC), decreasing = TRUE)
	sorter <- 1:ncol(FC) # temporary
	loadings <- array(cbind(PP, RS, PS, RP, FC), dim = c(dim(PP), 5),
				dimnames = list(rownames(Lambda), NULL, 
				c("PP", "RS", "PS", "RP", "FC")))

        correlations <- array(cbind(Phi, Psi, D), dim = c(rep(ncol(D), 2), 3),
                                dimnames = list(NULL, NULL, c("PF", "RF", "PR")))

	signs <- ifelse(colSums(loadings[,,1]) >= 0, 1, -1)
	signs[signs == 0] <- 1
	if(any(signs != 1)) {
		loadings <- sweep(loadings, 2, signs, FUN = "*")
		loadings[,,5] <- FC
	        correlations  <- sweep(sweep(correlations, 1, signs, FUN = "*"),
         	                                           2, signs, FUN = "*")
		Phi <- Phi * tcrossprod(signs)
	}

	restrictions <- FAobject@restrictions
	restrictions@beta@x <- PP
	## FIXME: think about freedom
	Phi <- new("parameter.cormat", x = Phi, free = lower.tri(Phi),
			num_free = sum(lower.tri(Phi)), invalid = 0.0)
	new_restrictions <- new("restrictions.1storder.EFA", 
				factors = restrictions@factors,
				nvars = restrictions@nvars, dof = restrictions@dof,
				Domains = restrictions@Domains, model = "EFA", 
				discrepancy = restrictions@discrepancy, Phi = Phi, 
				beta = restrictions@beta, Omega = restrictions@Omega, 
				criteria = restrictions@criteria, Tcriteria = criteria,
				Lambda = Lambda, orthogonal = orthogonal, Tmat = Tmat, 
				free = restrictions@free)
	scores <- FAobject@call$scores
	if(is.null(scores)) scores <- "none"
 	scores  <- FAiR_scores(scores, manifest, PP, cormat(Phi), uniquenesses(FAobject))
	trans_mats <- array(cbind(rotmat_primary, rotmat_reference, Tmat), 
				dim = c(rep(ncol(D), 2), 3), dimnames = list(NULL, NULL,
				c("primary", "reference", "T")))

	correlations <- correlations[sorter,sorter,]
	attributes(correlations)$orthogonal <- orthogonal	
	new_FAobject <- new("FA.EFA", loadings = loadings[,sorter,], 
				correlations = correlations,
				uniquenesses = FAobject@uniquenesses,
				scale = FAobject@scale, restrictions = new_restrictions,
				Jacobian = FAobject@Jacobian, vcov = FAobject@vcov,
				scores = scores, manifest = FAobject@manifest,
				optimization = FAobject@optimization, 
				call = FAobject@call, seeds = FAobject@seeds,
				Lambda = Lambda, trans_mats = trans_mats, rotated = TRUE)

	colnames(new_FAobject@loadings) <-  colnames(new_FAobject@correlations)  <-
	rownames(new_FAobject@correlations) <- colnames(new_FAobject@trans_mats) <-
	colnames(new_FAobject@trans_mats) <- paste("F", 1:ncol(D), sep = "")

	return(new_FAobject)
}

## BHHH covariance estimator for models where the data are MVN, not currently used
FAiR_BHHH <- 
function(par, restrictions, manifest) {
	stop("this function should not have been called")
	free <- par != 0
	free_sum <- sum(free)
	if(restrictions@discrepancy != "MLE") {
		warning("BHHH can only be calculated under maximum likelihood")
		return(matrix(NA_real_, free_sum, free_sum))
	}
	if(manifest@n.obs < sum(free)) {
		warning("BHHH can only be calculated when the number of observations",
			"is greater than or equal to the number of free parameters")
		return(matrix(NA_real_, free_sum, free_sum))
	}
	if(!is(manifest, "manifest.data")) {
		warning("BHHH cannot be calculated because the raw data were not",
			" passed to Factanal()")
		return(matrix(NA_real_, free_sum, free_sum))
	}
	restrictions@model <- "CFA"
	z <- sweep(manifest@zz, 2, sqrt(diag(manifest@cov)), FUN = "*")
	helper2 <- bfgs_helpS4(par, object = restrictions, done = TRUE, 
			       S = manifest@cor, lower = sqrt(.Machine$double.eps))
	G <- t(apply(z, 1, FUN = function(x) {
			gr_fitS4(par, object = restrictions, helper2, 
				  S = tcrossprod(x), lower = sqrt(.Machine$double.eps))
		}))
	G <- G[,free] * sqrt(0.5)
	GG_inv <- try(chol2inv(chol(crossprod(G))))
	if(is.matrix(GG_inv)) return(GG_inv)
	else {
		warning("variance-covariance appears not to be positive definite")
		GG_inv <- try(solve(crossprod(G)))
		if(is.matrix(GG_inv)) return(GG_inv)
		else return(matrix(NA_real_, free_sum, free_sum))
	}
}

## calculates the gradient of the ML discrepancy function wrt the parameters
FAiR_gradient_MLE <-
function(par, restrictions, manifest, helper, lower) { # Joreskog (1971) version

	# Check bounds
	if(any(par < restrictions@Domains[,1]) | any(par > restrictions@Domains[,2])) {
		do.call(return, args = list(x = rep(.Machine$double.eps, length(par))),
			envir = parent.frame())
	}

	SEFA <- restrictions@model == "SEFA"
	if(SEFA) { # squash some coefficients and pretend it's a CFA
		small <- par[helper$squashed]
		par[helper$squashed] <- 0
	}

	# Fill in model
	model <- restrictions2model(par, restrictions, manifest, lower, FALSE)
	criteria <- model$criteria
	if(any(criteria != -1)) { # violated a constraint
		do.call(return, args = list(x = rep(.Machine$double.eps, length(par))),
			envir = parent.frame())
	}
	restrictions <- model$restrictions

	# Get all lexical criteria
	fits <- FAiR_lexical_driver(restrictions, manifest, lower = lower) 
	marker <- which(fits != -1)[1]
	if(marker < length(fits)) { # violated a constraint
		do.call(return, args = list(x = rep(.Machine$double.eps, length(par))),
			envir = parent.frame())
	}

	# Make reproduced covariance matrix
	R <- fitted(restrictions, reduced = TRUE)
	diag(R) <- 1
	sds <- restrictions@Omega@x
	scale_mat <- tcrossprod(sds)
	C <- R * scale_mat
	C_inverse <- chol2inv(chol(C))

	# Make "middle" matrix, called Omega by Joreskog
	S <- manifest@cov
	middle <- C_inverse %*% (C - S) %*% C_inverse

	# Calculate gradient wrt manifest standard deviations (diagonal of B in Joreskog)
	dF_dscale <- diag(2 * middle %*% (sds * R)) * sds

	beta   <- coef(restrictions)
	Phi  <- cormat(restrictions)
	betaPhi2 <- beta %*% Phi * 2

	# Calculate gradient wrt primary pattern matrix
	middle_scaled <- middle * scale_mat
	middle_scaled_diag <- diag(middle_scaled)
	dF_dbeta <- middle_scaled %*% betaPhi2
	dF_dbeta <- dF_dbeta - middle_scaled_diag * betaPhi2 # due to embeddedness

	# Calculate gradient wrt correlation matrix among primary factors
	if(!is(restrictions, "restrictions.orthonormal")) {
		dF_dPhi <- 2 * t(beta) %*% middle_scaled %*% beta
		diag(dF_dPhi) <- 0 # diag(Phi) = 1 and hence are not free parameters
	}
	else dF_dPhi <- NULL
	out <- list(dF_dPhi = dF_dPhi, dF_dbeta = dF_dbeta, dF_dscale = dF_dscale)
	return(out)
}

## calculates gradient of the ADF discrepancy function wrt the parameters
FAiR_gradient_QD <-
function(par, restrictions, manifest, helper, lower) {

	# Check bounds
	if(any(par < restrictions@Domains[,1]) | any(par > restrictions@Domains[,2])) {
		do.call(return, args = list(x = rep(.Machine$double.eps, length(par))),
			envir = parent.frame())
	}

	if(helper$marker == length(helper$fits)) {
		model <- restrictions2model(par, restrictions, manifest, lower, TRUE)
		criteria <- model$criteria
		if(any(criteria != -1)) return(rep(.Machine$double.eps, length(par)))
		restrictions <- model$restrictions
	
		fits <- FAiR_lexical_driver(restrictions, manifest, lower = lower) 
		marker <- which(fits != -1)[1]
		if(marker < length(fits)) return(rep(.Machine$double.eps, length(par)))
	
		l <- length(restrictions@criteria)
		middle <- formals(restrictions@criteria[[l]])$middle
		s <- formals(restrictions@criteria[[l]])$s

		if(length(s) == length(middle)) { # DWLS
			dc_dpar <- FAiR_deriv_matrix(par, restrictions, manifest, 
							lower, FALSE)
			C <- fitted(restrictions, reduced = TRUE)
			diag(C) <- 1
			C <- C * tcrossprod(restrictions@Omega@x)
	
			c <- C[lower.tri(C, FALSE)]
			x <- s - c
			length_x <- length(x)

			dF_dc <- -2 * middle * x
		}
		else {
			dc_dpar <- FAiR_deriv_matrix(par, restrictions, manifest, 
							lower, TRUE)
			C <- fitted(restrictions, reduced = TRUE)
			diag(C) <- 1
			C <- C * tcrossprod(restrictions@Omega@x)
	
			c <- C[lower.tri(C, TRUE)]
			x <- s - c
			length_x <- length(x)
			dF_dc <- .C(FAiR_QD_grad, grad = numeric(length_x), x = x, 
					middle = middle, length = length_x, 
					holder = numeric(length_x), DUP = FALSE, 
					NAOK = TRUE)$grad
		}
		gradient <- dF_dc %*% dc_dpar
		if(restrictions@model == "SEFA") {
			gradient[helper$squashed] <- 2*par[helper$squashed] * !helper$done
		}

	}
	else { # find numeric gradient of whatever the marginal criterion is
		gradient <- FAiR_numeric_gradient(par, restrictions, manifest, 
							helper, lower)
	}
	return(gradient)
}

## used to be used in create_start() methods; not used currently
FAiR_beta_under_orthogonality <- 
function(x, fixed, communalities) {
	beta <- fixed
	beta[is.na(fixed)] <- x
	beta <- sweep(beta, 1, sqrt(rowSums(beta^2) / communalities), "/")
	signs <- ifelse(colSums(beta) >= 0, 1, -1)
	beta <- sweep(beta, 2, signs, FUN = "*")
	return(beta[is.na(fixed)])
}

## this critical function is the driver for SEFA and CFA; called by fit_S4() default
FAiR_lexical_driver <- 
function(restrictions, manifest, lower) {
	out  <- rep(0, 2 + length(restrictions@criteria))
	MLE  <- restrictions@discrepancy == "MLE"

	R  <- fitted(restrictions, reduced = TRUE)   # make reduced cormat
	h2 <- diag(R)
	h2_max <- max(h2)
	if(h2_max > 1) {                             # check for Heywood cases
		out[1] <- h2_max
		return(out)
	}

	out[1]  <- -1
	diag(R) <-  1
	C <- R * tcrossprod(restrictions@Omega@x)   # scale to covariance matrix

	if(MLE) {
		ev <- eigen(C, symmetric = TRUE)
		smallest <- ev$values[ncol(C)]
		if(smallest < lower) {              # check positive definiteness
			out[2] <- -smallest
			return(out)
		}
		else out[2] <- -1.0

		log_det_C <- sum(log(ev$values))    # calculate log determinant
		C_inv <- crossprod(t(ev$vectors) / sqrt(ev$values)) # invert C
	}
	else {
		ev <- eigen(C, TRUE, TRUE)
		smallest <- ev$values[ncol(C)]
		if(smallest < lower) {              # check positive definiteness
			out[2] <- -smallest
			return(out)
		}
		else out[2] <- -1.0
	}

	# calculate all remaining lexical criteria
	if(l <- length(restrictions@criteria) - 1) {
		beta <- coef(restrictions)
		Phi <- cormat(restrictions)
		FC <- beta * (beta %*% Phi)
		for(i in 1:l) {
			environment(restrictions@criteria[[i]]) <- environment()
			out[2L + i] <- restrictions@criteria[[i]]()
		}
	}
	environment(restrictions@criteria[[l+1]]) <- environment()
	out[length(out)] <- restrictions@criteria[[l+1]]()
	return(out)
}

## estimates the coefficents for a 3x3 correlation matrix, not used now
FAiR_triads <- 
function(Phi) {
	cors <- Phi[lower.tri(Phi)]
	triads <- prod(cors) / cors^2
	Lambda <- (sqrt(triads) * ifelse(cors >= 0, 1, -1))[c(2,3,1)]
	return(Lambda)
}

## estimates the coefficents for a 4x4 correlation matrix, not used now
FAiR_tetrads <- 
function(Phi) {
	A <- matrix(NA, nrow = 4, ncol = 3)
	A[1,] <- c(Phi[1,2] * Phi[1,3] / Phi[2,3],
		   Phi[1,2] * Phi[1,4] / Phi[2,4],
		   Phi[1,3] * Phi[1,4] / Phi[3,4])
	A[2,] <- c(Phi[2,1] * Phi[2,3] / Phi[1,3],
		   Phi[2,1] * Phi[2,4] / Phi[1,4],
		   Phi[2,3] * Phi[2,4] / Phi[3,4])
	A[3,] <- c(Phi[3,1] * Phi[3,2] / Phi[1,2],
		   Phi[3,1] * Phi[3,4] / Phi[1,4],
		   Phi[3,2] * Phi[3,4] / Phi[2,4])
	A[4,] <- c(Phi[4,1] * Phi[4,2] / Phi[1,2],
		   Phi[4,1] * Phi[4,3] / Phi[1,3],
		   Phi[4,2] * Phi[4,3] / Phi[2,3])
	## What do I do with A? See which column fits best? Look at Cureton & D'Agostino
	return(A)
}

## numeric gradient
FAiR_numeric_gradient <- 
function(par, restrictions, manifest, helper, lower) {
	gradient.env <- new.env()
	assign("par", par, envir = gradient.env)
	assign("restrictions", restrictions, envir = gradient.env)
	assign("helper", helper, envir = gradient.env)
	assign("manifest", manifest, envir = gradient.env)
	assign("lower", lower, envir = gradient.env)
	gradient <- as.double(attr(numericDeriv(quote(bfgs_fitS4(par, restrictions,
			manifest, helper, lower)), 
			theta = c("par"), gradient.env), "gradient"))
	return(gradient)
}

## makes (inverse) "duplicator" matrix, see Browne (1974)
FAiR_duplicator <- 
function(C, diag = TRUE, inverse = TRUE) { # NOTE: C must be symmetric (duh)
	vec  <- as.vector(C)
	vecs <- vec[c(lower.tri(C, diag))]

	Kp <- matrix(0, nrow = length(vec), ncol = length(vecs))
	intmat  <- matrix(0, nrow = nrow(C), ncol = ncol(C))
	lt <- lower.tri(C,diag)
	intmat[lt] <- 1:length(vecs)
	t_intmat <- t(intmat)
	mark <- 1
	for(j in 1:ncol(C)) for(i in 1:nrow(C)) {
		if(i >= j) Kp[mark,  intmat[mark]] <- 1
		else       Kp[mark,t_intmat[mark]] <- 1
		mark <- mark + 1
	}
	Kp_prime <- round(solve(crossprod(Kp)) %*% t(Kp), 2) # can be only 1, 0, and 0.5
	if(!inverse) return(Kp)
	return(Kp_prime)
}

## makes Gamma_N from Browne (1974)
FAiR_Browne1974 <- 
function(C) {
	Kp_prime <- FAiR_duplicator(C)
	C_chol <- chol(C)
	Gamma_N <- 2 * tcrossprod(Kp_prime %*% t(kronecker(C_chol, C_chol)))
	return(Gamma_N)
}

## outline for an EM algorithm (not complete or used)
FAiR_EM <- 
function() {
# 	stopifnot(require(Design))
	stop("write this function")

	notconverged <- TRUE
	while(notconverged) { # break if a constraint is not met
		## E step
	
		## Get unconstrained beta
	
		## Transform to constrained beta
	
		## Get Phi or Delta and Xi
	
		## Get uniquenesses

		## Check convergence

		## Move temp stuff to real stuff
	}
}

## construct various factor scores
FAiR_scores <- 
function(scores, manifest, beta, Phi, uniquenesses) {
	if(scores == "none") return(matrix(NA_real_, 0, ncol(Phi)))
	else if(!FAiR_is.manifest.data(manifest)) {
		warning("could not calculate factor scores because the raw data on the",
			" outcome variables was not supplied")
		return(matrix(NA_real_, 1, 1))
	}
	zz <- manifest@X
	zz <- sweep(colMeans(zz), 2, manifest@center, "-")
	zz <- sweep(colMeans(zz), 2, manifest@sds,    "/")
	Theta2inv <- diag(1/uniquenesses)
	S <- cormat(manifest)
	S_inv <- chol2inv(chol(S))
	if(scores == "regression") {
		B <- solve(S, beta) %*% Phi
	}
	else if(scores == "Bartlett") { # check if same
		B <- Theta2inv %*% beta %*% solve(t(beta) %*% Theta2inv %*% beta)
	}
	else if(scores == "Thurstone") {
		B <- S_inv %*% beta %*% Phi
	}
	else if(scores == "Ledermann") {
		B <- Theta2inv %*% beta %*% solve( t(beta) %*% Theta2inv %*% beta +
			chol2inv(chol(Phi)) )
	}
	else if(scores == "Anderson-Rubin") {
		B <- Theta2inv %*% beta %*% diag(diag( t(beta) %*% Theta2inv %*% S %*% 
						Theta2inv %*% beta )^(-1/2))
	}
	else if(scores == "McDonald") {
		N <- t(chol(Phi)) # check transposition
		B <- Theta2inv %*% beta %*% N %*% diag(diag(t(N) %*% t(beta) %*% 
			Theta2inv %*% S %*% Theta2inv %*% beta %*% N)^(-1/2)) %*% t(N)
	}
	else if(scores == "Krinjen") {
		Phi_half <- diag(diag(Phi)^(1/2))
		B <- S_inv %*% beta %*% Phi_half %*% diag(diag(Phi_half %*%
			t(beta) %*% S_inv %*% beta %*% Phi_half)^(-1/2)) %*% Phi_half
	}
	else if(scores == "Takeuchi") {
		B <- S_inv %*% beta %*% diag(diag(t(beta) %*% S_inv %*% beta)^(-1/2))
	}
	else if(scores == "Harman") {
		B <- beta %*% chol2inv(chol(crossprod(beta)))
	}
	sc <- zz %*% B
	return(sc)
}

## check whether model "is" identified using a theorem in Reiersol (1950)
FAiR_indeterminator_SEFA <-
function(factors, zeros) {
	if(any(zeros[[1]] < (factors[1] - 1) )) {
		stop("a SEFA model with ", factors[1], " first-order factors",
			" should have at least ", factors[1], " zeros",
			" per first-order factor")
	}
	else if( any(zeros[[1]] < factors[1]) ) {
		warning("a SEFA model with ", factors[1], " first-order factors",
			" should have at least ", factors[1], " zeros",
			" per first-order factor\nThis model depends on inequality",
			" constraints for identification")
	}
	else if(factors[2] <= 1) return(FALSE)
	else if(any(zeros[[2]] < (factors[2] - 1) )) {
		stop("a SEFA model with ", factors[2], " second-order factors",
			" should have at least ", factors[2], " zeros",
			" per second-order factor")
	}
	else if(any(zeros[[2]] < (factors[2] - 1) )) {
		warning("a SEFA model with ", factors[2], " second-order factors",
			" should have at least ",factors[1], " zeros",
			" per second-order factor\nThe second level of this model",
			" depends on inequality constraints for identification")
	}
	return(FALSE)
}

## check whether model is free of rotational indeterminancy using a theorem in Howe (1955)
FAiR_indeterminator_CFA <- 
function(factors, zeros, nonzeros, fixed, fixed2) {
	if(any(zeros[[1]] < (factors[1] - 1) )) {
		stop("a CFA model with ", factors[1], " first-order factors",
			" requires at least ", factors[1] - 1, " specific zeros",
			" per first-order factor")
	}
	else if(any(zeros[[1]] > (nrow(fixed) - 2))) {
		warning("you have less than three nonzero coefficients for at least one",
			"factor at level 1")
	}
	else if(factors[2] == 0) {
		fixed[fixed != 0] <- NA_real_
		uniques <- nrow(unique(t(fixed)))
		if(uniques != ncol(fixed)) {
			stop("no two factors can have all their fixed zeros",
			      " in the same rows at level 1")
		}
	}
	else if(any(zeros[[2]] < (factors[2] - 1) )) {
		stop("a CFA model with ", factors[2], " second-order factors",
			" requires at least ", factors[2] - 1, " specific zeros",
			" per second-order factor")
	}
	else if(any(zeros[[2]] > (nrow(fixed) - 2))) {
		warning("you have less than three nonzero coefficients for at least one",
			"factor at level 2")
	}
	else if(factors[2] > 1) {
		fixed2[fixed2 != 0] <- NA_real_
		uniques <- nrow(unique(t(fixed2)))
		if(uniques != ncol(fixed2)) {
			stop("no two second-order factors can have all their",
			      " fixed zeros in the same rows")
		}
	}
	return(FALSE)
}

## a driver for the show() methods for restrictions
FAiR_show_coef <-
function(object, level) {
	hc <- 1.5
	if(level == 1) {
# 		fix_args <- formals(object@beta@mapping_rule)
		coefs  <- loadings(object)
		select <- object@beta@select
		any_fixed <- matrix(!object@beta@free, nrow(coefs), ncol(coefs))
		dimnames(any_fixed) <- dimnames(coefs)
	}
	else {
# 		fix_args <- formals(object@Delta@mapping_rule)
		coefs  <- loadings(object, level = 2)
		select <- object@Delta@select
		any_fixed <- matrix(!object@Delta@free, nrow(coefs), ncol(coefs))
		dimnames(any_fixed) <- dimnames(coefs)
	}
	if(any(any_fixed)) {
		cat("Fixed coefficients at level ", level, "\n")
		mark <- apply(any_fixed, 1, any)
		fixed <- any_fixed * NA
		fixed[any_fixed] <- coefs[any_fixed]
		print(fixed[mark,,drop = FALSE])
	}
	else cat("No fixed coefficients at level ", level, "\n")

	Domains_coefs <- object@Domains[select,,drop = FALSE]
	any_bounds <- apply(Domains_coefs, 1, FUN = 
			function(x) {
				out <- (x[1] != -hc) | (x[2] != hc)
				return(out)
			})

	if(any(any_bounds)) {
		cat("\nBounded coefficients at level ", level, "\n")
		print(Domains_coefs[any_bounds,,drop = FALSE])
	}
	else cat("\nAll coefficients at level ", level, "in [", -hc, ",", hc, "]\n")
	return(NULL)
}

## a driver for the show() methods for SEFA models
FAiR_show_sefa <- 
function(coef, level) {
	coef_args <- formals(coef@mapping_rule)
	cat("\nZeros per factor at level ", level, "\n")
	if(!is.numeric(coef_args$zeros)) coef_args$zeros <- rep(ncol(coef@x),ncol(coef@x))
	mat <- matrix(coef_args$zeros, nrow = 1)
	rownames(mat) <- "zeros"
	colnames(mat) <- paste("F", 1:ncol(mat), sep = "")
	print(mat)

	if(is.na(coef_args$row_complexity[1])) {
		if(length(coef_args$row_complexity) == 1) {
			row_complexity <- coef_args$row_complexity
			if(is.na(row_complexity)) row_complexity <- ncol(mat)
			cat("All outcomes are of maximum complexity ", 
			    row_complexity, "\n")
		}
		else cat("Outcomes have various complexities\n")
	}
	else {
		factors <- ncol(mat)
		cat("\nStipulations on zeros at level ", level, "\n",
		if(coef_args$communality) paste("\tAt least one zero per factor on a",
						"high communality variable\n") else NULL,
		if(coef_args$quasi_Yates) "\tEncourage cohyperplanarity\n" else NULL,
		if(coef_args$weak_Thurstone)"\tWeak simple structure\n" else NULL,
		if(coef_args$Butler) "\tUnifactorial basis\n" else NULL,
		if(coef_args$viral) paste("\t", 0.5 * factors * (factors - 1),
				  "outcomes each of complexity", factors - 2) else NULL,
		"\n")
	}
	return(NULL)
}

## a driver for the show() methods for restrictions
FAiR_show_constraints <- 
function(criteria, EFA = FALSE) {
	CRITERIA <- c(FAiR_constraints_2nd(), FAiR_constraints_1st(), FAiR_constraints())
	CRITERIA <- sapply(CRITERIA, FUN = function(x) x[[2]])
	l <- length(criteria)
	if(l > 1) {
		cat("\nLexical criteria:\n")
		for(i in 1:(l - 1)) {
			cat(CRITERIA[names(criteria)[i]], "\n")
		}
	}
	if(EFA) cat("Analytic criterion: ",   names(criteria)[l], "\n")
	else    cat("Discrepancy function: ", names(criteria)[l], "\n")
	return(NULL)
}

## the start of a function to (hopefully) pass info to SAGE, not used right now
FAiR_make_symbolic_C <- 
function(beta, Phi, Delta = NULL, Xi = NULL) {
	stopifnot(ncol(beta) == nrow(Phi))
	Z <- matrix("", nrow(beta), ncol(Phi))
	C <- matrix("", nrow(beta), nrow(beta))
	for(i in 1:nrow(beta)) for(j in 1:ncol(beta)) {
		string <- as.character(NULL)
		for(k in 1:ncol(Phi)) {
			ministring <- paste(paste("beta", i, k, sep = ""), "*",
					if(k > j) paste("Phi",  k, j, sep = "") else
						  paste("Phi",  j, k, sep = ""))
			string <- c(string, ministring)
		}
		Z[i,j] <- paste("(", paste(string, collapse = " + "), ")", sep = "")
	}

	for(i in 1:nrow(Z)) for(j in 1:i) {
		string <- as.character(NULL)
		for(k in 1:ncol(Phi)) {
			ministring <- paste(Z[i,k], "*", paste("beta", j, k, sep = ""))
			string <- c(string, ministring)
		}
		C[i,j] <- paste(paste("Sigma", i, j, sep = ""), "== (", 
				paste(string, collapse = " + "), ")")
		if(i == j) {
			for(k in 1:ncol(Phi)) {
				token <- paste("Phi", k, k, sep = "")
				C[i,j] <- sub(token, "0", C[i,j], fixed = TRUE)
			}
		}
	}
	diag(C) <- paste(diag(C), paste("Theta", 1:nrow(C), 1:ncol(C), sep = ""), 
			sep = " + ")
	vars <- c(paste("beta", row(beta), col(beta), sep = ""), 
		paste("Phi",    row(Phi)[lower.tri(Phi, TRUE)], 
				col(Phi)[lower.tri(Phi, TRUE)], sep = ""),
		paste("Theta", 1:nrow(beta), 1:nrow(beta), sep = ""),
		paste("Sigma",  row(C)[lower.tri(C, TRUE)], 
				col(C)[lower.tri(C, TRUE)], sep = ""))
	vars2 <- format(paste(vars, " = var('", vars, "');", sep = ""))
	C <- C[lower.tri(C, TRUE)]
	eqs <- format(paste("eq", 1:length(C), " = ", C, "\n", sep = ""))
	out <- list(vars2 = vars2, eqs = eqs, 
		math = format(paste("eq", 1:length(eqs), ",", sep = "")),
		vars = format(paste(vars, ", ", sep = "")))
	return(out)
}

## estimate 4th order moments unbiasedly, see Browne(1982), not used
FAiR_ADF_unbiased <- function(x, unbiased = TRUE, robust = TRUE, alpha = 3/4) {
	p_star <- 0.5 * ncol(x) * (ncol(x) + 1)
	Gamma <- matrix(0, p_star, p_star)
	lowers <- which(lower.tri(Gamma, TRUE))
	n <- nrow(x)
	n_1 <- n - 1
	count <- 1
	if(robust) {
		covmat <- CovMcd(x, alpha = alpha)
		x <- sweep(x, 2, covmat$center)
	}
	else { 
		covmat <- cov.wt(x)
		x <- sweep(x, 2, colMeans(x))
	}
	for(i in 1:ncol(x)) for(j in i:ncol(x)) {
		s_ij <- covmat$cov[i,j]
		for(k in i:ncol(x)) for(l in k:ncol(x)) {
			if(k == i && l < j) next
			s_kl <- covmat$cov[k,l]
			z <- x[,i] * x[,j] * x[,k] * x[,l]
			if(robust) s_ijkl <- huberM(z, weights = covmat$mcd.wt)$mu
			else       s_ijkl <- mean(z)
			if(unbiased) {
				s_ik <- covmat$cov[i,k]
				s_jl <- covmat$cov[j,l]
				s_il <- covmat$cov[i,l]
				s_jk <- covmat$cov[j,k]
				Gamma[lowers[count]] <- n * (n_1) * (s_ijkl -
				  (n_1 / n)^2 * s_ij * s_kl) /
				  ((n - 2) * (n - 3)) -
				  n_1^2 * (s_ik * s_jl + s_il * s_jk -
				  2 * s_ij * s_kl / n_1) / ( n *
				  (n - 2) * (n - 3) )
			}
			Gamma[lowers[count]] <- s_ijkl - s_ij * s_kl
			count <- count + 1
		}
	}
	Gamma <- Gamma + t(Gamma)
	diag(Gamma) <- diag(Gamma) / 2
	if(nrow(Gamma) > nrow(x)) {
		warning("the ADF estimator can only be inverted when the number of",
			" observations is greater than the number of unique covariances")
	}
	return(Gamma)
}

## estimate sample fourth-order central moments via maximum likelihood
FAiR_ADF <- 
function(x) {
	z <- sweep(x, 2, colMeans(x), "-")
	mark <- lower.tri(tcrossprod(z[1,]), TRUE)
	foo <- function(x) tcrossprod(x)[mark]
	z <- t(apply(z, 1, foo))
	Gamma <- crossprod(Matrix(sweep(z, 2, colMeans(z)) / sqrt(nrow(z))))
	Gamma <- as(Gamma, "dppMatrix")
	if(nrow(Gamma) > nrow(z)) {
		warning("the ADF estimator can only be inverted when the number of",
			" observations is greater than the number of unique covariances")
	}
	return(Gamma)
}

## estimate sample fourth-order moments robustly, not used
FAiR_ADF_robust <-
function(zz, alpha = 3/4, ...) {
	mark <- lower.tri(tcrossprod(zz[1,]), TRUE)
	foo <- function(x) tcrossprod(x)[mark]
	zz  <- sweep(zz, 2, colMeans(zz), "-")
	zz  <- t(apply(zz, 1, foo))
	mcd <- CovMcd(zz, alpha = alpha, ...)
	Gamma <- mcd@cov
	if(nrow(Gamma) > nrow(zz)) {
		warning("the ADF estimator can only be inverted when the number of",
			" observations is greater than the number of unique covariances")
	}
	return(Gamma)
}

## numerically estimate dC_dpar
FAiR_deriv_matrix <- 
function(par, restrictions, manifest, lower, diag = TRUE) {
	myenv <- new.env()
	assign("par", par, envir = myenv)
	assign("restrictions", restrictions, envir = myenv)
	assign("manifest", manifest, envir = myenv)
	assign("lower", lower, envir = myenv)
	fn <- function(par, restrictions, manifest, lower) {
		model <- restrictions2model(par, restrictions, manifest, lower, TRUE)
		if(is.null(model$restrictions)) { # bad
			warning("Boundary solution encountered")
			temp_res <- restrictions
			temp_res@Omega@x[temp_res@scale@free] <-
				exp(par[temp_res@scale@select])
			temp_res@beta@beta[temp_res@beta@free] <-
				par[temp_res@beta@select]
			if(is(temp_res, "restrictions.1storder")) {
				Phi <- temp_res@Phi@x
				Phi[temp_res@Phi@free] <- par[temp_res@Phi@select]
				Phi <- Phi + t(Phi)
				temp_res@Phi@x <- Phi
			}
			else if(is(temp_res, "restrictions.general")) {
				temp_res@Delta@Delta[temp_res@Delta@free] <- 
					par[temp_res@Delta@select]
				Phi <- tcrossprod(temp_res@Delta@Delta)
				diag(Phi) <- 1
				temp_res@Phi@x <- Phi
			}
			else if(is(temp_res, "restrictions.2ndorder")) {
				Xi <- temp_res@Xi@x
				Xi[temp_res@Xi@free] <- par[temp_res@Xi@select]
				Xi <- Xi + t(Xi)
				Delta <- temp_res@Delta@Delta
				Delta[temp_res@Delta@free] <- par[temp_res@Delta@select]
				Phi <- crossprod(chol(Xi) %*% t(Delta))
				diag(Phi) <- 1
				temp_res@Xi@x <- Xi
				temp_res@Delta@Delta <- Delta
				temp_res@Phi@x <- Phi
			}
			model$restrictions <- temp_res
		}
		C <- fitted(model$restrictions, reduced = TRUE)
		diag(C) <- 1
		C <- C * tcrossprod(model$restrictions@Omega@x)
		c <- C[lower.tri(C, diag)]
		return(c)
	}
	assign("fn", fn, envir = myenv)
	gradient <- numericDeriv(quote(fn(par, restrictions, manifest, lower)), 
				c("par"), myenv)
	return(attributes(gradient)$gradient)
}

FAiR_mve <-  ## FIXME, need to put this inside Factanal()
function(beta) {
	combos <- pi ## remove this line
	X <- rbind(diag(ncol(beta)), beta)
	criterion  <- .C("exhaustive_mve", number = as.integer(0), 
			combos = combos,	
			lengthcombos = length(combos), 
			Xdata = as.double(X), 
			Xrow = nrow(X), 
			Xcol =	ncol(X),
			DUP = FALSE, PACKAGE = "FAiR")
	return(-1 + criterion$number / ncol(beta))
}

## analytic derivative of unique elements of C wrt beta (usually)
# FAiR_deriv_C_wrt_coef <- 
# function(coef, cormat, diag = TRUE) { ## keep thinking about how to speed this up
# 	betaPhi <- coef %*% cormat
# 	if(diag) {
# 		gradients <- matrix(0,  nrow = 0.5 * nrow(coef) * (nrow(coef) + 1), 
# 					ncol = length(coef))
# 	}
# 	else gradients <- matrix(0, nrow = 0.5 * nrow(coef) * (nrow(coef) - 1),
# 				    ncol = length(coef))
# 	mark <- 1
# 	for(j in 1:ncol(coef)) for(i in 1:nrow(coef)) {
# 		temp <- matrix(0, nrow(coef), nrow(coef))
# 		temp[,i] <- betaPhi[,j]
# 		temp[i,] <- temp[i,] + betaPhi[,j]
# 		gradients[,mark] <- temp[lower.tri(temp, diag)]
# 		mark <- mark + 1
# 	}
# 	return(gradients)
# }

## analytic derivative of unique elements of C wrt Phi (usually)
# FAiR_deriv_C_wrt_cormat <-
# function(coef, cormat, D = NULL) {
# 	if(is.null(D)) D <- FAiR_duplicator(cormat, inverse = FALSE)
# 	return(kronecker(coef, coef) %*% D)
# }

## check whether x is of class FA or inherits from class FA
FAiR_is.FA <-
function(x) {
	is(x, "FA")
}

## sanity check on equality restrictions
FAiR_check_equalities <-
function(equalities, level) {
	if(length(equalities) == 0) return(list())
	if(!is.list(equalities)) {
		text <- " must be a list of equality_restriction objects"
		if(level == 1) stop("'equalities_1'", text)
		else           stop("'equalities_2'", text)
	}
	else if(!all(sapply(equalities, is, class2 = "equality_restriction"))) {
		begin <- "all elements of "
		end <- " must be equality_restriction objects"
		if(level == 1) stop(begin, 'equalities_1', end)
		else           stop(begin, 'equalities_2', end)
	}
	return(equalities)
}

## driver for show() methods for objects that inherit from class "manifest"
FAiR_show_manifest <-
function(manifest) {
	S <- cormat(manifest)
	N <- manifest@n.obs
	n <- ncol(S)
	cat("Number of observations: ", N, "\n")
	cat("Number of manifest variables: ", n, "\n")
	cat("Proportion of positive correlations: ", mean( S[lower.tri(S)] > 0), "\n")

	Bartlett <- -determinant(S)$modulus * (N - 1 - (2*n+5) / 6 )
	p.val <- pchisq(Bartlett, df = 0.5 * n * (n - 1), lower.tail = FALSE)
	cat("p-value for null hypothesis that manifest variables are uncorrelated: ",
		p.val, "\n")

	AICM <- cov2cor(chol2inv(chol(S)))
	Bartlett <- -determinant(AICM)$modulus * (N - 1 - (2*n+5) / 6 )
	p.val <- pchisq(Bartlett, df = 0.5 * n * (n - 1), lower.tail = FALSE)
	cat("p-value for null hypothesis that the anti-images are uncorrelated: ",
		p.val, "\n")

	MSA <- (sum(S^2) - n) / ( sum(S^2) - n + sum(AICM^2) - n )
	cat("Kaiser-Meyer-Okin Measure of Sampling Adequacy: ", MSA, "\n")
	MSA1 <- colSums(S^2) - 1
	MSA2 <- colSums(AICM^2) - 1
	MSAj <- MSA1 / (MSA1 + MSA2)
	cat("Kaiser-Meyer-Okin Measure of Homogeneity of Each Manifest Variable\n")
	print(as.matrix(MSAj))

	cat("\nEigenvalues of sample correlation matrix\n")
	print(eigen(S, TRUE, TRUE)$value)
}

## arithmetic mean to geometric mean ratio for use in high communality mapping rule
FAiR_AM_GM_ratio <-
function(x) {
	y <- 0
	z <- 1
	for(i in 1:length(x)) {
		if(x[i] < 0) return(NA_real_)
		y <- y + x[i]
		z <- z * x[i]
	}
	y <-    y / length(x)
	z <- z^(1 / length(x))
	return(y/z)
}

## this function shrinks a covariance estimate a la Dey and Srinivasan (1982)
FAiR_shrink <-
function(x, n.obs) {
	if(is.na(n.obs)) {
		stop("'n.obs' must be specified for the shrinkage estimator")
	}
	p <- ncol(x)
	ev <- eigen(x, symmetric = TRUE)

	constant  <- n.obs / (n.obs + p + 1 - 2 * (1:p))
	ev$values <- ev$values * constant
	Sigma <- crossprod(t(ev$vectors) * sqrt(ev$values))
	rownames(Sigma) <- colnames(Sigma) <- rownames(x)
	if(all(diag(x) == 1)) {
		warning("shrinkage may lead to odd behavior when a correlation matrix",
			" is passed")
		Sigma <- cov2cor(Sigma)
	}
	return(Sigma)
}

## hopefully will never be called :) keep it on one line just in case
FAiR_oops <- function(envir) { browser(); ls(envir = envir); stop("oops")}

## see enable_asserts.sh and disable_asserts.sh in FAiR/R (keep as one line!!!)
FAiR_assert <- function(..., warn=TRUE) if(warn) FAiR_warnifnot(...) else stopifnot(...)

## checks whether maximum likelihood was used
FAiR_is.ML <-
function(object) {
	if(FAiR_is.FA(object)) discrepancy <- object@restrictions@discrepancy
	else if(is(object, "restrictions")) discrepancy <- object@discrepancy
	else stop("object must be of class 'FA' or class 'restrictions'")
	return(discrepancy == "MLE")
}

## checks whether raw data are available
FAiR_is.manifest.data <-
function(object) {
	if(FAiR_is.FA(object)) return(is(object@manifest, "manifest.data"))
	else if(is(object, "manifest")) return(is(object, "manifest.data"))
	else stop("object must be of class 'FA' or class 'manifest'")
}

## estimate S with MAR & normality assumptions
FAiR_mlest <-
function(z, wt, bootstrap, shrink, ...) {
	mle <- try(mlest(z))
	if(!is.list(mle)) {
		if(bootstrap >= 0) stop("variance estimator could not be calculated")
		else if(bootstrap < 0) return(list(sigmahat = NA_real_))
	}
	if(is.null(wt)) wt <- rep(1/nrow(z), nrow(z))

	if(shrink) covmat <- FAiR_shrink(mle$sigmahat, nrow(z))
	else       covmat <- mle$sigmahat

	if(bootstrap < 0) return(covmat[lower.tri(covmat, TRUE)])
	else if(bootstrap > 0) {
		holder <- matrix(NA_real_, nrow = bootstrap, 
				ncol = 0.5 * nrow(covmat) * (nrow(covmat) + 1) )
		
		for(i in 1:bootstrap) {
			invalid <- TRUE
			while(invalid) {
				z_rows <- sample(1:nrow(z), nrow(z), replace = TRUE, 
						prob = wt)
				mle_bs <- FAiR_mlest(z[z_rows,], NA, -1, shrink)
				if(is.numeric(mle_bs)) invalid <- FALSE
			}
			holder[i,] <- mle_bs
			if( (i %% 100) == 0 ) print( paste(i * 100 / nrow(holder),
						    "% bootstraps completed") )
		}
		acov <- FAiR_bootstrap2dppMatrix(holder, nrow(z)) # might be wrong
	}
	else {
		acov <- FAiR_Browne1974(covmat) # probably wrong with missing data
		acov <- as(acov, "dppMatrix")
	}
	colnames(covmat) <- rownames(covmat) <- colnames(z)
	manifest <- new("manifest.data", cov = covmat, cor = cov2cor(covmat),
			center = mle$muhat, sds = sqrt(diag(covmat)),
			n.obs = nrow(z), how = "mar", wt = wt,
			acov = acov, X = z, diag = TRUE)
	return(manifest)
}

## robust estimate of S
FAiR_Mcd <-
function(z, wt, bootstrap, shrink, ...) {
	covmat <- CovMcd(z, ...)
	if(is.null(wt)) wt <- rep(1/nrow(z), nrow(z))

	if(shrink) {
		if(any(is.na(covmat@cov))) {
			if(bootstrap >= 0) stop("MCD estimator could not be calculated")
			else if(bootstrap < 0) return(covmat@cov)
		}
		covmat@cov <- FAiR_shrink(covmat@cov, sum(covmat@wt))
	}

	if(bootstrap < 0) return(covmat@cov[lower.tri(covmat@cov, TRUE)])
	else if(bootstrap > 0) {
		holder <- matrix(NA_real_, nrow = bootstrap, 
				ncol = 0.5 * nrow(covmat@cov) * (nrow(covmat@cov) + 1) )
		for(i in 1:bootstrap) {
			invalid <- TRUE
			while(invalid) {
				z_rows <- sample(1:nrow(z), nrow(z), replace = TRUE, 
						prob = wt)
				mcd <- FAiR_Mcd(z[z_rows,], NA, -1, shrink)
				if(all(is.finite(mcd))) invalid <- FALSE
			} 
			holder[i,] <- mcd
			if( (i %% 100) == 0 ) print( paste(i * 100 / nrow(holder),
						    "% bootstraps completed") )
		}
		acov <- FAiR_bootstrap2dppMatrix(holder, nrow(z)) # might be wrong
	}
	else    acov <- FAiR_ADF(z[as.logical(covmat@wt),])
	manifest <- new("manifest.data.mcd", cov = covmat@cov,
			cor = cov2cor(covmat@cov), 
			sds = sqrt(diag(covmat@cov)), center = covmat@center,
			wt = covmat@wt, how = "mcd",
			acov = acov, X = z, mcd = covmat,
			n.obs = as.integer(sum(covmat@wt)), diag = TRUE)
	return(manifest)
}

## optimally shrunk covariance estimate
FAiR_cov.shrink <-
function(z, wt, bootstrap, ...) {
	if(is.null(wt)) wt <- rep(1/nrow(z), nrow(z))
	covmat <- cov.shrink(z, w = wt)

	if(bootstrap < 0) return(covmat[lower.tri(covmat, TRUE)])
	else if(bootstrap > 0) {
		holder <- matrix(NA_real_, nrow = bootstrap, 
				ncol = 0.5 * nrow(covmat) * (nrow(covmat) + 1) )
		wts <- rep(1/nrow(z), nrow(z))
		for(i in 1:bootstrap) {
			z_rows <- sample(1:nrow(z), nrow(z), replace = TRUE, prob = wt)
			holder[i,] <- FAiR_cov.shrink(z[z_rows,], wts, -1)
			if( (i %% 100) == 0 ) print( paste(i * 100 / nrow(holder),
						    "% bootstraps completed") )
		}
		acov <- FAiR_bootstrap2dppMatrix(holder, nrow(z)) # might be wrong
	}
	else {
		acov <- Diagonal(n = 0, x = NA_real_)
# 		acov <- as(acov, "dspMatrix")
		warning("cannot estimate uncertainty without bootstrapping when",
			" how = 'lambda'")
	}
	means <- apply(z, 2, weighted.mean, w = wt)
	covmat <- matrix(c(covmat), nrow(covmat), ncol(covmat))
	colnames(covmat) <- rownames(covmat) <- colnames(z)
	manifest <- new("manifest.data", cov = covmat, cor = cov2cor(covmat),
			sds = sqrt(diag(covmat)), center = means,
			how = "corpcor", n.obs = nrow(z), wt = wt,
			X = z, acov = acov, diag = TRUE)
	return(manifest)
}

## S based on Spearman correlations
FAiR_cor.rank <-
function(z, wt, bootstrap, ...) {
	if(is.null(wt)) wt <- rep(1/nrow(z), nrow(z))
	cormat <- cor(z, method = "spearman")
	if(bootstrap > 0) {
		holder <- matrix(NA_real_, nrow = bootstrap, 
				ncol = 0.5 * ncol(z) * (ncol(z) - 1) )
		mark <- lower.tri(cormat, diag = FALSE)
		for(i in 1:bootstrap) {
			z_rows <- sample(1:nrow(z), nrow(z), replace = TRUE, prob = wt)
			holder[i,] <- cor(z[z_rows,], method = "spearman")[mark]
			if( (i %% 100) == 0 ) print( paste(i * 100 / nrow(holder),
						    "% bootstraps completed") )
		}
		acov <- FAiR_bootstrap2dppMatrix(holder, nrow(z))
	}
	else {
		acov <- Diagonal(n = 0, x = NA_real_)
		warning("currently cannot estimate uncertainty without bootstrapping",
			" when how = 'ranks'")
	}
	means <- rep(0, nrow(cormat))
	sds <- rep(1, nrow(cormat))
	manifest <- new("manifest.data.ranks", cov = cormat, cor = cormat,
			sds = sds, center = means, how = "ranks", n.obs = nrow(z),
			X = as.matrix(z), acov = acov, wt = wt, diag = FALSE)
	return(manifest)
}

## estimate heterogenous correlations
FAiR_hetcor <-
function(z, wt, bootstrap, shrink, ...) {
	if(is.null(wt)) wt <- rep(1/nrow(z), nrow(z))
	cormat <- hetcor(z)
	SEs <- cormat$std.errors
	cormat <- cormat$correlations

	if(shrink) { # shrink like in cor.shrink from library(corpcor)	
		lambda <- sum(SEs^2) / 
			( sum(SEs^2) - ncol(z) )
		if(lambda >= 1) {
			stop("shrinkage parameter = 1, which implies that",
				" all shrunk correlations will be zero")
		}
		cormat <- cormat * (1 - lambda)
		diag(cormat) <- 1
	}

	if(bootstrap > 0) {
		mark <- lower.tri(cormat, diag = FALSE)
		holder <- matrix(NA_real_, nrow = bootstrap, 
				ncol = 0.5 * ncol(z) * (ncol(z) - 1) )
		for(i in 1:ncol(z)) if(is.numeric(z[,i])) z[,i] <- rank(z[,i], "keep")
		for(i in 1:bootstrap) {
			z_rows <- sample(1:nrow(z), nrow(z), replace = TRUE, prob = wt)
			holder[i,] <- 2 * sin(cor(z[z_rows,,drop = FALSE]) * pi/6)[mark]
			if( (i %% 100) == 0 ) print( paste(i * 100 / nrow(holder),
						    "% bootstraps completed") )
		}
		holder <- sweep(holder, 2, cormat[mark] / colMeans(holder), "*")
		acov <- FAiR_bootstrap2dppMatrix(holder, nrow(z))
	}
	else acov <- Diagonal(x = nrow(z) * SEs[lower.tri(SEs)]^2)

	manifest <- new("manifest.data.ordinal", cov = cormat,
			cor = cormat, sds = rep(1,ncol(z)), center = rep(0, ncol(z)),
			how = "hetcor", wt = wt, n.obs = nrow(z),
			X = as.matrix(z), acov = acov, diag = FALSE)
	return(manifest)
}

## makes an object of class "manifest.data" from all numeric data
FAiR_make_manifest.data_numeric <-
function(z, wt = NULL, how = "mcd", bootstrap = 0, shrink = FALSE, seed = 12345, ...) {
	if(ncol(z) < 3) stop("must have at least 3 manifest variables")

	if(!is.null(seed) && (bootstrap > 0 | how == "mcd")) set.seed(seed)
	if(any(is.na(z))) {     # estimate S with MAR & normality assumptions
		if(!require(mvnmle)) stop("the 'mvnmle' package",
			" must be installed when there are missing data")
		manifest <- FAiR_mlest(z, wt, bootstrap, shrink)
		how <- "mar"
	}
	else how <- match.arg(how, c("mcd", "lambda", "unbiased", "mle", "ranks"))

	if(how == "mcd") {    # robust estimate of S
		manifest <- FAiR_Mcd(z, wt, bootstrap, shrink, ...)
	}
	else if(how == "ranks") {  # Spearman correlations
		if(shrink) warning("shrink = TRUE not implemented when how = 'ranks'")
		manifest <- FAiR_cor.rank(z, wt, bootstrap, ...)
	}
	else if(how == "lambda") { # optimally shrunk estimate of S
		if(!require(corpcor)) {
			stop("the suggested package 'corpcor' must be installed")
		}
		if(shrink) {
			warning("shrink = TRUE redundant and ignored when how = 'lambda'")
		}
		manifest <- FAiR_cov.shrink(z, wt, bootstrap, ...)
	}
	else if(how != "mar") { # unbiased or ML estimate of S
		if(bootstrap > 0) warning("bootstrap ignored in this case")
		if(is.null(wt)) wt <- rep(1/nrow(z), nrow(z))
		covmat <- cov.wt(z, cor = FALSE, wt = wt, 
				method = ifelse(how == "unbiased", "unbiased", "ML"))
		if(shrink) covmat$cov <- FAiR_shrink(covmat$cov, nrow(z))
		sds <- sqrt(diag(covmat$cov))
		means <- apply(z, 2, weighted.mean, w = wt)
		acov <- FAiR_ADF(z)
		manifest <- new("manifest.data", cov = covmat$cov,
				cor = cov2cor(covmat$cov), sds = sds, center = means,
				how = how, n.obs = nrow(z), X = z, acov = acov, wt = wt,
				diag = TRUE)
	}

	if(is.null(dimnames(manifest@cov))) {
		warning("no variable names present, arbitrary ones added")
		rownames(manifest@cov) <- colnames(manifest@cov) <- names(manifest@sds) <-
			paste("Y", 1:ncol(manifest@cov), sep = "_")
	}
	else if(length(dimnames(manifest@cov)) == 1) {
		dimnames(manifest@cov) <- c(dimnames(manifest@cov),
					    dimnames(manifest@cov))
	}
	names(manifest@sds) <- rownames(manifest@cov)
	return(manifest)
}

## makes an object of class "manifest.data" involving some ordinal data
FAiR_make_manifest.data_ordinal <-
function(z, wt = NULL, how = "hetcor", bootstrap = 0, shrink = FALSE, ...) { 

	if(how == "ranks") {  # Spearman correlations
		if(shrink) stop("shrink is not implemented when how = 'ranks'")
		for(i in 1:ncol(z)) {
			if(is.factor(z[,i])) z[,i] <- as.integer(z[,i])
			else if(!is.numeric(z[,i])) {
				stop("variable", i, "is neither a factor nor numeric;",
					" please respecify")
			}
		}
		if(!is.null(wt)) wt <- rep(1/nrow(z), rep(z))
		manifest <- FAiR_cor.rank(z, wt, bootstrap, ...)
	}
	else { # bivariate normality assumptions
		if(!require(polycor)) {
			stop("the 'polycor' package must be installed",
			" when some variables are ordinal and how != 'ranks'")
		}
		if(bootstrap) {
			if(any(apply(z, 2, is.numeric))) {
				stop("bootstrapping is not yet implemented when there is",
					" a mix of continuous and ordinal data")
			}
			if(shrink) {
				warning("the behavior when shrink = TRUE, how = 'ranks'",
					"and bootstrap > 0 is unorthodox\n",
					"see ?make_manifest")
			}
		}
		for(i in 1:ncol(z)) {
			if(is.factor(z[,i])) z[,i] <- as.integer(z[,i])
			else if(!is.numeric(z[,i])) stop("variable", i, "is neither a
						factor nor numeric; please respecify")
		}
		if(!is.null(wt)) wt <- rep(1/nrow(z), rep(z))
		manifest <- FAiR_hetcor(z, wt, bootstrap, shrink, ...)
	}
	if(is.null(dimnames(manifest@cov))) {
		warning("no variable names present, arbitrary ones added")
		rownames(manifest@cov) <- colnames(manifest@cov) <-
			paste("Y", 1:ncol(manifest@cov), sep = "_")
	}
	else if(length(dimnames(manifest@cov)) == 1) {
		dimnames(manifest@cov) <- c(dimnames(manifest@cov),
					    dimnames(manifest@cov))
	}
	names(manifest@sds) <- rownames(manifest@cov)
	return(manifest)
}

## make restrictions.factanal or restrictions.orthonormal object; used in make_restricions
FAiR_res_EFA <-
function(factors, manifest, beta, discrepancy, auto) {
	S <- cormat(manifest)
	p <- nrow(S)
	dof <- as.integer(0.5 * ((p - factors[1])^2 - p - factors[1]))
# 	Phi <- diag(factors[1])
# 	attributes(Phi)$ev <- 1
	if(discrepancy  == "MLE") {
		if(auto) algorithm <- 2
		else algorithm <- FAiR_get_algorithm()

		if(algorithm < 3) { # Lawley-Maxwell algorithm
			Domains <- cbind(0, rep(1, p))
			rownames(Domains) <- rownames(S)
			restrictions <- new("restrictions.factanal", 
					factors = as.integer(factors), nvars = p,
					dof = dof, Domains = Domains, free = rep(TRUE, p),
					model = "EFA", discrepancy = "MLE",
					fast = algorithm == 1)
		}
		else {  # CFA with zeros in upper triangle and orthogonal factors
			top <- diag(factors)
			top[lower.tri(top, diag = TRUE)] <- NA_real_
			fixed <- rbind(top, matrix(NA_real_, nrow = p - factors[1],
							     ncol = factors[1]) )
			rownames(fixed) <- rownames(S)
			Domains <- cbind(-1, rep(1, sum(is.na(fixed))))
			Domains[1,1] <- 0
			select <- c(rep(TRUE, nrow(Domains)), rep(FALSE, p))
			Domains <- rbind(Domains, cbind(-18, log(2 * manifest@sds)))
			scale <- new("parameter.scale", x = rep(NA_real_, p), 
				select = !select, free = rep(TRUE, p), num_free = p,
				invalid = 0.0)
			beta <- new("parameter.coef", x = fixed, free = is.na(fixed),
					num_free = as.integer(sum(is.na(fixed))),
					select = select, invalid = 0.0)
			foo <- list(FAiR_discrepancy(discrepancy))
			if(discrepancy == "SHK") {
				criteria <- try(FAiR_constantanator(manifest, "SHK", foo))
				if(!is.list(criteria)) {
					discrepancy <- "ELLIPTICAL"
					foo <- list(FAiR_discrepancy(discrepancy))
					criteria <- FAiR_constantanator(manifest,
									discrepancy, foo)
				}
			}
			else criteria <- FAiR_constantanator(manifest, discrepancy, foo)
			names(criteria)[length(criteria)] <- discrepancy
			restrictions <- new("restrictions.orthonormal",
					factors = as.integer(factors), 
					nvars = nrow(Domains), dof = dof, model = "EFA", 
					Domains = Domains, discrepancy = discrepancy,
					beta = beta, scale = scale, criteria = criteria,
					free = c(beta@free, rep(TRUE, p)))
		}
	}
	else { # YWLS or ADF -> CFA with zeros in upper triange and orthogonality
		top <- diag(factors)
		top[lower.tri(top, diag = TRUE)] <- NA_real_
		fixed <- rbind(top, matrix(NA_real_, nrow = p - factors, 
						     ncol = factors) )
		rownames(fixed) <- rownames(S)
		Domains <- cbind(-1, rep(1, sum(is.na(fixed))))
		Domains[1,1] <- sqrt(.Machine$double.eps)
		select <- c(rep(TRUE, nrow(Domains)), rep(FALSE, p))
		Domains <- rbind(Domains, cbind(-18, log(2 * manifest@sds)))
		scale <- new("parameter.scale", x = rep(NA_real_, p), 
				select = !select, free = rep(TRUE, p), num_free = p,
				invalid = 0.0)
		beta <- new("parameter.coef", x = fixed, free = is.na(fixed),
				num_free = as.integer(sum(is.na(fixed))), select = select,
				invalid = 0.0)
		foo <- list(FAiR_discrepancy(discrepancy))
		if(discrepancy == "SHK") {
			criteria <- try(FAiR_constantanator(manifest, "SHK", foo))
			if(!is.list(criteria)) {
				discrepancy <- "ELLIPTICAL"
				foo <- list(FAiR_discrepancy(discrepancy))
				criteria <- FAiR_constantanator(manifest,
								discrepancy, foo)
			}
		}
		else criteria <- FAiR_constantanator(manifest, discrepancy, foo)
		names(criteria)[length(criteria)] <- discrepancy
		restrictions <- new("restrictions.orthonormal",
				factors = as.integer(factors), nvars = nrow(Domains), 
				dof = dof, model = "EFA", Domains = Domains, 
				discrepancy = discrepancy, beta = beta, scale = scale, 
				criteria = criteria, free = c(beta@free, rep(TRUE, p)))
	}
	return(restrictions)
}

## skeleton of the correctly-sized "fixed" matrices, used in make_restrictions()
FAiR_make_fixed <-
function(fixed, S, factors, levels) {
	if(is.null(fixed)) { # make skeleton fixed matrices for mixed SEFA & CFA
		fixed <- matrix(NA_real_, nrow = nrow(S), ncol = factors[1])
		rownames(fixed) <- rownames(S)
		colnames(fixed) <- paste("Factor_", 1:factors[1], sep = "")

		if(factors[2] > 0) { # skeleton for second level
			fixed2 <- matrix(NA_real_, factors[1], factors[2])
			rownames(fixed2) <- paste("1_Factor_", 1:factors[1], 
							sep = "")
			colnames(fixed2) <- paste("2_Factor_", 1:factors[2], 
							sep = "")
		}
		else fixed2 <- NULL
		FIXED <- list(fixed, fixed2)
	}
	else if(is.list(fixed)) {
		if(levels == 1) {
			stop("you passed 'fixed' as a list but specified a",
				"single-equation model, please respecify")
		}
		fixed2 <- fixed[[2]]
		fixed  <- fixed[[1]]

		if(is.null(fixed)) fixed <- matrix(NA_real_, nrow(S), factors[1])
		if(is.null(fixed2)) {
			stop("if 'fixed' is a list, it must have two elements")
		}

		if(nrow(fixed) != nrow(S)) {
			stop("first element of 'fixed' must have as many rows as",
				" there are manifest variables")
		}
		if(ncol(fixed) != factors[1]) {
			stop("first element of 'fixed' must have as many columns as",
				" are factors at level 1")
		}
		if(is.null(rownames(fixed))) rownames(fixed) <- rownames(S)

		if(nrow(fixed2) != factors[1]) {
			stop("second element of 'fixed' must have as many rows as",
				"factors at level 1")
		}
		else if(ncol(fixed2) != factors[2]) {
			stop("second element of 'fixed' must have as many columns as",
				"are factors are level 2")
		}
		if(is.null(rownames(fixed2))) {
			rownames(fixed2) <- paste("1_Factor_", 1:factors[1], sep = "")
		}
		if(is.null(colnames(fixed2))) {
			colnames(fixed2) <- paste("2_Factor_", 1:factors[2], sep = "")
		}
		FIXED <- list(fixed, fixed2)
	}
	else if(is.matrix(fixed)) {
		if(nrow(fixed) != nrow(S)) {
			stop("'fixed' must have as many rows as there are manifest",
				" variables")
		}
		if(ncol(fixed) != factors[1]) {
			stop("'fixed' must have as many columns as there are",
				" factors at level 1")
		}
		if(is.null(rownames(fixed))) rownames(fixed) <- rownames(S)

		if(factors[2] > 0) { # skeleton for second level
			fixed2 <- matrix(NA_real_, factors[1], factors[2])
			rownames(fixed2) <- paste("1_Factor_", 1:factors[1], sep = "")
			colnames(fixed2) <- paste("2_Factor_", 1:factors[2], sep = "")
		}
		else fixed2 <- NULL
		FIXED <- list(fixed, fixed2)
	}
	else stop("'fixed' must be a matrix, a list of two matrices, or NULL")
	return(FIXED)
}

## skeleton of the correctly-sized "fixed" matrices, used in make_restrictions()
FAiR_make_fixed <-
function(fixed, S, factors, level) {
	if(is.null(fixed)) { # make skeleton fixed matrices for mixed SEFA & CFA
		if(level == 1) {
			fixed <- matrix(NA_real_, nrow = nrow(S), ncol = factors[1])
			rownames(fixed) <- rownames(S)
			colnames(fixed) <- paste("F", 1:factors[1], sep = "")
		}
		else { # skeleton for second level
			fixed <- matrix(NA_real_, factors[1], factors[2])
			rownames(fixed) <- paste("F", 1:factors[1], sep = "")
							
			colnames(fixed) <- paste("G", 1:factors[2], sep = "")
		}
		return(fixed)
	}
	else if(!is.matrix(fixed)) {
		stop('fixed_', level, " must be NULL or a numeric matrix")
	}
	else if(!is.numeric(fixed)) {
		stop('fixed_', level, " must be NULL or a numeric matrix")
	}
	else if(level == 1) {
		if(!all(dim(fixed) == c(nrow(S), factors[1]))) {
			stop("the dimensions of 'fixed' must be ", nrow(S), " by ",
				factors[1])
		}
		if(is.null(rownames(fixed))) rownames(fixed) <- rownames(S)
		if(is.null(colnames(fixed))) {
			colnames(fixed) <- paste("F", 1:factors[1], sep = "")
		}
		return(fixed)
	}
	else if(level == 2) {
		if(!all(dim(fixed) == factors)) {
			stop("the dimensions of 'fixed' must be ", factors[1], " by ",
				factors[2])
		}
		if(is.null(rownames(fixed))) {
			rownames(fixed) <- paste("F", 1:factors[1], sep = "")
		}
		if(is.null(colnames(fixed))) {
			colnames(fixed) <- paste("G", 1:factors[2], sep = "")
		}
		return(fixed)
	}
}

## calculates the variance-covariance matrix of the estimates and z-statistics
FAiR_uncertainty <-
function(restrictions, manifest, factors, sorter, sorter_2nd = NULL, analytic) {
	discrepancy <- restrictions@discrepancy
	if(is.na(manifest@n.obs)) {
		n_star <- 0.5 * nrow(manifest@cov) * (nrow(manifest@cov) + 1)
		derivs <- matrix(NA_real_, nrow = n_star, ncol = restrictions@nvars)
		vcov <- matrix(NA_real_, restrictions@nvars, restrictions@nvars)
		warning("could not calculate measures of uncertainty because the",
			" number of observations is unknown")
		return(list(Jacobian = derivs, vcov = vcov))
	}
	else if(discrepancy == "YWLS") {
		n_star <- 0.5 * nrow(manifest@cov) * (nrow(manifest@cov) + 1)
		derivs <- matrix(NA_real_, nrow = n_star, ncol = restrictions@nvars)
		vcov <- matrix(NA_real_, restrictions@nvars, restrictions@nvars)
		warning("cannot calculate measures of uncertainty when Yates",
			" weighted least squares criterion is used")
		return(list(Jacobian = derivs, vcov = vcov))
	}

	R <- fitted(restrictions, reduced = FALSE)
	scale <- restrictions@Omega@x
	C <- R * tcrossprod(scale)
	par <- FAiR_restrictions2par(restrictions)
	lower <- sqrt(.Machine$double.eps)
	if(FAiR_is.ML(restrictions)) diag <- TRUE
	else diag <- manifest@diag
	derivs <- FAiR_deriv_matrix(par, restrictions, manifest, lower, diag)
	scaleR <- scale * R
	Jacobian_scale <- matrix(NA_real_, nrow = nrow(derivs), ncol = ncol(R))
	mark <- lower.tri(R, TRUE)
	for(i in 1:ncol(R)) {
		J <- matrix(0L, ncol(R), ncol(R))
		J[i,i] <- 1L
		scaleRJ <- scaleR %*% J
		Jacobian_scale[,i] <- scale[i] * (scaleRJ + t(scaleRJ))[mark]
	}
	Jacobian_scale <- Jacobian_scale[,restrictions@Omega@free]
	derivs[,(ncol(derivs) - ncol(Jacobian_scale) + 1):ncol(derivs)] <- Jacobian_scale
	mark <- (par != 0) | restrictions@Omega@select
	derivs <- derivs[,mark]
	parnames <- names(par[mark])

	if(discrepancy == "MLE" | (FAiR_is.manifest.data(manifest) && 
					is(manifest@acov, "diagonalMatrix"))) {
		Gamma_N <- FAiR_Browne1974(C)

		W <- chol2inv(chol(Gamma_N))
		W_chol <- chol(W)
		dWd_inv <- chol2inv(chol((crossprod(W_chol %*% derivs))))
		if(FAiR_is.manifest.data(manifest) && 
		       !is(manifest@acov, "diagonalMatrix") && 
			manifest@n.obs > nrow(manifest@acov)) {
			Gamma_chol <- chol(manifest@acov)
			middle <- chol(crossprod(Gamma_chol %*% W %*% derivs))
			vcov <- crossprod(middle %*% dWd_inv)
			vcov <- as.matrix(vcov) / manifest@n.obs
		}
		else    vcov <- dWd_inv / manifest@n.obs
	}
	else if(FAiR_is.QD(restrictions)) {
		foo <- restrictions@criteria[[length(restrictions@criteria)]]
		middle <- formals(foo)$middle
		t_W_chol <- matrix(0, nrow(derivs), nrow(derivs))
		t_W_chol[lower.tri(t_W_chol, TRUE)] <- middle
		vcov <- chol2inv(chol(tcrossprod(t(derivs) %*% t_W_chol)))
		vcov <- vcov / manifest@n.obs
	}
	else {
		n_star <- 0.5 * nrow(manifest@cov) * (nrow(manifest@cov) + 1)
		derivs <- matrix(NA_real_, nrow = n_star, ncol = restrictions@nvars)
		vcov <- matrix(NA_real_, restrictions@nvars, restrictions@nvars)
		warning("could not calculate measures of uncertainty because the",
			" discrepancy function is not recognized")
		return(list(Jacobian = derivs, vcov = vcov))
	}
	rownames(vcov) <- colnames(vcov) <- parnames
	return(list(Jacobian = derivs, vcov = vcov))
}

## make bootstrapped covariance matrix of covariances into a dppMatrix
FAiR_bootstrap2dppMatrix <-
function(x, n.obs) {
	Gamma <- crossprod(Matrix(sweep(x, 2, colMeans(x)) / n.obs))
	Gamma <- as(Gamma, "dppMatrix")
	return(Gamma)
}

## make some stuff that is common to almost all objects that inherit from class FA
FAiR_restrictions2FA <-
function(restrictions, manifest, scores) {
	# extract the basics
	factors <- restrictions@factors
	Phi  <- cormat(restrictions)
	beta <-   coef(restrictions)
	uniquenesses <- uniquenesses(restrictions)

	# construct primary structure and factor contributions
	Pi <- beta %*% Phi
	FC <- beta * Pi
	names(uniquenesses) <- rownames(beta)
	sorter <- order(colSums(FC), decreasing = TRUE)
	sorter <- 1:ncol(beta) # temporary
	# construct reference structure and reference pattern
	Phi_inv <- chol2inv(chol(Phi))
	D <- 1/sqrt(diag(Phi_inv))
	Psi <- Phi_inv * tcrossprod(D)
	Upsilon <- sweep(beta, 2, D, FUN = "*")
	RP <- sweep(Pi, 2, D, FUN = "/")

	# bundle everything up
	loadings <- array(cbind(beta, Upsilon, Pi, RP, FC),
			dim = c(nrow(beta), factors[1], 5),
			dimnames = list(rownames(beta), NULL, 
					c("PP", "RS", "PS", "RP", "FC")))
	loadings <- loadings[,sorter,,drop = FALSE]
# 	if(ncol(Phi) > 1) {
		correlations <- array(cbind(Phi, Psi, diag(D)),
				dim = c(factors[1], factors[1], 3), 
				dimnames = list(NULL, NULL, c("PF", "RF", "PR")))
		correlations <- correlations[sorter,sorter,,drop = FALSE]
		attributes(correlations)$orthogonal <- restrictions@model == "EFA"
# 	}
# 	else {
# 		correlations <- array(1, c(1,1,3), dimnames = list(NULL, NULL, 
# 					c("PF", "RF", "PR")))
# 	}
	rownames(correlations) <- colnames(correlations) <- colnames(loadings) <- 
							paste("F", 1:factors[1], sep = "")
	
	scores <- FAiR_scores(scores, manifest, beta, Phi, uniquenesses)
	out <- list(loadings = loadings, correlations = correlations, sorter = sorter,
			uniquenesses = uniquenesses, scores = scores)
			
	return(out)
}

## function checks whether analytic gradients are possible
FAiR_analytic_safe <-
function(restrictions) {
	if(is(restrictions, "restrictions.independent"))             return(TRUE)
	if(restrictions@discrepancy == "IRLS")                       return(FALSE)
	if(FAiR_is.EFA(restrictions))                                return(TRUE)
	if(is(restrictions, "restrictions.general")) {
		if(is(restrictions@Delta, "parameter.coef.nl"))      return(FALSE)
		if(is(restrictions@Delta, "parameter.coef.SEFA.nl")) return(FALSE)
	}
	if(is(restrictions@beta, "parameter.coef.nl"))               return(FALSE)
	if(is(restrictions@beta, "parameter.coef.SEFA.nl"))          return(FALSE)
	return(TRUE)
}

## fills in formals() of discrepancy functions
FAiR_constantanator <-
function(manifest, discrepancy, criteria) {
	S <- model.matrix(manifest, standardized = FALSE)
	     if(discrepancy == "YWLS") return(criteria)
	else if(discrepancy == "MLE") {
		formals(criteria[[length(criteria)]])$constant <- 
							determinant(S)$modulus + nrow(S)
		return(criteria)
	}
	else if(discrepancy == "IRLS") {
		formals(criteria[[length(criteria)]])$s <- S[lower.tri(S, TRUE)]
		return(criteria)
	}

	     if(discrepancy == "ELLIPTICAL")   Gamma <- FAiR_elliptical(manifest)
	else if(discrepancy == "HK")           Gamma <- FAiR_HK(manifest, shrink = FALSE)
	else if(discrepancy == "SHK")          Gamma <- FAiR_HK(manifest, shrink = TRUE )
	else                                   Gamma <- manifest@acov

	if(eigen(Gamma, TRUE, TRUE)$values[ncol(Gamma)] < sqrt(.Machine$double.eps)) {
		warning("4th order moment matrix is not positive definite and",
			" something bad will probably happen")
	}
	diag <- manifest@diag
	mark <- lower.tri(S, diag)
	s    <- S[mark]
	if(is(Gamma, "diagonalMatrix")) {
		middle <- Gamma@x #/ manifest@n.obs # why were you dividing?
		middle <- middle^(-0.5)
	}
	else {
		middle <- solve(Gamma)
		middle <- as.matrix(t(chol(middle)))
		middle <- middle[lower.tri(middle, TRUE)]
	}
	formals(criteria[[length(criteria)]])$diag <- diag
	formals(criteria[[length(criteria)]])$middle <- middle
	formals(criteria[[length(criteria)]])$s <- s
	return(criteria)
}

## convert to 4th order moments under an elliptical distribution
FAiR_elliptical <- function(manifest) {
	S <- model.matrix(manifest, standardized = FALSE)
	p <- ncol(S)
	p_star <- 0.5 * p * (p + 1)
	X  <- sweep(manifest@X, 2, manifest@center)
	S_inv <- chol2inv(chol(S))
	k1 <- mean(apply(X, 1, FUN = function(x) (x %*% S_inv %*% x)^2)) / ( p * (p + 2) )

	Gamma  <- matrix(0, p_star, p_star)
	lowers <- which(lower.tri(Gamma, TRUE))
	count  <- 1
	for(i in 1:ncol(X)) for(j in i:ncol(X)) {
		s_ij <- S[i,j]
		for(k in i:ncol(X)) {
			s_ik <- S[i,k]
			s_jk <- S[j,k]
			for(l in k:ncol(X)) {
				if(k == i && l < j) next
				s_kl <- S[k,l]
				s_jl <- S[j,l]
				s_il <- S[i,l]
				Gamma[lowers[count]] <- k1 * (s_ij * s_kl  + s_ik * s_jl +
							      s_il * s_jk) - s_ij * s_kl
				count <- count + 1
			}
		}
	}
	Gamma <- Gamma +  t(Gamma)
	diag(Gamma) <- diag(Gamma) / 2
	return(Gamma)
}

## convert to 4th order moments under an HK distribution
FAiR_HK <- function(manifest, shrink = FALSE) {
	S <- model.matrix(manifest, standardized = FALSE)
	p <- ncol(S)
	p_star <- 0.5 * p * (p + 1)
	N <- manifest@n.obs
	X <- sweep(manifest@X, 2, manifest@center, FUN = "-")
	K <- apply(X, 2, FUN = function(x) sqrt(N * sum(x^4) / (3 * crossprod(x)^2)))
	if(shrink) {
		var_k <- 6 * N * (N - 1) / ( (N - 2) * (N + 1) * (N + 3) ) *
			 4 *   (N^2 - 1) / ( (N - 3) * (N + 5) * 9)
		deriv <- 0.5 / K
		var_k <- sum(var_k * deriv^2)
# 		KK <- K^2
		med_K <- median(K)
		lambda <- var_k / sum( (K - med_K)^2 )
		if(lambda > 1) {
			stop("shrinkage parameter > 1.0, switching automatically to ",
				"discrepancy = 'ELLIPTICAL'")
		}
		else cat("shrinkage parameter is ", lambda, "\n")
		K <- lambda * median(K) + (1 - lambda) * K
	}
	Gamma  <- matrix(0, p_star, p_star)
	lowers <- which(lower.tri(Gamma, TRUE))
	count  <- 1
	for(i in 1:ncol(X)) for(j in i:ncol(X)) {
		s_ij <- S[i,j]
		a_ij <- (K[i] + K[j]) / 2
		for(k in i:ncol(X)) {
			s_jk <- S[j,k]
			s_ik <- S[i,k]
			a_jk <- (K[j] + K[k]) / 2
			a_ik <- (K[i] + K[k]) / 2
			for(l in k:ncol(X)) {
				if(k == i && l < j) next
				s_kl <- S[k,l]
				s_jl <- S[j,l]
				s_il <- S[i,l]
				a_kl <- (K[k] + K[l]) / 2
				a_jl <- (K[j] + K[l]) / 2
				a_il <- (K[i] + K[l]) / 2
				Gamma[lowers[count]] <- a_ij * a_kl * s_ij * s_kl +
							a_ik * a_jl * s_ik * s_jl +
							a_il * a_jk * s_il * s_jk -
							s_ij * s_kl
				count <- count + 1
			}
		}
	}
	Gamma <- Gamma +  t(Gamma)
	diag(Gamma) <- diag(Gamma) / 2
	return(Gamma)
}

## checks for EFAness
FAiR_is.EFA <-
function(object) {
	if(FAiR_is.FA(object)) return(is(object, "FA.EFA"))
	else if(is(object, "restrictions")) return(object@model == "EFA")
	else stop("object must be of class 'FA' or class 'restrictions'")
}

## checks for orthogonal factors
FAiR_is.orthogonal <-
function(object) {
	if(!FAiR_is.FA(object)) stop("'object' must be of class 'FA'")
	else if(FAiR_is.EFA(object)) return(!object@rotated)
	else return(FALSE)
}

## constructor for the criteria to be used by Rotate()
FAiR_make_criteria <- 
function(FAobject, criteria, methodArgs) {
	CRITERIA <- c(FAiR_constraints(), FAiR_constraints_1st())
	CRITERIA[["cohyperplanarity_1st"]] <- CRITERIA[["dist_cols_1st"]] <-
		CRITERIA[["volume_1st"]] <- NULL
	if(l <- length(criteria)) {
		CRITERIA_names <- names(CRITERIA)
		l <- length(criteria)
		if(l > 1) for(i in 1:(l-1)) {
			if(is.character(criteria[[i]])) {
				foo <- match.arg(criteria[[i]], CRITERIA_names)
				criteria[[i]] <- CRITERIA[[foo]][[1]]
				names(criteria)[i] <- foo
			}
			else if(!is.function(criteria[[i]])) {
				stop("each element of 'criteria' should be a function or",
					" a character string naming a function\nit is",
					" usually best to leave 'criteria' unspecified")
			}
			else names(criteria)[i] <- paste("fun", i, sep = "_")
		}
		if(is.character(criteria[[l]])) {
			analytics <- FAiR_analytic()
			analytics_names <- names(analytics)
			foo <- criteria[[l]]
			if(foo %in% analytics_names) {
				criteria[[l]] <- analytics[[foo]][[1]]
			}
			names(criteria)[l] <- foo
		}
		else if(!is.function(criteria[[l]])) {
			stop("each element of 'criteria' should be a function or",
				" a character string naming a function\nit is",
				" usually best to leave 'criteria' unspecified")
		}
		else names(criteria)[l] <- paste("fun", l, sep = "_")

		if(names(criteria)[1] != "no_factor_collapse") {
			criteria <- c(CRITERIA[["no_factor_collapse"]][[1]], criteria)
			names(criteria)[1] <- "no_factor_collapse"
		}
	}
	else criteria <- FAiR_make_criteria_GUI()
	criteria <- FAiR_fill_methodArgs_Rotate(FAobject, criteria, methodArgs)
	return(criteria)
}

## sanity check for target matrix
FAiR_check_target <-
function(Lambda, target) {
	if(!is.matrix(target)) {
		stop("'target' must be a matrix")
	}
	else if(!is.numeric(target)) {
		stop("'target' must be a numeric matrix")
	}
	else if(!identical(dim(Lambda), dim(target))) {
		stop("'target' must have the same dimensions as the loadings matrix")
	}
	else if(any(is.na(target))) {
		stop("'target' must be fully specified and not contain any NAs")
	}
	return(TRUE)
}

## sanity check for partially-specified target matrix
FAiR_check_pst <-
function(Lambda, target) {
	if(!is.matrix(target)) {
		stop("'target' must be a matrix")
	}
	else if(!is.numeric(target)) {
		stop("'target' must be a numeric matrix")
	}
	else if(!identical(dim(Lambda), dim(target))) {
		stop("'target' must have the same dimensions as the loadings matrix")
	}
	return(TRUE)
}

## make confidence intervals from simulations
FAiR_draws2CI <-
function(draws, bounds) {
	out <- list()
	for(i in 1:length(draws)) {
		arr <- array(NA_real_, c(dim(draws[[i]])[1:2], 2))
		rownames(arr) <- rownames(draws[[i]])
		colnames(arr) <- colnames(draws[[i]])
		dimnames(arr)[[3]] <- paste(bounds, "%")
		arr[,,1] <- apply(draws[[i]], 1:2, quantile, probs = bounds[1])
		arr[,,2] <- apply(draws[[i]], 1:2, quantile, probs = bounds[2])
		out[[i]] <- arr
	}
	names(out) <- names(draws)
	return(out)
}

## objective function used by optim() during rotation in restrictions2draws method
FAiR_optim_Rotate <-
function(par, A, criteria) {
	fits <- FAiR_Rotate(par, A, criteria)
	if(which(fits != -1)[1] < length(fits)) return(NA)
	else return(fits[length(fits)])
}

## rotate quickly with help from library(GPArotation)
FAiR_quick_Rotate <-
function(restrictions, Target) {
	if(restrictions@orthogonal) {
		GPA <- targetT(restrictions@Lambda, restrictions@Tmat, Target)
		out <- GPForth(restrictions@Lambda, Tmat = GPA$Th, method = method)$Th
		out <- FAiR_make_Tmat(opt$par)
		return(out)
	}
	GPA <- try(targetQ(restrictions@Lambda, restrictions@Tmat, Target), silent = TRUE)
	if(!is.list(GPA)) GPA <- list(Th = restrictions@Tmat)
	if(is.character(method <- restrictions@Tcriteria[[1]])) { # from GPArotation
		out <- GPFoblq(restrictions@Lambda, Tmat = GPA$Th, method = method)$Th
		out <- FAiR_make_Tmat(opt$par)
		return(out)
	}
	par  <- FAiR_Tmat2par(GPA$Th)
	fits <- FAiR_Rotate(par, restrictions@Lambda, restrictions@Tcriteria)
	if(which(fits != -1)[1] < length(fits)) return(NA_real_)
	opt <- try(optim(par, fn = FAiR_optim_Rotate, method = "Nelder-Mead", 
			A = restrictions@Lambda, criteria = restrictions@Tcriteria),
			silent = TRUE)
	if(!is.list(opt)) return(NA_real_)
	out <- FAiR_make_Tmat(opt$par)
	return(out)
}

## checks whether ADF discrepancy function (or special case) was used
FAiR_is.QD <-
function(object) {
	if(FAiR_is.FA(object)) object <- object@restrictions
	if(is(object, "restrictions")) {
		return(object@discrepancy %in% c("ADF", "DWLS", "HK", "SHK",
						"ELLIPTICAL"))
	}
	else stop("'object' must inherit from class 'FA' or from class 'restrictions'")
}

## add fuzz to a criterion
FAiR_uniquify <-
function(x, bigger = TRUE) {
	k <- .Machine$double.eps
	noise <- runif(1, min = k, max = k * 10)
	if(bigger) x <- x + noise
	else       x <- x - noise
	return(x)
}

## convert from a restrictions object to a vector of parameters
FAiR_restrictions2par <-
function(restrictions) {
	if(is(restrictions, "restrictions.factanal")) {
		stop("this should not have happened")
	}
	beta <- coef(restrictions)
	par <- c(beta, log(restrictions@Omega@x))

	if(is(restrictions, "restrictions.2ndorder")) {
		Xi <-  cormat(restrictions, level = 2)
		Delta <- loadings(restrictions, level = 2)
		par <- c(Xi[lower.tri(Xi)], Delta, par) 
	}
	else if(is(restrictions, "restrictions.general")) {
		Delta <- loadings(restrictions, level = 2)
		par <- c(Delta, par)
	}
	else if(is(restrictions, "restrictions.1storder")) {
		Phi <- cormat(restrictions)
		par <- c(Phi[lower.tri(Phi)], par)
	}
	par <- par[restrictions@free]
	names(par) <- rownames(restrictions@Domains)
	return(par)
}

## make inverse of Gamma matrix
FAiR_make_W <-
function(FAobject) {
	if(!FAiR_is.FA(FAobject)) {
		stop("'FAobject' must be of class 'FA'")
	}
	W <- array(0, dim = dim(FAobject@manifest@acov))
	l <- length(FAobject@restrictions@criteria)
	middle <- formals(FAobject@restrictions@criteria[[l]])$middle
	if(is(FAobject@manifest@acov, "diagonalMatrix")) diag(W) <- middle^2
	else {
		W[lower.tri(W, TRUE)] <- middle
		W <- tcrossprod(W)
	}
	return(W)
}

## Do every postestimation command with every option
FAiR_stress_test <-
function(FAobject) {
	if(require(Rgraphviz)) {
		plot(FAobject)
		dev.off()
	}
	show(FAobject)
	summary(FAobject)
# 	plot(summary(FAobject)) # why is this commented out?
# 	dev.off()
	summary(FAobject, standardized = FALSE)
	summary(FAobject, nsim = 1001)
	summary(FAobject, nsim = 1001, standardized = FALSE)
	draws <- simulate(FAobject, nsim = 1001)
	draws <- simulate(FAobject, nsim = 1001, standardized = FALSE)
	fitted(FAobject, reduced = TRUE, standardized = TRUE)
	fitted(FAobject, reduced = FALSE, standardized = TRUE)
	fitted(FAobject, reduced = TRUE, standardized = FALSE)
	fitted(FAobject, reduced = FALSE, standardized = FALSE)
	if(!FAiR_is.QD(FAobject)) influence(FAobject)
	model.matrix(FAobject, standardized = FALSE)
	model.matrix(FAobject, standardized = TRUE)
	pairs(FAobject)
	dev.off()
	residuals(FAobject, standardized = TRUE)
	residuals(FAobject, standardized = FALSE)
	rstandard(FAobject)
	if(!FAiR_is.QD(FAobject)) weights(FAobject)
	coef(FAobject)
	PF <- cormat(FAobject, matrix = "PF")
	RF <- cormat(FAobject, matrix = "RF")
	D <- cormat(FAobject, matrix = "PR")
	stopifnot(all.equal(1/sqrt(diag(solve(RF))), diag(D)))
	PP <- loadings(FAobject, matrix = "PP")
	RS <- loadings(FAobject, matrix = "RS")
	stopifnot(all.equal(PP %*% D, RS[,1:ncol(RS)]))
	PS <- loadings(FAobject, matrix = "PS")
	RP <- loadings(FAobject, matrix = "RP")
	stopifnot(all.equal(PS[,1:ncol(PS)], RP %*% D))
	FC <- loadings(FAobject, matrix = "FC")
	stopifnot(all.equal(FC, RS * RP))
	stopifnot(all.equal(FC, PS * PP))
	if(is(FAobject, "FA.2ndorder")) {
		PP <- loadings(FAobject, matrix = "PP", level = 2)
		PS <- loadings(FAobject, matrix = "PS", level = 2)
		RP <- loadings(FAobject, matrix = "RP", level = 2)
		RS <- loadings(FAobject, matrix = "RS", level = 2)
		FC <- loadings(FAobject, matrix = "FC", level = 2)
		stopifnot(all.equal(FC, RP * RS))
		stopifnot(all.equal(FC, PP * PS))
		stopifnot(all(uniquenesses(FAobject, level = 2) >= 0))
	}
	else if(is(FAobject, "FA.general")) {
		PP <- loadings(FAobject, matrix = "PP", level = 2)
		PS <- loadings(FAobject, matrix = "PS", level = 2)
		stopifnot(all.equal(PP, PS))
		RP <- loadings(FAobject, matrix = "RP", level = 2)
		RS <- loadings(FAobject, matrix = "RS", level = 2)
		stopifnot(all.equal(RP, RS))
		FC <- loadings(FAobject, matrix = "FC", level = 2)
		stopifnot(all.equal(FC, RP * RS))
		stopifnot(all.equal(FC, PP * PS))
		stopifnot(all(uniquenesses(FAobject, level = 2) >= 0))
	}
	uniquenesses(FAobject, standardized = FALSE)
	uniquenesses(FAobject, standardized = TRUE)
	stopifnot(all(sqrt(diag(vcov(FAobject))) > 0))
	logLik(FAobject)
	BIC(FAobject)
	CI <- confint(FAobject, level = .90)
	if(!FAiR_is.EFA(FAobject)) profile(FAobject, plot.it = FALSE)
	screeplot(FAobject)
	dev.off()
	model_comparison(FAobject)
	return(TRUE)
}

## make restrictions.independent
FAiR_make_independent <-
function(restrictions) {
	new_nvars <- as.integer(restrictions@Omega@num_free)
	new_dof <- restrictions@dof + restrictions@nvars - new_nvars
	new_criteria <- restrictions@criteria[length(restrictions@criteria)]
	new_Omega <- restrictions@Omega
	old_select <- new_Omega@select
	new_Omega@select <- new_Omega@free
	new_restrictions <- new("restrictions.independent", factors = c(0L,0L),
				nvars = new_nvars, dof = new_dof, Omega = new_Omega,
				Domains = restrictions@Domains[old_select,],
				model = "CFA", discrepancy = restrictions@discrepancy, 
				free = restrictions@Omega@free, criteria = new_criteria)
	return(new_restrictions)
}

## internal function for make_manifest() with covmat specified
FAiR_make_manifest_list <-
function(covmat, shrink) { # covmat is a *list*
	how <- "unknown"
	if(is.null(covmat$cov)) {
		stop("if 'covmat' is a list, it must contain an element",
			" named 'cov' that is a covariance matrix")
	}
	else if(!isTRUE(all.equal(covmat$cov, t(covmat$cov)))) {
		stop("covariance matrix is not symmetric")
	}
	else if(nrow(covmat$cov) < 3) {
		stop("order of covariance matrix must be at least 3")
	}

	warning("it is strongly preferable to pass the raw data to make_manifest()")

	if(is.null(covmat$n.obs)) {
		warning("you will not be able to calculate",
		" measures of uncertainty unless you specify 'n.obs'")
	}

	if(is.null(colnames(covmat$cov))) {
		if(is.null(rownames(covmat$cov))) {
			warning("no variable names in 'covmat', made some up")
			rownames(covmat$cov) <- colnames(covmat$cov) <- 
				paste("V", 1:ncol(covmat$cov), sep = "")
		}
		else colnames(covmat$cov) <- rownames(covmat$cov)
	}
	else rownames(covmat$cov) <- colnames(covmat$cov)

	if(!is.null(covmat$sds)) {
		if(any(diag(covmat$cov) != 1L)) {
			stop("if 'sds' is specified, the 'cov' element of the covariance",
				" list must have 1.0 for each of its diagonal elements")
		}
		else if(length(sds <- covmat$sds) != ncol(covmat$cov)) {
			stop("the length of 'sds' must be the same as the number of",
				" manifest variables")
		}
		else if(any(sds <= 0)) {
			stop("each element of 'sds' must be positive")
		}
		covmat$cov <- covmat$cov * tcrossprod(sds)
	}
	else sds <- sqrt(diag(covmat$cov))

	diag <- any(sds != 1L)


	if(!is.null(covmat$W)) {
		if(!is.matrix(W <- covmat$W)) {
			stop("if 'W' is supplied it must be a positive definite matrix")
		}
		else if(nrow(W) != ncol(W)) {
			stop("if 'W' is supplied it must be a positive definite matrix")
		}
		else if(!isTRUE(all.equal(W, t(W)))) {
			stop("'W' is not symmetric")
		}

		acov <- try(chol2inv(chol(W)))
		if(!is.matrix(acov)) {
			stop("'W' is not positive definite")
		}
		else if(all(c(!acov[lower.tri(acov)], !acov[upper.tri(acov)]))) {
			acov <- Diagonal(x = diag(acov))
		}
		else acov <- as(acov, "dppMatrix")

		if(shrink) {
			warning("setting 'shrink = TRUE' is inappropriate when passing",
				" a user-specified weight matrix (W)")
		}

		manifest <- new("manifest.basic.userW", cov = covmat$cov, 
				cor = cov2cor(covmat$cov), sds = sds, 
				center = as.numeric(covmat$center), 
				n.obs = as.integer(covmat$n.obs), how = how,
				diag = diag, acov = acov)
		return(manifest)
	}
	if(shrink) covmat$cov <- FAiR_shrink(covmat$cov, covmat$n.obs)
	manifest <- new("manifest.basic", cov = covmat$cov,
			cor = cov2cor(covmat$cov), sds = sds,
			center = as.numeric(covmat$center),
			n.obs = as.integer(covmat$n.obs), how = how, diag = diag)
	return(manifest)
}

## copy an object of class restrictions
FAiR_copy_restrictions <-
function(restrictions) {
	sn <- slotNames(restrictions)
	expr <- paste("new('", class(restrictions), "', ", 
			paste(sn, " = ", "restrictions@", sn, sep = "", collapse = ", "), 
			")", sep = "")
	return(eval(parse(text = expr)))
}

## helper for profile() at level 1
FAiR_profile <-
function(fitted, delta, number, plot.it, ...) {
	beta <- coef(fitted@restrictions)
	notfree <-  !fitted@restrictions@beta@free | !beta

	width <- seq(from = -delta, to = delta, length = number)
	outlist <- list()
	count <- 1
	parnames <- rownames(beta)
	if(plot.it) {
		ask <- par("ask")
		par("ask" = TRUE)
		on.exit(par("ask" = ask))
		cat("Factors may arbitrarily be plotted in a different order than they",
			" appear in summary()\n")
	}
	for(p in 1:ncol(notfree)) for(j in 1:nrow(notfree)) if(notfree[j,p]) {
		y <- x <- rep(NA_real_, length(width))
		val <- beta[j,p]
		for(i in 1:length(x)) {
			x[i] <- fitted@restrictions@beta@x[j,p] <- val + width[i]
			y[i] <- deviance(fitted)
		}
		outlist[[count]] <- list(x = x, y = y)
		names(outlist)[count] <- paste(parnames[j], p, sep = "_")
		fitted@restrictions@beta@x[j,p] <- val
		if(plot.it) {
			plot(x, y, type = "l", ylab = "Value of Discrepancy Function",
				xlab = paste("Factor", p, "for", parnames[j]), ...)
			abline(v = val, col = "gray", lty = "dotted")
		}
		count <- count + 1
	}
	return(outlist)
}

## establish Domains for a cormat
FAiR_clean_cormat <-
function(cormat, level = 1) {
	if(!is(cormat, "parameter.cormat")) {
		stop("'cormat' must be of class 'parameter.cormat'")
	}
	if(!length(cormat@Domains)) {
		upper <- lower <- cormat(cormat)
		lower[cormat@free] <- -1 + .Machine$double.eps
		upper[cormat@free] <-  1 - .Machine$double.eps
		Domains <- array(cbind(lower, upper), c(dim(lower), 2))
		cormat@Domains <- Domains
	}
	if(is.null(rownames(cormat@x))) {
		letter <- ifelse(level == 1, "F", "G")
		r <- ncol(cormat@x)
		rownames(cormat@x) <- colnames(cormat@x) <- paste(letter, 1:r, sep = "")
	}
	rownames(cormat@Domains) <- rownames(cormat@free) <- rownames(cormat@x)
	colnames(cormat@Domains) <- rownames(cormat@free) <- colnames(cormat@x)
	return(cormat)
}

## establish Domains for a coef
FAiR_clean_coef <-
function(coef, manifest) {
	if(!is(coef, "parameter.coef")) {
		stop("'coef' must be of class 'parameter.coef'")
	}
	if(!length(coef@Domains)) {
		upper <- lower <- coef(coef)
		factors <- ncol(lower)
		lower[coef@free] <- if(factors[1] == 1) -1.0 else -1.5
		upper[coef@free] <- if(factors[1] == 1)  1.0 else  1.5
		coef@Domains <- array(cbind(lower, upper), c(dim(lower), 2))
	}

	if(is.null(manifest) && is.null(rownames(coef@x))) {
		rownames(coef@x) <- paste("F", 1:nrow(coef@x), sep = "")
	}
	else if(is.null(rownames(coef@x))) {
		rownames(coef@x) <- rownames(cormat(manifest))
	}
	rownames(coef@Domains) <- rownames(coef@free) <- rownames(coef@x)
	colnames(coef@Domains) <- colnames(coef@free) <- colnames(coef@x)
	return(coef)
}

## establish Domains for a scale
FAiR_clean_Omega <-
function(Omega, beta, sds) {
	if(!is(Omega, "parameter.scale")) {
		stop("'Omega' must be of class parameter.scale")
	}
	else if(!is(beta, "parameter.coef")) {
		stop("'beta' must inherit from class parameter.coef")
	}
	else if(length(Omega@Domains)) return(Omega)

	beta <- coef(beta)
	foo <- function(x) any(x == 1) && sum(x == 1) == 1 && sum(x == 0) == length(x) - 1
	scale_fixed <- apply(beta, 1, foo)
	scale_fixed[is.na(scale_fixed)] <- FALSE
	if(any(scale_fixed)) {
		warning("It appears as if you are trying to execute the ",
			"'LISREL trick' of\nfixing one coefficient to 1.0 and",
			" the rest of the coefficients in that row to 0.0.\n",
			" This trick is allowed but very discouraged in FAiR.\n",
			"Instead of fixing one coefficient to 1.0, it is usually",
			" better to designate\nit as the best indicator of a ",
			"latent factor that is measured with error.\n",
			"See the options for inequality restrictions.")
	}
	scale_free <- !scale_fixed
	Domains <- array(cbind(-18, log(2 * sds[scale_free])), c(nrow(beta), 1, 2))
	rownames(Domains) <- paste(names(sds)[scale_free], "logsd", sep = "_")
	Omega@Domains <- Domains
	names(Omega@x) <- names(Omega@free) <- rownames(Domains)
	return(Omega)
}

FAiR_default_discrepancy <-
function(discrepancy, manifest) {
	if(discrepancy[1] == "default") {
		if(is(manifest, "manifest.data")) {
			if(is(manifest@acov, "diagonalMatrix")) {
				discrepancy <- "DWLS"
			}
			else    discrepancy <- "ADF"
		}
		else if(is(manifest, "manifest.userW")) {
			if(is(manifest@acov, "diagonalMatrix")) {
				discrepancy <- "DWLS"
			}
			else    discrepancy <- "ADF"
		}
		else discrepancy <- "MLE"
	}
	else discrepancy <- match.arg(toupper(discrepancy), c("MLE", "IRLS",
					"ADF", "ELLIPTICAL", "HK", "SHK","YWLS"))
	return(discrepancy)
}
