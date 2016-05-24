#     This file is part of FAiR, a program to conduct Factor Analysis in R
#     Copyright 2008 Benjamin King Goodrich
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

## This file contains (most of) the functions that are used as lexical criteria

## NOTE: This file is meant to be read with 90 columns and 8 space tabs

## discrepancy functions used by Factanal()
FAiR_discrepancy <-
function(fn = NULL) {

	## Placeholders to shut up R CMD check
	S <- NULL    # placeholder for sample correlation matrix
	C <- NULL    # placeholder for model-reproduced correlation matrix
	llik <- NULL # placeholder for log-likelihood
	log_det_C <- NULL # placeholder for log(det(C))
	manifest <- NULL  # placeholder for manifest
	C_inv <- NULL     # placeholder for solve(C)
	restrictions <- NULL # placeholder for restrictions
	h2 <- NULL # placeholder for communalities
	R  <- NULL # placeholder for correlation matrix

	# log-likelihood discrepancy function
	FAiR_discrepancy_MLE <- 
	function(constant) {
		llik  <- log_det_C + as.vector(crossprod( as.vector(manifest@cov), 
							  as.vector(C_inv)))
		return(llik - constant)
	}

	# iteratively reweighted least squares discrepancy function for use with hetcor()
# 	FAiR_discrepancy_IRLS <-
# 	function(s) { ## FIXME, check transposition
# 		C_inv_chol <- chol(C_inv)
# 		W <- tcrossprod(manifest@D %*% t(kronecker(C_inv_chol, C_inv_chol)))
# 		W <- crossprod((kronecker(C_inv_chol, C_inv_chol)) %*% manifest@D)
# 		c <- C[lower.tri(C, TRUE)]
# 		x <- s - c
# 		middle <- t(chol(W))[lower.tri(W, TRUE)]
# 		misfit <- .C(FAiR_QD_sum, sum = numeric(1), x = x, middle = middle, 
# 				length = length(x), DUP = FALSE, NAOK = TRUE)$sum
# 		return(misfit)
# 	}

	# Yates' (1987) weighted least squares discrepancy 
	FAiR_discrepancy_YWLS <- # this isn't a discrepancy because w can be a zero matrix
	function() {
		diag(R) <- h2
# 		w <- 1 - sweep(sweep(C^2, 1, communalities, FUN = "/", FALSE),
# 					  2, communalities, FUN = "/", FALSE)
		w <- 1 - (R^2) / tcrossprod(h2)
		return(0.5 * sum(w * (manifest@cor - R)^2))
	}

	# Asymptotic Distribution Free discrepancy function
	FAiR_discrepancy_QD <- 
	function(middle, s, diag = TRUE) {
		c <- C[lower.tri(C, diag)]
		x <- s - c
		misfit <- .C(FAiR_QD_sum, sum = numeric(1), x = x, middle = middle, 
				length = length(x), DUP = FALSE, NAOK = TRUE)$sum
		return(misfit)
	}

	# Asymptotic Distribution Free discrepancy function with diagonal W
	FAiR_discrepancy_DWLS <-
	function(middle, s) {
		c <- C[lower.tri(C)]
		x <- s - c
		misfit <- c(crossprod(middle * x))
		return(misfit)
	}

	discrepancy <- list(MLE = FAiR_discrepancy_MLE, #IRLS = FAiR_discrepancy_IRLS,
			ELLIPTICAL = FAiR_discrepancy_QD,
			YWLS = FAiR_discrepancy_YWLS, ADF = FAiR_discrepancy_QD,
			DWLS = FAiR_discrepancy_DWLS,
			HK = FAiR_discrepancy_QD, SHK = FAiR_discrepancy_QD)

	if(is.null(fn)) return(discrepancy)
	else            return(discrepancy[[fn]])
}

## constraints used by Factanal() in two-level models
FAiR_constraints_2nd <- 
function(fn = NULL, onecol = FALSE) {

	## Placeholders to shut up R CMD check
	restrictions <- NULL # placehold for object that inherits from class restrictions

	FAiR_criterion_ranks_2nd <- # imposes order constraints on FCs
	function(user_ranks, MARGIN) {
		Delta <- loadings(restrictions, level = 2)
		if(ncol(Delta) > 1) {
			FC <- Delta * (Delta %*% cormat(restrictions, level = 2))
		}
		else    FC <- Delta^2
		ranks <- apply(-FC, MARGIN = MARGIN, FUN = rank, ties.method = "min")
		if(MARGIN == 1) ranks <- t(ranks)
		diff <- ranks - user_ranks
		return(-mean(diff <= 0))
	}

	FAiR_criterion_indicators_2nd <- # designate best 1st order factors
	function(indicators_2nd) {
		Delta <- loadings(restrictions, level = 2)
		if(ncol(Delta) > 1) {
			FC <- Delta * (Delta %*% cormat(restrictions, level = 2))
		}
		else    FC <- Delta^2
		ranks <- apply(-FC, MARGIN = 2, FUN = rank, ties.method = "min")
		out <- 0
		denom <- ncol(ranks)
		for(i in 1:ncol(ranks)) {
			mark <- indicators_2nd[i]
			if(is.na(mark)) denom <- denom - 1
			else            out <- out + ranks[mark,i]
		}
		out <- out / denom - 2
		return(out)
	}

	FAiR_criterion_dist_cols_2nd <- # require a minimum binary distance
	function(cutpoint_2nd) {
		coefs <- t(!loadings(restrictions, level = 2))
		dist_min <- min(dist(coefs, method = "binary"))
		if(dist_min >= cutpoint_2nd) return(-1.0)
		else return(-dist_min)
	}

	FAiR_criterion_evRF_2nd <- # effective variance of the 2nd level reference factors
	function() {               # must be greater than the effective variance of the
                                   # 1st level primary factors
		Xi  <- cormat(restrictions, level = 2)
		Phi <- cormat(restrictions, level = 1)
		ev_RF <- det(cov2cor(chol2inv(chol(Xi))))^(1/ncol(Xi))
		ev_Phi <- det(Phi)^(1/ncol(Phi))
		ev_diff <- ev_Phi - ev_RF
		if(ev_diff <= 0) return(-1.0)
		else return(ev_diff)
	}

	FAiR_criterion_evPF_2nd <- # effective variance of the 2nd level primary factors
	function() {               # must be greater than the effective variance of the
                                   # 1st level primary factors
		Xi  <- cormat(restrictions, level = 2)
		Phi <- cormat(restrictions, level = 1)
		ev_Xi <- det(Xi)^(1/ncol(Xi))
		ev_Phi <- det(Phi)^(1/ncol(Phi))
		ev_diff <- ev_Phi - ev_Xi
		if(ev_diff <= 0) return(-1.0)
		else return(ev_diff)
	}

	FAiR_criterion_h2_over_FC_2nd <- # communality must be greater than all factor
	function() {			 # contributions at level 2
		Delta <- loadings(restrictions, level = 2)
		FC <- Delta * (Delta %*% cormat(restrictions, level = 2))
		communalities <- rowSums(FC)
		return(-mean(FC <= communalities + .Machine$double.eps))
	}

	FAiR_criterion_distinguishability_2nd <- # prohibits suppressor variables among
	function() {				 # the purest 2nd-order factors
		Delta <- loadings(restrictions, level = 2)
		FC <- Delta * (Delta %*% cormat(restrictions, level = 2))
		maxs <- apply(FC, 2, which.max)
		dups <- duplicated(maxs)
		if(any(dups)) return(mean(dups))
		return(-mean(FC[maxs,] + .Machine$double.eps >= 0))
	}

	FAiR_criterion_mve_FC_2nd <- # forces FC to be inside the minimum volume ellipsoid
	function() {                 # when stacked on top of an identity matrix
		Delta <- loadings(restrictions, level = 2)
		FC <- Delta * (Delta %*% cormat(restrictions, level = 2))
		return(FAiR_mve(FC))
	} # disabled

	FAiR_criterion_no_neg_suppressors_2nd <- # no negative suppressors at level 2
	function(FC_threshold_2nd = -0.01) {
		Delta <- loadings(restrictions, level = 2)
		FC <- Delta * (Delta %*% cormat(restrictions, level = 2))
		return(-mean(FC + .Machine$double.eps >= FC_threshold_2nd))
	}

	FAiR_criterion_gv_2nd <- # primary factors must have less generalized variance
	function() {		 # than do the reference factors at level 2
		Xi_chol <- chol(cormat(restrictions, level = 2))
		D <- 1/sqrt(diag(chol2inv(Xi_chol)))
		det_Xi <- prod(diag(Xi_chol))^2
		gv_diff <- det_Xi - prod(D)
		if(gv_diff <= 0) return(-1.0)
		else             return(gv_diff)
	}

	FAiR_criterion_mve_coef_2nd <- # forces Delta to be in the min. volume ellipsoid
	function() {                   # when stacked on top of an identity matrix
		return(FAiR_mve(loadings(restrictions, level = 2)))
	} # disabled

	FAiR_criterion_block_2nd <- # prohibits some coefficients from being squashed
	function(blockers_2nd) {
		coefs <- as.logical(loadings(restrictions, level = 2))
		out   <- -mean(coefs == blockers_2nd, na.rm = TRUE)
		return(out)
	}

	FAiR_criterion_cohyperplanarity_2nd <- # forces tests within hyperplanes to have
	function() {                           # more effective variance than the primary
					       # factors at level 1
		Delta <- loadings(restrictions, level = 2)
		Phi   <-   cormat(restrictions, level = 1)
		ev_Phi <- det(Phi)^(1/ncol(Phi))
		out  <- rep(NA, ncol(Delta))
		for(i in 1:ncol(Delta)) {
			mark <- Delta[,i] == 0
			out[i] <- (det(Phi[mark,mark])^(1/sum(mark)) >= ev_Phi)
		}
		return(-mean(out))
	}

	constraint <- list(ranks_rows_2nd = list(FAiR_criterion_ranks_2nd,
			paste("Factor contributions at level 2 must respect row-wise",
				"rank order constraints")),

			ranks_cols_2nd = list(FAiR_criterion_ranks_2nd,
			paste("Factor contributions at level 2 must respect column-wise",
				"rank order constraints")),

			indicators_2nd = list(FAiR_criterion_indicators_2nd,
			paste("Designate which factor at level 1 is the best indicator",
				"of each factor at level 2")),

			dist_cols_2nd  = list(FAiR_criterion_dist_cols_2nd,
			paste("Minimum binary distance among the columns of the",
				"logically negated loadings matrix at level 2")),

			evRF_2nd = list(FAiR_criterion_evRF_2nd,
			paste("Reference factors at level 2 must have more effective",
			"variance than do the primary factors at level 1") ),

			   evPF_2nd = list(FAiR_criterion_evPF_2nd,
			paste("Primary factors at level 2 must have more effective",
			"variance than do the primary factors at level 1") ),

		     h2_over_FC_2nd = list(FAiR_criterion_h2_over_FC_2nd,
			paste("Communalities at level 2 must exceed all the factor",
			"contributions that constitute them") ),

	     distinguishability_2nd = list(FAiR_criterion_distinguishability_2nd,
			paste("The best indicator for each factor at level 2 must have",
				"non-negative factor contributions") ),

# 			 mve_FC_2nd = list(FAiR_criterion_mve_FC_2nd,
# 			paste("Factor contribtuion matrix at level 2 must be inside the",
# 				"minimum volume ellipsoid")),

	     no_neg_suppressors_2nd = list(FAiR_criterion_no_neg_suppressors_2nd,
			"No negative suppressor variables at level 2"),

			     gv_2nd = list(FAiR_criterion_gv_2nd,
			paste("Reference factors at level 2 must have more generalized",
			"variance than do the primary factors at level 2") ),

			  block_2nd = list(FAiR_criterion_block_2nd,
			c("Prohibit some coefficients at level 2 from being zero") ),

		cohyperplanarity_2nd= list(FAiR_criterion_cohyperplanarity_2nd,
			paste("Outcomes in hyperplanes at level 2 must have more",
			"effective variance than do the primary factors at level 1") )

# 		       mve_coef_2nd = list(FAiR_criterion_mve_coef_2nd,
# 			paste("Primary pattern matrix at level 2 must be inside the",
# 				"minimum volume ellipsoid")),

 		)

	if(onecol) {
		constraint$ranks_rows_2nd <- constraint$evRF_2nd <-
		constraint$dist_cols_2nd  <- constraint$evPF_2nd <-
		constraint$h2_over_FC_2nd <- constraint$no_neg_suppressors_2nd <-
		constraint$gv_2nd <- constraint$distinguishability_2nd <- 
		constraint$block_2nd <- constraint$cohyperplanarity_2nd <- NULL
	}
	if(is.null(fn)) return(constraint)
	else            return(constraint[[fn]])
}

## constraints used by Factanal() at level 1 and used by Rotate()
FAiR_constraints_1st <- 
function(fn = NULL, onecol = FALSE) {

	## Placeholders to shut up R CMD check
	A <- NULL    # placeholder for preliminary factor pattern 
	rotmat_primary <- # placeholder for rotation matrix for the primary pattern
	Phi <- NULL  # placeholder for primary factor intercorrelation matrix at level 1
	Phi_inv <- NULL # likewise for its inverse
	beta <- NULL # placeholder for primary pattern matrix at level 1
	Pi <- NULL   # placeholder for primary structure matrix at level 1
	FC <- NULL   # placeholder for factor contribution matrix at level 1
	C <- NULL    # placeholder for model-reproduced covariance matrix
	R <- NULL    # placeholder for model-reproduced correlation matrix
	Upsilon <- NULL # placeholder for reference structure matrix at level 1
	MLE <- NULL  # placeholder for whether MLE
	log_det_C <- NULL # placeholder for log of determinant of C
	ev <- NULL   # placeholder for eigenstuff
	restrictions <- NULL # placeholder for restrictions

	FAiR_criterion_ranks_1st <- # imposing ordering constraints on FCs
	function(user_ranks, MARGIN, EFA = FALSE) {
		ranks <- apply(-FC, MARGIN = MARGIN, FUN = rank, ties.method = "min")
		if(MARGIN == 1) ranks <- t(ranks)
		diff <- ranks - user_ranks
		return(-mean(diff <= 0))
	}

	FAiR_criterion_indicators_1st <- # designate best manifest variables
	function(indicators, EFA = FALSE) {
		ranks <- apply(-FC, MARGIN = 2, FUN = rank, ties.method = "min")
		out <- 0
		denom <- ncol(ranks)
		for(i in 1:ncol(ranks)) {
			mark <- indicators[i]
			if(is.na(mark)) denom <- denom - 1
			else            out <- out + ranks[mark,i]
		}
		out <- out / denom - 2
		return(out)
	}

	FAiR_criterion_evRF_1st <-   # effective variance of the 1st level reference
	function(threshold, EFA = FALSE) { # factors must be greater than the
                # effective variance of the battery as a whole (in common factor space)
		if(EFA) ev_Psi <- det(cov2cor(Phi_inv))^(1/ncol(Phi_inv))
		else {
			ev_Psi <- det(cov2cor(chol2inv(chol(Phi))))^(1/ncol(Phi))
			threshold <- det(R)^(1/ncol(R))
		}
		ev_diff <- threshold - ev_Psi
		if(ev_diff <= 0) return(-1.0)
		else return(ev_diff)
	}
	
	FAiR_criterion_evPF_1st <-   # effective variance of the 1st level primary
	function(threshold, EFA = FALSE) { # factors must be greater than the
                # effective variance of the battery as a whole (in common factor space)
		if(!EFA) threshold <- det(R)^(1/ncol(R))
		ev_Phi <- det(Phi)^(1/ncol(Phi))
		ev_diff <- threshold - ev_Phi
		if(ev_diff <= 0) return(-1.0)
		else return(ev_diff)
	}
	
	FAiR_criterion_h2_over_FC_1st <- # communality must be greater than all factor
	function(EFA = FALSE) {          # contributions at level 1
		if(EFA) h2 <- rowSums(FC)
		return(-mean(FC <= h2 + .Machine$double.eps))
	}

	FAiR_criterion_mve_FC_1st <- # forces FC to be inside the minimum volume ellipsoid
	function(EFA = FALSE) {      # when stacked on top of an identity matrix
		return(FAiR_mve(FC))
	} # disabled

	FAiR_criterion_distinguishability_1st <- # prohibits suppressor variables among
	function(EFA = FALSE) {                  # the purest 1st-order factors
		maxs <- apply(FC, 2, which.max)
		dups <- duplicated(maxs)
		if(any(dups)) return(mean(dups))
		return( -mean(FC[maxs,] + .Machine$double.eps >= 0))
	}
	
	FAiR_criterion_no_neg_suppressors_1st <- # no negative suppressors at level 1
	function(FC_threshold = -0.01, EFA = FALSE) {
		return(-mean(FC + .Machine$double.eps >= FC_threshold))
	}
	
	FAiR_criterion_gv_1st <- # primary factors must have less generalized variance
	function(EFA = FALSE) {  # than do the reference factors at level 1
		Phi_chol <- chol(Phi)
		D <- 1/sqrt(diag(chol2inv(Phi_chol)))
		det_Phi <- prod(diag(Phi_chol))^2
		gv_diff <-  det_Phi - prod(D)
		if(gv_diff <= 0) return(c(det1 = -1.0))
		else             return(c(det1 = gv_diff))
	}
	
	FAiR_criterion_mve_coef_1st <- # forces beta to be in the minimum volume ellipsoid
	function(EFA = FALSE) {        # when stacked on top of an identity matrix
		if(EFA) beta <- t(t(Upsilon) / sqrt(diag(Phi_inv)))
		return(FAiR_mve(beta))
	} # disabled

	## NOT used by Rotate()
	FAiR_criterion_cohyperplanarity_1st <- # forces tests within hyperplanes to have
	function() {                           # more effective variance than the battery
					       # as a whole in common factor space
		ev_R <- det(R)^(1/ncol(R))
		out  <- rep(NA, ncol(beta))
		for(i in 1:ncol(beta)) {
			mark <- beta[,i] == 0
			out[i] <- (det(R[mark,mark])^(1/sum(mark)) >= ev_R)
		}
		return(-mean(out))
	}

	FAiR_criterion_dist_cols_1st <- # require a minimum binary distance
	function(cutpoint, EFA = FALSE) {
		if(EFA) coefs <- t(!Upsilon)
		else    coefs <- t(!beta)
		dist_min <- min(dist(coefs, method = "binary"))
		if(dist_min >= cutpoint) return(-1.0)
		else return(-dist_min)
	}

	FAiR_criterion_volume_1st <- # requires det(R) >= det(cormat(manifest))
	function(det_S) {
		diff <- det_S - det(R)
		if(diff <= 0) return(-1)
		else          return(diff)
	}

	FAiR_criterion_block_1st <- # prohibits some coefficients from being squashed
	function(blockers) {
		coefs <- as.logical(beta)
		out   <- -mean(coefs == blockers, na.rm = TRUE)
		return(out)
	}

	constraint <- list(ranks_rows_1st = list(FAiR_criterion_ranks_1st,
			paste("Factor contributions at level 1 must respect row-wise",
				"rank order constraints")),

			ranks_cols_1st = list(FAiR_criterion_ranks_1st,
			paste("Factor contributions at level 1 must respect column-wise",
				"rank order constraints")),

			indicators_1st = list(FAiR_criterion_indicators_1st,
			paste("Designate which outcome variable is the best indicator",
				"of each factor at level 1")),

			evRF_1st = list(FAiR_criterion_evRF_1st,
			paste("Reference factors at level 1 must have more effective",
			"variance than does the battery as a whole") ),

			   evPF_1st = list(FAiR_criterion_evPF_1st,
			paste("Primary factors at level 1 must have more effective",
			"variance than does the battery as a whole") ),

		     h2_over_FC_1st = list(FAiR_criterion_h2_over_FC_1st,
			paste("Communalities at level 1 must exceed all the factor",
			"contributions that constitute them") ),

	     distinguishability_1st = list(FAiR_criterion_distinguishability_1st,
			paste("The best indicator of each factor at level 1 must have",
				"non-negative factor contributions") ),

# 			 mve_FC_1st = list(FAiR_criterion_mve_FC_1st,
# 			paste("Factor contribtuion matrix at level 1 must be inside the",
# 				"minimum volume ellipsoid")),

	     no_neg_suppressors_1st = list(FAiR_criterion_no_neg_suppressors_1st,
			"No negative suppressor variables at level 1"),

			     gv_1st = list(FAiR_criterion_gv_1st,
			paste("Reference factors at level 1 must have more generalized",
			"variance than do the primary factors at level 1") ),

# 		       mve_coef_1st = list(FAiR_criterion_mve_coef_1st,
# 			paste("Primary pattern matrix at level 1 must be inside the",
# 				"minimum volume ellipsoid")),

	       cohyperplanarity_1st = list(FAiR_criterion_cohyperplanarity_1st,
			paste("Outcomes in hyperplanes at level 1 must have more",
			"effective variance than does the battery as a whole") ),

			dist_cols_1st  = list(FAiR_criterion_dist_cols_1st,
			paste("Minimum binary distance among the columns of the",
				"logically negated loadings matrix at level 1")),

			 volume_1st = list(FAiR_criterion_volume_1st,
			paste("The generalized variance of the reproduced correlation",
			"matrix must be greater than that of the sample correlation",
			"matrix") ),

			  block_1st =  list(FAiR_criterion_block_1st,
			c("Prohibit some coefficients at level 1 from being zero") ) )


	if(onecol) {
		constraint$ranks_rows_1st <- constraint$evRF_1st <- constraint$evPF_1st <-
		constraint$h2_over_FC_1st <- constraint$no_neg_suppressors_1st <-
		constraint$gv_1st <- constraint$distinguishability_1st <-
		constraint$cohyperplanarity_1st <- constraint$dist_cols_1st <- 
		constraint$block_1st <- NULL
	}

	if(is.null(fn)) return(constraint)
	else            return(constraint[[fn]])

}

## analytic transformation criteria used by Rotate()
FAiR_analytic <- 
function(fn = NULL) {

	## Placeholders to shut up R CMD check
	Upsilon <- NULL  # placeholder for reference structure matrix
	Phi <- NULL      # placeholder for correlation matrix among primary factors
	A <- NULL	 # placeholder for preliminary primary pattern matrix
	rotmat_primary <- NULL # placeholder for transformation matrix

	FAiR_criterion_phi <- # (log of) Thurstone's criterion
	function(c = 1) {
		return(c(logphi = log(sum(exp(rowSums(log(Upsilon^2)/c))))))
	}

	FAiR_criterion_varphi <- # BG's generalization of Thurstone's criterion
	function(weights = NULL) {
		# Each row of the reference structure matrix is sorted by magnitude
		sorted <- t(apply(Upsilon^2, 1, sort))
	
		# Start with a plain Thurstone's criterion
		varphi <- phi <- sum(exp(rowSums(log(sorted))))
	
		# Exclude a column, recalculate Thurstone's criterion, weight, cumulate
		if(weights[1] < 0) for(i in 1:(ncol(sorted) - 1)) {
			varphi <- varphi + sqrt(max(sorted[,i+1])) * # dynamic weights
				sum(exp(rowSums(log(sorted[, -c(1:i), drop = FALSE]))))
		}
		else for(i in 1:(ncol(sorted) - 1)) {
			varphi <- varphi + weights[i] * 
				sum(exp(rowSums(log(sorted[, -c(1:i), drop = FALSE]))))
		}
		return(c(logvarphi = log(varphi)))
	}

	FAiR_criterion_minimaximin <-   # Another way of operationalizing
	function() {			# Thurstonian simple structure
		return(c(logmmm = log(max(apply(Upsilon^2, 1, min)))))
	}

	FAiR_criterion_LS <- # Loading Simplicity Index by Lorenzo-Seva (2004) p.50-1
	function(eps = .Machine$double.eps, scale = 10,
		E = (1 / ncol(Upsilon) + eps)^(scale / ncol(Upsilon))) {
		C <- diag(colSums(Upsilon^2))
# 		H <- diag(diag(sweep(Upsilon, 2, 1/diag(C), "*", FALSE) %*% t(Upsilon)))
		t_Upsilon <- t(Upsilon)
		H <- diag(diag(t(t_Upsilon / diag(C)) %*% t_Upsilon))
# 		B <- sweep(sweep(Upsilon, 1, diag(H)^(-.5), "*", FALSE), 
# 					  2, diag(C)^(-.5), "*", FALSE)
		B <- t(t(Upsilon / sqrt(diag(H))) / sqrt(diag(C)))
		w <- mean( (B^2 + eps)^(scale * B^2) )
	
		LS <- (w - E) / (1 - E)
		if(is.nan(LS)) LS <- 1  # absolutely perfect clustering -> NaN
		return(-LS)
	}

	FAiR_criterion_SIC <- # Stochastic Information Criterion by Preacher (2006)
	function() {
		return(NULL)
	} # disabled

	FAiR_criterion_GPA <-
	function(matrix, method, methodArgs) {
		     if(matrix == "PP") L <- A %*% rotmat_primary 
		else if(matrix == "RS") L <- Upsilon
		else if(matrix == "FC") {
						beta <- A %*% rotmat_primary
						L <- beta * (beta %*% Phi)
					}


		methodArgs$L <- L
		return(do.call(paste("vgQ", method, sep = "."), args = methodArgs))
	}

	analytic <- list(phi = list(FAiR_criterion_phi,
				"phi: criterion in Thurstone (1935)"),

		      varphi = list(FAiR_criterion_varphi,
				"varphi: a generalization of phi"),

		 minimaximin = list(FAiR_criterion_minimaximin,
				"minimaximin: another simiple structure criterion"),

			  LS = list(FAiR_criterion_LS,
			"Loading Simplicity Index: criterion by Lorenzo-Seva (2003)"),

# 			 SIC = list(FAiR_criterion_SIC,
# 			"Stochastic Information Complexity: by Preacher (2006)"),

			 GPA = list(FAiR_criterion_GPA,
			"Some oblique criterion from the GPArotation package"))

	if(is.null(fn)) return(analytic)
	else            return(analytic[[fn]])
}

## some constraints used by Rotate(); there are more in constraints_1st()
FAiR_constraints <-
function(fn = NULL) {

	## Placeholders to shut up R CMD check
	Phi <- NULL  # placeholder for primary factor intercorrelation matrix at level 1
	beta <- NULL # placeholder for primary pattern matrix at level 1
	C <- NULL    # placeholder for model-reproduced correlation matrix
	Upsilon <- NULL # placeholder for reference structure matrix at level 1

	FAiR_criterion_no_factor_collapse <- # forces effective variance of the primary
	function(Phi, nfc_threshold = 0.25) {# factor intercorrelation matrix to exceed a
					     # threshold (0 -> factor collapse)
		effective_variance <- det(Phi)^(1/ncol(Phi))
		if(effective_variance >= nfc_threshold) return(-1.0)
		else return(-effective_variance)
	}

	FAiR_criterion_limit_correlations <- # curtail pairwise correlations between
	function(lower, upper, EFA = TRUE) { # primary factors
		# Can be used to make inter-factor correlations positive (i.e. lower = 0)
		# or to prevent factor fission (upper << 1). Note that no_factor_collapse
		# and possibly no_suppressors can be used to accomplish a similar goal
		cors <- Phi[lower.tri(Phi)]
		toolow  <- cors < lower
		toohigh <- cors > upper
		inadmissable <- toolow | toohigh
		if(any(inadmissable)) {
			if(any(toolow))  x <- (cors[toolow]  - lower)^2
			else x <- 0
	
			if(any(toohigh)) y <- (cors[toohigh] - upper)^2
			else y <- 0
	
			return(sqrt(max(c(x,y))))
		}
		else return(-1.0)
	}

	FAiR_criterion_positive_manifold <- # forces reference structure correlations to
	function(pm_threshold = -0.1, EFA = TRUE) { # exceed a threshold
		pm <- mean(Upsilon + .Machine$double.eps >= pm_threshold)
		if(pm < 1) return(-min(Upsilon))
		else       return(-1.0)
	}

	constraint <- list(no_factor_collapse = list(FAiR_criterion_no_factor_collapse,
						"Factor collapse must be avoided"),

			   limit_correlations = list(FAiR_criterion_limit_correlations,
					"Bounds on correlations among primary factors"),

			    positive_manifold = list(FAiR_criterion_positive_manifold,
						"Positive manifold"))

	if(is.null(fn)) return(constraint)
	else            return(constraint[[fn]])

}
