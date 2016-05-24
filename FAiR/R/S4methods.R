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

## This file defines S4 methods

## NOTE: This file is meant to be read with 90 columns with 8 space tabs

setMethod("make_manifest", signature(x = "missing", data = "missing", covmat = "hetcor"),
function(covmat, shrink = FALSE) {

	if(!is.null(covmat$W)) W <- covmat$W
	else if(!is.null(covmat$std.errors)) {
		W <- diag(covmat$std.errors[lower.tri(covmat$std.errors)]^(-2))
	}
	else    W <- NULL
	covmat <- list(cov = covmat$correlations, n.obs = covmat$n, W = W)
	return(FAiR_make_manifest_list(covmat, shrink))
})

setMethod("make_manifest", signature(x = "missing", data = "missing", covmat = "list"),
function(covmat, n.obs = NA_integer_, shrink = FALSE) {

	if(is.null(covmat$n.obs)) covmat$n.obs <- n.obs

	return(FAiR_make_manifest_list(covmat, shrink))
})

setMethod("make_manifest", signature(x = "missing", data = "missing", covmat = "matrix"),
function(covmat, n.obs = NA_integer_, shrink = FALSE, sds = NULL) {

	covmat_list <- list(cov = covmat, n.obs = n.obs, sds = sds)
	return(make_manifest(covmat = covmat_list, shrink = shrink))
})

setMethod("make_manifest", signature(x = "missing", data = "missing", covmat = "CovMcd"),
function(covmat) {
	if(is.null(covmat@X)) {
		stop("covmat must have a slot named 'X' with the data matrix in it")
	}
	else acov <- FAiR_ADF(covmat@X[as.logical(covmat@wt),])

	manifest <- new("manifest.data.mcd", cov = covmat@cov,
			cor = cov2cor(covmat@cov), 
			sds = sqrt(diag(covmat@cov)), center = covmat@center,
			wt = covmat@wt, how = "mcd",
			acov = acov, X = covmat@X, mcd = covmat,
			n.obs = as.integer(sum(covmat@wt)), diag = TRUE)
	return(manifest)
})

setMethod("make_manifest", signature("data.frame", data = "missing", covmat = "missing"),
function(x, subset, shrink = FALSE, 
	bootstrap = 0, how = "default", seed = 12345, wt = NULL, ...) {
	if(!missing(subset)) x <- x[subset,,drop = FALSE]
	if(all(sapply(x, is.numeric))) {
		x <- as.matrix(x)
		return(make_manifest(x, bootstrap = bootstrap, how = how, seed = seed, 
					wt = wt, ...))
	}
	return(FAiR_make_manifest.data_ordinal(z = x, shrink = shrink, 
		bootstrap = bootstrap, how = how, seed = seed, wt = wt, ...))

})

setMethod("make_manifest", signature("missing", data = "data.frame", covmat = "missing"),
function(data, subset, shrink = FALSE,
	bootstrap = 0, how = "default", seed = 12345, wt = NULL, ...) {
	if(!missing(subset)) data <- data[subset,,drop = FALSE]
	if(all(sapply(data, is.numeric))) {
		x <- as.matrix(data)
		return(make_manifest(x, bootstrap = bootstrap, how = how, seed = seed, 
					wt = wt, ...))
	}
	return(FAiR_make_manifest.data_ordinal(z = data, shrink = shrink,
		bootstrap = bootstrap, how = how, seed = seed, wt = wt, ...))
})

setMethod("make_manifest", signature("missing", data = "matrix", covmat = "missing"),
function(data, subset, shrink = FALSE, 
	bootstrap = 0, how = "default", seed = 12345, wt = NULL, ...) {
	if(!missing(subset)) data <- data[subset,,drop = FALSE]
	return(make_manifest(x = data, shrink = shrink, 
		bootstrap = bootstrap, how = how, seed = seed, wt = wt, ...))
})

setMethod("make_manifest", signature(x = "matrix", data = "missing", covmat = "missing"),
function(x, subset, shrink = FALSE, 
	bootstrap = 0, how = "default", seed = 12345, wt = NULL, ...) {
	if(!missing(subset)) x <- x[subset,,drop = FALSE]
	return(FAiR_make_manifest.data_numeric(z = x, shrink = shrink,
		bootstrap = bootstrap, how = how, seed = seed, wt = wt, ...))
})

setMethod("make_manifest", signature("formula", data = "data.frame", covmat = "missing"),
function(x, data, subset, shrink = FALSE, na.action = "na.pass",
	bootstrap = 0, how = "default", seed = 12345, wt = NULL, ...) {
	z <- get_all_vars(x, data)
	if(all(sapply(z, is.numeric))) {
		z <- FAiR_parse_formula(x, data, na.action)
		if(!missing(subset)) {
			z <- z[subset,,drop = FALSE]
		}
		return(make_manifest(x = z, shrink = shrink,
			bootstrap = bootstrap, how = how, seed = seed, wt = wt, ...))
	}
	else {
		if(!missing(subset)) z <- Z[subset,,drop = FALSE]
		return(FAiR_make_manifest.data_ordinal(z = z, shrink = shrink,
			bootstrap = bootstrap, how = how, seed = seed, wt = wt, ...))
	}
})

## a basic constructor for objects that inherit from class "restrictions"
setMethod("make_restrictions", signature(manifest = "manifest.basic", 
					    Omega = "missing",
					    beta  = "missing",
					      Phi = "missing",
					    Delta = "missing",
					       Xi = "missing"), def =
function(manifest, factors = NULL, model = c("SEFA", "EFA", "CFA"), 
	discrepancy = "default", nl_1 = NULL, nl_2 = NULL) {

	S <- cormat(manifest)
	n <- nrow(S)

	model <- match.arg(toupper(model), c("SEFA", "EFA", "CFA"))

	## Deal with the number of factors
	if(is.null(factors)) {
		levels <- 1
		text <- "How many factors should be extracted at level 1?"
		factors <- as.integer(FAiR_get_number(text, from = 1, to = floor(n / 2),
						value = ifelse(n >= 5, 2, 1), by = 1))
		if(factors <= 2)        factors <- c(factors, 0L)
		else if(model == "EFA") factors <- c(factors, 0L)
		else {
			text <- paste("Would you like to estimate a simultaneous",
					"second-order model?")
			second <- FAiR_yesno(text)
			levels <- levels + second

			if(levels > 1) {
				text <- "How many factors should be extracted at level 2?"
				if(  factors < 5) factors <- c(factors, 1L)
				else factors[2] <- as.integer(FAiR_get_number(text, 
							from = 1, to = floor(factors/2),
							by = 1, value = 1))
			}
			else factors <- c(factors, 0L)
		}
	}
	else if(length(factors) > 2) {
		stop("length of factors must be at most two")
	}
	else if(length(factors) == 2) {
		if(!isTRUE(all(factors == as.integer(factors)))) {
			stop("if specified, 'factors' must contain integers only")
		}
		else if(factors[1] == 0) factors <- c(0L, 0L)
		else if(factors[1] <= factors[2]) {
			stop("the number of first-order factors must (greatly) exceed",
				" the number of second-order factors")
		}
		else if(factors[1] <= 2 && factors[2] >= 1) {
			stop("there must be at least three first-order factors to",
				"estimate a two-level model")
		}
		else factors <- as.integer(factors)
		
		levels <- 1 + (factors[2] > 0)
	}
	else {
		if(!isTRUE(factors == as.integer(factors))) {
			stop("if specified, 'factors' must be an integer")
		}
		factors <- c(as.integer(factors), 0L)
		levels <- 1
	}

	if(factors[1] == 0) {
		x <- rep(NA_real_, n)
		free <- rep(TRUE, n)
		Omega <- new("parameter.scale", x = x, free = free, num_free = n)
		restrictions <- make_restrictions(manifest, Omega, discrepancy)
		return(restrictions)
	}

	if(factors[1] == 1 && model == "SEFA") {
		stop("SEFA with one factor does not make sense, please respecify")
	}


	if(model == "EFA") {
		if(!is.null(nl_1)) stop("providing 'nl_1' is incompatable with EFA")
		if(!is.null(nl_2)) stop("providing 'nl_2' is incompatable with EFA")

		discrepancy <- FAiR_default_discrepancy(discrepancy, manifest)

		if(discrepancy != "MLE") {
			x <- matrix(NA_real_, nrow = n, ncol = factors[1])
			x[upper.tri(x)] <- 0
			rownames(x) <- rownames(S)
			free <- is.na(x)
			beta <- new("parameter.coef", x = x, free = free, 
					num_free = sum(free), invalid = 0.0)
			return(make_restrictions(manifest = manifest, beta = beta,
						discrepancy = discrepancy))
		}

		dof <- as.integer(0.5 * ((n - factors[1])^2 - n - factors[1]))
		Domains <- cbind(0, rep(1, n))
		rownames(Domains) <- rownames(S)
		restrictions <- new("restrictions.factanal", factors = factors, nvars = n,
				dof = dof, Domains = Domains, free = rep(TRUE, n),
				model = "EFA", discrepancy = "MLE")
		return(restrictions)
	}

	choices <- FAiR_menunator(factors, model == "SEFA")

	## Possibly fix some cells of the primary pattern matrix (at level 1)
	x <- matrix(NA_real_, nrow = n, ncol = factors[1])
	rownames(x) <- rownames(S)
	colnames(x) <- paste("F", 1:factors[1], sep = "")
	if(model == "CFA")                x <- FAiR_peg_coefficients(x, level = 1)
	else if(is.function(nl_1))        x <- FAiR_peg_coefficients(x, level = 1)
	else if(choices[[1]]["peg_coef"]) x <- FAiR_peg_coefficients(x, level = 1)

	## Possibly impose equality restrictions (at level 1)
	if(choices[[1]]["equalities"]) equalities <- FAiR_equality_restrictions(x, 1)
	else equalities <- list()
	if(l <- length(equalities)) for(i in 1:l) x[equalities[[i]]@fixed] <- Inf

	## Possibly impose bounds on coefficients (at level 1)
	Domains <- FAiR_bounds_coef(x, 1, choices[[1]]["bounds_coef"])
	bound <- Domains[,,1] == Domains[,,2]
	x[bound] <- Domains[,,1][bound]

	free <- is.na(x)
	if(model == "CFA") {
		beta <- if(is.function(nl_1)) new("parameter.coef.nl", x = x, free = free,
			num_free = sum(free), equalities = equalities, 
			nonlinearities = nl_1, invalid = 0.0, Domains = Domains) else
			new("parameter.coef", x = x, free = free, invalid = 0.0,
			num_free = sum(free), equalities = equalities, Domains = Domains)
	}
	else { # SEFA
		beta <- if(is.function(nl_1)) new("parameter.coef.SEFA.nl", x = x, 
				free = free, invalid = 0.0, num_free = sum(free), 
				equalities = equalities, nonlinearities = nl_1, 
				mapping_rule = mapping_rule, rankcheck = "reiersol",
				Domains = Domains) else
			new("parameter.coef.SEFA", x = x, free = free, invalid = 0.0,
				num_free = sum(free), equalities = equalities,
				mapping_rule = mapping_rule, rankcheck = "reiersol",
				Domains = Domains)
		mr    <- choices[[1]]["mapping_rule"]
		zeros <- choices[[1]]["zeros"]
		formals(beta@mapping_rule) <- FAiR_get_mapping_rule_args(x, 1, mr, zeros)
		if(any(formals(beta@mapping_rule)$zeros < factors[1])) {
			beta@rankcheck <- "howe"
		}
	}

	if(levels == 1) {
		## Make correlation matrix among first-order factors
		x  <- diag(factors[1])
		Domains <- FAiR_bounds_cormat(x, 1, choices[[1]]["bounds_cormat"])
		bound <- Domains[,,1] == Domains[,,2]
		x[bound] <- Domains[,,1][bound]
		free <- !bound
		Phi  <- new("parameter.cormat", x = x, free = free, 
				num_free = sum(free), invalid = 0.0, Domains = Domains)

		## Set up inequality restrictions
		criteria <- FAiR_inequalities(  choices[[1]]["inequalities"], FALSE,
						factors, model == "SEFA" )

		return(make_restrictions(manifest = manifest, beta = beta, Phi = Phi,
					discrepancy = discrepancy, criteria = criteria))
	}

	if(factors[2] > 1) {
		## Make correlation matrix among second-order factors
		x <- diag(factors[2])
		Domains <- FAiR_bounds_cormat(x, 2, choices[[2]]["bounds_cormat"])
		bound <- Domains[,,1] == Domains[,,2]
		x[bound] <- Domains[,,1][bound]
		free <- !bound
		Xi <- new("parameter.cormat", x = x, free = free, num_free = sum(free),
				Domains = Domains, invalid = 0.0)
	}

	## Possibly fix some cells of the primary pattern matrix at level 2
	x <- matrix(NA_real_, nrow = factors[1], ncol = factors[2])
	rownames(x) <- paste("F",  1:factors[1], sep = "")
	colnames(x) <- paste("G",  1:factors[2], sep = "")

	if(model == "CFA" && factors[2] > 1 ) x <- FAiR_peg_coefficients(x, level = 2)
	else if(is.function(nl_2))            x <- FAiR_peg_coefficients(x, level = 2)
	else if(choices[[2]]["peg_coef"])     x <- FAiR_peg_coefficients(x, level = 2)

	## Possibly impose equality restrictions at level 2
	if(choices[[2]]["equalities"]) equalities <- FAiR_equality_restrictions(x, 2)
	else equalities <- list()
	if(l <- length(equalities)) for(i in 1:l) x[equalities[[i]]@fixed] <- Inf

	## Possibly impose bounds on coefficients (at level 2)
	Domains <- FAiR_bounds_coef(x, 2, choices[[2]]["bounds_coef"])
	bound <- Domains[,,1] == Domains[,,2]
	x[bound] <- Domains[,,1][bound]

	free <- is.na(x)

	## Set up inequality restrictions
	criteria <- FAiR_inequalities(  choices[[1]]["inequalities"],
					choices[[2]]["inequalities"], factors,
					model == "SEFA" )

	if(factors[2] == 1) { # general second-order factor
		Delta <- if(is.function(nl_2)) new("parameter.coef.nl", x = x, free=free,
				num_free = sum(free), invalid = 0.0, Domains = Domains,
				equalities = equalities, nonlinearities = nl_2) else
			new("parameter.coef", x = x, free = free, num_free = sum(free),
				equalities = equalities, invalid = 0.0, Domains = Domains)

		return(make_restrictions(manifest = manifest, beta = beta, Delta = Delta,
					discrepancy = discrepancy, criteria = criteria))
	}
	else if(model == "CFA") {
		Delta <- if(is.function(nl_2)) new("parameter.coef.nl", x = x, free=free,
				num_free = sum(free), invalid = 0.0, Domains = Domains,
				equalities = equalities, nonlinearities = nl_2) else
			new("parameter.coef", x = x, free = free, num_free = sum(free),
				equalities = equalities, invalid = 0.0, Domains = Domains)
	}
	else { # SEFA
		Delta <- if(is.function(nl_2)) new("parameter.coef.SEFA.nl", x = x, 
			free = free, invalid = 0.0, num_free = sum(free),
			equalities = equalities, nonlinearities = nl_2, Domains = Domains,
			mapping_rule = mapping_rule, rankcheck = "reiersol") else
			new("parameter.coef.SEFA", x = x, free = free, invalid = 0.0,
			num_free = sum(free), equalities = equalities, Domains = Domains,
			mapping_rule = mapping_rule, rankcheck = "reiersol")

		mr    <- choices[[2]]["mapping_rule"]
		zeros <- choices[[2]]["zeros"]
		formals(Delta@mapping_rule) <- FAiR_get_mapping_rule_args(x, 2, mr, zeros)
		if(any(formals(Delta@mapping_rule)$zeros < factors[1])) {
			Delta@rankcheck <- "howe"
		}

	}

	return(make_restrictions(manifest = manifest, beta = beta, Delta = Delta, Xi = Xi,
				discrepancy = discrepancy, criteria = criteria))
})

## make restrictions.independent
setMethod("make_restrictions", signature(manifest = "manifest.basic", 
					    Omega = "parameter.scale",
					    beta  = "missing",
					      Phi = "missing",
					    Delta = "missing",
					       Xi = "missing"), def =
function(manifest, Omega, discrepancy = "default") {

	factors <- c(0L, 0L)
	if(!length(Omega@Domains)) {
		Domains <- cbind(-18, log(2 * manfiest@sds))
		rownames(Domains) <- paste(rownames(cormat(manifest))[Omega@free],
						"logsd", sep = "_")
		colnames(Domains) <- c("lower", "upper")
		Omega@Domains <- Domains
	}
	else    Domains <- Omega@Domains

	discrepancy <- FAiR_default_discrepancy(discrepancy, manifest)

	n <- nrow(cormat(manifest))
	dof <- as.integer(0.5 * n * (n + 1) - nrow(Domains))
	restrictions <- new("restrictions.independent", factors = factors, dof = dof,
				Domains = Domains, nvars = nrow(Domains),
				Omega = Omega, criteria = criteria, model = "CFA",
				discrepancy = discrepancy, free = Omega@free)
	return(restrictions)
})

## make restrictions.orthonormal
setMethod("make_restrictions", signature(manifest = "manifest.basic", 
					    Omega = "missing",
					    beta  = "parameter.coef",
					      Phi = "missing",
					    Delta = "missing",
					       Xi = "missing"), def =
function(manifest, beta, discrepancy = "default") {

	S <- cormat(manifest)
	n <- nrow(S)
	factors <- ncol(coef(beta))
	dof <- as.integer(0.5 * ((n - factors)^2 - n - factors))

	beta <- FAiR_clean_coef(beta, manifest)

	if(!length(beta@Domains)) {
		dims <- dim(coef(beta))
		beta@Domains <- array(cbind(array(-1, dims), array(1, dims)), c(dims, 2))
	}
	Domains <- apply(beta@Domains, 3, FUN = function(x) x[beta@free])

	rn <- as.character(NULL)
	for(i in 1:factors) rn <- c(rn, paste(rownames(S)[i:n], i, sep = "_"))
	rn <- rn[beta@free]
	rownames(Domains) <- rn

	select <- c(rep(TRUE, nrow(Domains)), rep(FALSE, n))
	beta@select <- select

	Omega <- new("parameter.scale", x = manifest@sds, free = rep(TRUE, n), 
			num_free = n, select = !select, invalid = 0.0)
	Omega <- FAiR_clean_Omega(Omega, beta, manifest@sds)

	temp_Domains <- apply(Omega@Domains, 3, FUN = function(x) x[Omega@free])
	rownames(temp_Domains) <- names(Omega@x)[Omega@free]

	Domains <- rbind(Domains, temp_Domains)
	colnames(Domains) <- c("lower", "upper")

	discrepancy <- FAiR_default_discrepancy(discrepancy, manifest)
	factors <- c(factors, 0L)
	criteria <- FAiR_criterionator_extraction(list(), list(), discrepancy, 
							factors, manifest)

	restrictions <- new("restrictions.orthonormal",
				factors = factors, nvars = nrow(Domains), 
				dof = dof, model = "EFA", Domains = Domains, 
				discrepancy = discrepancy, beta = beta, Omega = Omega, 
				criteria = criteria, free = c(beta@free, Omega@free))
	return(restrictions)
})

## make restrictions.1storder with missing Omega
setMethod("make_restrictions", signature(manifest = "manifest.basic", 
					    Omega = "missing",
					    beta  = "parameter.coef",
					      Phi = "parameter.cormat",
					    Delta = "missing",
					       Xi = "missing"), def =
function(manifest, beta, Phi, discrepancy = "default", 
	criteria = list(), methodArgs = list()) {

	S <- cormat(manifest)
	n <- nrow(S)
	factors <- ncol(coef(beta))

	Omega <- new("parameter.scale", x = rep(NA_real_, n), free = rep(TRUE, n),
			num_free = n, invalid = 0.0)
	Omega <- FAiR_clean_Omega(Omega, beta, manifest@sds)

	return(make_restrictions(manifest = manifest, Omega = Omega, beta = beta,
				Phi = Phi, discrepancy = discrepancy,
				criteria = criteria, methodArgs = methodArgs))
})

## make restrictions.1storder with present Omega
setMethod("make_restrictions", signature(manifest = "manifest.basic", 
					    Omega = "parameter.scale",
					    beta  = "parameter.coef",
					      Phi = "parameter.cormat",
					    Delta = "missing",
					       Xi = "missing"), def =
function(manifest, Omega, beta, Phi, discrepancy = "default", 
	criteria = list(), methodArgs = list()) {

	coef <- coef(beta)
	factors <- c(ncol(coef), 0L)
	SEFA <- is(beta, "parameter.coef.SEFA")
	S <- cormat(manifest)

	Phi <- FAiR_clean_cormat(Phi)
	Domains <- apply(Phi@Domains, 3, FUN = function(x) x[Phi@free])
	if(length(Domains) == 0) Domains <- as.numeric(NULL)
	else {
		if(!is.matrix(Domains)) Domains <- matrix(Domains, nrow = 1)
		rn <- as.character(NULL)
		for(i in 1:(factors[1]-1)) for(j in (i+1):factors[1]) {
			rn <- c(rn, paste("F", j, " <-> ", "F", i, sep = ""))
		}
		rn <- rn[Phi@free[lower.tri(Phi@free)]]
		rownames(Domains) <- rn
	}

	beta <- FAiR_clean_coef(beta, manifest)
	temp_Domains <- apply(beta@Domains, 3, FUN = function(x) x[beta@free])

	rn <- as.character(NULL)
	for(i in 1:factors[1]) rn <- c(rn, paste(rownames(S), "_{", i, "}", sep = ""))
	rn <- rn[beta@free]
	rownames(temp_Domains) <- rn
	Domains <- rbind(Domains, temp_Domains)

	Omega <- FAiR_clean_Omega(Omega, beta, manifest@sds)
	temp_Domains <- apply(Omega@Domains, 3, FUN = function(x) x[Omega@free])
	rownames(temp_Domains) <- names(Omega@x)[Omega@free]

	Domains <- rbind(Domains, temp_Domains)
	colnames(Domains) <- c("lower", "upper")

	discrepancy <- FAiR_default_discrepancy(discrepancy, manifest)
	criteria <- FAiR_criterionator_extraction(criteria, methodArgs, 
							discrepancy, factors, manifest)

	n <- nrow(S)
	dof <- 0.5 * n * (n + 1) - nrow(Domains) + sum(Domains[,1] == Domains[,2])
	if(SEFA) dof <- dof + factors[1]^2
	dof <- as.integer(dof)

	free <- c(Phi@free[lower.tri(Phi@x)], beta@free, Omega@free)
	  Phi@select <- c(rep(TRUE,    Phi@num_free), 
			  rep(FALSE,  beta@num_free + Omega@num_free))
	 beta@select <- c(rep(FALSE,   Phi@num_free), rep(TRUE, beta@num_free),
			  rep(FALSE, Omega@num_free))
	Omega@select <- c(rep(FALSE,   Phi@num_free + beta@num_free),
			  rep(TRUE,  Omega@num_free))

	restrictions <- new("restrictions.1storder", factors = factors, 
				nvars = nrow(Domains), dof = dof, Domains = Domains,
				model = if(SEFA) "SEFA" else "CFA", free = free,
				criteria = criteria, discrepancy = discrepancy, Phi = Phi,
				beta = beta, Omega = Omega)
	return(restrictions)
})

## make restrictions.general with missing Omega
setMethod("make_restrictions", signature(manifest = "manifest.basic", 
					    Omega = "missing",
					    beta  = "parameter.coef",
					      Phi = "missing",
					    Delta = "parameter.coef",
					       Xi = "missing"), def =
function(manifest, beta, Delta, discrepancy = "default",
	criteria = list(), methodArgs = list()) {

	S <- cormat(manifest)
	n <- nrow(S)
	factors <- ncol(coef(beta))

	Omega <- new("parameter.scale", x = rep(NA_real_, n), free = rep(TRUE, n),
			num_free = n, invalid = 0.0)
	Omega <- FAiR_clean_Omega(Omega, beta, manifest@sds)
	return(make_restrictions(manifest = manifest, Omega = Omega, beta = beta,
				Delta = Delta, discrepancy = discrepancy,
				criteria = criteria, methodArgs = methodArgs))
})

## make restrictions.general with present Omega
setMethod("make_restrictions", signature(manifest = "manifest.basic", 
					    Omega = "parameter.scale",
					    beta  = "parameter.coef",
					      Phi = "missing",
					    Delta = "parameter.coef",
					       Xi = "missing"), def =
function(manifest, Omega, beta, Delta, discrepancy = "default",
	criteria = list(), methodArgs = list()) {

	coef <- coef(beta)
	factors <- c(ncol(coef), 1L)
	SEFA <- is(beta, "parameter.coef.SEFA")
	S <- cormat(manifest)

	Delta <- FAiR_clean_coef(Delta, NULL)
	Domains <- apply(Delta@Domains, 3, FUN = function(x) x[Delta@free])
	if(length(Domains) == 0) Domains <- as.numeric(NULL)
	else {
		if(!is.matrix(Domains)) Domains <- matrix(Domains, nrow = 1)
		rownames(Domains) <- paste("G -> F", 1:factors[1], sep = "")[Delta@free]
	}

	beta <- FAiR_clean_coef(beta, manifest)
	temp_Domains <- apply(beta@Domains, 3, FUN = function(x) x[beta@free])

	rn <- as.character(NULL)
	for(i in 1:factors[1]) rn <- c(rn, paste(rownames(S), "_{", i, "}", sep = ""))
	rn <- rn[beta@free]
	rownames(temp_Domains) <- rn
	Domains <- rbind(Domains, temp_Domains)

	Omega <- FAiR_clean_Omega(Omega, beta, manifest@sds)
	temp_Domains <- apply(Omega@Domains, 3, FUN = function(x) x[Omega@free])
	rownames(temp_Domains) <- names(Omega@x)[Omega@free]

	Domains <- rbind(Domains, temp_Domains)
	colnames(Domains) <- c("lower", "upper")

	discrepancy <- FAiR_default_discrepancy(discrepancy, manifest)
	criteria <- FAiR_criterionator_extraction(criteria, methodArgs, 
							discrepancy, factors, manifest)

	n <- nrow(coef)
	dof <- 0.5 * n * (n + 1) - nrow(Domains) + sum(Domains[,1] == Domains[,2])
	if(SEFA) dof <- dof + factors[1]^2
	dof <- as.integer(dof)

	free <- c(Delta@free, beta@free, Omega@free)
	Delta@select <- c(rep(TRUE,  Delta@num_free), 
			  rep(FALSE,  beta@num_free + Omega@num_free))
	 beta@select <- c(rep(FALSE, Delta@num_free), rep(TRUE, beta@num_free),
			  rep(FALSE, Omega@num_free))
	Omega@select <- c(rep(FALSE, Delta@num_free + beta@num_free),
			  rep(TRUE,  Omega@num_free))

	Phi <- diag(factors[1])
	Phi <- new("parameter.cormat", x = Phi, num_free = 0L, invalid = 0.0,
			free = matrix(FALSE, factors[1], factors[1]))

	restrictions <- new("restrictions.general", factors = factors, 
				nvars = nrow(Domains), dof = dof, Domains = Domains,
				model = if(SEFA) "SEFA" else "CFA",
				discrepancy = discrepancy, Delta = Delta,
				Phi = Phi, beta = beta, Omega = Omega, 
				criteria = criteria, free = free)
	return(restrictions)
})

## make restrictions.2ndorder with missing Omega
setMethod("make_restrictions", signature(manifest = "manifest.basic", 
					    Omega = "missing",
					    beta  = "parameter.coef",
					      Phi = "missing",
					    Delta = "parameter.coef",
					       Xi = "parameter.cormat"), def =
function(manifest, beta, Delta, Xi, discrepancy = "default", 
	criteria = list(), methodArgs = list()) {

	S <- cormat(manifest)
	n <- nrow(S)
	factors <- ncol(coef(beta))

	Omega <- new("parameter.scale", x = rep(NA_real_, n), free = rep(TRUE, n),
			num_free = n, invalid = 0.0)
	Omega <- FAiR_clean_Omega(Omega, beta, manifest@sds)
	return(make_restrictions(manifest = manifest, Omega = Omega, beta = beta,
				Delta = Delta, Xi = Xi, discrepancy = discrepancy,
				criteria = criteria, methodArgs = methodArgs))
})

## make restrictions.2ndorder with present Omega
setMethod("make_restrictions", signature(manifest = "manifest.basic", 
					    Omega = "parameter.scale",
					    beta  = "parameter.coef",
					      Phi = "missing",
					    Delta = "parameter.coef",
					       Xi = "parameter.cormat"), def =
function(manifest, Omega, beta, Delta, Xi, discrepancy = "default",
	criteria = list(), methodArgs = list()) {

	coef.1  <- coef(beta)
	coef.2  <- coef(Delta)
	factors <- c(ncol(coef.1), ncol(coef.2))
	SEFA <- is(beta,  "parameter.coef.SEFA") |
		is(Delta, "parameter.coef.SEFA")
	S <- cormat(manifest)

	Xi <- FAiR_clean_cormat(Xi, level = 2)
	Domains <- apply(Xi@Domains, 3, FUN = function(x) x[Xi@free])
	if(length(Domains) == 0) Domains <- as.numeric(NULL)
	else {
		if(!is.matrix(Domains)) Domains <- matrix(Domains, nrow = 1)
	
		rn <- as.character(NULL)
		for(i in 1:(factors[2]-1)) for(j in (i+1):factors[2]) {
			rn <- c(rn, paste("G", j, " <-> ", "G", i, sep = ""))
		}
		rn <- rn[Xi@free[lower.tri(Xi@free)]]
		rownames(Domains) <- rn
	}

	Delta <- FAiR_clean_coef(Delta, NULL)
	temp_Domains <- apply(Delta@Domains, 3, FUN = function(x) x[Delta@free])
	rn <- as.character(NULL)
	for(i in 1:factors[2]) rn <- c(rn, paste(paste("F", 1:factors[1], sep = ""),
						"_{", i, "}", sep = ""))
	rn <- rn[Delta@free]
	rownames(temp_Domains) <- rn
	Domains <- rbind(Domains, temp_Domains)

	beta <- FAiR_clean_coef(beta, manifest)
	temp_Domains <- apply(beta@Domains, 3, FUN = function(x) x[beta@free])

	rn <- as.character(NULL)
	for(i in 1:factors[1]) rn <- c(rn, paste(rownames(S), "_{", i, "}", sep = ""))
	rn <- rn[beta@free]
	rownames(temp_Domains) <- rn
	Domains <- rbind(Domains, temp_Domains)

	Omega <- FAiR_clean_Omega(Omega, beta, manifest@sds)
	temp_Domains <- apply(Omega@Domains, 3, FUN = function(x) x[Omega@free])
	rownames(temp_Domains) <- names(Omega@x)[Omega@free]

	Domains <- rbind(Domains, temp_Domains)
	colnames(Domains) <- c("lower", "upper")

	discrepancy <- FAiR_default_discrepancy(discrepancy, manifest)
	criteria <- FAiR_criterionator_extraction(criteria, methodArgs, 
							discrepancy, factors, manifest)

	n <- nrow(coef.1)
	dof <- 0.5 * n * (n + 1) - nrow(Domains) + sum(Domains[,1] == Domains[,2])
	if(SEFA) dof <- dof + factors[1]^2 + factors[2]^2
	dof <- as.integer(dof)

	free <- c(Xi@free[lower.tri(Xi@x)], Delta@free, beta@free, Omega@free)
	   Xi@select <- c(rep(TRUE,     Xi@num_free), rep(FALSE, Delta@num_free + 
				      beta@num_free + Omega@num_free))
	Delta@select <- c(rep(FALSE,    Xi@num_free), rep(TRUE,  Delta@num_free),
			  rep(FALSE,  beta@num_free + Omega@num_free))
	 beta@select <- c(rep(FALSE,    Xi@num_free + Delta@num_free), 
			  rep(TRUE,   beta@num_free), rep(FALSE, Omega@num_free))
	Omega@select <- c(rep(FALSE,    Xi@num_free + Delta@num_free + beta@num_free),
			  rep(TRUE,  Omega@num_free))

	Phi <- diag(factors[1])
	Phi <- new("parameter.cormat", x = Phi, num_free = 0L, invalid = 0.0,
			free = matrix(FALSE, factors[1], factors[1]))

	restrictions <- new("restrictions.2ndorder", factors = factors, 
				nvars = nrow(Domains), Domains = Domains, dof = dof,
				model = if(SEFA) "SEFA" else "CFA", free = free,
				discrepancy = discrepancy, Xi = Xi, Delta = Delta,
				Phi = Phi, beta = beta, Omega = Omega,
				criteria = criteria)
	return(restrictions)
})

## make_parameter methods
setMethod("make_parameter", signature(object = "parameter.scale"),  definition = 
function(par, object) {
	sds <- object@x
	sds[object@free] <- exp(par[object@select])
	FAiR_set_slot(object, "x", value = sds)
	return(object)
})

setMethod("make_parameter", signature(object = "parameter.cormat"), definition = 
function(par, object, lower) {
	cormat <- diag(rep(0.5, ncol(object@x)))
	cormat[object@free] <- par[object@select]
	fixed <- lower.tri(cormat) & !object@free
	cormat[fixed] <- object@x[fixed]
	cormat <- cormat + t(cormat)
	cormat <- FAiR_nearPD(cormat, posd.tol = lower) # fudge if not pd
	check <- attributes(cormat)$ev
	if(check < lower) { # not sufficiently positive definite -> bailout from env
		FAiR_set_slot(object, "invalid", value = -check)
	}
	else    FAiR_set_slot(object, "invalid", value = 0.0)
	FAiR_set_slot(object, "x", value = cormat)
	return(object)
})

setMethod("make_parameter", signature(object = "parameter.coef"), definition = 
function(par, object, ...) {
	coefs <- object@x
	coefs[object@free] <- par[object@select]               # fill free cells

	if(l <- length(object@equalities)) for(i in 1:l) {     # deal with equalities
		 coefs[object@equalities[[i]]@fixed] <- coefs[object@equalities[[i]]@free]
	}

	signs <- ifelse(colSums(coefs) >= 0, 1, -1)
	if( (check <- -mean(signs)) != -1 ) {                  # check positiveness
		check <- FAiR_uniquify(check)
		FAiR_set_slot(object, "invalid", value = check)
	}
	else    FAiR_set_slot(object, "invalid", value = 0.0)
	FAiR_set_slot(object, "x", value = coefs)
	return(object)
})

setMethod("make_parameter", signature(object = "parameter.coef.nl"), definition = 
function(par, object, ...) {
	coefs <- object@x
	coefs[object@free] <- par[object@select]               # fill free cells

	if(l <- length(object@equalities)) for(i in 1:l) {  # deal with equalities
		 coefs[object@equalities[[i]]@fixed] <- coefs[object@equalities[[i]]@free]
	}

	coefs <- object@nonlinearities(x = coefs)              # nonlinear restrictions

	signs <- ifelse(colSums(coefs) >= 0, 1, -1)
	if( (check <- -mean(signs)) != -1 ) {                  # check positiveness
		check <- FAiR_uniquify(check)
		FAiR_set_slot(object, "invalid", value = check)
	}
	else    FAiR_set_slot(object, "invalid", value = 0.0)

	FAiR_set_slot(object, "x", value = coefs)
	return(object)
})

setMethod("make_parameter", signature(object = "parameter.coef.SEFA"), definition = 
function(par, object, cormat, mapping_rule) {
	coefs <- object@x
	coefs[object@free] <- par[object@select]               # fill free cells

	if(l <- length(object@equalities)) for(i in 1:l) {     # deal with equalities
		 coefs[object@equalities[[i]]@fixed] <- coefs[object@equalities[[i]]@free]
	}

	if(mapping_rule) coefs <- object@mapping_rule(coefs, cormat)
	FAiR_set_slot(object, "x", value = coefs)

	if(object@rankcheck == "reiersol") check <- FAiR_check_Reiersol(coefs)
	else                               check <- FAiR_check_Howe(coefs)
	if(check != -1) { # fail submatrix rank check
		FAiR_set_slot(object, "invalid", value = check)
		return(object)
	}

	signs <- ifelse(colSums(coefs) >= 0, 1, -1)
	if( (check <- -mean(signs)) != -1 ) {                  # check positiveness
		check <- FAiR_uniquify(check)
		FAiR_set_slot(object, "invalid", value = check)
	}
	else    FAiR_set_slot(object, "invalid", value = 0.0)
	return(object)
})

setMethod("make_parameter", signature(object = "parameter.coef.SEFA.nl"), definition = 
function(par, object, cormat, mapping_rule) {
	coefs <- object@x
	coefs[object@free] <- par[object@select]               # fill free cells

	if(l <- length(object@equalities)) for(i in 1:l) {  # deal with equalities
		 coefs[object@equalities[[i]]@fixed] <- coefs[object@equalities[[i]]@free]
	}

	coefs <- object@nonlinearities(x = coefs)              # nonlinear restrictions
	if(mapping_rule) coefs <- object@mapping_rule(coefs, cormat)
	FAiR_set_slot(object, "x", value = coefs)

	if(object@rankcheck == "reiersol") check <- FAiR_check_Reiersol(coefs)
	else                               check <- FAiR_check_Howe(coefs)
	if(check != -1) {                                      # fail submatrix rank check
		FAiR_set_slot(object, "invalid", value = check)
		return(object)
	}

	signs <- ifelse(colSums(coefs) >= 0, 1, -1)
	if( (check <- -mean(signs)) != -1 ) {                  # check positiveness
		check <- FAiR_uniquify(check)
		FAiR_set_slot(object, "invalid", value = check)
	}
	else    FAiR_set_slot(object, "invalid", value = 0.0)
	return(object)
})

## Default methods that will be called if there are no specific methods available.

# lexical fits
setMethod("fitS4", signature(restrictions = "restrictions", 
				 manifest = "manifest.basic"), definition = 
function(par, restrictions, manifest, lower, mapping_rule) {
	model <- restrictions2model(par, restrictions, manifest, lower, 
					mapping_rule = mapping_rule)
	criteria <- model$criteria
	if(any(criteria != -1)) { # fails a constraint
		return(c(criteria, rep(0, 2 + length(restrictions@criteria))))
	}
	restrictions <- model$restrictions
	return(c(criteria, FAiR_lexical_driver(restrictions, manifest, lower)))
})

# scalar fit (discrepancy)
setMethod("bfgs_fitS4", signature(  restrictions = "restrictions",
					manifest = "manifest.basic"), 
definition = function(par, restrictions, manifest, helper, lower) {

	# Check bounds
	if(any(par < restrictions@Domains[,1]) | any(par > restrictions@Domains[,2])) {
		return(.Machine$double.xmax)
	}

	SEFA <- restrictions@model == "SEFA"
	if(helper$marker == length(helper$fits)) { # candidate was admissable at outset
		if(SEFA) { # squash some coefficiets and pretend it's a CFA
			small <- par[helper$squashed]
			par[helper$squashed] <- 0
			out <- crossprod(small)
		}
		else out <- 0

		# get *all* fits
		fits <- fitS4(par, restrictions, manifest, lower,
				mapping_rule = helper$done)
		if(helper$done) return(fits[length(fits)]) # ignores some constraints

		marker <- which(fits != -1)[1]
		if(marker < length(helper$fits)) {         # fails a constraint
			return(.Machine$double.xmax)
		}
		else    return(out + fits[marker])
	}
	else {  # optimize with respect to the marginal criterion
		fits <- fitS4(par, restrictions, manifest, lower, mapping_rule = TRUE)
		marker <- which(fits != -1)[1]
		if(marker > helper$marker) return(-.Machine$double.xmax)
		return(fits[marker])
	}
})

# numeric gradient at scalar fit (discrepancy)
setMethod("gr_fitS4", signature(    restrictions = "restrictions", 
					manifest = "manifest.basic"), definition = 
function(par, restrictions, manifest, helper, lower) {
	gradient <- FAiR_numeric_gradient(par, restrictions, manifest, helper, lower)
	return(gradient)
})

# helper
setMethod("bfgs_helpS4", signature( restrictions = "restrictions", 
					manifest = "manifest.basic"), definition = 
function(initial, restrictions, manifest, done, lower) {
	SEFA <- restrictions@model == "SEFA"
	model <- restrictions2model(initial, restrictions, manifest, lower, 
					mapping_rule = TRUE)
	criteria <- model$criteria
	if(any(criteria != -1)) { # fails some constraint
		marker <- which(criteria != -1)[1]
		fits <- c(criteria, rep(0, length(restrictions@criteria)))
		return(list(fits = fits, marker = marker, done = done))
	}

	# calculate initial state of genetic individual
	restrictions <- model$restrictions
	fits <- c(criteria, FAiR_lexical_driver(restrictions, manifest, lower = lower))
	marker <- which(fits != -1)[1]
	if( (marker < length(fits)) | !SEFA ) {
		return(list(fits = fits, marker = marker, done = done))
	}

	# in case of SEFA, convert to CFA
	if(is(restrictions, "restrictions.2ndorder")) {
		# Mark which elements of Delta and beta got squashed to zero
		squashed_Delta <- restrictions@Delta@select
		squashed_Delta[squashed_Delta] <- c(!restrictions@Delta@x[
						     restrictions@Delta@free])
	
		squashed_beta <- restrictions@beta@select
		squashed_beta[squashed_beta] <- c(!restrictions@beta@x[
						   restrictions@beta@free])
	
		squashed <- squashed_Delta | squashed_beta
	}
	else {
		# Mark which elements of beta got squashed to zero
		squashed <- restrictions@beta@select
		squashed[squashed] <- c(!restrictions@beta@x[restrictions@beta@free])
	}
	return(list(fits = fits, marker = marker, squashed = squashed, done = done))
})

# create random starting values
setMethod("create_start", signature(restrictions = "restrictions",
					manifest = "manifest.basic"), definition =
function(number, start, restrictions, manifest, reject = FALSE) {
	out <- apply(restrictions@Domains, 1, FUN = function(x) {
			return(runif(number, min = x[1], max = x[2]))
		})
	return(out)
})

# create a half-decent FA object 
# this method *really* should be (and is) superceded by a specific one
setMethod("create_FAobject", signature(restrictions = "restrictions",
					   manifest = "manifest.basic"), definition = 
function(restrictions, manifest, opt, call, scores, lower, analytic) {
	fits <- opt$value
	marker <- which(fits != -1)[1]
	if(marker < length(fits)) warning("constraints did not bind")

	par <- opt$par
	model <- restrictions2model(par, restrictions, manifest, lower, 
					mapping_rule = TRUE)
	restrictions <- model$restrictions # has all the parameter estimates

	# create a bunch of skeletons
	S <- cormat(manifest)
	p <- nrow(S)
	factors <- restrictions@factors
	factornames <- paste("F", 1:factors[1], sep = "")

	loadings <- array(NA_real_, dim = c(p, factors[1], 5), dimnames = c(
				rownames(S), factornames,c("PP", "PS", "RP", "RS", "FC")))

	correlations <- array(NA_real_, dim = c(factors[1], factors[1], 3), dimnames = c(
				factornames, factornames, c("PF", "RF", "PR")))

	attributes(correlations)$orthogonal <- FALSE

	n_star <- 0.5 * nrow(S) * (nrow(S) + 1)
	derivs <- matrix(NA_real_, nrow = n_star, ncol = restrictions@nvars)
	vcov   <- matrix(NA_real_, restrictions@nvars, restrictions@nvars)

	scores <- matrix(NA_real_, nrow = if(is.null(manifest@n.obs)) 0 else
				manifest@n.obs, ncol = factors[1])
	
	seeds <- matrix(NA_integer_, nrow = 1, ncol = 2)
	FAobject <- new("FA",   loadings = loadings, correlations = correlations,
			uniquenesses = uniquenesses, call = call,
			scale = restrictions@Omega@x, restrictions = restrictions, 
			Jacobian = derivs, vcov = vcov, scores = scores, 
			manifest = manifest, optimization = list(extraction = opt))
	warning("signature of restrictions object not recognized",
		" but the estimates are buried somewhere in the restrictions slot")
	return(FAobject)
})

## Begin methods for restrictions.independent
setMethod("restrictions2model", signature(  restrictions = "restrictions.independent",
						manifest = "manifest.basic"), definition =
function(par, restrictions, manifest, lower, mapping_rule = FALSE) {
	# Pull in estimated manifest standard deviations
	FAiR_set_slot(restrictions, "Omega", 1, make_parameter(par, restrictions@Omega))
	return(list(criteria = -1, restrictions = restrictions))
})

setMethod("fitS4", signature(restrictions = "restrictions.independent", 
				manifest  = "manifest.basic"), definition = 
function(par, restrictions, manifest, lower, mapping_rule = FALSE) {
	model <- restrictions2model(par, restrictions, manifest, lower)
	C <- fitted(model$restrictions, standardized = FALSE)
	environment(restrictions@criteria[[1]]) <- environment()
	return(c(model$criteria, restrictions@criteria[[1]]()))
})

setMethod("create_FAobject", signature(restrictions = "restrictions.independent",
					   manifest = "manifest.basic"), definition = 
function(restrictions, manifest, opt, call, scores, lower, analytic) {
	fits <- opt$value

	par <- opt$par
	model <- restrictions2model(par, restrictions, manifest, lower)
	restrictions <- model$restrictions

	S <- cormat(manifest)
	sorter <- NA # temporary
	uncertainty <- try(FAiR_uncertainty(restrictions, manifest, factors, 
						sorter, NULL, analytic), TRUE)
	if(!is.list(uncertainty)) {
		warning("there was a problem calculating the variance-covariance matrix",
			" of the parameters")
		cols <- restrictions@nvars
		uncertainty <- list(vcov = matrix(NA, nrow = cols, ncol = cols),
				Jacobian = matrix(NA_real_, nrow = nrow(S), ncol = cols))
	}

	seeds <- matrix(NA_integer_, nrow = 1, ncol = 2)

	loadings <- array(0, dim = c(nrow(S), 0, 5), dimnames = 
			list(rownames(S), NULL, c("PP", "RS", "PS", "RP", "FC")))
	correlations <- array(0, dim = c(0, 0, 3), 
				dimnames = list(NULL, NULL, c("PF", "RF", "PR")))


	FAobject <- new("FA", loadings = loadings, correlations = correlations,
			uniquenesses = rep(1, ncol(S)), scale = restrictions@Omega@x,
			restrictions = restrictions, Jacobian = uncertainty$Jacobian, 
			vcov = uncertainty$vcov, scores = matrix(NA_real_, 0, 0),
			call = call, manifest = manifest, 
			optimization = list(extraction = opt))
	return(FAobject)
})
## End methods for restrictions.independent

## Note: methods for restrictions.factanal are in GPLv3_or_later.R

## Begin methods for restrictions.orthonormal, which is an EFA model with orthonormal
## factors and exact zeros in the upper triangle of the primary pattern matrix
## Note: fit_S4, bfgs_fitS4, and bfgs_helpS4 are covered by the defaults
setMethod("restrictions2model", signature(restrictions = "restrictions.orthonormal", 
				 	      manifest = "manifest.basic"), definition = 
function(par, restrictions, manifest, lower, mapping_rule = FALSE) {
	# Pull in estimated manifest standard deviations
	FAiR_set_slot(restrictions, "Omega", 1, make_parameter(par, restrictions@Omega))

	# fill in the primary pattern matrix
	FAiR_set_slot(restrictions, "beta",  1, make_parameter(par, restrictions@beta))

	if(restrictions@beta@invalid) {
		return(list(criteria = restrictions@beta@invalid, 
			restrictions = restrictions))
	}
	else    return(list(criteria = -1, restrictions = restrictions))
})

setMethod("gr_fitS4", signature(    restrictions = "restrictions.orthonormal", 
					manifest = "manifest.basic"), definition = 
function(par, restrictions, manifest, helper, lower) {
	if(helper$marker == length(helper$fits)) {
		gradient <- FAiR_gradient_MLE(par, restrictions, manifest, 
						helper, lower)
		gradient <- c(gradient$dF_dbeta[restrictions@beta@free],
			      gradient$dF_dscale)
	}
	else { # find numeric gradient of whatever the marginal criterion is
		gradient <- FAiR_numeric_gradient(par, restrictions, manifest, 
							helper, lower)
	}
	return(gradient)
})

setMethod("create_start", signature(restrictions = "restrictions.orthonormal", 
					manifest = "manifest.basic"), definition =
function(number, start, restrictions, manifest, reject = TRUE) {
	beta_fun <- function() {
			beta <- restrictions@beta@x
			beta[restrictions@beta@free] <- runif(restrictions@beta@num_free,
								min = -1, max = 1)
			signs <- sign(colSums(beta))
			beta  <- t(t(beta) * signs)
			communalities <- rowSums(beta^2)
			beta[communalities > 1,] <- beta[communalities > 1,] /
						sqrt(communalities[communalities > 1])
			return(beta[restrictions@beta@free])
		}
	beta  <- t(replicate(number, beta_fun()))
	draws <- matrix(rnorm(number * nrow(manifest@cov), mean = log(manifest@sds),
				sd = 0.5), nrow = number)
	out <- cbind(beta, draws)
	out <- t(apply(out, 1, FUN = function(x) {
			too_small <- x < restrictions@Domains[,1]
			x[too_small]  <- restrictions@Domains[too_small,1]
			too_big <- x > restrictions@Domains[,2]
			x[too_big]  <- restrictions@Domains[too_big,2]
			return(x)
		}))
	return(out)
})

setMethod("create_FAobject", signature(restrictions = "restrictions.orthonormal",
				 	   manifest = "manifest.basic"), definition =
function(restrictions, manifest, opt, call, scores, lower, analytic) {
	fits <- opt$value
	marker <- which(fits != -1)[1]
	if(marker < length(fits)) warning("constraints did not bind")

	par <- opt$par
	model <- restrictions2model(par, restrictions, manifest, lower)
	restrictions <- model$restrictions

	S <- cormat(manifest)
	
	stuff <- FAiR_restrictions2FA(restrictions, manifest, scores)
	sorter <- stuff$sorter
	factors <- restrictions@factors[1]
	sorter <- 1:factors # temporary
	uncertainty <- try(FAiR_uncertainty(restrictions, manifest, factors, 
						sorter, NULL, analytic), FALSE)
	if(!is.list(uncertainty)) {
		warning("there was a problem calculating the variance-covariance matrix",
			" of the parameters")
		cols <- restrictions@nvars
		uncertainty <- list(vcov = matrix(NA, nrow = cols, ncol = cols),
				Jacobian = matrix(NA_real_, nrow = nrow(S), ncol = cols))
	}

	seeds <- matrix(NA_integer_, nrow = 1, ncol = 2)
	trans_mats <- array(diag(factors), c(factors, factors, 3), 
			dimnames = list(NULL, NULL, c("primary", "reference", "T")))

	FAobject <- new("FA.EFA", loadings = stuff$loadings, 
			correlations = stuff$correlations, trans_mats = trans_mats,
			uniquenesses = stuff$uniquenesses, scale = restrictions@Omega@x,
			restrictions = restrictions, Jacobian = uncertainty$Jacobian, 
			vcov = uncertainty$vcov, scores = stuff$scores, call = call,
			manifest = manifest, optimization = list(extraction = opt), 
			Lambda = coef(model$restrictions), rotated = FALSE)
	return(FAobject)
})

setMethod("restrictions2draws", signature(restrictions = "restrictions.orthonormal",
	"manifest.basic"), def = function(restrictions, manifest, vcov, nsim,
	covariances, standardized) {

	Domains <- restrictions@Domains
	S <- model.matrix(manifest, standardized = !covariances)
	vcov_chol <- chol(vcov)
	beta <- coef(restrictions)
	FC <- beta^2
	sorter <- order(colSums(FC), decreasing = TRUE)
	sorter <- 1:ncol(beta) # temporary
	free <- restrictions@beta@free
	mu <- c(beta[free], log(restrictions@Omega@x))

	if(covariances) { # setup to return reproduced covariances
		holder <- array(NA_real_, c(dim(S), nsim))
		rownames(holder) <- colnames(holder) <- rownames(S)
	}
	else {            # setup to return parameter estimates
		holder_beta   <- array(NA_real_, c(dim(beta), nsim))
		holder_Theta2 <- holder_Omega <- array(NA_real_, c(nrow(beta), 1, nsim))
		rownames(holder_beta) <- rownames(holder_Theta2) <- 
					 rownames(holder_Omega)  <- rownames(S)
		colnames(holder_beta) <- paste("F", 1:ncol(beta), sep = "")
		colnames(holder_Theta2) <- "uniqueness"
		colnames(holder_Omega)  <- "sds"
	}

	lower <- sqrt(.Machine$double.eps)
	while(nsim) {  # rejection sampling

		# draw from multivariate normal distribution
		rnorms <- rnorm(length(mu))
		rnorms <- rnorms %*% vcov_chol + mu

		# check if the draw is legal
		if(any(rnorms < Domains[,1] | rnorms > Domains[,2])) next
		model <- restrictions2model(rnorms, restrictions, manifest, lower)
		if(any(model$criteria != -1)) next
		beta <- coef(model$restrictions)

		scale <- model$restrictions@Omega@x
		if(covariances) {
			# make reproduced covariance matrix
			C <- tcrossprod(beta)
			diag(C) <- 1

			if(!standardized) C <- C * tcrossprod(scale)
			holder[,,nsim] <- C
		}
		else {
			uniquenesses <- uniquenesses(model$restrictions)
			if(!standardized) {
				beta <- beta * scale
				uniquenesses <- uniquenesses * scale^2
			}

			# stick in holders
			holder_beta[,,nsim] <- beta
			holder_Theta2[,,nsim] <- uniquenesses
			holder_Omega[,,nsim] <- scale
		}
		nsim <- nsim - 1
		if( (nsim %% 100) == 0 ) print(paste(nsim, "simulations remaining"))
	}
	if(covariances) holder <- list(covariances = holder)
	else holder <- list(beta = holder_beta[,sorter,], 
					Theta2 = holder_Theta2, Omega = holder_Omega)
	return(holder)
})

setMethod("restrictions2RAM", "restrictions.orthonormal", definition = 
function(restrictions) {
	sds <- restrictions@Omega@x
	# make first column (e.g. X -> Y)
	Path <- as.character(NULL)
	loadings <- sweep(coef(restrictions), 1, sds, "*")
	for(i in 1:ncol(loadings)) { # e.g. F1 -> Y1
		Path <- c(Path, paste("F", i, " -> ", rownames(loadings), sep = ""))
	}
	for(i in 1:ncol(loadings)) { # e.g. F1 <-> F1
		Path <- c(Path, paste("F", i, " <-> ", "F", i, sep = ""))
	}
	# tack on the uniquenesses   # e.g. Y1 <-> Y1
	Path <- c(Path, paste(rownames(loadings), "<->", rownames(loadings)))

	# make second column (parameter names)
	Parameter <- as.character(NULL)
	for(i in 1:ncol(loadings)) { # loadings
		Parameter <- c(Parameter, paste("beta_", 1:nrow(loadings), i, sep = ""))
	}
	# exclude unfree parameters in upper triangle of beta
	Parameter[!restrictions@beta@free] <- NA_character_
	Parameter <- c(Parameter, rep(NA_character_, ncol(loadings)))
	Parameter <- c(Parameter, paste("theta_", 1:nrow(loadings), sep = ""))

	# make third column (start value)
	StartValue <- as.character(c(loadings, rep(1, ncol(loadings)), 
					uniquenesses(restrictions) * sds^2))
	RAM <- cbind(Path, Parameter, StartValue)
	class(RAM) <- "mod"
	return(RAM)
})
## End methods for restrictions.orthonormal

## Begin methods for restrictions.1storder, which is a SEFA or CFA model (or transformed
## EFA model) with correlated primary factors but no second-order factors
## Note: fit_S4, bfgs_fitS4, and bfgs_helpS4 are covered by the defaults
setMethod("restrictions2model", signature(restrictions = "restrictions.1storder", 
					      manifest = "manifest.basic"), definition = 
function(par, restrictions, manifest, lower, mapping_rule) {
	# Pull in estimated manifest standard deviations
	FAiR_set_slot(restrictions, "Omega", 1, make_parameter(par, restrictions@Omega))

	# Make correlation matrix among primary factors
	FAiR_set_slot(restrictions, "Phi", TRUE, 
			make_parameter(par, restrictions@Phi, lower))

	# Make primary pattern matrix
	FAiR_set_slot(restrictions, "beta", 1, make_parameter(par, restrictions@beta,
						restrictions@Phi@x, mapping_rule))

	criteria <- rep(-1, 2)
	if(restrictions@Phi@invalid)  criteria[1] <- restrictions@Phi@invalid
	if(restrictions@beta@invalid) criteria[2] <- restrictions@beta@invalid
	return(list(criteria = criteria, restrictions = restrictions))
})

setMethod("gr_fitS4", signature(restrictions = "restrictions.1storder", 
				    manifest = "manifest.basic"), 
definition = function(par, restrictions, manifest, helper, lower) {
	SEFA <- restrictions@model == "SEFA"
	if(helper$marker == length(helper$fits)) {
		gradient <- FAiR_gradient_MLE(par, restrictions, manifest, helper, lower)

		# sum derivatives when there are equality restrictions on coefficients
		if(l <- length(restrictions@beta@equalities)) for(i in 1:l) {
			gradient$dF_dbeta[restrictions@beta@equalities[[i]]@free] <- 
			gradient$dF_dbeta[restrictions@beta@equalities[[i]]@free] + 
		    sum(gradient$dF_dbeta[restrictions@beta@equalities[[i]]@fixed])
		}

		gradient <- c(  gradient$dF_dPhi[restrictions@Phi@free] * 2,
				gradient$dF_dbeta[restrictions@beta@free],
				gradient$dF_dscale)
		if(SEFA) { # quadratic loss around zero for small coefficients if !done
			gradient[helper$squashed] <- 2*par[helper$squashed] * !helper$done
		}
	}
	else { # find numeric gradient of whatever the marginal criterion is
		gradient <- FAiR_numeric_gradient(par, restrictions, manifest, 
							helper, lower)
	}
	return(gradient)
})

setMethod("create_start", signature(restrictions = "restrictions.1storder", 
				manifest = "manifest.basic"), definition =
function(number, start, restrictions, manifest, reject = TRUE) {
	if(restrictions@model == "CFA") reject <- FALSE
	if(!reject) {
		out <- apply(restrictions@Domains, 1, FUN = function(x) {
				return(runif(number, min = x[1], max = x[2]))
			})
		return(out)
	}
	factors <- restrictions@factors[1]
	Lambda  <- FAiR_get_loadings(1 - c(start), manifest@cor, factors)
	if(require(GPArotation)) {
		efa <- GPForth(Lambda, method = "tandemII")
		Lambda <- loadings(efa)
	}
	else Lambda <- Lambda %*% t(solve(FAiR_Landahl(factors)))
	signs <- ifelse(colSums(Lambda) >= 0, 1, -1)
	if(any(signs != 1)) Lambda <- sweep(Lambda, 2, signs, FUN = "*")
	out  <- matrix(NA_real_, nrow = number, ncol = restrictions@nvars)
	counter <- 0
	while(number) {
		counter <- counter + 1
		Tmat <- FAiR_make_Tmat(runif(factors^2, min = -1))
		beta  <- try(Lambda %*% t(solve(Tmat)))
		if(!is.matrix(beta)) next
		signs <- ifelse(colSums(beta) >= 0, 1, -1)
		if(any(signs != 1)) {
			beta <- t(t(beta) * signs)
			Phi_temp <- crossprod(Tmat)
			Phi_temp <- Phi_temp * tcrossprod(signs)
		}
		else Phi_temp <- crossprod(Tmat)
		Phi_temp <- FAiR_nearPD(Phi_temp)
		FC <- beta * (beta %*% Phi_temp)
		sorter <- order(colSums(FC), decreasing = TRUE)
		beta <- beta[,sorter]
		Phi_temp <- Phi_temp[sorter, sorter]
		scale <- rnorm(nrow(beta), log(manifest@sds), sd = .5)
		par <- c(Phi_temp[restrictions@Phi@free], beta[restrictions@beta@free],
			scale[restrictions@Omega@free])
		par <-  ifelse( par <   restrictions@Domains[,1], 
					restrictions@Domains[,1],
			ifelse( par >   restrictions@Domains[,2], 
					restrictions@Domains[,2], par))
		fits <- fitS4(par, restrictions, manifest, 
				    lower = sqrt(.Machine$double.eps), TRUE)
		if(which(fits != -1)[1] != length(fits)) next
		out[number,] <- par
		if( (number %% 100) == 0 ) { 
			print(paste(number, "starting values left to create"))
		}
		number <- number - 1
	}
	cat("The realized probability of satisfying the specified restrictions is:", 
		nrow(out) / counter)
	Phi <- diag(factors[1])
	par <- c(Phi[lower.tri(Phi)], Lambda, log(manifest@sds))
	par <- par[restrictions@free]
	out[1,] <- par
	return(out)
})

setMethod("create_FAobject", signature(restrictions = "restrictions.1storder",
					   manifest = "manifest.basic"), definition =
function(restrictions, manifest, opt, call, scores, lower, analytic) {
	par <- opt$par
	model <- restrictions2model(par, restrictions, manifest, lower, 
					mapping_rule = TRUE)
	restrictions <- model$restrictions

	SEFA <- restrictions@model == "SEFA"
	S <- cormat(manifest)

	factors <- restrictions@factors[1]
	fits <- opt$value
	marker <- which(fits != -1)[1]
	if(marker < length(fits)) warning("constraints did not bind")

	stuff  <- FAiR_restrictions2FA(restrictions, manifest, scores)
	sorter <- stuff$sorter
	sorter <- 1:factors # temporary
	uncertainty <- try(FAiR_uncertainty(restrictions, manifest, factors, 
						sorter, NULL, analytic), FALSE)
	if(!is.list(uncertainty)) {
		warning("there was a problem calculating the variance-covariance matrix",
			" of the parameters")
		cols <- restrictions@nvars
		uncertainty <- list(vcov = matrix(NA, nrow = cols, ncol = cols),
				Jacobian = matrix(NA_real_, nrow = nrow(S), ncol = cols))
	}
	
	seeds <- matrix(NA_integer_, nrow = 1, ncol = 2)
	if(SEFA) { # recalculate degrees of freedom correctly
		help <- bfgs_helpS4(par, restrictions, manifest, done = TRUE, lower)
		restrictions@Domains <- restrictions@Domains[!help$squashed,]
		restrictions@beta@select  <- restrictions@beta@select[!help$squashed]
		restrictions@Phi@select   <- restrictions@Phi@select[!help$squashed]
		restrictions@Omega@select <- restrictions@Omega@select[!help$squashed]
		n <- nrow(S)
		dof <- 0.5 * n * (n + 1) - nrow(restrictions@Domains)
		restrictions@dof <- as.integer(dof)
		restrictions@nvars <- nrow(restrictions@Domains)
		restrictions@beta@squashed <- !coef(restrictions) & restrictions@beta@free
		dimnames(restrictions@beta@squashed) <- dimnames(restrictions@beta@x)
		mark <- lower.tri(restrictions@Phi@x)
		restrictions@free <- c(restrictions@Phi@free[mark], 
					 restrictions@beta@free & 
					!restrictions@beta@squashed,
					 restrictions@Omega@free)
	}
	FAobject <- new("FA", loadings = stuff$loadings, 
			correlations = stuff$correlations,
			uniquenesses = stuff$uniquenesses, scale = restrictions@Omega@x,
			restrictions = restrictions, Jacobian = uncertainty$Jacobian, 
			vcov = uncertainty$vcov, scores = stuff$scores, call = call,
			manifest = manifest, optimization = list(extraction = opt))
	return(FAobject)
})

setMethod("restrictions2draws", signature(restrictions = "restrictions.1storder",
	"manifest.basic"), def = function(restrictions, manifest, vcov, nsim, 
	covariances, standardized) {
	S <- model.matrix(manifest, standardized = !covariances)
	Domains <- restrictions@Domains
	vcov_chol <- chol(vcov)
	Phi  <- cormat(restrictions)
	restrictions@Phi@x <- diag(rep(0.5, ncol(Phi)))
	beta <- coef(restrictions)
	FC <- beta * (beta %*% Phi)
	sorter <- order(colSums(FC), decreasing = TRUE)
	sorter <- 1:ncol(beta) # temporary
	if(restrictions@model == "SEFA") {
# 		free <- restrictions@beta@free
# 		zero <- c(beta == 0)
# 		squashed <- zero[free]
# 		really_free <- free & !zero
# 		restrictions@beta@free <- really_free
# 		temp_Domains <- Domains[restrictions@beta@select,]
# 		temp_Domains <- temp_Domains[!squashed,]
# 		mark <- 0.5 * ncol(Phi) * (ncol(Phi) - 1)
# 		Domains <- rbind(Domains[1:mark,], temp_Domains, 
# 				 Domains[restrictions@Omega@select,])
# 		select <- c(rep(TRUE, sum(really_free)), 
# 			    rep(FALSE, sum(restrictions@Omega@select)))
# 		restrictions@beta@select  <- c(rep(FALSE, mark),  select)
# 		restrictions@Omega@select <- c(rep(FALSE, mark), !select)
# 		restrictions@model <- "CFA"
		really_free <- restrictions@beta@free & !restrictions@beta@squashed
		restrictions@beta@free <- really_free
		mu <- c(Phi[restrictions@Phi@free], beta[really_free],
			log(restrictions@Omega@x[restrictions@Omega@free]))
	}
	else {  # CFA
		mu <- c(Phi[restrictions@Phi@free], beta[restrictions@beta@free], 
			log(restrictions@Omega@x[restrictions@Omega@free]))
	}

	if(covariances) {
		holder <- array(NA_real_, c(dim(S), nsim))
		rownames(holder) <- colnames(holder) <- rownames(S)
	}
	else {
		holder_Phi <- array(NA_real_, c(dim(Phi), nsim))
		mark <- 0.5 * ncol(Phi) * (ncol(Phi) - 1)
		rownames(holder_Phi) <- colnames(holder_Phi) <- 
							paste("F", 1:ncol(Phi), sep = "")
		holder_beta <- array(NA_real_, c(dim(beta), nsim))
		holder_Theta2 <- holder_Omega <- array(NA_real_, c(nrow(beta), 1, nsim))
		rownames(holder_beta) <- rownames(holder_Theta2) <- 
					 rownames(holder_Omega)  <- rownames(S)
		colnames(holder_beta) <- paste("F", 1:ncol(beta), sep = "")
		colnames(holder_Theta2) <- "uniqueness"
		colnames(holder_Omega)  <- "sds"
	}
	lower <- sqrt(.Machine$double.eps)

	while(nsim) {
		rnorms <- rnorm(length(mu))
		rnorms <- c(rnorms %*% vcov_chol + mu)

		if(any(rnorms < Domains[,1] | rnorms > Domains[,2])) next
		model <- restrictions2model(rnorms, restrictions, manifest, lower, 
						mapping_rule = FALSE)
		if(any(model$criteria != -1)) next
		Phi   <- cormat(model$restrictions)
		beta  <- coef(  model$restrictions)
		scale <- model$restrictions@Omega@x
		
		criteria <- FAiR_lexical_driver(model$restrictions, manifest, lower)
		if(which(criteria != -1)[1] != length(criteria)) next

		C <- fitted(model$restrictions, reduced = TRUE)
		diag(C) <- 1
		if(covariances) {
			if(!standardized) C <- C * tcrossprod(scale)
			holder[,,nsim] <- C
		}
		else {
			uniquenesses <- uniquenesses(model$restrictions)
			if(!standardized) {
				beta <- beta * scale
				uniquenesses <- uniquenesses * scale^2
			}
			holder_beta[,,nsim]   <- beta
			holder_Theta2[,,nsim] <- uniquenesses
			holder_Omega[,,nsim]  <- scale
			holder_Phi[,,nsim]    <- Phi
		}
		nsim <- nsim - 1
		if( (nsim %% 100) == 0 ) print(paste(nsim, "simulations remaining"))
	}
	if(covariances) holder <- list(covariances = holder)
	else holder <- list(Phi = holder_Phi[sorter,sorter,], Omega = holder_Omega, 
				beta = holder_beta[,sorter,], Theta2 = holder_Theta2)
	return(holder)
})

setMethod("restrictions2RAM", "restrictions.1storder", definition = 
function(restrictions) {
	sds <- restrictions@Omega@x
	Path <- as.character(NULL)
	loadings <- sweep(coef(restrictions), 1, sds, "*")
	for(i in 1:ncol(loadings)) {
		Path <- c(Path, paste("F", i, " -> ", rownames(loadings), sep = ""))
	}
	for(i in 1:ncol(loadings)) for(j in i:ncol(loadings)) {
		Path <- c(Path, paste("F", j, " <-> ", "F", i, sep = ""))
	}
	Path <- c(Path, paste(rownames(loadings), "<->", rownames(loadings)))

	Parameter <- as.character(NULL)
	for(i in 1:ncol(loadings)) {
		Parameter <- c(Parameter, paste("beta_", 1:nrow(loadings), i, sep = ""))
	}
	Parameter[!restrictions@beta@free] <- NA_character_
	Parameter[loadings == 0] <- NA_character_
	for(i in 1:ncol(loadings)) for(j in i:ncol(loadings)) {
		if( (i == j) | !restrictions@Phi@free[j,i] ) {
			Parameter <- c(Parameter, NA_character_)
		}
		else    Parameter <- c(Parameter, paste("phi_", j, i, sep = ""))
	}
	Parameter <- c(Parameter, paste("theta_", 1:nrow(loadings), sep = ""))

	Phi <- cormat(restrictions)
	diag(Phi) <- 1
	StartValue <- as.character(c(loadings, Phi[lower.tri(Phi, TRUE)], 
					uniquenesses(restrictions) * sds^2))
	RAM <- cbind(Path, Parameter, StartValue)
	class(RAM) <- "mod"
	return(RAM)
})
## End methods for restrictions.1storder

## Begin methods for restrictions.1storder.EFA
setMethod("restrictions2model", signature(restrictions = "restrictions.1storder.EFA", 
				 	      manifest = "manifest.basic"), definition = 
function(par, restrictions, manifest, lower, mapping_rule = FALSE) {

	# Check Domains
	if(any(par < restrictions@Domains[,1]) | any(par > restrictions@Domains[,2])) {
		return(list(criteria = rep(0, 2)))
	}

	# Pull in estimated manifest standard deviations
	FAiR_set_slot(restrictions, "Omega", 1, make_parameter(par, restrictions@Omega))

	# fill in the preliminary primary pattern matrix
	restrictions@Lambda[restrictions@beta@free] <- par[1:restrictions@beta@num_free]

	if(any(colSums(restrictions@Lambda) < 0)) {
		return(list(criteria = c(-1, 0)))
	}

	return(list(criteria = rep(-1, 2), restrictions = restrictions))
})

setMethod("restrictions2draws", signature(restrictions = "restrictions.1storder.EFA",
	"manifest.basic"), def = function(restrictions, manifest, vcov, nsim,
	covariances, standardized) {

	Domains <- restrictions@Domains
	S <- model.matrix(manifest, standardized = !covariances)
	vcov_chol <- chol(vcov)
	Phi  <- cormat(restrictions)
	beta <- coef(restrictions)
	FC <- beta * (beta %*% Phi)
	sorter <- order(colSums(FC), decreasing = TRUE)
	sorter <- 1:ncol(beta) # temporary
	free <- restrictions@beta@free
	Lambda <- restrictions@Lambda
	mu <- c(Lambda[free], log(restrictions@Omega@x))

	if(covariances) { # setup to return reproduced covariances
		holder <- array(NA_real_, c(dim(S), nsim))
		rownames(holder) <- colnames(holder) <- rownames(S)
	}
	else {            # setup to return parameter estimates
		if(!require(GPArotation)) {
			stop("the suggested GPArotation package must be ",
				"installed to obtain confidence intervals for ",
				"rotated loadings")
		}
		Target <- beta
		holder_Phi <- array(NA_real_, c(dim(Phi), nsim))
		rownames(holder_Phi) <- colnames(holder_Phi) <- 
						paste("F", 1:ncol(Phi), sep = "")
		holder_beta   <- array(NA_real_, c(dim(beta), nsim))
		holder_Theta2 <- holder_Omega <- array(NA_real_, c(nrow(beta), 1, nsim))
		rownames(holder_beta) <- rownames(holder_Theta2) <- 
					 rownames(holder_Omega)  <- rownames(S)
		colnames(holder_beta) <- paste("F", 1:ncol(beta), sep = "")
		colnames(holder_Theta2) <- "uniqueness"
		colnames(holder_Omega)  <- "sds"
	}

	lower <- sqrt(.Machine$double.eps)
	while(nsim) {  # rejection sampling

		# draw from multivariate normal distribution
		rnorms <- rnorm(length(mu))
		rnorms <- rnorms %*% vcov_chol + mu

		# check if the draw is legal
		model <- restrictions2model(rnorms, restrictions, manifest, lower, 
						mapping_rule = FALSE)
		if(any(model$criteria != -1)) next

		scale <- model$restrictions@Omega@x
		if(covariances) {
			# make reproduced covariance matrix
			C <- tcrossprod(model$restrictions@Lambda)
			diag(C) <- 1

			if(!standardized) C <- C * tcrossprod(scale)
			holder[,,nsim] <- C
		}
		else {  # rotate parameters
			Tmat2 <- FAiR_quick_Rotate(model$restrictions, Target)
			if(any(is.na(Tmat2))) next
			Phi  <- crossprod(Tmat2)
			beta <- Lambda %*% Tmat2 %*% chol2inv(chol(Phi))
			holder_Phi[,,nsim] <- Phi

			uniquenesses <- uniquenesses(model$restrictions)
			if(!standardized) {
				beta <- beta * scale
				uniquenesses <- uniquenesses * scale^2
			}

			# stick in holders
			holder_beta[,,nsim] <- beta
			holder_Theta2[,,nsim] <- uniquenesses
			holder_Omega[,,nsim] <- scale
		}
		nsim <- nsim - 1
		if( (nsim %% 100) == 0 ) print(paste(nsim, "simulations remaining"))
	}
	if(covariances) holder <- list(covariances = holder)
	else holder <- list(Phi = holder_Phi[sorter,sorter,], 
				beta = holder_beta[,sorter,],
				Theta2 = holder_Theta2, Omega = holder_Omega)
	return(holder)
})
## End methods for restrictions.1storder.EFA

## Begin methods for restrictions.general, which is a SEFA or CFA model with one 
## second-order factor
## Note: fit_S4, bfgs_fitS4, and bfgs_helpS4 are covered by the defaults
setMethod("restrictions2model", signature(restrictions = "restrictions.general", 
					      manifest = "manifest.basic"), definition = 
function(par, restrictions, manifest, lower, mapping_rule) {
	# Pull in estimated manifest standard deviations
	FAiR_set_slot(restrictions, "Omega", 1, make_parameter(par, restrictions@Omega))

	# Make primary pattern matrix at level 2 (general factor)
	FAiR_set_slot(restrictions, "Delta", 1, make_parameter(par, restrictions@Delta))

	criteria <- rep(-1, 3)
	if(restrictions@Delta@invalid) criteria[1] <- restrictions@Delta@invalid

	# Make correlation matrix among primary factors at level 1
	Phi <- tcrossprod(restrictions@Delta@x)
	diag(Phi) <- 1
	check <- eigen(Phi, TRUE, TRUE)$values[ncol(Phi)]
	if(check < lower) criteria[2] <- -check
	restrictions@Phi@x <- Phi

	# Make primary pattern matrix at level 1
	FAiR_set_slot(restrictions, "beta", 1, make_parameter(par, restrictions@beta,
							restrictions@Phi@x, mapping_rule))
	if(restrictions@beta@invalid) criteria[3] <- restrictions@beta@invalid

	return(list(criteria = criteria, restrictions = restrictions))
})

setMethod("gr_fitS4", signature(    restrictions = "restrictions.general", 
					manifest = "manifest.basic"), definition = 
function(par, restrictions, manifest, helper, lower) {
	SEFA <- restrictions@model == "SEFA"
	if(helper$marker == length(helper$fits)) {
		gradient <- FAiR_gradient_MLE(par, restrictions, manifest, helper, lower)

		# Calculate derivatives with respect to second-order coefficients
		dF_dPhi <- gradient$dF_dPhi
		dF_dDelta <- dF_dPhi %*% par[restrictions@Delta@select] * 2

		# sum derivatives when there are equality restrictions at level 2
		if(l <- length(restrictions@Delta@equalities)) for(i in 1:l) {
			gradient$dF_dDelta[restrictions@Delta@equalities[[i]]@free] <- 
			gradient$dF_dDelta[restrictions@Delta@equalities[[i]]@free] + 
		    sum(gradient$dF_dDelta[restrictions@Delta@equalities[[i]]@fixed])
		}

		# sum derivatives when there are equality restrictions at level 1
		if(l <- length(restrictions@beta@equalities)) for(i in 1:l) {
			gradient$dF_dbeta[restrictions@beta@equalities[[i]]@free] <- 
			gradient$dF_dbeta[restrictions@beta@equalities[[i]]@free] + 
		    sum(gradient$dF_dbeta[restrictions@beta@equalities[[i]]@fixed])
		}
		dF_dDelta <- dF_dDelta[restrictions@Delta@free]
		gradient <- c(dF_dDelta, gradient$dF_dbeta[restrictions@beta@free], 
				gradient$dF_dscale)

		if(SEFA) { # quadratic loss around zero for small coefficients if !done
			gradient[helper$squashed] <- 2*par[helper$squashed] * !helper$done
		}
	}
	else { # find numeric gradient of whatever the marginal criterion is
		gradient <- FAiR_numeric_gradient(par, restrictions, manifest, 
							helper, lower)
	}
	return(gradient)
})

setMethod("create_start", signature(restrictions = "restrictions.general", 
					manifest = "manifest.basic"), definition =
function(number, start, restrictions, manifest, reject = TRUE) {
	factors <- restrictions@factors[1]
	if(restrictions@model == "CFA") reject <- FALSE
	if(!reject) {
		out <- apply(restrictions@Domains, 1, FUN = function(x) {
				return(runif(number, min = max(x[1], -0.75), 
						     max = min(x[2],  0.75)))
			})
		return(out)
	}
	Lambda  <- FAiR_get_loadings(1 - c(start), manifest@cor, factors)
	if(require(GPArotation)) {
		efa <- GPForth(Lambda, method = "tandemII")
		Lambda <- loadings(efa)
	}
	else Lambda <- Lambda %*% t(solve(FAiR_Landahl(factors)))
	signs <- ifelse(colSums(Lambda) >= 0, 1, -1)
	if(any(signs != 1)) Lambda <- sweep(Lambda, 2, signs, FUN = "*")
	out  <- matrix(NA_real_, nrow = number, ncol = restrictions@nvars)
	counter <- 0
	while(number) {
		counter <- counter + 1
		Delta <- as.matrix(runif(factors, min = restrictions@Domains[1:factors,1],
					 max = restrictions@Domains[1:factors,2]))
		if(sum(Delta) < 0) Delta <- Delta * -1
		Phi <- tcrossprod(Delta)
		diag(Phi) <- 1
		Tmat <- try(chol(Phi), silent = TRUE)
		if(!is.matrix(Tmat)) next
		beta <- Lambda %*% t(chol2inv(Tmat))
		signs <- ifelse(colSums(beta) >= 0, 1, -1)
		if(any(signs != 1)) beta <- t(t(beta) * signs)
		FC <- beta * (beta %*% Phi)
		sorter <- order(colSums(FC), decreasing = TRUE)
		beta <- beta[,sorter]
		Delta <- Delta[sorter,]
		scale <- rnorm(nrow(beta), log(manifest@sds), sd = .5)
		par <- c(Delta[restrictions@Delta@free], beta[restrictions@beta@free], 
			 scale[restrictions@Omega@free])
		par <- ifelse(  par <   restrictions@Domains[,1], 
					restrictions@Domains[,1],
			ifelse( par >   restrictions@Domains[,2], 
					restrictions@Domains[,2], par))
		fits <- fitS4(par, restrictions, manifest, 
				lower = sqrt(.Machine$double.eps), mapping_rule = TRUE)
		if(which(fits != -1)[1] != length(fits)) next
		out[number,] <- par
		if( (number %% 100) == 0 ) { 
			print(paste(number, "starting values left to create"))
		}
		number <- number - 1
	}
	cat("The realized probability of satisfying the specified restrictions is:", 
		nrow(out) / counter)
	par <- c(rep(0, restrictions@factors[1]), Lambda, log(manifest@sds))
	par <- par[restrictions@free]
	out[1,] <- par
	return(out)
})

setMethod("create_FAobject", signature(restrictions = "restrictions.general", 
					   manifest = "manifest.basic"), definition =
function(restrictions, manifest, opt, call, scores, lower, analytic) {
	par <- opt$par
	model <- restrictions2model(par, restrictions, manifest, lower, 
					mapping_rule = TRUE)
	restrictions <- model$restrictions
	SEFA <- restrictions@model == "SEFA"
	S <- cormat(manifest)

	factors <- restrictions@factors[1]
	fits <- opt$value
	marker <- which(fits != -1)[1]
	if(marker < length(fits)) warning("constraints did not bind")

	Delta <- loadings(restrictions, level = 2)
	uniquenesses_2nd <- 1 - c(Delta^2)
	names(uniquenesses_2nd) <- rownames(Delta)
	loadings_2nd <- Delta

	stuff <- FAiR_restrictions2FA(restrictions, manifest, scores)
	sorter <- stuff$sorter
	sorter <- 1:nrow(Delta) # temporary
	uncertainty <- try(FAiR_uncertainty(restrictions, manifest, factors, 
						sorter, NULL, analytic), TRUE)
	if(!is.list(uncertainty)) {
		warning("there was a problem calculating the variance-covariance matrix",
			" of the parameters")
		cols <- restrictions@nvars
		uncertainty <- list(vcov = matrix(NA, nrow = cols, ncol = cols),
				Jacobian = matrix(NA_real_, nrow = nrow(S), ncol = cols))
	}
	
	seeds <- matrix(NA_integer_, nrow = 1, ncol = 2)
	loadings_2nd <- loadings_2nd[sorter,, drop = TRUE]
	loadings_2nd <- array(loadings_2nd, dim = c(length(loadings_2nd), 1, 5))
	loadings_2nd[,,5] <- loadings_2nd[,,5]^2
	dimnames(loadings_2nd) <- list(rownames(Delta), colnames(Delta),
					c("PP", "RS", "PS", "RP", "FC"))
	if(SEFA) { # recalculate degrees of freedom correctly
		help <- bfgs_helpS4(par, restrictions, manifest, done = TRUE, lower)
		restrictions@Domains <- restrictions@Domains[!help$squashed,]
		restrictions@beta@select  <- restrictions@beta@select[!help$squashed]
		restrictions@Delta@select <- restrictions@Delta@select[!help$squashed]
		restrictions@Omega@select <- restrictions@Omega@select[!help$squashed]
		n <- nrow(S)
		dof <- 0.5 * n * (n + 1) - nrow(restrictions@Domains)
		restrictions@dof <- as.integer(dof)
		restrictions@nvars <- nrow(restrictions@Domains)
		restrictions@beta@squashed <- !coef(restrictions) & restrictions@beta@free
		dimnames(restrictions@beta@squashed) <- dimnames(restrictions@beta@x)
		restrictions@free <- c(  restrictions@Delta@free, 
					 restrictions@beta@free & 
					!restrictions@beta@squashed,
					 restrictions@Omega@free )
	}

	FAobject <- new("FA.general", loadings = stuff$loadings, 
			loadings_2nd = loadings_2nd, correlations = stuff$correlations,
			uniquenesses = stuff$uniquenesses, call = call, 
			scale = restrictions@Omega@x, restrictions = restrictions,
			uniquenesses_2nd = uniquenesses_2nd, vcov = uncertainty$vcov, 
			Jacobian = uncertainty$Jacobian, scores = stuff$scores, 
			manifest = manifest, optimization = list(extraction = opt))
	return(FAobject)
})

setMethod("restrictions2draws", signature(restrictions = "restrictions.general",
	"manifest.basic"), def = function(restrictions, manifest, vcov, nsim, 
	covariances, standardized) {
	S <- model.matrix(manifest, standardized = !covariances)
	Domains <- restrictions@Domains
	vcov_chol <- chol(vcov)
	Delta <- loadings(restrictions, level = 2)
	Phi <- cormat(restrictions)
	beta  <- loadings(restrictions, level = 1)
	FC <- beta * (beta %*% Phi)
	sorter <- order(colSums(FC), decreasing = TRUE)
	sorter <- 1:ncol(beta) # temporary
	if(restrictions@model == "SEFA") {
# 		free <- restrictions@beta@free
# 		zero <- c(beta == 0)
# 		squashed <- zero[free]
# 		really_free <- free & !zero
# 		restrictions@beta@free <- really_free
# 		temp_Domains <- Domains[restrictions@beta@select,]
# 		temp_Domains <- temp_Domains[!squashed,]
# 		Delta_free <- restrictions@Delta@free
# 		Domains <- rbind(Domains[restrictions@Delta@select,], temp_Domains, 
# 				 Domains[restrictions@Omega@select,])
# 		select <- c(rep(FALSE, sum(Delta_free)), rep(TRUE, sum(really_free)), 
# 			    rep(FALSE, nrow(S)))
# 		restrictions@beta@select <- select
# 		select <- !select
# 		select[1:sum(Delta_free)] <- FALSE
# 		restrictions@Omega@select <- select
# 		restrictions@model <- "CFA"
		really_free <- restrictions@beta@free & !restrictions@beta@squashed
		restrictions@beta@free <- really_free
		mu <- c(Delta[restrictions@Delta@free], beta[really_free], 
			log(restrictions@Omega@x[restrictions@Omega@free]))
	}
	else {
		mu <- c(Delta[restrictions@Delta@free], beta[restrictions@beta@free], 
			log(restrictions@Omega@x[restrictions@Omega@free]))
	}
	
	if(covariances) {
		holder <- array(NA_real_, c(dim(S), nsim))
		rownames(holder) <- colnames(holder) <- rownames(S)
	}
	else {
		holder_Delta <- array(NA_real_, c(dim(Delta), nsim))
		rownames(holder_Delta) <- rownames(Delta)
		colnames(holder_Delta) <- "G1"
		holder_Phi <- array(NA_real_, c(nrow(Delta), nrow(Delta), nsim))
		rownames(holder_Phi) <- colnames(holder_Phi) <- 
							paste("F", 1:ncol(Phi), sep = "")
		holder_beta <- array(NA_real_, c(dim(beta), nsim))
		holder_Theta2 <- holder_Omega <- array(NA_real_, c(nrow(beta), 1, nsim))
		rownames(holder_beta) <- rownames(holder_Theta2) <- 
					 rownames(holder_Omega)  <- rownames(S)
		colnames(holder_beta) <- paste("F", 1:ncol(beta), sep = "")
		colnames(holder_Theta2) <- "uniqueness"
		colnames(holder_Omega)  <- "sds"
	}
	lower <- sqrt(.Machine$double.eps)
	
	while(nsim) {
		rnorms <- rnorm(length(mu))
		rnorms <- rnorms %*% vcov_chol + mu

		if(any(rnorms < Domains[,1] | rnorms > Domains[,2])) next
		model <- restrictions2model(rnorms, restrictions, manifest, lower, 
						mapping_rule = FALSE)
		if(any(model$criteria != -1)) next
		Delta <- loadings(model$restrictions, level = 2)
		beta  <- loadings(model$restrictions)
		scale <- model$restrictions@Omega@x
		criteria <- FAiR_lexical_driver(model$restrictions, manifest, lower)
		if(which(criteria != -1)[1] != length(criteria)) next

		C <- fitted(model$restrictions, reduced = TRUE)
		diag(C) <- 1
		if(covariances) {
			if(!standardized) C <- C * tcrossprod(scale)
			holder[,,nsim] <- C
		}
		else {
			uniquenesses <- uniquenesses(model$restrictions)
			if(!standardized) {
				beta <- beta * scale
				uniquenesses <- uniquenesses * scale^2
			}
			holder_beta[,,nsim]   <- beta
			holder_Theta2[,,nsim] <- uniquenesses
			holder_Omega[,,nsim]  <- scale
			holder_Delta[,,nsim]  <- Delta
			holder_Phi[,,nsim]    <- cormat(model$restrictions)
		}
		nsim <- nsim - 1
		if( (nsim %% 100) == 0 ) print(paste(nsim, "simulations remaining"))
	}
	if(covariances) holder <- list(covariances = holder)
	else holder <- list(Delta = holder_Delta[sorter,,,drop = FALSE], 
				Phi = holder_Phi[sorter,sorter,], Omega = holder_Omega,
				beta = holder_beta[,sorter,,drop = FALSE], 
				Theta2 = holder_Theta2)
	return(holder)
})

setMethod("restrictions2RAM", "restrictions.general", definition = 
function(restrictions) {
	sds <- restrictions@Omega@x
	Path <- as.character(NULL)
	Delta <- loadings(restrictions, level = 2)
	scale.factor <- sqrt(Delta^2 + 1)
	Delta <- sweep(Delta, 1, scale.factor, "*")
	loadings <- sweep(loadings(restrictions), 1, sds, "*")
	loadings <- sweep(loadings, 2, scale.factor, "/")
	for(i in 1:ncol(loadings)) {
		Path <- c(Path, paste("F", i, " -> ", rownames(loadings), sep = ""))
	}
	for(i in 1:ncol(loadings)) {
		Path <- c(Path, paste("F", ncol(loadings) + 1, " -> ", "F", i, sep = ""))
	}
	for(i in 1:(ncol(loadings) + 1)) {
		Path <- c(Path, paste("F", i, " <-> ", "F", i, sep = ""))
	}
	Path <- c(Path, paste(rownames(loadings), "<->", rownames(loadings)))

	Parameter <- as.character(NULL)
	for(i in 1:ncol(loadings)) {
		Parameter <- c(Parameter, paste("beta_", 1:nrow(loadings), i, sep = ""))
	}
	Parameter[!restrictions@beta@free] <- NA_character_
	Parameter[loadings == 0] <- NA_character_
	for(i in 1:ncol(loadings)) {
		Parameter <- c(Parameter, if(restrictions@Delta@free[i]) 
				paste("Delta_", i, sep = "") else NA_character_)
	}
	Parameter <- c(Parameter, rep(NA_character_, ncol(loadings) + 1))
	Parameter <- c(Parameter, paste("theta_", 1:nrow(loadings), sep = ""))

	StartValue <- as.character(c(loadings, Delta, rep(1, ncol(loadings) + 1),
					uniquenesses(restrictions) * sds^2))
	RAM <- cbind(Path, Parameter, StartValue)
	warning("RAM version of the model has different normalizations")
	class(RAM) <- "mod"
	return(RAM)
})
## End methods for restrictions.general

## begin methods for restrictions.2ndorder, which is a SEFA or CFA model with multiple
## second-order factors
setMethod("restrictions2model", signature(restrictions = "restrictions.2ndorder", 
					      manifest = "manifest.basic"), definition = 
function(par, restrictions, manifest, lower, mapping_rule) {
	SEFA <- restrictions@model == "SEFA"
	criteria <- rep(-1, 5)

	# Pull in estimated manifest standard deviations
	FAiR_set_slot(restrictions, "Omega", 1, make_parameter(par, restrictions@Omega))

	# Make correlation matrix among primary factors at level 2
	FAiR_set_slot(restrictions, "Xi", TRUE, 
			make_parameter(par, restrictions@Xi, lower))
	if(restrictions@Xi@invalid) criteria[1] <- restrictions@Xi@invalid

	# Make primary pattern matrix at level 2 (general factor)
	FAiR_set_slot(restrictions, "Delta", 1, make_parameter(par, restrictions@Delta,
							restrictions@Xi@x, mapping_rule))
	if(restrictions@Delta@invalid) criteria[2] <- restrictions@Delta@invalid

	# Make correlation matrix among primary factors at level 1
	Phi   <- crossprod(chol(restrictions@Xi@x) %*% t(restrictions@Delta@x))
	check <- max(diag(Phi))
	if(check > 1) criteria[3] <- check # Heywood case
	diag(Phi) <- 1
	
	check <- eigen(Phi, TRUE, TRUE)$values[ncol(Phi)]
	if(check < lower) criteria[4] <- -check # not posdef
	restrictions@Phi@x <- Phi

	# Make primary pattern matrix at level 1
	FAiR_set_slot(restrictions, "beta", 1, make_parameter(par, restrictions@beta,
							restrictions@Phi@x, mapping_rule))
	if(restrictions@beta@invalid) criteria[5] <- restrictions@beta@invalid

	return(list(criteria = criteria, restrictions = restrictions))
})

setMethod("gr_fitS4", signature(restrictions = "restrictions.2ndorder",
				    manifest = "manifest.basic"), definition = 
function(par, restrictions, manifest, helper, lower) {
	SEFA <- restrictions@model == "SEFA"
	if(helper$marker == length(helper$fits)) {
		gradient <- FAiR_gradient_MLE(par, restrictions, manifest, helper, lower)

		# Get matrices for second level
		model <- restrictions2model(par, restrictions, manifest, lower, 
						mapping_rule = FALSE)
		Delta <- loadings(model$restrictions, level = 2)
		Xi    <- cormat(model$restrictions, level = 2)

		# Calculate gradient for second-order coefficients
		dF_dPhi <- gradient$dF_dPhi
		dF_dDelta <- dF_dPhi %*% Delta %*% Xi * 2

		# Calculate gradient for free elements of second-order cormat
		dF_dXi <- t(Delta) %*% dF_dPhi %*% Delta * 2

		# sum derivatives when there are equality restrictions on Delta
		if(l <- length(restrictions@Delta@equalities)) for(i in 1:l) {
			gradient$dF_dDelta[restrictions@Delta@equalities[[i]]@free] <- 
			gradient$dF_dDelta[restrictions@Delta@equalities[[i]]@free] + 
		    sum(gradient$dF_dDelta[restrictions@Delta@equalities[[i]]@fixed])
		}

		# sum derivatives when there are equality restrictions on beta
		if(l <- length(restrictions@beta@equalities)) for(i in 1:l) {
			gradient$dF_dbeta[restrictions@beta@equalities[[i]]@free] <- 
			gradient$dF_dbeta[restrictions@beta@equalities[[i]]@free] + 
		    sum(gradient$dF_dbeta[restrictions@beta@equalities[[i]]@fixed])
		}

		gradient <- c(dF_dXi[lower.tri(Xi)], dF_dDelta[restrictions@Delta@free], 
				gradient$dF_dbeta[restrictions@beta@free],
				gradient$dF_dscale)
		if(SEFA) { # quadratic loss around zero for small coefficients if !done
			gradient[helper$squashed] <- 2*par[helper$squashed] * !helper$done
		}
	}
	else { # find numeric gradient of whatever the marginal criterion is
		gradient <- FAiR_numeric_gradient(par, restrictions, manifest, 
							helper, lower)
	}
	return(gradient)
})

setMethod("create_start", signature(restrictions = "restrictions.2ndorder", 
					manifest = "manifest.basic"), definition =
function( number, start, restrictions, manifest, reject = TRUE) {
	factors  <- restrictions@factors[1]
	factors2 <- restrictions@factors[2]
	if(restrictions@model == "CFA") reject <- FALSE
	if(!reject) {
		out <- apply(restrictions@Domains, 1, FUN = function(x) {
				return(runif(number, min = max(x[1], -0.75), 
						     max = min(x[2],  0.75)))
			})
		return(out)
	}
	Lambda  <- FAiR_get_loadings(1 - c(start), cormat(manifest), factors)
	if(require(GPArotation)) {
		efa <- GPForth(Lambda, method = "tandemII")
		Lambda <- loadings(efa)
	}
	else Lambda <- Lambda %*% t(solve(FAiR_Landahl(factors)))
	signs <- ifelse(colSums(Lambda) >= 0, 1, -1)
	if(any(signs != 1)) Lambda <- sweep(Lambda, 2, signs, FUN = "*")
	out  <- matrix(NA_real_, nrow = number, ncol = restrictions@nvars)
	mark1 <- 0.5 * factors2 * (factors2 - 1)
	mark2 <- mark1 + sum(restrictions@Delta@free)
	counter <- 0
	while(number) {
		counter <- counter + 1
		draws <- runif(mark2, min = restrictions@Domains[1:mark2, 1],
				      max = restrictions@Domains[1:mark2, 2])
		draws <- draws / 2
		Xi <- diag(factors2)
		Xi[lower.tri(Xi)] <- draws[1:mark1]
		Xi <- Xi + t(Xi)
		diag(Xi) <- 1
		Xi <- FAiR_nearPD(Xi)
		Delta <- restrictions@Delta@Delta
		Delta[restrictions@Delta@free] <- draws[(mark1 + 1):mark2]
		signs <- ifelse(colSums(Delta) >= 0, 1, -1)
		if(any(signs != 1)) Delta <- t(t(Delta) * signs)
		FC <- Delta * (Delta %*% Xi)
		sorter2 <- order(colSums(FC), decreasing = TRUE)
		Phi <- try(crossprod(chol(Xi) %*% t(Delta)), TRUE)
		if(!is.matrix(Phi)) next
		if(any(diag(Phi) > 1)) next
		diag(Phi) <- 1
		Tmat <- try(chol(Phi), silent = TRUE)
		if(!is.matrix(Tmat)) next
		beta <- Lambda %*% t(chol2inv(Tmat))
		signs <- ifelse(colSums(beta) >= 0, 1, -1)
		if(any(signs != 1)) beta <- t(t(beta) * signs)
		FC <- beta * (beta %*% Phi)
		sorter <- order(colSums(FC), decreasing = TRUE)
		beta <- beta[,sorter]
		Delta <- Delta[sorter,sorter2]
		Xi <- Xi[sorter2, sorter2]
		scale <- rnorm(nrow(beta), log(manifest@sds), sd = .1)
		par <- c(Xi[lower.tri(Xi)], Delta[restrictions@Delta@free], 
			beta[restrictions@beta@free], scale[restrictions@Omega@free])
		par <-  ifelse(par <    restrictions@Domains[,1],
					restrictions@Domains[,1],
			ifelse(par >    restrictions@Domains[,2],
					restrictions@Domains[,2], par))
		fits <- fitS4(par, restrictions, manifest, 
				lower = sqrt(.Machine$double.eps), mapping_rule = TRUE)
		if(which(fits != -1)[1] != length(fits)) next
		out[number,] <- par
		if( (number %% 100) == 0 ) {
			print(paste(number, "starting values left to create"))
		}
		number <- number - 1
	}
	cat("The realized probability of satisfying the specified restrictions is:", 
		nrow(out) / counter)
	Xi <- diag(restrictions@factors[2])
	par <- c(Xi[lower.tri(Xi)], rep(0, prod(restrictions@factors)), Lambda, 
		log(manifest@sds))
	par <- par[restrictions@free]
	out[1,] <- par
	return(out)
})

setMethod("create_FAobject", signature(restrictions = "restrictions.2ndorder", 
					   manifest = "manifest.basic"), definition =
function(restrictions, manifest, opt, call, scores, lower, analytic) {
	par <- opt$par
	model <- restrictions2model(par, restrictions, manifest, lower, 
					mapping_rule = TRUE)
	restrictions <- model$restrictions
	S <- cormat(manifest)
	factors  <- restrictions@factors[1]
	factors2 <- restrictions@factors[2]
	SEFA <- restrictions@model == "SEFA"
	fits <- opt$value
	marker <- which(fits != -1)[1]
	if(marker < length(fits)) warning("constraints did not bind")

	Xi    <-   cormat(restrictions, level = 2)
	Delta <- loadings(restrictions, level = 2)
	Phi   <-   cormat(restrictions, level = 1)
	beta  <- loadings(restrictions, level = 1)

	# Make loadings_2nd and correlations_2nd
	Pi2 <- Delta %*% Xi
	FC2 <- Pi2 * Delta
	sorter_2nd <- order(colSums(FC2), decreasing = TRUE)
	sorter_2nd <- 1:ncol(Delta) # temporary
	uniquenesses_2nd <- 1 - rowSums(FC2)

	Xi_inv <- chol2inv(chol(Xi))
	D2 <- 1/sqrt(diag(Xi_inv))
	Psi2 <- Xi_inv * tcrossprod(D2)
	Upsilon2 <- sweep(Delta, 2, D2, FUN = "*")
	RP2 <- sweep(Pi2, 2, D2, FUN = "/")

	loadings_2nd <- array(cbind(Delta, Upsilon2, Pi2, RP2, FC2),
				dim = c(factors, factors2, 5),
				dimnames = list(NULL, NULL, 
					c("PP", "RS", "PS", "RP", "FC")))

	correlations_2nd <- array(cbind(Xi, Psi2, diag(D2)),
			dim = c(factors2, factors2, 3), 
			dimnames = list(NULL, NULL, c("PF", "RF", "PR")))

	stuff <- FAiR_restrictions2FA(restrictions, manifest, scores)
	sorter <- stuff$sorter
	sorter <- 1:ncol(beta) # temporary
	uncertainty <- try(FAiR_uncertainty(restrictions, manifest, c(factors, factors2),
						sorter, sorter_2nd, analytic), TRUE)
	if(!is.list(uncertainty)) {
		warning("there was a problem calculating the variance-covariance matrix",
			" of the parameters")
		cols <- restrictions@nvars
		uncertainty <- list(vcov = matrix(NA, nrow = cols, ncol = cols),
				Jacobian = matrix(NA_real_, nrow = nrow(S), ncol = cols))
	}

	loadings_2nd <- loadings_2nd[sorter,sorter_2nd,]
	correlations_2nd <- correlations_2nd[sorter_2nd,sorter_2nd,]
	
	seeds <- matrix(NA_integer_, nrow = 1, ncol = 2)
	if(SEFA) { # recalculate degrees of freedom correctly
		help <- bfgs_helpS4(par, restrictions, manifest, done = TRUE, lower)
		restrictions@Domains <- restrictions@Domains[!help$squashed,]
		restrictions@Xi@select    <- restrictions@Xi@select[!help$squashed]
		restrictions@beta@select  <- restrictions@beta@select[!help$squashed]
		restrictions@Delta@select <- restrictions@Delta@select[!help$squashed]
		restrictions@Omega@select <- restrictions@Omega@select[!help$squashed]
		n <- nrow(S)
		dof <- 0.5 * n * (n + 1) - nrow(restrictions@Domains)
		restrictions@dof <- as.integer(dof)
		restrictions@nvars <- nrow(restrictions@Domains)
		if(is(restrictions@beta, "parameter.coef.SEFA")) {
			restrictions@beta@squashed <- !beta & restrictions@beta@free
			dimnames(restrictions@beta@squashed) <- 
				dimnames(restrictions@beta@x)
		}
		if(is(restrictions@Delta, "parameter.coef.SEFA")) {
			restrictions@Delta@squashed <- !Delta & restrictions@Delta@free
			dimnames(restrictions@Delta@squashed) <- 
				dimnames(restrictions@Delta@x)
		}
		restrictions@free <- c(restrictions@Xi@free[lower.tri(restrictions@Xi@x)],
					 restrictions@Delta@free &
					!restrictions@Delta@squashed,
					 restrictions@beta@free & 
					!restrictions@beta@squashed,
					 restrictions@Omega@free)
	}

	FAobject <- new("FA.2ndorder", loadings = stuff$loadings, call = call, 
			loadings_2nd = loadings_2nd, correlations = stuff$correlations,
			correlations_2nd = correlations_2nd, restrictions = restrictions,
			uniquenesses = stuff$uniquenesses, vcov = uncertainty$vcov,
			scale = restrictions@Omega@x, Jacobian = uncertainty$Jacobian,
			uniquenesses_2nd = uniquenesses_2nd, scores = stuff$scores, 
			manifest = manifest, optimization = list(extraction = opt))
	return(FAobject)
})

setMethod("restrictions2draws", signature(restrictions = "restrictions.2ndorder",
	"manifest.basic"), def = function(restrictions, manifest, vcov, nsim, 
	covariances, standardized) {
	S <- model.matrix(manifest, standardized = !covariances)
	Domains <- restrictions@Domains
	vcov_chol <- chol(vcov)
	Xi    <- cormat(restrictions, level = 2)
	Phi   <- cormat(restrictions, level = 1)
	Delta <- loadings(restrictions, level = 2)
	beta  <- loadings(restrictions, level = 1)
	FC <- beta * (beta %*% Phi)
	sorter <- order(colSums(FC), decreasing = TRUE)
	sorter <- 1:ncol(beta) # temporary
	FC2 <- Delta * (Delta %*% Xi)
	sorter2 <- order(colSums(FC2), decreasing = TRUE)
	sorter2 <- 1:ncol(Delta) # temporary
	if(restrictions@model == "SEFA") {
# 		beta_free <- restrictions@beta@free
# 		beta_zero <- c(beta == 0)
# 		beta_squashed <- beta_zero[beta_free]
# 		beta_really_free <- beta_free & !beta_zero
# 		restrictions@beta@free <- beta_really_free
# 		temp_beta_Domains <- Domains[restrictions@beta@select,]
# 		temp_beta_Domains <- temp_beta_Domains[!beta_squashed,]
# 		Delta_free <- restrictions@Delta@free
# 		Delta_zero <- c(Delta == 0)
# 		Delta_squashed <- Delta_zero[Delta_free]
# 		temp_Delta_Domains <- Domains[restrictions@Delta@select,]
# 		temp_Delta_Domains <- temp_Delta_Domains[!Delta_squashed,]
# 		Delta_really_free <- Delta_free & !Delta_zero
# 		restrictions@Delta@free <- Delta_really_free
# 		mark <- 0.5 * ncol(Xi) * (ncol(Xi) - 1)
# 		Domains <- rbind(Domains[1:mark,], temp_Delta_Domains,
# 				temp_beta_Domains, Domains[restrictions@Omega@select,])
# 
# 		select <- c(rep(FALSE, mark), rep(TRUE, sum(Delta_really_free)),
# 			    rep(FALSE, sum(beta_really_free)), rep(FALSE, nrow(S)))
# 		restrictions@Delta@select <- select
# 		select <- c(rep(FALSE, mark + sum(Delta_really_free)), 
# 			    rep(TRUE, sum(beta_really_free)), rep(FALSE, nrow(S)))
# 		restrictions@beta@select <- select
# 		select <- c(rep(FALSE, mark + sum(Delta_really_free) +
# 					      sum( beta_really_free)), rep(TRUE, nrow(S)))
# 		restrictions@Omega@select <- select
# 		restrictions@model <- "CFA"
		Delta_free <- restrictions@Delta@free & !restrictions@Delta@squashed
		restrictions@Delta@free <- Delta_free
		 beta_free <- restrictions@beta@free  & !restrictions@beta@squashed
		restrictions@beta@free  <- beta_free
		mu <- c(Xi[restrictions@Xi@free], Delta[Delta_free], beta[beta_free],
			log(restrictions@Omega@x[restrictions@Omega@free]))
	}
	else {
		mu <- c(Xi[restrictions@Xi@free], Delta[restrictions@Delta@free], 
			beta[restrictions@beta@free], 
			log(restrictions@Omega@x[restrictions@Omega@free]))
	}
	if(covariances) {
		holder <- array(NA_real_, c(dim(S), nsim))
		rownames(holder) <- colnames(holder) <- rownames(S)
	}
	else {
		holder_Xi <- array(NA_real_, c(dim(Xi), nsim))
		rownames(holder_Xi) <- colnames(Xi) <- paste("G", 1:ncol(Xi), sep = "")
		holder_Delta <- array(NA_real_, c(dim(Delta), nsim))
		rownames(holder_Delta) <- paste("F", 1:nrow(Delta), sep = "")
		colnames(holder_Delta) <- paste("G", 1:ncol(Delta), sep = "")
		holder_Phi <- array(NA_real_, c(dim(Phi), nsim))
		rownames(holder_Phi) <- colnames(holder_Phi) <- 
							paste("F", 1:ncol(Phi), sep = "")
		holder_beta <- array(NA_real_, c(dim(beta), nsim))
		holder_Theta2 <- holder_Omega <- array(NA_real_, c(nrow(beta), 1, nsim))
		rownames(holder_beta) <- rownames(holder_Theta2) <- 
					 rownames(holder_Omega)  <- rownames(S)
		colnames(holder_beta) <- paste("F", 1:ncol(beta), sep = "")
		colnames(holder_Theta2) <- "uniqueness"
		colnames(holder_Omega)  <- "sds"
	}
	lower <- sqrt(.Machine$double.eps)
	
	while(nsim) {
		rnorms <- rnorm(length(mu))
		rnorms <- rnorms %*% vcov_chol + mu
		if(any(rnorms < Domains[,1] | rnorms > Domains[,2])) next
		model <- restrictions2model(rnorms, restrictions, manifest, lower, 
						mapping_rule = FALSE)
		if(any(model$criteria != -1)) next
		Xi     <- cormat(model$restrictions, level = 2)
		Delta  <-   loadings(model$restrictions, level = 2)
		Phi    <- cormat(model$restrictions, level = 1)
		beta   <-   loadings(model$restrictions, level = 1)
		scale  <- model$restrictions@Omega@x
		criteria <- FAiR_lexical_driver(model$restrictions, manifest,lower)
		if(which(criteria != -1)[1] != length(criteria)) next

		C <- fitted(model$restrictions, reduced = TRUE)
		diag(C) <- 1
		if(covariances) {
			if(!standardized) C <- C * tcrossprod(scale)
			holder[,,nsim] <- C
		}
		else {
			uniquenesses <- uniquenesses(model$restrictions)
			if(!standardized) {
				beta <- beta * scale
				uniquenesses <- uniquenesses * scale^2
			}
			holder_beta[,,nsim]   <- beta
			holder_Theta2[,,nsim] <- uniquenesses
			holder_Omega[,,nsim]  <- scale
			holder_Delta[,,nsim]  <- Delta
			holder_Xi[,,nsim]  <- Xi
			holder_Phi[,,nsim] <- Phi
		}
		nsim <- nsim - 1
		if( (nsim %% 100) == 0 ) print(paste(nsim, "simulations remaining"))
	}
	if(covariances) holder <- list(covariances = holder)
	else holder <- list(Xi = holder_Xi[sorter2, sorter2,], 
				Delta = holder_Delta[sorter,sorter2,],
				Phi = holder_Phi[sorter,sorter,],
				beta = holder_beta[,sorter,], Theta2 = holder_Theta2,
				Omega = holder_Omega)
	return(holder)
})

setMethod("restrictions2RAM", "restrictions.2ndorder", definition = 
function(restrictions) {
	sds <- restrictions@Omega@x
	Path <- as.character(NULL)
	Delta <- loadings(restrictions, level = 2)
	Xi <- cormat(restrictions, level = 2)
	FC <- Delta * (Delta %*% Xi)
	scale.factor <- sqrt(rowSums(FC) + 1)
	Delta <- sweep(Delta, 1, scale.factor, "*")
	loadings <- sweep(loadings(restrictions), 1, sds, "*")
	loadings <- sweep(loadings, 2, scale.factor, "/")
	for(i in 1:ncol(loadings)) {
		Path <- c(Path, paste("F", i, " -> ", rownames(loadings), sep = ""))
	}
	Path <- c(Path, paste("F", 1:ncol(loadings), " <-> F", 1:ncol(loadings),sep = ""))
	for(j in 1:ncol(Delta)) for(i in 1:ncol(loadings)) {
		Path <- c(Path, paste("F", ncol(loadings) + j, " -> ", "F", i, sep = ""))
	}
	for(i in 1:ncol(Delta)) for(j in i:ncol(Delta)) {
		Path <- c(Path, paste(  "F", ncol(loadings) + i, " <-> ",
					"F", ncol(loadings) + j, sep = "" ))
	}
	Path <- c(Path, paste(rownames(loadings), "<->", rownames(loadings)))

	Parameter <- as.character(NULL)
	for(i in 1:ncol(loadings)) {
		Parameter <- c(Parameter, paste("beta_", 1:nrow(loadings), i, sep = ""))
	}
	Parameter[!restrictions@beta@free] <- NA_character_
	Parameter[loadings == 0]   <- NA_character_
	Parameter <- c(Parameter, rep(NA_character_, ncol(loadings)))
	temp_Parameter <- as.character(NULL)
	for(i in 1:ncol(Delta)) {
		temp_Parameter <- c(temp_Parameter, paste("Delta_", 1:ncol(loadings), 
								i, sep = ""))
	}
	temp_Parameter[!restrictions@Delta@free] <- NA_character_
	temp_Parameter[Delta == 0] <- NA_character_
	Parameter <- c(Parameter, temp_Parameter)
	for(i in 1:ncol(Delta)) for(j in i:ncol(Delta)) {
		Parameter <- c(Parameter, paste("Xi_", j, i, sep = ""))
	}
	Parameter <- c(Parameter, paste("theta_", 1:nrow(loadings), sep = ""))

	StartValue <- as.character(c(loadings, rep(1, ncol(loadings)),
					Delta, Xi[lower.tri(Xi, TRUE)], 
					uniquenesses(restrictions) * sds^2))
	RAM <- cbind(Path, Parameter, StartValue)
	warning("RAM version of the model has different normalizations")
	class(RAM) <- "mod"
	return(RAM)
})
## End methods for restrictions.2ndorder

## The rest are various common methods; several are established in the stats4 library

# these loadings methods get the primary pattern matrix
setMethod("loadings", "restrictions", 
function(x, standardized = TRUE) {
	if(length(standardized) != 1) stop("'standardized' must be TRUE or FALSE")
	if(!is.logical(standardized)) stop("'standardized' must be TRUE or FALSE")

	beta <- coef(x)
	if(!standardized) beta <- beta * x@Omega@x
	class(beta) <- "loadings"
	return(beta)
})

setMethod("loadings", "restrictions.factanal", 
function(x, standardized = TRUE) {
	stop("loadings is not defined for objects of class 'restrictions.factanal'")
})

setMethod("loadings", "restrictions.general", # inherited by restrictions.2ndorder
function(x, standardized = TRUE, level = 1) {
	if(length(standardized) != 1) stop("'standardized' must be TRUE or FALSE")
	if(!is.logical(standardized)) stop("'standardized' must be TRUE or FALSE")
	if(length(level) != 1) stop("'level' must be either 1 or 2")
	if(!(level %in% 1:2))  stop("'level' must be either 1 or 2")
	if(!standardized && level == 2) {
		stop("requesting 'standardized = FALSE' with 'level = 2' is inconsistent")
	}

	if(level == 1) {
		mat <- x@beta@x
		if(!standardized) beta <- beta * x@Omega@x
	}
	else    mat <- x@Delta@x
	class(mat)  <- "loadings"
	return(mat)
})

# these methods get the "loadings", not necessarily the primary pattern or at level 1
setMethod("loadings", "FA", 
function(x, matrix = "PP", standardized = TRUE) {
	matrix <- toupper(matrix)
	if(length(matrix) != 1) {
		stop("'matrix' must be exactly one of 'PP', 'RS', 'PS', 'RP', or 'FC'")
	}
	if(!(matrix %in% c("PP", "RS", "PS", "RP", "FC"))) {
		stop("'matrix' must be exactly one of 'PP', 'RS', 'PS', 'RP', or 'FC'")
	}
	if(length(standardized) != 1) stop("'standardized' must be TRUE or FALSE")
	if(!is.logical(standardized)) stop("'standardized' must be TRUE or FALSE")

	mat <- as.matrix(x@loadings[,,matrix])
	if(!standardized) mat <- mat * x@scale
	class(mat) <- "loadings"
	if(!attributes(x@correlations)$orthogonal) attributes(mat)$covariance <- TRUE
	return(mat)
})

setMethod("loadings", "FA.general", # inherited by FA.2ndorder
function(x, matrix = "PP", standardized = TRUE, level = 1) {
	matrix <- toupper(matrix)
	if(length(matrix) != 1) {
		stop("'matrix' must be exactly one of 'PP', 'RS', 'PS', 'RP', or 'FC'")
	}
	if(!(matrix %in% c("PP", "RS", "PS", "RP", "FC"))) {
		stop("'matrix' must be exactly one of 'PP', 'RS', 'PS', 'RP', or 'FC'")
	}
	if(length(standardized) != 1) stop("'standardized' must be TRUE or FALSE")
	if(!is.logical(standardized)) stop("'standardized' must be TRUE or FALSE")
	if(length(level) != 1) stop("'level' must be either 1 or 2")
	if(!(level %in% 1:2))  stop("'level' must be either 1 or 2")
	if(!standardized && level == 2) {
		stop("requesting 'standardized = FALSE' with 'level = 2' is inconsistent")
	}

	if(level == 1) {
		mat <- as.matrix(x@loadings[,,matrix])
		if(!standardized) mat <- mat * x@scale
	}
	else    mat <- as.matrix(x@loadings_2nd[,,matrix])
	class(mat)  <- "loadings"
	if(level == 1) attributes(mat)$covariance <- TRUE
	return(mat)
})

# these coef methods get the primary pattern matrix at level 1
setMethod("coef", "parameter.coef",
function(object) {
	return(object@x)
})

setMethod("coef", "restrictions", 
function(object) {
	return(object@beta@x)
})

setMethod("coef", "restrictions.independent",
function(object) {
	return(as.numeric(NULL))
})

setMethod("coef", "restrictions.factanal", 
function(object) {
	stop("coef is not defined for objects of class 'restrictions.factanal'")
})

setMethod("coef", "FA",
function(object) {
	return(object@loadings[,,"PP"])
})

# these cormat methods get the correlation matrix among primary factors
setMethod("cormat", "parameter.cormat",
function(object) {
	return(object@x)
})

setMethod("cormat", "restrictions", 
function(object) {
	return(object@Phi@x)
})

setMethod("cormat", "restrictions.factanal", 
function(object) {
	return(diag(object@factors[1]))
})

setMethod("cormat", "restrictions.orthonormal", 
function(object) {
	return(diag(object@factors[1]))
})

setMethod("cormat", "restrictions.2ndorder", 
function(object, level = 1) {
	if(level == 1)      return(object@Phi@x)
	else if(level == 2) return(object@Xi@x)
})

# these cormat methods get the correlation matrix among "factors" but not necessarily
# the primary factors and not necessarily at level 1
setMethod("cormat", "FA", 
function(object, matrix = "PF") {
	matrix <- toupper(matrix)
	if(length(matrix) != 1) {
		stop("'matrix' must be exactly one of 'PF', 'RF', 'PR'")
	}
	if(!(matrix %in% c("PF", "RF", "PR"))) {
		stop("'matrix' must be exactly one of 'PF', 'RF', 'PR'")
	}

	return(object@correlations[,,toupper(matrix)])
})

setMethod("cormat", "FA.2ndorder", 
function(object, matrix = "PF", level = 1) {
	matrix <- toupper(matrix)
	if(length(matrix) != 1) {
		stop("'matrix' must be exactly one of 'PF', 'RF', or 'PR'")
	}
	if(!(matrix %in% c("PF", "RF", "PR"))) {
		stop("'matrix' must be exactly one of 'PF', 'RF', or 'PR'")
	}
	if(length(level) != 1) stop("'level' must be either 1 or 2")
	if(!(level %in% 1:2))  stop("'level' must be either 1 or 2")

	if(level == 1)      return(object@correlations[,,matrix])
	else if(level == 2) return(object@correlations_2nd[,,matrix])
})

# this cormat method just gets the correlation matrix among manifest variables
setMethod("cormat", "manifest.basic",
function(object) {
	return(object@cor)
})

# these uniquenesses methods get the uniquenesses
setMethod("uniquenesses", "restrictions", 
function(object, standardized = TRUE) {
	beta <- coef(object)
	Phi  <- cormat(object)
	uniquenesses <- 1 - rowSums(beta * (beta %*% Phi))
	if(!standardized) uniquenesses <- uniquenesses * restrictions@Omega@x^2
	return(uniquenesses)
})

setMethod("uniquenesses", "restrictions.factanal",
function(object, ...) {
	stop("uniquenesses method is not defined for restrictions.factanal")
})

setMethod("uniquenesses", "FA", 
function(object, standardized = TRUE) {
	if(length(standardized) != 1) stop("'standardized' must be TRUE or FALSE")
	if(!is.logical(standardized)) stop("'standardized' must be TRUE or FALSE")

	out <- object@uniquenesses
	if(!standardized) out <- out * object@scale^2
	return(out)
})

setMethod("uniquenesses", "FA.general", 
function(object, standardized = TRUE, level = 1) {
	if(length(standardized) != 1) stop("'standardized' must be TRUE or FALSE")
	if(!is.logical(standardized)) stop("'standardized' must be TRUE or FALSE")

	if(length(level) != 1) stop("'level' must be either 1 or 2")
	if(!(level %in% 1:2))  stop("'level' must be either 1 or 2")

	if(!standardized && level == 2) {
		stop("unstandardized uniquenesses are not possible at level 2")
	}

	if(level == 1) {
		out <- object@uniquenesses
		if(!standardized) out <- out * object@scale^2
	}
	else    out <- object@uniquenesses_2nd
	return(out)
})

# get estimated variance-covariance matrix of the estimates
setMethod("vcov", "FA", function(object) object@vcov) 

# get log-likelihood
setMethod("logLik", "FA",
function (object) {
	S <- model.matrix(object, standardized = FALSE)
	C <- fitted(object, reduced = FALSE, standardized = FALSE)
	p <- ncol(S)
	npars <- 0.5 * p * (p + 1) - df.residual(object)
	n.obs <- object@manifest@n.obs

	if(!FAiR_is.ML(object)) {
		warning("log-likelihood can only be calculated for MLEs")
		val <- NA_real_
	}
	C_inv <- chol2inv(chol(C))	
	val <- -0.5 * (n.obs - 1) * (log(det(C)) + crossprod(c(S), c(C_inv)))
	attr(val, "df")   <- npars
	attr(val, "nobs") <- n.obs
	class(val)        <- "logLik"
	return(val)
})

# get Bayesian Information Criterion (no 2 * pi correction)
# setMethod("BIC", signature(object = "FA"),
# function (object) { 
#         BIC(logLik(object))
# })

# get CIs from FAobject
setMethod("confint", signature(object = "FA"), definition = 
function (object, parm, level = 0.95, nsim = 1001, seed = NULL) {
	if(level <= 0 | level >= 1) stop("'level' must be strictly between 0 and 1")
        draws <- FA2draws(object, nsim = nsim, seed = seed, covariances = FALSE)
	if(!missing(parm)) warning("'parm' ignored")
	bounds <- c( (1 - level) / 2, 1 - (1 - level) / 2)
	CIs <- FAiR_draws2CI(draws, bounds)
	return(CIs)
})

# get CIs from summary.FA
setMethod("confint", signature(object = "summary.FA"), definition = 
function (object, parm, level = 0.95){ # get CIs
	if(level <= 0 | level >= 1) stop("'level' must be strictly between 0 and 1")
	if(!missing(parm)) warning("'parm' ignored")
	bounds <- c( (1 - level) / 2, 1 - (1 - level) / 2)
	CIs <- FAiR_draws2CI(object@draws, bounds)
	return(CIs)
})

# get profiles of discrepancy function
setMethod("profile", signature(fitted = "FA"), definition = 
function (fitted, delta = 0.05, number = 100, plot.it = TRUE, ...) {
	if(FAiR_is.EFA(fitted)) {
		stop("profile is not currently defined for EFA solutions")
	}
	outlist <- FAiR_profile(fitted, delta, number, plot.it, ...)
	return(invisible(outlist))
})

setMethod("profile", signature(fitted = "FA.general"), definition = 
function (fitted, delta = 0.05, number = 100, plot.it = TRUE, ...) {
	if(FAiR_is.EFA(fitted)) {
		stop("profile is not currently defined for EFA solutions")
	}
	Delta <- fitted@restrictions@Delta@x
	notfree <- !fitted@restrictions@Delta@free
	if(any(notfree)) {
		width <- seq(from = -delta, to = delta, length = number)
		outlist <- list()
		count <- 1
		parnames <- rownames(Delta)
		if(plot.it) {
			ask <- par("ask")
			par("ask" = TRUE)
			on.exit(par("ask" = ask))
			cat("Factors may arbitrarily be plotted in a different order",
				" than they appear in summary()\n")
		}
		for(j in 1:nrow(notfree)) if(notfree[j,1]) {
			y <- x <- rep(NA_real_, length(width))
			val <- Delta[j,1]
			for(i in 1:length(x)) {
				x[i] <- Delta[j,1] <- val + width[i]
				Phi <- tcrossprod(Delta)
				diag(Phi) <- 1
				fitted@restrictions@Phi@x <- Phi
				y[i] <- deviance(fitted)
			}
			outlist[[count]] <- list(x = x, y = y)
			names(outlist)[count] <- paste(parnames[j], p, sep = "_")
			Delta[j,1] <- val
			if(plot.it) {
				plot(x, y, type = "l", 
					ylab = "Value of Discrepancy Function",
					xlab = paste("Factor", p, "for", parnames[j]),
							"at level 2")
				abline(v = val, col = "gray", lty = "dotted")
			}
			count <- count + 1
		}
		Phi <- tcrossprod(Delta)
		diag(Phi) <- 1
		fitted@restrictions@Phi@x <- Phi
	}
	else outlist <- list()
	outlist <- c(outlist, FAiR_profile(fitted, delta, number, plot.it, ...))
	return(invisible(outlist))
})

setMethod("profile", signature(fitted = "FA.2ndorder"), definition = 
function (fitted, delta = 0.05, number = 100, plot.it = TRUE, ...) {
	if(FAiR_is.EFA(fitted)) {
		stop("profile is not currently defined for EFA solutions")
	}
	Xi_chol <- chol(cormat(fitted, level = 2))
	Delta <- fitted@restrictions@Delta@x
	notfree <-  !fitted@restrictions@Delta@free | !Delta

	width <- seq(from = -delta, to = delta, length = number)
	outlist <- list()
	count <- 1
	parnames <- rownames(Delta)
	if(plot.it) {
		ask <- par("ask")
		par("ask" = TRUE)
		on.exit(par("ask" = ask))
		cat("Factors may arbitrarily be plotted in a different order than they",
			" appear in summary()\n")
	}
	for(p in 1:ncol(notfree)) for(j in 1:nrow(notfree)) if(notfree[j,p]) {
		y <- x <- rep(NA_real_, length(width))
		val <- Delta[j,p]
		for(i in 1:length(x)) {
			x[i] <- Delta[j,p] <- val + width[i]
			Phi  <- crossprod(Xi_chol %*% t(Delta))
			diag(Phi) <- 1
			fitted@restrictions@Phi@x <- Phi
			y[i] <- deviance(fitted)
		}
		outlist[[count]] <- list(x = x, y = y)
		names(outlist)[count] <- paste(parnames[j], p, sep = "_")
		Delta[j,p] <- val
		Phi  <- crossprod(Xi_chol %*% t(Delta))
		diag(Phi) <- 1
		fitted@restrictions@Phi@x <- Phi
		if(plot.it) {
			plot(x, y, type = "l", ylab = "Value of Discrepancy Function",
				xlab = paste("Factor", p, "for", parnames[j]), ...)
			abline(v = val, col = "gray", lty = "dotted")
		}
		count <- count + 1
	}
	outlist <- c(outlist, FAiR_profile(fitted, delta, number, plot.it, ...))
	return(invisible(outlist))
})

# show (basically print) method for FA
setMethod("show", "FA", 
function(object) {
	cat("\nCall:\n")
	print(object@call)
	n.obs <- object@manifest@n.obs
	cat("\nNumber of observations: ", n.obs, "\n")
	discrepancy <- deviance(object)
	cat("\nDiscrepancy: ", discrepancy, "\n")
	show(object@restrictions)
	cat("\n")
})

# summary method for FA
setMethod("summary", signature("FA"), def =
function(object, standardized = TRUE, conf.level = 0.95, 
	nsim = 0, seed = NULL, order = 1:ncol(loadings(object)), 
	polarity = rep(1L, ncol(loadings(object))), ...) {

	if(is.null(order)) {
		FC <- loadings(object, matrix = "FC")
		order <- order(colSums(FC), decreasing = TRUE)
	}
	if(conf.level <= 0 | conf.level >= 1) {
		stop("'conf.level' must be strictly between 0 and 1")
	}
	if(nsim  <  0)              stop(" 'nsim' must be nonnegative")
	else if(nsim > 0 && any(is.na(vcov(object)))) {
		warning("simulations could not be conducted because the",
			" variance-covariance matrix was not calculated")
		nsim <- 0
	}
	factors <- ncol(loadings(object))
	if(!all(order %in% 1:factors)) {
		stop("'order' must be a vector of unique integers between 1 and ",
			factors)
	}
	else if(length(order) != factors) {
		stop("'order' must be a vector of unique integers whose length is equal",
			" to the number of factors (", factors, ")")
	}
	else if(!all(abs(polarity) == 1)) {
		stop("all elements of polarity must be +1 or -1")
	}
	else if(length(polarity) != factors) {
		stop("'polarity' must be a vector of +1 or -1 whose length is ",
			"equal to the number of factors (", factors, ")")
	}
	call <- object@call
	sds  <- object@manifest@sds

	if(nsim) {
		draws <- FA2draws(object, nsim, seed, covariances = FALSE,
					standardized = standardized)
		out <- new("summary.FA", restrictions = object@restrictions, call = call,
			draws = draws, order = as.integer(order), conf.level = conf.level,
			orthogonal = FAiR_is.orthogonal(object),
			standardized = standardized, polarity = as.integer(polarity))
	}
	else {
		out <- new("summary.FA", restrictions = object@restrictions,
			call = call, draws = list(), order = as.integer(order), 
			orthogonal = FAiR_is.orthogonal(object), conf.level = conf.level,
			standardized = standardized, polarity = as.integer(polarity))
	}
	return(out)
})

# this prints the summary to the screen
setMethod("show", "summary.FA", function(object) {
	cat("\nCall:\n")
	print(object@call)
	restrictions <- object@restrictions
	if(object@orthogonal) {
		uniquenesses <- uniquenesses(restrictions)
		if(!object@standardized) {
			uniquenesses <- uniquenesses * restrictions@Omega@x
		}
		mat <- matrix(uniquenesses, ncol = 1)
		rownames(mat) <- rownames(restrictions@beta@x)
		if(length(object@draws)) {
			mat <- cbind(mat, apply(object@draws$Theta2, 1:2, sd))
			colnames(mat) <- c("uniqueness", "std.err")
		}
		else colnames(mat) <- "uniqueness"
		print(round(mat, 3))
		return(NULL)
	}

	beta <- loadings(restrictions)
	uniquenesses <- uniquenesses(restrictions)
	Phi <- cormat(restrictions)
	FC <- beta * (beta %*% Phi)
	sorter <- order(colSums(FC), decreasing = TRUE)
	sorter <- 1:ncol(beta) # temporary
	beta <- beta[,sorter,drop = FALSE]
	beta <- sweep(beta[,object@order, drop = FALSE], 2, object@polarity, "*")
	Phi  <- Phi[ sorter,sorter]
	Phi  <- Phi[object@order, object@order] * tcrossprod(object@polarity)
	if(!object@standardized) {
		beta <- beta * restrictions@Omega@x
		uniquenesses <- uniquenesses * restrictions@Omega@x^2
	}

	mat <- cbind(beta, 0, uniquenesses)
	colnames(mat) <- c(paste("F", 1:ncol(beta), sep = ""), "", "Uniqueness")
	varnames <- rownames(beta)

	if(ncol(beta) > 1) {
		mat <- rbind( mat, 0, cbind(Phi, 0, 0) )
		varnames <- c(varnames, "", paste("F", 1:ncol(beta), sep = ""))
	}
	if(is(restrictions, "restrictions.2ndorder")) {
		Delta <- loadings(restrictions, level = 2)[sorter,]
		Delta <- sweep(Delta[object@order,], 1, object@polarity, "*")
		mat <- rbind( mat, 0, cbind(t(Delta), 0, 0) )
		Xi <- cormat(restrictions, level = 2)
		mat <- rbind(mat, 0, cbind(Xi, matrix(0, nrow(Xi), ncol(mat) - ncol(Xi))))
		varnames <- c(varnames, "", paste("G", 1:ncol(Xi), " -> F", sep = ""), "",
				paste("G", 1:ncol(Xi), sep = "") )
	}
	else if(is(restrictions, "restrictions.general")) {
		Delta <- loadings(restrictions, level = 2)[sorter,,drop = FALSE]
		Delta <- sweep(Delta[object@order,,drop = FALSE], 1, object@polarity, "*")
		mat <- rbind(mat, 0, cbind(t(Delta), 0, 0) )
		varnames <- c(varnames, "", "2nd order factors")
	}
	rownames(mat) <- varnames

	cat("\nPoint estimates (blanks, if any, are exact zeros):\n")
	FAiR_print.loadings(mat)

	if(length(object@draws) == 0) return(NULL)

	CI_list <- confint(object, level = object@conf.level)
	pegged  <- sapply(CI_list, FUN = function(x) x[,,1] != x[,,2])

	## Upper bounds
	mat_CI <- mat
	CI_beta <- CI_list$beta[,,2] * pegged$beta
	CI_beta <- sweep(CI_beta[,object@order], 2, object@polarity, "*")
	mat_CI[1:nrow(beta), 1:ncol(beta)] <- CI_beta
	mat_CI[1:nrow(beta), ncol(mat)]    <- CI_list$Theta2[,,2] * pegged$Theta2
	if(ncol(beta) > 1) {
		mark_start <- nrow(beta) + 2
		mark_end   <- mark_start + nrow(Phi) - 1
		diag(pegged$Phi) <- TRUE
		CI_Phi <- CI_list$Phi[,,2] * pegged$Phi
		CI_Phi <- CI_Phi[object@order,object@order] * tcrossprod(object@polarity)
		mat_CI[mark_start:mark_end, 1:ncol(Phi)] <- CI_Phi
	}
	if(is(restrictions, "restrictions.2ndorder")) {
		mark_start <- mark_end + 2
		mark_end   <- mark_start + nrow(Delta) - 1
		CI_Delta <- CI_list$Delta[,,2] * pegged$Delta
		CI_Delta <- sweep(CI_Delta[object@order,], 1, object@polarity, "*")
		mat_CI[mark_start:mark_end, 1:ncol(Delta)] <- CI_Delta
		mark_start <- mark_end + 2
		mark_end   <- mark_start + nrow(Xi) - 1
		mat_CI[mark_start:mark_end, 1:ncol(Xi)] <- CI_list$Xi[,,2] * pegged$Xi
	}
	else if(is(restrictions, "restrictions.general")) {
		mark_start <- mark_end + 2
		mark_end   <- mark_start
		CI_Delta <- c(CI_list$Delta[,,2] * pegged$Delta)
		CI_Delta <- CI_Delta[object@order] * object@polarity
		mat_CI[mark_start:mark_end, 1:nrow(Delta)] <- CI_Delta
	}

	cat("\nUpper confidence bounds (blanks, if any, are restricted)\n")
	FAiR_print.loadings(mat_CI)

	## Lower bounds
	CI_beta <- CI_list$beta[,,1] * pegged$beta
	CI_beta <- sweep(CI_beta[,object@order], 2, object@polarity, "*")
	mat_CI[1:nrow(beta), 1:ncol(beta)] <- CI_beta
	mat_CI[1:nrow(beta), ncol(mat)]    <- CI_list$Theta2[,,1] * pegged$Theta2
	if(ncol(beta) > 1) {
		mark_start <- nrow(beta) + 2
		mark_end   <- mark_start + ncol(Phi) - 1
		CI_Phi <- CI_list$Phi[,,1] * pegged$Phi
		CI_Phi <- CI_Phi[object@order,object@order] * tcrossprod(object@polarity)
		mat_CI[mark_start:mark_end, 1:ncol(Phi)] <- CI_Phi
	}
	if(is(restrictions, "restrictions.2ndorder")) {
		mark_start <- mark_end + 2
		mark_end   <- mark_start + nrow(Delta) - 1
		CI_Delta <- CI_list$Delta[,,1] * pegged$Delta
		CI_Delta <- sweep(CI_Delta[object@order,], 1, object@polarity, "*")
		mat_CI[mark_start:mark_end, 1:ncol(Delta)] <-  CI_Delta
		mark_start <- mark_end + 2
		mark_end   <- mark_start + nrow(Xi) - 1
		mat_CI[mark_start:mark_end, 1:ncol(Xi)] <- CI_list$Xi[,,1] * pegged$Xi
	}
	else if(is(restrictions, "restrictions.general")) {
		mark_start <- mark_end + 2
		mark_end   <- mark_start
		CI_Delta <- c(CI_list$Delta[,,1]   * pegged$Delta)
		CI_Delta <- CI_Delta[object@order] * object@polarity
		mat_CI[mark_start:mark_end, 1:nrow(Delta)] <- CI_Delta
	}

	cat("\nLower confidence bounds (blanks, if any, are restricted)\n")
	FAiR_print.loadings(mat_CI)

	return(NULL)
})

# plot method for FA
setMethod("plot", signature(x = "summary.FA"), definition = 
function (x, y = NULL, ...) { 
	heatmap(t(coef(x@restrictions)))
})

# show() methods for parameter.*
setMethod("show", signature(object = "parameter.coef"), definition = 
function(object) {
	onefactor <- ncol(coef(object)) == 1
	bound <- ifelse(onefactor, 1, 1.5)
	if(length(object@Domains)) {
		if(any( ( (lowers <- object@Domains[,,1]) != -bound)[object@free] ) | 
		   any( ( (uppers <- object@Domains[,,2]) !=  bound)[object@free] ) ) {
			bounds <- cbind(lowers, 0, uppers)
			cat("Lower and upper bounds on coefficients\n")
			FAiR_print.loadings(bounds)
		}
		else cat("All coefficients on the [", -bound, ",", bound, "] interval\n")
	}
	else cat("domains on coefficients have not been specified yet\n")

	if(any(!object@free)) {
		fixed <- coef(object) + .Machine$double.eps
		fixed[object@free] <- 0
		cat("Fixed coefficients\n")
		FAiR_print.loadings(fixed)
	}

	if(l <- length(object@equalities)) {
		equal <- matrix("", nrow = nrow(coef(object)), ncol = ncol(coef(object)))
		for(i in 1:l) {
			equal[object@equalities[[i]]@fixed] <- letters[i]
			equal[object@equalities[[i]]@free]  <- letters[i]
		}
		cat("Equality restrictions (same letters indicates equality)\n")
		print(equal)
	}

	if(is(object, "parameter.coef.SEFA")) {
		coef_args <- formals(object@mapping_rule)
		cat("\nZeros per factor\n")
		if(!is.numeric(coef_args$zeros)) {
			coef_args$zeros <- rep(ncol(object@x), ncol(object@x))
		}
		mat <- matrix(coef_args$zeros, nrow = 1)
		rownames(mat) <- "zeros"
		colnames(mat) <- LETTERS[1:ncol(coef(object))]
		print(mat)

		if(!identical(names(coef_args), names(formals(mapping_rule)))) {
			cat("Customized mapping rule\n")
			return(invisible(NULL))
		}

		if(!is.na(coef_args$row_complexity[1])) {
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
			cat("Mapping rule: ")
			if(coef_args$communality) {
				cat("At least one zero per factor on a",
					"high communality variable\n")
			}
			else if(coef_args$quasi_Yates) {
				cat("Encourage cohyperplanarity\n")
			}
			else if(coef_args$weak_Thurstone) {
				cat("Weak simple structure\n")
			}
			else if(coef_args$Butler) {
				cat("Unifactorial basis\n")
			}
			else if(coef_args$viral) 
				cat(0.5 * factors * (factors - 1),
					"outcomes each of complexity", factors - 2, "\n")
			else cat("default\n")
		}
	}
})

setMethod("show", signature(object = "parameter.cormat"), definition = 
function(object) {
	cormat <- cormat(object)

	if(!length(object@Domains)) {
		cat("Domains of the factor intercorrelations have not yet been specified")
		cat("\n")
		return(invisible(NULL))
	}
	else if(all(object@Domains[,,1][object@free] <= (-1 + .Machine$double.eps)) &
		all(object@Domains[,,2][object@free] >= (-1 - .Machine$double.eps))) {
		cat("All free factor intercorrelations are on the [-1,1] interval\n")
	}
	else {
		Domains <- object@Domains[,,1]
		Domains[upper.tri(Domains)] <- t(object@Domains[,,2])[upper.tri(Domains)]
		diag(Domains) <- 1
		cat(" Upper bounds on factor intercorrelations are in the upper triangle",
		    "\n",
		     "lower bounds on factor intercorrelations are in the lower triangle",
		    "\n")
		print(Domains, digits = 3)
	}

	if(any(!(object@free[lower.tri(object@free)]))) {
		fixed <- cormat
		fixed[object@free] <- 0
		fixed[upper.tri(fixed)] <- 0
		diag(fixed) <- 1
		cat("Fixed elements of the factor correlation matrix\n")
		FAiR_print.loadings(fixed)
	}
})

# All these show() methods for restrictions.* give information about the restrictions
setMethod("show", signature(object = "restrictions"), 
	definition = function (object) str(object))

setMethod("show", signature(object = "restrictions.factanal"), definition = 
function (object) {
	cat("\nExploratory factor anaylsis\n")
	cat(object@factors[1], "factors\n")
	cat(df.residual(object), "degrees of freedom\n")
})

setMethod("show", signature(object = "restrictions.orthonormal"), definition = 
function (object) {
	cat("\nExploratory factor anaylsis\n")
	cat(object@factors[1], "factors\n")
	cat(df.residual(object), "degrees of freedom\n")
})

setMethod("show", signature(object = "restrictions.1storder"), definition = 
function (object) {
	factors <- object@factors[1]
	if(object@model == "SEFA") {
		cat("\nSemi-exploratory factor analysis with ", factors, "factors\n")
	}
	else    cat("\nConfirmatory factor analysis with ", factors, "factors\n")

	show(object@Phi)
	cat("\n")
	show(object@beta)
	cat("\n")
	FAiR_show_constraints(object@criteria)
	cat("\n", object@dof, "degrees of freedom\n")
})

setMethod("show", signature(object = "restrictions.1storder.EFA"), definition = 
function(object) {
	cat("\nExploratory factor anaylsis\n")
	cat(object@factors[1], "factors\n")
	FAiR_show_constraints(object@Tcriteria, EFA = TRUE)
	cat(df.residual(object), "degrees of freedom\n")
})

setMethod("show", signature(object = "restrictions.general"), definition = 
function (object) {
	factors <- object@factors[1]
	if(object@model == "SEFA") {
		cat("\nSemi-exploratory factor anaylsis with ", factors, "first-order ",
			"factors and one second-order factor\n")
	}
	else {
		cat("\nConfirmatory factor anaylsis with ", factors, "first-order ",
			"factors and one second-order factor\n")
	}
	cat("Level two\n")
	show(object@Delta)
	cat("\nLevel one\n")
	show(object@beta)
	cat("\n")
	FAiR_show_constraints(object@criteria)
	cat("\n", object@dof, "degrees of freedom\n")
})

setMethod("show", signature(object = "restrictions.2ndorder"), definition = 
function (object) {
	factors <- object@factors
	if(object@model == "SEFA") {
		cat("\nSemi-exploratory factor anaylsis with ", factors[1], "first-order",
			" factors and ", factors[2], "second-order factor\n")
	}
	else {
		cat("\nConfirmatory factor anaylsis with ", factors[1], "first-order ",
			"factors and ", factors[2], "second-order factor\n")
	}
	cat("Level two\n")
	show(object@Xi)
	cat("\n")
	show(object@Delta)
	cat("\nLevel one\n")
	show(object@beta)
	cat("\n")
	FAiR_show_constraints(object@criteria)
	cat("\n", object@dof, "degrees of freedom\n")
})

# these show() methods for manifest.* give information about the sample covariance matrix
setMethod("show", "manifest.basic", definition = 
function(object) {
	FAiR_show_manifest(object)
})

setMethod("show", "manifest.data", definition = 
function(object) {
	FAiR_show_manifest(object)

	if(require(mvnormtest)) {
		if(!any(is.na(object@X))) print(mshapiro.test(object@X))
		else cat("mshapiro.test cannot handle NAs in data")
	}
	else cat("for more tests, install the mvnormtest package\n")
	if(require(energy)) { 
		if(!any(is.na(object@X))) print(mvnorm.etest(object@X))
		else cat("mvnorm.etest cannot handle NAs in data")
	}
	else cat("for more tests, install the energy package\n")
})

setMethod("show", "manifest.data.ordinal", definition = 
function(object) {
	FAiR_show_manifest(object)
})

setMethod("show", "manifest.data.ranks", definition = 
function(object) {
	FAiR_show_manifest(object)
})

setMethod("show", "equality_restriction", definition = 
function(object) {
	mat <- matrix("", nrow = object@dims[1], ncol = object@dims[2])
	if(object@level == 1) {
		if(length(object@rownames) > 1) rownames(mat) <- object@rownames
		colnames(mat) <- paste("F", 1:ncol(mat), sep = "")
	}
	else {
		rownames(mat) <- paste("F", 1:nrow(mat), sep = "")
		colnames(mat) <- paste("G", 1:ncol(mat), sep = "")
	}

	mat[object@fixed] <- mat[object@free] <- "equal"
	cat("Equality restrictions at level ", object@level, ":\n")
	print(mat)
})

# plot cormat
setMethod("plot", signature(x = "manifest.basic", y = "missing"),
function(x, ...) { 
	heatmap(cormat(x), ...)
})

# plot mcd stuff and cormat
setMethod("plot", signature(x = "manifest.data.mcd", y = "missing"),
function(x, ...) {
	plot(x@mcd, ask = TRUE, ...)
	heatmap(cormat(x), ...)
})

# advanced scree plots on C (not S)
setMethod("screeplot", "FA",
function(x, ...) { # keep ... for when Gilles adds them to plotnScree
	if(!require(nFactors)) {
		stop("screeplot requires that the nFactors package to be installed")
	}
	C <- fitted(x, reduced = TRUE, standardized = TRUE)
	if(nrow(C) < 5) {
		stop("analyzing eigenvalues is fairly meaningless when the number of",
			" manifest variables is four or fewer")
	}
	eigenvalues <- eigen(C, TRUE, TRUE)$values
	variables <- length(eigenvalues)
	nsubjects <- x@manifest@n.obs
	if(is.numeric(nsubjects)) {
		aparallel <- parallel(var = variables, subject = nsubjects)$eigen$qevpea
		results   <- nScree(eig = eigenvalues, aparallel = aparallel)
	}
	else    results   <- nScree(eig = eigenvalues)
	plotnScree(results, main = paste("Non Graphical Solutions to Scree Test\n",
			"Based on the Reduced Correlation Matrix", sep = ""))
	return(invisible(results))
})

# advanced scree plots on S (not C)
setMethod("screeplot", "manifest.basic",
function(x, ...) {
	if(!require(nFactors)) {
		stop("screeplot requires that the nFactors package to be installed")
	}
	S <- cormat(x)
	if(nrow(S) < 5) {
		stop("analyzing eigenvalues is fairly meaningless when the number of",
			" manifest variables is four or fewer")
	}
	eigenvalues <- eigen(S, TRUE, TRUE)$values
	variables <- length(eigenvalues)
	nsubjects <- x@n.obs
	if(is.numeric(nsubjects)) {
		aparallel <- parallel(var = variables, subject = nsubjects)$eigen$qevpea
		results   <- nScree(eig = eigenvalues, aparallel = aparallel)
	}
	else    results   <- nScree(eig = eigenvalues)
	plotnScree(results)
	return(invisible(results))
})

# advanced scree plots on S (not C) and Gamma
setMethod("screeplot", "manifest.data",
function(x, ...) {
	if(!require(nFactors)) {
		stop("screeplot requires that the nFactors package to be installed")
	}
	S <- cormat(x)
	if(nrow(S) < 5) {
		stop("analyzing eigenvalues is fairly meaningless when the number of",
			" manifest variables is four or fewer")
	}
	eigenvalues <- eigen(S, TRUE, TRUE)$values
	variables <- length(eigenvalues)
	nsubjects <- x@n.obs
	aparallel <- parallel(var = variables, subject = nsubjects)$eigen$qevpea
	results   <- nScree(eig = eigenvalues, aparallel = aparallel)
	title <- "Non Graphical Solutions to Scree Test on Manifest Correlations"
	plotnScree(results, main = title)
	if(is(x@acov, "diagonalMatrix")) return(invisible(results))

	Gamma <- cov2cor(x@acov)
	eigenvalues <- eigen(Gamma, TRUE, TRUE)$values
	variables <- length(eigenvalues)
	aparallel <- parallel(var = variables, subject = nsubjects)$eigen$qevpea
	results_Gamma <- nScree(eig = eigenvalues, aparallel = aparallel)
	title <- "Non Graphical Solutions to Scree Test on 4th-order Correlations"
	plotnScree(results_Gamma, main = title)
	return(invisible(results))
})

# predict outcomes
setMethod("predict", "FA", def = 
function(object) {
	if(!FAiR_is.manifest.data(object)) {
		stop("cannot predict outcomes because the raw data are not available")
	}
	if(is.na(object@scores[1,1])) {
		warning("factor scores were not requested in the call to Factanal();",
			" using Anderson-Rubin scores")
		scores <- FAiR_scores("Anderson-Rubin", object@manifest, coef(object),
				      cormat(object), uniquenesses(object))
	}
	else scores <- object@scores
	predictions <- scores %*% coef(object)
# 	if(is(object@manifest, "manifest.data.ranks")) {
# 		predictions <- apply(predictions, 2, rank)
# 	}
	return(predictions)
})

# these plot methods produce DAGs by post-hoc hacking of FAiR_DAG() which is basically
# just Bill Revelle's fa.graph() function
setMethod("plot", signature(x = "FA", y = "missing"),
function (x, out.file = NULL, ...) {
	graph <- FAiR_DAG(x, ...)
	if(!FAiR_is.orthogonal(x)) {
		Phi <- cormat(x)
		factor_names <- paste("F", 1:ncol(Phi), sep = "")
		for(j in 1:(ncol(Phi) - 1)) for(i in (j+1):(ncol(Phi))) {
			graph$x <- addEdge(factor_names[j], factor_names[i], graph$x, 1)
			temp_label <- round(Phi[j, i], 1)
			names(temp_label) <- paste(factor_names[j], "~", factor_names[i],
							sep = "")
			graph$edgeAttrs$label <- c(graph$edgeAttrs$label, temp_label)

			graph$x <- addEdge(factor_names[i], factor_names[j], graph$x, 1)
			temp_label <- round(Phi[i, j], 1)
			names(temp_label) <- paste(factor_names[i], "~", factor_names[j],
							sep = "")
			graph$edgeAttrs$label <- c(graph$edgeAttrs$label, temp_label)

		}
	}
	if (!is.null(out.file)) {
		names(graph[[1]]) <- "graph"
		graph$filename <- out.file
		do.call(toDot, args = graph)
	}
	else {
		graph$recipEdges <- "combined"
		do.call(plot,  args = graph)
	}
	invisible(graph[[1]])
})

setMethod("plot", signature(x = "FA.general", y = "missing"),
function (x, out.file = NULL, ...) {
	graph  <- FAiR_DAG(x, ...)
	graph2 <- FAiR_DAG(x, level = 2, ...)
	graph$nodeAttrs$shape <- c(graph2$nodeAttrs$shape, graph$nodeAttrs$shape)
	graph$nodeAttrs$rank  <- c(graph2$nodeAttrs$rank,  graph$nodeAttrs$rank)
	graph$edgeAttrs$label <- c(graph2$edgeAttrs$label, graph$edgeAttrs$label)
	graph[[1]] <- join(graph2[[1]], graph[[1]])
	if (!is.null(out.file)) {
		names(graph[[1]]) <- "graph"
		graph$filename <- out.file
		do.call(toDot, args = graph)
	}
	else    do.call(plot,  args = graph)
	invisible(graph[[1]])
})

setMethod("plot", signature(x = "FA.2ndorder", y = "missing"),
function (x, out.file = NULL, ...) {
	graph  <- FAiR_DAG(x, ...)
	graph2 <- FAiR_DAG(x, level = 2, ...)
	Xi <- cormat(x, level = 2)
	factor_names <- paste("G", 1:ncol(Xi), sep = "")
	for(j in 1:(ncol(Xi) - 1)) for(i in (j+1):(ncol(Xi))) {
		graph2$x <- addEdge(factor_names[j], factor_names[i], graph2$x, 1)
		temp_label <- round(Xi[j, i], 1)
		names(temp_label) <- paste(factor_names[j], "~", factor_names[i],
						sep = "")
		graph2$edgeAttrs$label <- c(graph2$edgeAttrs$label, temp_label)

		graph2$x <- addEdge(factor_names[i], factor_names[j], graph2$x, 1)
		temp_label <- round(Xi[i, j], 1)
		names(temp_label) <- paste(factor_names[i], "~", factor_names[j],
						sep = "")
		graph2$edgeAttrs$label <- c(graph2$edgeAttrs$label, temp_label)
	}

	graph$nodeAttrs$shape <- c(graph2$nodeAttrs$shape, graph$nodeAttrs$shape)
	graph$nodeAttrs$rank  <- c(graph2$nodeAttrs$rank,  graph$nodeAttrs$rank)
	graph$edgeAttrs$label <- c(graph2$edgeAttrs$label, graph$edgeAttrs$label)
	graph[[1]] <- join(graph2[[1]], graph[[1]])
	if (!is.null(out.file)) {
		names(graph[[1]]) <- "graph"
		graph$filename <- out.file
		do.call(toDot, args = graph)
	}
	else    do.call(plot,  args = graph)
	invisible(graph[[1]])
})
