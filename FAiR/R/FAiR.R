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

## This file defines original functions that are meant to be accessed by users directly.

## NOTE: This file is intended to be read with 90 columns and 8 space tabs

## Factanal() estimates all factor analysis models
Factanal <-
function(manifest, restrictions, scores = "none", seeds = 12345, 
	lower = sqrt(.Machine$double.eps), analytic = TRUE, 
	reject = TRUE, NelderMead = TRUE, impatient = FALSE, ...) {
## Arguments
#  manifest: an object that inherits from manifest class; see make_manifest
#  restrictions: an object that inherits from restrictions class; see make_restrictions
#  scores: character indicating whether / how to calculate factor scores
#  seeds: PRNG seeds for the unif.seed and int.seed arguments of genoud()
#  lower: lower bound on uniquenesses in factanal()-style EFA models and the lower bound
#         on eigenvalues when determining whether a matrix is computationally posdef
#  analytic: use analytic gradients, etc.?
#  reject: reject starting values that fail one or more constraints?
# NelderMead: use method = "Nelder-Mead" in a call to optim() to polish the solution?
# impatient: logical indicating whether to skip the slow parts
#    ...: arguments that get passed to genoud()

	## Preliminaries
	if(!is(manifest, "manifest")) {
		stop("'manifest' must inherit from class 'manifest'")
	}
	if(!is(restrictions, "restrictions")) {
		stop("'restrictions' must inherit from class 'restrictions'")
	}
	if(lower < 0) {
		stop("'lower' must be nonnegative")
	}
        kall <- match.call()

	S <- cormat(manifest)
	factors  <- restrictions@factors
	scores   <-  match.arg(scores, c("none", "regression", "Bartlett", "Thurstone",
					"Ledermann", "Anderson-Rubin", "McDonald",
					"Krinjen", "Takeuchi", "Harman"))

	if(analytic) {
		safe <- FAiR_analytic_safe(restrictions)
		if(!safe) {
			analytic <- FALSE
			warning("'analytic' coerced to FALSE because analytic gradients",
				" are not possible in this case")
		}
	}

	## copy restrictions to give it new memory
	restrictions_copy <- FAiR_copy_restrictions(restrictions)

	# shortcut when the user wants the factanal() behavior
	if(is(restrictions_copy, "restrictions.factanal") && impatient) {
		opt <- factanal(covmat = cormat(manifest), scores = "none",
				factors = factors[1],
				rotation = "none", control = list(lower = lower))
		opt$par <- opt$uniquenesses
		FAobject <- create_FAobject(restrictions_copy, manifest, opt, 
						kall, scores, lower, analytic)
		FAobject@seeds <- matrix(NA_integer_, ncol = 2,
					dimnames = list("extraction", 
							c("unif.seed", "int.seed")))
		return(FAobject)
	} # else proceed to use genoud()

	## Prepare to call genoud() via the well-known model.frame() trick
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name("genoud")
        mc[names(formals(Factanal))] <- NULL

	# Arguments for genoud() that are *logically required* by Factanal()
	mc$nvars    <- restrictions_copy@nvars
	mc$max      <- FALSE
	mc$hessian  <- FALSE
	mc$lexical  <- ifelse(is(restrictions_copy, "restrictions.factanal"), FALSE, TRUE)
	mc$Domains  <- restrictions_copy@Domains
	mc$default.domains <- NULL
	mc$data.type.int   <- FALSE
	if(any(restrictions_copy@Domains[,1] == restrictions_copy@Domains[,2])) {
		mc$boundary.enforcement <- 2
	}
	if(restrictions@model == "SEFA") mc$fn <- function(par) { # lexical fit function
				fitS4(par, restrictions_copy, manifest, lower, TRUE)
			}
	else mc$fn <- function(par) { # lexical fit function
			fitS4(par, restrictions_copy, manifest, lower, FALSE)
	}

	if(is(restrictions_copy, "restrictions.factanal")) mc$BFGSfn <- NULL
	else mc$BFGSfn <- function(par, helper = NA) { # scalar fit function (continuous)
			bfgs_fitS4(par, restrictions_copy, manifest, helper, lower)
		}
	mc$BFGShelp <- function(initial, done = FALSE) { # helper for BFGSfn and gr
			bfgs_helpS4(initial, restrictions_copy, manifest, done, lower)
		}
	if(restrictions_copy@discrepancy == "MLE") {
		mc$gr <- function(par, helper) { # gradient for ML
			gr_fitS4(par, restrictions_copy, manifest, helper, lower)
		}
	}
	else if(FAiR_is.QD(restrictions_copy)) {
		mc$gr <- function(par, helper) { # gradient for QD
				FAiR_gradient_QD(par, restrictions_copy, manifest, 
						helper, lower)
			}
	}
	else          mc$gr <- NULL # numeric gradient

	if(!analytic) mc$gr <- NULL # overwrite whatever was previous
	
	# Workaround to get replicatability from genoud()
	if(any(!is.null(mc$unif.seed) | !is.null(mc$int.seed))) {
		warning("Use the seeds argument to Factanal() instead of the unif.seed",
			"and int.seed arguments to genoud(). Using 12345 as the seed.")
	}
	if(!is.null(seeds)) {
		mc$unif.seed <- seeds[1]
		mc$int.seed  <- if(length(seeds) == 1) seeds else seeds[2]
	}
	else { # use genoud() defaults
		mc$unif.seed <- 812821
		mc$int.seed  <- 53058
	}

	# "Default" arguments for genoud() that can be superceded if explicitly specified
	if(is.null(mc$boundary.enforcement)) mc$boundary.enforcement <- 1
	if(is.null(mc$pop.size))             mc$pop.size <- 1000
	if(is.null(mc$MemoryMatrix))         mc$MemoryMatrix <- FALSE
	if(is.null(mc$print.level))          mc$print.level <- 1
	if(is.null(mc$P9mix))                mc$P9mix <- 1
	if(is.null(mc$BFGSburnin))           mc$BFGSburnin <- -1
	if(is.null(mc$max.generations))      mc$max.generations <- 1000
	if(is.null(mc$project.path))         mc$project.path <- paste(tempfile(), 
								"Factanal.txt", sep = "")

	# Deal with starting values
	if(FAiR_is.FA(SV <- eval(mc$starting.values))) { # use old result to start
		par <- c(coef(SV), log(SV@scale))
		if(is(restrictions_copy, "restrictions.factanal")) {
			par <- uniquenesses(SV)
		}
		else if(is(restrictions_copy, "restrictions.independent")) {
			par <- log(SV@scale)[restrictions_copy@free]
		}
		else if(is(restrictions_copy, "restrictions.orthonormal")) {
			Phi <- diag(factors[1])
			par <- c(Phi[lower.tri(Phi)], par)
			par <- par[restrictions_copy@free]
		}
		else if(is(restrictions_copy, "restrictions.2ndorder")) {
			if(is(SV, "FA.2ndorder")) {
				Xi <- cormat(SV, level = 2)
				Delta <- loadings(SV, level = 2)
				par <- c(Xi[lower.tri(Xi)], Delta, par)
			}
			else {
				Xi <- diag(factors[2])
				par <- c(Xi[lower.tri(Xi)], rep(0, prod(factors)), par)
			}
			par <- par[restrictions_copy@free]
		}
		else if(is(restrictions_copy, "restrictions.general")) {
			if(is(SV, "FA.general")) {
				par <- c(loadings(SV, level = 2), par)
			}
			else    par <- c(rep(0, factors[1]), par)
			par  <- par[restrictions_copy@free]
		}
		else if(is(restrictions_copy, "restrictions.1storder")) {
			Phi <- cormat(SV)
			par <- c(Phi[lower.tri(Phi)], par)
			par <- par[restrictions_copy@free]
		}
		else if(length(SV@optimization$extraction$par) == mc$nvars) {
			par <- SV@optimization$extraction$par
		}
		else {
			stop("it does not seem possible to extract starting values from ",
				"the object passed to 'starting.values'\n",
				"Please supply alternate starting values or leave ",
				"'starting.values' unspecified")
		}
		mc$starting.values <- par
	}
	else if(is(restrictions_copy, "restrictions.independent")) { # independence model
		if(is.null(SV)) {
			mc$starting.values <- log(manifest@sds)[restrictions_copy@free]
		}
	}

	if(is.null(mc$starting.values)) { # usual case of no starting values
		cat("Generating good starting values, have patience ...\n")
		flush.console()
		pop.size <- eval(mc$pop.size)
		if(impatient) {
			efa <- factanal(covmat = S, factors = factors[1],
					rotation = "none")
			start <- 1 - efa$uniquenesses
		}
		else { # Use approximate PACE to get initial communality estimates
			start <- FAiR_PACE_by_RGENOUD(S, factors[1], seeds = seeds)
			start <- as.matrix(start)
		}

		# Call method for creating start values for all parameters
		if(!is.null(seeds)) set.seed(mc$unif.seed)
		mc$starting.values <- create_start(pop.size, start, restrictions_copy,
							manifest, reject)
	}
	else if(any(is.na(eval(mc$starting.values)))) mc$starting.values <- NULL
	else if(length(eval(mc$starting.values)) == ncol(S)) { # starting communalities
		pop.size <- eval(mc$pop.size)
		# Call method for creating start values for all parameters
		start <- as.matrix(eval(mc$starting.values))
		if(ncol(start) > nrow(start)) start <- t(start)
		if(!is.null(seeds)) set.seed(mc$unif.seed)
		mc$starting.values <- create_start(pop.size, start, restrictions_copy,
							manifest, reject)
	}
	else if(length(eval(mc$starting.values)) == mc$nvars) {
		if(any(eval(mc$starting.values) < restrictions_copy@Domains[,1])) {
			stop("some starting values are below their lower bounds")
		}
		if(any(eval(mc$starting.values) > restrictions_copy@Domains[,2])) {
			stop("some starting values are above their upper bounds")
		}
	}
	else if(!is.matrix(eval(mc$starting.values))) {
		stop("'starting.values' must either be an object of class FA, a numeric",
			" vector of length ", mc$nvars, " or a numeric matrix")
	}

	## Estimate model
	opt <- eval(mc) # does all the real estimation work

	if(NelderMead && !is(restrictions_copy, "restrictions.factanal") &&
			 !is(restrictions_copy, "restrictions.independent")) {
		help  <- bfgs_helpS4(opt$par, restrictions_copy, manifest, FALSE, lower)
		foo <- function(par) {
			bfgs_fitS4(par, restrictions_copy, manifest, help, lower)
		}
		par <- opt$par
		if(restrictions_copy@model == "SEFA") par[help$squashed] <- 0
		optNM <- optim(par = par, fn = foo, method = "Nelder-Mead")
		if(optNM$value < opt$value[help$marker]) {
			opt$par <- optNM$par
			opt$value <- fitS4(opt$par, restrictions_copy, 
						manifest, lower, TRUE)
			opt$counts <- optNM$counts
			opt$convergence <- optNM$convergence
			opt$message <- optNM$message
			if(optNM$convergence == 10) {
				warning("Nelder-Mead solution is probably on a boundary")
				print(paste("Convergence code:", opt$convergence, 
						"Convergence message:", opt$message))
				cat("\n")
			}
		}
		else if(optNM$value == opt$value[help$marker]) {
			cat("Nelder-Mead resulted in no improvement; ",
				"convergence presumably achieved\n")
		}
		else cat("Nelder-Mead optimizer went backwards; thus ignored\n")

		if(restrictions_copy@model == "SEFA") par[help$squashed] <- 0
	}

	# Call method to postprocess opt and bake object of correct class
	FAobject <- create_FAobject(restrictions_copy, manifest, opt, 
					kall, scores, lower, analytic)
	FAobject@seeds <- matrix(c(mc$unif.seed, mc$int.seed), ncol = 2,
				dimnames = list("extraction", c("unif.seed", "int.seed")))
	return(FAobject)
}

## Rotate() finds an oblique transformations of the factors following EFA extraction
Rotate <-
function(FAobject, criteria = list(), methodArgs = list(), 
	normalize = rep(1, nrow(loadings(FAobject))), seeds = 12345, 
	NelderMead = TRUE, ...) {
## Arguments
# FAobject: an object of class FA; see Factanal()
# criteria: a list of functions or a list of character strings with names of functions or
#           leave it unspecified to get help from the GUI
# methodArgs: a list of arguments to be used by various criteria
# normalize: a logical, or a function, or a vector for row-normalization
#    seeds: PRNG seeds for the unif.seed and int.seed arguments of genoud()
# NelderMead: use method = "Nelder-Mead" in a call to optim() to polish the solution?
#      ...: More arguments, which are passed to genoud()

	if(!FAiR_is.EFA(FAobject)) {
 		stop("Rotate() only works for exploratory factor analysis")
	}
	if(!FAiR_is.orthogonal(FAobject)) {
		stop("'FAobject' is already rotated")
	}

	## Establish rotation criteria
	criteria <- FAiR_make_criteria(FAobject, criteria, methodArgs)

	## Prepare to call genoud via the well-known model.frame() trick
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name("genoud")
        mc[names(formals(Rotate))] <- NULL

	# Arguments for genoud() that are *logically required* by Rotate()
	mc$nvars    <- ncol(loadings(FAobject))^2
	mc$max      <- FALSE
	mc$hessian  <- FALSE
	mc$lexical  <- TRUE
	mc$Domains  <- NULL
	mc$default.domains <- 1
	mc$data.type.int  <- FALSE

	A <- loadings(FAobject)
	A <- sweep(A, 1, NormalizingWeight(A, normalize), "/")

	mc$fn       <- function(par) FAiR_Rotate(par, A, criteria)
	mc$BFGShelp <- function(initial, done = FALSE) {
				fits <- FAiR_Rotate(initial, A, criteria)
				marker <- which(fits != -1)[1]
				return(list(fits = fits, marker = marker))
			}
	mc$BFGSfn <- function(par, helper) {
				fits <- FAiR_Rotate(par, A, criteria)
				marker <- which(fits != -1)[1]
				if(marker > helper$marker) return(-.Machine$double.xmax)
				if(marker < helper$marker) return( .Machine$double.xmax)
				                           return(fits[marker])
			}
 	mc$gr <- NULL  # use numeric gradients even if an analytic gradient is known

	# Workaround to get replicatability from genoud()
	if(any(!is.null(mc$unif.seed) | !is.null(mc$int.seed))) {
		warning("Use the seeds argument to Factanal() instead of the unif.seed",
			"and int.seed arguments to genoud(). Using 12345 as the seed.")
	}
	if(!is.null(seeds)) {
		mc$unif.seed <- seeds[1]
		mc$int.seed  <- if(length(seeds) == 1) seeds else seeds[2]
	}
	else { # use genoud() default seeds
		mc$unif.seed <- 812821
		mc$int.seed  <- 53058
	}

	# "Default" arguments for genoud() that can be superceded if explicitly specified
        if(is.null(mc$pop.size))        mc$pop.size <- formals(genoud)$pop.size
	if(is.null(mc$print.level))     mc$print.level <- 1
	if(is.null(mc$max.generations)) mc$max.generations <- 1000
        if(is.null(mc$boundary.enforcement)) mc$boundary.enforcement <- 1
	if(is.null(mc$MemoryMatrix))    mc$MemoryMatrix <- FALSE
	if(is.null(mc$P9mix))           mc$P9mix <- 1
	if(is.null(mc$BFGSburnin))      mc$BFGSburnin <- 5
	if(is.null(mc$project.path))    mc$project.path <- paste(tempfile(), 
								"Rotate.txt", sep = "")
	if(is.null(mc$starting.values)) { # start with identity matrix
		Tmat <- diag(ncol(A))
		mc$starting.values <- c(Tmat[-ncol(A),], Tmat[ncol(A),])
	}
	else if(is.matrix(eval(mc$starting.values))) {
		if(eval(nrow(mc$starting.values) == mc$nvars)) {
			mc$starting.values <- t(mc$starting.values)
		}
		else stop("starting.values has the wrong number of columns")
		
		bar <- function(x) FAiR_Tmat2par(matrix(x, ncol(A), ncol(A)))
		mc$starting.values <- apply(eval(mc$starting.values), 1, bar)
	}
	else if(length(eval(mc$starting.values)) == mc$nvars) {
		baz <- function(x) FAiR_Tmat2par(matrix(x, ncol(A), ncol(A)))
		mc$starting.values <- baz(eval(mc$starting.values))
	}
	else stop("'starting.values' not in the proper format")

	## Choose T
	opt <- eval(mc) # does all the real work
	if(NelderMead) {
		help  <- mc$BFGShelp(opt$par, TRUE)
		foo <- function(par, helper = help) {
			fits <- FAiR_Rotate(par, A, criteria)
			marker <- which(fits != -1)[1]
			if(marker > helper$marker) return(-.Machine$double.xmax)
			if(marker < helper$marker) return( .Machine$double.xmax)
							return(fits[marker])
		}
		par <- opt$par
		optNM <- optim(par = par, fn = foo, method = "Nelder-Mead")
		if(optNM$value < opt$value[help$marker]) {
			opt$par <- optNM$par
			opt$value <- FAiR_Rotate(opt$par, A, criteria)
			opt$counts <- optNM$counts
			opt$convergence <- optNM$convergence
			opt$message <- optNM$message
			if(optNM$convergence == 10) {
				warning("Nelder-Mead solution is probably on a boundary")
				print(paste("Convergence code:", opt$convergence, 
						"Convergence message:", opt$message))
				cat("\n")
			}
		}
		else if(optNM$value == opt$value[help$marker]) {
			cat("Nelder-Mead resulted in no improvement; ",
				"convergence presumably achieved\n")
		}
		else {
			cat("Nelder-Mead optimizer went backwards; thus ignored\n")
		}
	}
	FAobject <- FAiR_opt2FAobject(opt, FAobject, criteria) # bake FA object
	FAobject@seeds <- rbind(FAobject@seeds, "rotation" = c(mc$unif.seed, mc$int.seed))
	return(FAobject)
}

## read.cefa() tries to read files produced by the Comprehensive Exploratory Factor
## Analysis (CEFA) package but is probably pretty fragile
read.cefa <- read.CEFA <- 
function(file) {
	lines <- readLines(file)

	n.obs <- as.integer(strsplit(lines[[1]], split = " ")[[1]][1])

	datatype <- grep("([D|d]ata", lines)
	if(length(datatype) == 0) datatype <- lines[2]
	datatype <- as.integer(strsplit(lines[datatype[1]], split = "")[[1]][1])

	if(datatype == 3) {
		stop("importing a factor loading matrix is not supportable", 
			" please import the raw data or the covariance matrix")
	}

	varnames_marker <- grep("([V|v]ariable", lines)
	if(length(varnames_marker) == 0) varnames_marker <- 4
	else varnames_marker <- varnames_marker[1]

	varnames <- as.integer(strsplit(lines[varnames_marker], split = "")[[1]][1])
	if(varnames == 1) {
		varnames <- as.character(NULL)
		varnames_marker <- varnames_marker + 1
		while(varnames_marker) {
			if(length(grep("[a-z]", lines[varnames_marker], 
					ignore.case = TRUE)) > 0) {
					
					tempnames <- strsplit(lines[varnames_marker],
							split = " ")
					if(length(grep("^[\\s]*[0-9]",tempnames[[1]][1])))
						break
					varnames <- c(varnames, unlist(sapply(tempnames,
							FUN = function(x) x[x!=""])))
					varnames_marker <- varnames_marker + 1
				}
			else if(length(grep("[0-9]", lines[varnames_marker])) > 0) break
			else {
				varnames <- NULL
				warning("there was problem assigning variable names")
				break
			}
		}
	}
	else if(varnames == 0) varnames <- NULL
	else {
		varnames <- NULL
		warning("could not determine whether there are variable names")
	}

	lines <- lines[-c(1:6)]
	mark <- grep("[a-z]", lines, ignore.case = TRUE)
	if(length(mark) > 0) lines <- lines[-mark]
	
	mark <- grep("[\\?|\\+|\\-]", lines)
	if(length(mark) > 0) lines <- lines[-mark]

	lines <- lines[grep("[^0123456789]",  lines)]

	if(datatype == 1) { # covariance matrix
		lines <- strsplit(lines, split = " ")
		lines <- lapply(lines, FUN = function(x) as.numeric(x[x!=""]))
		out <- matrix(0, nrow = length(lines), 
				 ncol = length(lines[[length(lines)]]))
		for(i in 1:nrow(out)) out[i,1:length(lines[[i]])] <- lines[[i]]
		out <- out + t(out)
		diag(out) <- diag(out) / 2
		rownames(out) <- colnames(out) <- varnames
		covlist <- list(cov = out, n.obs = n.obs)
		return(covlist)
	}
	else if(datatype == 2) { # numeric data
		out <- matrix(unlist(sapply(strsplit(lines, split = " "), 
					FUN = function(x) {
						as.numeric(x[x!=""])
					})), nrow = length(lines), byrow = TRUE)
		out <- as.data.frame(out)
		colnames(out) <- varnames
		return(out)
	}
	else if(datatype == 4) { # ordinal data
		out <- matrix(unlist(sapply(strsplit(lines, split = " "), 
					FUN = function(x) {
						as.numeric(x[x!=""])
					})), nrow = length(lines), byrow = TRUE)	
		out <- as.data.frame(out)
		for(i in 1:ncol(out)) out[,i] <- factor(out[,i], ordered = TRUE)
		colnames(out) <- varnames
		return(out)
	}
	else {
		stop("datatype not recognized, should be 1,2, or 4")
	}
}

## write.cefa() tries to write a file to be read by the Comprehensive Exploratory Factor
## Analysis (CEFA) package but is probably pretty fragile
write.cefa <- write.CEFA <-
function(file, FAobject) {
	if(!FAiR_is.FA(FAobject)) {
		stop("'FAobject' must be produced by Factanal()")
	}
	if(!FAiR_is.orthogonal(FAobject)) {
		stop("CEFA presumes the loading matrix has an orthogonal basis")
	}
	Lambda <- loadings(FAobject)
	cat(FAobject@manifest@n.obs, nrow(Lambda), "\t", "(Ncases, p)\n", file = file)
	cat("3", "\t", "(DataType)\n", file = file, append = TRUE)
	cat("0", "0", "\t", "(Random Seed)\n", file = file, append = TRUE)
	cat("1", "\t", "(Varnames)\n", file = file, append = TRUE)
	cat(rownames(Lambda), "\n", file = file, append = TRUE)
	cat("0", "\t", "(Facnames)\n", file = file, append = TRUE)
	cat("0", "\t", "(Tarmat)\n\n", file = file, append = TRUE)
	write.table(Lambda, file = file, append = TRUE, row.names = FALSE, 
							col.names = FALSE)
	cat("\n\nRemarks:\n", 
		"This matrix was written by a little-tested function in FAiR.\n",
		"You may need to edit it for taste or conformance with CEFA.\n",
		"In particular, you may need to change varnames to 0 and delete ",
		"the variable names for as yet unknown reasons.\n",
		file = file, append = TRUE)
	cat("File written to", deparse(substitute(file)), "\n", file = "")
}

## compares a model to a Null model or to another estimated model (not implemented yet)
model_comparison <-
function(..., correction = c("swain", "bartlett", "none"), conf.level = .9, nsim = 1001){
	
	## need to add how to arg list
	models <- list(...)
	correction <- match.arg(correction)
	if(conf.level <= 0 | conf.level >= 1) {
		stop("conf.level must strictly be between 0 and 1")
	}

	if(length(models) == 0) {
		stop("must pass an object of class 'FA' through the ...")
	}
	else if(!all(sapply(models, FUN = FAiR_is.FA))) {
		stop("all objects passed through ... must inherit from class 'FA'")
	}

	correction <- match.arg(correction)
	if(length(models) == 1) {
		FAobject <- models[[1]]
		out <- list(restrictions = FAobject@restrictions)
		if(FAiR_is.QD(FAobject)) {
			T_ADF <- FAiR_test_YB98( FAobject)
			T_1 <- T_ADF$statistic
			T_Browne <- FAiR_test_Browne84(FAobject, how = "ADF")
			null_model <- FAiR_independence_model(FAobject)
			T_0 <- FAiR_test_YB98(null_model)$statistic
# 			draws <- simulate(FAobject, nsim = 1001, standardized = FALSE)
# 			T_np <- FAiR_nonparametric(draws, FAobject)
			out$exact_fit <- list(T_ADF = T_ADF, T_Browne = T_Browne)
		}
		else if(FAobject@restrictions@discrepancy == "YWLS") {
			T_1 <- deviance(FAobject) # warning?
			T_0 <- NA_real_
		}

		if(FAiR_is.ML(FAobject)) {
			T_ML <- FAiR_test_exact_fit(FAobject, correction)
			T_1  <- T_ML$statistic
			null_model <- FAiR_independence_model(FAobject)
			T_0 <- FAiR_test_exact_fit(null_model, correction)$statistic
# 			draws <- simulate(FAobject, nsim = nsim, standardized = FALSE)
# 			T_np <- FAiR_nonparametric(draws, FAobject)
			out$exact_fit <- list(T_ML = T_ML)
			out$infocriteria <- FAiR_infocriteria(FAobject, null_model)
		}

		RMSEA <- FAiR_RMSEA(   FAobject, T_1, conf.level)
		gamma <- FAiR_Steiger( FAobject, RMSEA)

		out$close_fit <- list(RMSEA = RMSEA, gamma = gamma)

		GFI      <- FAiR_GFI( FAobject)
		AGFI     <- FAiR_AGFI(FAobject, GFI)
		McDonald <- FAiR_McDonald(FAobject, T_1)
		SRMR     <- FAiR_SRMR(FAobject)

		TLI  <- FAiR_TLI( FAobject, T_1, T_0)
		CFI  <- FAiR_CFI( FAobject, T_1, T_0)
		NFI  <- FAiR_NFI( T_1, T_0)
		NNFI <- FAiR_NNFI(FAobject, T_1, T_0)

		out$fit_indices <- list(GFI = GFI, AGFI = AGFI, McDonald = McDonald,
					SRMR = SRMR, TLI = TLI, CFI = CFI, NFI = NFI,
					NNFI = NNFI)
	}
	else {
		man <- models[[1]]@manifest
		if(!all(sapply(models, FUN = function(x) identical(x@manifest, man)))) {
			stop("all 'manifest' slots for the objects passed through the",
				" ... must be identical to conduct model comparisons")
		}
		out <- lapply(models, FUN = model_comparison, correction = correction,
				conf.level = conf.level, nsim = nsim)
	}
	return(out)
}

## This function compares two nested models
paired_comparison <-
function(M_0, M_1) {
	if(!FAiR_is.FA(M_0)) {
		stop("M_0 must inherit from class 'FA'")
	}
	else if(!FAiR_is.FA(M_1)) {
		stop("M_1 must inherit from class 'FA'")
	}
	else if(!identical(M_0@manifest, M_1@manifest)) {
		stop("M_0 and M_1 must utilize the same manifest object")
	}
	warning("paired comparison test does not verify whether the M_1 is nested ",
		"within M_0")
	return(FAiR_paired(M_0, M_1, how = "normal"))
}

## pass model info to Mathomatic
restrictions2Mathomatic <- restrictions2mathomatic <-
function(object, file = "") {
	if(FAiR_is.FA(object)) object <- object@restrictions
	else if(!is(object, "restrictions")) {
		stop("'object' must be of class 'FA' or class 'restrictions'")
	}
	beta <-   coef(object)
	Phi  <- cormat(object)
	Z <- matrix("", nrow(beta), ncol(Phi))
	C <- matrix("", nrow(beta), nrow(beta))
	for(i in 1:nrow(beta)) for(j in 1:ncol(beta)) {
		string <- as.character(NULL)
		for(k in 1:ncol(Phi)) {
			ministring <- paste(paste("beta_R", i, "C", k, sep = ""), "*",
					if(k > j) paste("Phi_R", k, "C", j, sep = "") else
					if(k==j)  paste("1") else
						  paste("Phi_R", j, "C", k, sep = ""))
			string <- c(string, ministring)
		}
		Z[i,j] <- paste("(", paste(string, collapse = " + "), ")", sep = "")
	}

	for(i in 1:nrow(Z)) for(j in 1:i) {
		string <- as.character(NULL)
		for(k in 1:ncol(Phi)) {
			ministring <- paste(Z[i,k], "*", 
					paste("beta_R", j, "C", k, sep=""))
			string <- c(string, ministring)
		}
		C[i,j] <- paste(paste("\n\nSigma_R", i, "C", j, sep = ""), "= (", 
				paste(string, collapse = " + "), ")")
	}

	C <- C[lower.tri(C, FALSE)]
	C <- paste(C, "; equation #", 1:length(C))
	cat(C, file = file)
	if(file != "") cat("\nfile written to ", file)
	cat("\nIn principle, you can solve for the free parameters in terms of Sigma_*")
	cat("\nSee ?eliminate, ?solve, and ?simplify in Mathomatic\n")
	cat("\nYou have to add equations for the constraints on your own")
	cat("\nFor example, enter (in Mathomatic) beta_R1C2 = 0 to signify that ")
	cat("\nRow 1, Column 2 of the beta matrix is constrained to be zero")
	cat("\nGood luck!\n")
}

## convert an object produced by GPFoblq(loadings(FAobject)) back to class 'FA'
GPA2FA <-
function(GPAobject, FAobject) {
	if(class(GPAobject) != "GPArotation") {
		stop("'GPAobject' must be of class 'GPArotation'") 
	}
	else if(! "Th" %in% names(GPAobject)) {
		stop("'GPAobject' needs to have an element named 'Th'")
	}
	else if(!is.matrix(GPAobject$Th)) {
		stop(" the 'Th' element of 'GPAobject' must be a matrix")
	}
	else if(!isTRUE(all.equal(colSums(GPAobject$Th^2), rep(1, ncol(GPAobject$Th))))) {
		stop(" the 'Th' element of 'GPAobject' must have normal columns")
	}
	else if(! "orthogonal" %in% names(GPAobject)) {
		stop("'GPAobjet' needs to have an element named 'orthogonal'")
	}

	if(!is(FAobject, "FA.EFA")) stop("FAobject must be produced by Factanal()")
	else if(FAobject@rotated) {
		stop("'FAobject' is already rotated; pass the unrotated version")
	}

	FAobject <- FAiR_opt2FAobject(opt = GPAobject, FAobject, GPA = TRUE)
	return(FAobject)
}

## converts to RAM version of the model like that used in library(sem)
FA2RAM <-
function(FAobject) {
	if(!FAiR_is.FA(FAobject)) {
		stop("FAobject must be produced by Factanal()")
	}
	return(restrictions2RAM(FAobject@restrictions))
}
