# FREGAT (c) 2016

check.ini <- function(formula, phenodata, kin) {

	if (missing(formula)) stop("'formula' is missing, with no default")
	if (missing(phenodata)) stop("'phenodata' is missing, with no default")
#	if (missing(kin)) stop("'kin' is missing, with no default")

	H2est <- ifelse (!is.null(kin), TRUE, FALSE)

	if (!class(phenodata) %in% c('matrix', 'data.frame')) stop("'phenodata' should be a data.frame or matrix")

	n <- dim(phenodata)[1]
	if (H2est) {
		if (dim(kin)[1] != dim(kin)[2] | dim(kin)[1] != n) stop("Dimensions of 'kin' and 'phenodata' do not match")
	}
#return(class(formula))
#browser()
	if (is(try(as.formula(formula), silent = TRUE), "try-error")) {
#return(as.character(eval(formula)))
		formula <- phenodata[, as.character(formula)]
#		formula <- eval(formula,as.data.frame(phenodata))
	} else {
		if (!is(try(as.formula(formula), silent = TRUE), "try-error")) { formula <- as.formula(formula) }
	}

	return(list(n = n, formula = formula, H2est = H2est))
}

check.geno <- function(genodata, regions, n, ...) {
	if (grepl('vcf', genodata)) {
		if (requireNamespace("seqminer", quietly = TRUE)) {
			if (!local.file.exists(genodata)) stop(paste("Input file '", genodata,"' does not exist", sep = ''))
			if (hasIndex(genodata)) {
				if (!'geneFile' %in% names(match.call())) stop("'geneFile' should be specified to use VCF file as input")
				geneFile <- eval(match.call()$geneFile)
				if (class(geneFile) == 'character' & length(geneFile) == 1) stopifnot(local.file.exists(geneFile))
				if (is.null(regions)) stop("'regions' should be specified to use VCF file as input")
			}
			if (getRversion() >= "3.2.0") {
				invisible(capture.output(vcfn <- dim(seqminer::readVCFToMatrixByRange(genodata, '1:0-0', '')[[1]])[2], type = 'output'))
			} else { vcfn <- dim(seqminer::readVCFToMatrixByRange(genodata, '1:0-0', '')[[1]])[2] }
			if (vcfn != n) stop("Dimensions of 'phenodata' and 'genodata' do not match")
			annoType <- ifelse('annoType' %in% names(match.call()), eval(match.call()$annoType), '')
			return(list(gtype = 3, geneFile = geneFile, annoType = annoType))
		} else { stop(paste("Please install 'seqminer' package to process VCF file", sep = '')) }
	} else {
		tmp <- read.plink.names(genodata, ...)
		for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])
		if (length(idnames) != n) stop("Dimensions of 'phenodata' and 'genodata' do not match")
		return(list(gtype = 2, k = length(snvnames), snvnames = snvnames, idnames = idnames, bed = bed))
	}
}

check.cores <- function(ncores, nreg) {
	if (ncores < 1) ncores <- 1
	if (ncores > 1) {
		ncores <- min(ncores,nreg)
		if (ncores > detectCores()) {
			ncores <- detectCores()
			warning("Number of CPUs available: ", ncores)
		}
	}
	if (ncores > 1) {
		if (!(requireNamespace("foreach", quietly = TRUE) & requireNamespace("doParallel", quietly = TRUE))) {
			warning(paste("Please install packages 'foreach' and 'doParallel' to enable parallel calculations"))
			ncores <- 1
		}
	}
	return(ncores)
}

check.method <- function(method) {
	if (method == 'Davies') method <- 'davies'
	if (method == 'Kuonen') method <- 'kuonen'
	method <- match.arg(method, c('davies', 'kuonen'))
	return(method)
}

check.regions <- function(k, regions, sliding.window) {#, snvnames) {
	if (!is.null(regions)) {
		rtype <- dim(regions)[2]
		if (rtype == 1) {
			if (length(regions) != k) {
				if (!k) rtype <- 3 else stop("Dimensions of regions and genodata do not match")
			}
			l <- unique(regions)
		} else {
			if (rtype > 2) {
				regions <- regions[,1:2]
				rtype <- 2
			}
			l <- unique(regions[, 2])
#			if (is.null(snvnames)) stop("colnames not found in 'genodata'")
		}
		nreg <- length(l)
		return(list(l = l, nreg = nreg, rtype = rtype))
	} else {
		if (k < sliding.window[1]) sliding.window[1] <- k
		if (sliding.window[1] < 1) stop("Window size should be >= 1")
		if (sliding.window[2] < 1) stop("Step should be >= 1")
		nreg <- floor((k - sliding.window[1]) / sliding.window[2] + 1)
		p1 <- sapply(1:nreg, function(i) (i - 1) * sliding.window[2] + 1)
		p2 <- p1 + sliding.window[1] - 1
		l <- paste(p1, p2, sep = "_")
		return(list(l = l, nreg = nreg, p1 = p1, p2 = p2))
	}
}

check.weights <- function(weights, k, beta.par) {
	if (is.null(weights)) {
		if (beta.par[1] <= 0 | beta.par[2] <= 0) stop("beta.par should be > 0")
#		fweights <- function(maf, beta.par) ifelse(maf > 0, dbeta(maf, beta.par[1], beta.par[2]), 0)
		fweights <- function(maf, a = as.numeric(beta.par[1]), b = as.numeric(beta.par[2])) ifelse(maf > 0, dbeta(maf, a, b), 0)
	} else {
		if (is.function(weights)) {
			fweights <- weights
		} else {
			fweights <- NULL
			if (length(weights) != k & k) stop("Dimensions of weights and genodata do not match")
		}
	}
	return(fweights)
}

check.covariates <- function(n, formula, phenodata) {
	if (class(formula) == 'formula'){
		if (is(try(covariates <- as.matrix(model.frame(formula, phenodata)), silent = TRUE), "try-error")) stop("Check whether all covariates are given in 'phenodata'")
		X <- as.matrix(model.frame(formula, phenodata, na.action = NULL, drop.unused.levels = TRUE))
		measured.ids <- sapply(1:n, function(x) !any(is.na(X[x,])))
#	return(list(X = model.matrix(formula, phenodata, na.action = NULL, drop.unused.levels = TRUE), measured.ids = measured.ids))
		X <- cbind(as.numeric(X[as.logical(measured.ids), 1]), model.matrix(formula, phenodata, na.action = NULL, drop.unused.levels = TRUE))
#	return(list(X = X, measured.ids = measured.ids))
	} else {
		X <- as.matrix(formula)#;return(X)
		measured.ids <- !is.na(X)
		X <- as.matrix(cbind(X, 1))
		colnames(X)[2] <- '(Intercept)'
		X <- X[as.logical(measured.ids), ]
#	return(list(X = X, measured.ids = measured.ids))
	}
	return(list(X = X, measured.ids = measured.ids))
}

check.nullmod <- function(nullmod, X) {
	run.null <- 1
	if (!missing(nullmod)) {
		if (is(try(nullmod, silent = TRUE), "try-error")){
			warning("'nullmod' object not found. Running null model...")
		} else if (is(try(nullmod$total.var, silent = TRUE), "try-error") | is(try(nullmod$h2, silent = TRUE), "try-error") | is(try(nullmod$alpha, silent = TRUE), "try-error")){
			warning("'nullmod' should contain 'alpha', 'total.var' and 'h2' parameters. Running null model...")
		} else if (is.null(nullmod$total.var) | is.null(nullmod$h2) | is.null(nullmod$alpha)){
			warning("'nullmod' should contain 'alpha', 'total.var' and 'h2' parameters. Running null model...")
		} else {
			if (is.null(dim(nullmod$alpha))) {
				warning("nullmod$alpha table not found. Running null model...")
			} else if (dim(X)[2] != (dim(nullmod$alpha)[1] + 1)) {
				warning("Dimensions of nullmod$alpha and covariates do not match. Running null model...")
			} else { run.null <- 0 }
		}
	}
	return(run.null)
}

check.basis <- function(value, name, base = 'none', order) {
	if (base == 'bspline') {
		if (value < order) {
			value <- order
			warning(paste("bspline basis cannot be less than order, set to", value))
		}
	}
	if (value < 1) {
		value <- 1
		warning(paste(name, "cannot be less than 1, set to", value))
	}
	if (base == 'fourier' & (value + 1) %% 2 != 0) {
		if (ceiling(value) %% 2 != 0) { value <- ceiling(value)
		} else if (floor(value) %% 2 != 0) { value <- floor(value)
		} else value <- value - 1
		warning(paste("fourier basis should be an odd integer, set to", value))
	}
	if (value %% 1 != 0) {
		value <- round(value)
		warning(paste(name, "should be an integer number, set to", value))
	}
	return(value)
}

check.rho <- function(rho, kernel) {
	if (kernel != "linear.weighted") {
		stop("rho can only be used with linear.weighted kernel")
	}
	if (length(rho) == 1 & rho) {
		rho <- (0:10)/10
		#rho=c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2,0.5^2, 0.5, 1)
	} else {
		for (i in 1:length(rho)) {
			if (rho[i] < 0 || rho[i] > 1) {
				stop("rho should be >= 0 and <= 1")
			}
		}
	}
	return(rho)
}
