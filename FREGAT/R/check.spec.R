# FREGAT (c) 2016 Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

check.spec.famBT <- function() {

	fweights <- check.weights(weights, k, beta.par)
	if (!is.null(fweights)) weights <- NULL
	list(fweights = fweights, weights = weights)

}

check.spec.famSKAT <- function() {

	fweights <- check.weights(weights, k, beta.par)
	if (!is.null(fweights)) weights <- NULL
	method <- check.method(method)
	kernel <- match.arg(kernel, c('linear.weighted', 'quadratic', 'IBS', 'IBS.weighted', '2wayIX'))
	if (length(rho) > 1 | (length(rho) == 1 & rho)) {
		rhos <- check.rho(rho, kernel)
		rho <- TRUE
		return.variance.explained <- FALSE
		sapply(c('fweights', 'weights', 'method', 'kernel', 'rhos', 'rho', 'return.variance.explained'),
			function(x) get(x), simplify = FALSE, USE.NAMES = TRUE)
	} else {
		rho <- FALSE
		sapply(c('fweights', 'weights', 'method', 'kernel', 'rho'),
			function(x) get(x), simplify = FALSE, USE.NAMES = TRUE)
	}
}

check.spec.famFLM <- function() {

	stat <- match.arg(stat, c('F', 'Chisq', 'LRT'))
	BSF <- match.arg(BSF, c('fourier', 'bspline'))
	if (is.null(positions)) stop("'positions' are missing, with no default")
	if (length(positions) != k) stop("Dimensions of 'positions' and 'genodata' do not match")
 
	g <- TRUE

#	if (missing(GVF)) g <- FALSE
	if (is.null(GVF)) g <- FALSE
	if (is.logical(GVF)){
		if (!GVF) g <- FALSE else GVF <- 'bspline'
	}
	if (g) GVF <- match.arg(GVF, c('fourier', 'bspline'))

	kb0 <- kg0 <- order0 <- FALSE

	if (BSF == 'bspline') {
		order0 <- check.basis(order, 'order')
		kb0 <- check.basis(kb, 'basis', 'bspline', order0)
	} else { kb0 <- check.basis(kb, 'basis', 'fourier') }

	if (g) {
		if (GVF == 'bspline') {
			if (!order0) order0 <- check.basis(order, 'order')
			kg0 <- check.basis(kg, 'basis', 'bspline', order0)
		} else { kg0 <- check.basis(kg, 'basis', 'fourier') }
	}

	# base condition is kg >= kb

	if (g) {
		if (kg0 < kb0) {
			kb0 <- kg0
			warning (paste('kb cannot exceed kg, kb set to', kg0))
		}
		if (kg0 == kb0) {
			g <- FALSE
			if (GVF == 'bspline' & BSF == 'fourier') BSF <- 'bspline'
			if (GVF == 'fourier' & BSF == 'bspline') BSF <- 'fourier'
			warning ("GVF has no effect under kg = kb, the equivalent model with '", BSF, "' for BSF will be used")
		}
	}

	# if g then kg > kb
	# if kg == kb then !g

	if (BSF == 'bspline') betabasis0 = create.bspline.basis(norder = order0, nbasis = kb0)
	else betabasis0 <- create.fourier.basis(c(0, 1), nbasis = kb0)

	if (g) {
		if (GVF == 'bspline') { genobasis0 <- create.bspline.basis(norder = order0, nbasis = kg0)
		} else { genobasis0 <- create.fourier.basis(c(0, 1), nbasis = kg0) }
		J0 <- inprod(genobasis0, betabasis0)
		model0 <- paste(GVF, kg0, '-', BSF, kb0, sep = '')
	}
	else model0 <- paste(0, '-', BSF, kb0, sep = '')
	# all fourier bases are odd

	omit.linear.dependent <- TRUE

	if (g) sapply(c('stat', 'BSF', 'g', 'GVF', 'kb0', 'kg0', 'order0', 'model0', 'genobasis0', 'betabasis0', 'J0', 'omit.linear.dependent'),
			function(x) get(x), simplify = FALSE, USE.NAMES = TRUE)
	else sapply(c('stat', 'BSF', 'g', 'kb0', 'order0', 'model0', 'betabasis0', 'omit.linear.dependent'),
			function(x) get(x), simplify = FALSE, USE.NAMES = TRUE)

}
