# FREGAT (c) 2016 Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

'do.regional' <- function (formula, phenodata, genodata, kin = NULL, nullmod, #return.nullmod = TRUE,
regions = NULL, sliding.window = c(20, 10), mode = 'add', ncores = 1, return.time = FALSE, #return.sample.size = TRUE,
kernel = 'linear.weighted', beta.par = c(1, 25), weights = NULL, method = 'kuonen', acc = 1e-8, lim = 1e+6,
return.variance.explained = FALSE, reml = TRUE, positions = NULL, GVF = FALSE, BSF = 'fourier', kg = 30,
kb = 25, order = 4, stat = 'F', flip.genotypes = FALSE, impute.method = 'mean', rho = FALSE, test, ...) {

t0 <- proc.time()

############ COMMON CHECKS

if (missing(genodata)) stop("'genodata' is missing, with no default")

tmp <- check.ini(substitute(formula), phenodata, kin)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])

mode <- match.arg(mode, c('add', 'dom', 'rec'))
impute.method <- match.arg(impute.method, c('mean', 'blue'))

k <- FALSE
gtype <- 0

if (class(genodata) %in% c('gwaa.data', 'snp.data')) {
	if (requireNamespace("GenABEL", quietly = TRUE)) {
		if (class(genodata) == 'gwaa.data') genodata <- GenABEL::gtdata(genodata)
		gtype <- 1
	} else { stop(paste("'genodata' class '", class(genodata),"' cannot be processed without 'GenABEL' package installed", sep = '')) }
} else if (class(genodata) == 'character' & length(genodata) == 1) {
	tmp <- check.geno(genodata, regions, n, ...)
	for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])
} else if (class(genodata) == 'snpMatrix') {
	genodata <- as(genodata, 'numeric')
} else if (!class(genodata) %in% c('matrix', 'data.frame')) {
	if (is.numeric(genodata)) {
		genodata <- as.matrix(genodata)
	} else stop("Wrong 'genodata' class")
} else if (!is.numeric(genodata[1, 1])) {
	stop("'genodata' elements should be numeric")
}

if (gtype < 2) {
	if (dim(genodata)[1] != n) stop("Dimensions of 'phenodata' and 'genodata' do not match")
	k <- dim(genodata)[2]
}

X <- check.covariates(n, formula, phenodata)

measured.ids <- as.logical(X$measured.ids)
X <- as.matrix(X$X)
n1 <- sum(measured.ids)

if (n1 != n) {
	if (gtype < 2) genodata <- genodata[measured.ids, ]
	if (H2est) kin <- kin[measured.ids, measured.ids]
}

if (!is.null(regions)) regions <- as.matrix(regions)
tmp <- check.regions(k, regions, sliding.window)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])
if (!is.null(regions)) {
	if (rtype == 2 & gtype < 2) {
		snvnames <- colnames(genodata)
		if (is.null(snvnames)) stop("colnames not found in 'genodata'")
	}
}

ncores <- check.cores(ncores, nreg)

############ SPECIFIC CHECKS

fweights <- NULL # for genotypes() lazy estimate
omit.linear.dependent <- FALSE

if (test != 'MLR') { # no checks yet needed for MLR
	fun <- as.function(get(paste('check.spec', test, sep='.')))
	environment(fun) <- environment()
	tmp <- fun()
	for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])
}

############ NULL MODEL

t00 <- proc.time()

run.null <- check.nullmod(nullmod, X)

if (return.variance.explained & !reml) X0 <- X

if (run.null) nullmod <- suppressWarnings(NullMixedModel(X[, 1], X[, -1], kin * 2, H2est = H2est, ...))
#return(nullmod)

if (H2est) {
	SIGMA <- nullmod$total.var * (kin * 2 * nullmod$h2 + diag(n1) * (1 - nullmod$h2))
	SIGMAi <- chol2inv(chol(SIGMA))
	if (impute.method == 'blue') colInvOmega <- as.vector(colMeans(SIGMAi))
} else {
	SIGMA  <- nullmod$total.var * (diag(n1))
	SIGMAi <- (1. / nullmod$total.var) * (diag(n1))
	if (impute.method == 'blue') impute.method <- 'mean'
}

if (test == 'famFLM') {
	pheno <- X[, 1]
} else {
	pheno <- X[, 1] - X[, -1] %*% as.matrix(nullmod$alpha[match(colnames(X)[-1], rownames(nullmod$alpha)), 1])
}
X <- X[, -1]

ttt0 <- (proc.time() - t00)

############ PRELIMINARIES

if (test == 'famSKAT') {
	if (H2est) {
		SIG_res <- crossprod(SIGMAi, pheno)  # SIG_res = SIGMAi %*% residuals
		SiX <- crossprod(SIGMAi, X)  # SiX = SIGMAi %*% X
	} else {
		SIG_res <- pheno / nullmod$total.var
		SiX <- X / nullmod$total.var
	}
	if (kernel == 'linear.weighted') {
		SPS <- SIGMAi - SiX %*% solve(t(X) %*% SiX, t(SiX))
	} else {
		P11 <- diag(n1) - X %*% solve(t(X) %*% SiX, t(SiX))
		if (H2est) CholSigmaiP11 <- chol(SIGMAi) %*% P11 else CholSigmaiP11 <- sqrt(1./nullmod$total.var) *  P11
	}
}

if (test != 'famSKAT' | (test == 'famSKAT' & (rho | return.variance.explained))) {
	if (H2est) {
		if (test == 'famSKAT' & return.variance.explained) {
			beta.par <- c(0.5, 0.5)
			CholInvCo <- chol(SIGMAi) # Cov
		} else {
			Cori <- nullmod$total.var * SIGMAi
			CholInvCo <- chol(Cori) # Cor
		}
		covariate <- CholInvCo %*% X
		if (test != 'famFLM') {
			P11 <- diag(n1) - covariate %*% solve(t(covariate) %*% covariate, t(covariate))
			P11CholInvCo <- P11 %*% CholInvCo
		}
		if (!(test == 'famSKAT' & rho)) {
			pheno <- CholInvCo %*% pheno
		}
	} else {
		cf <- ifelse(test == 'famSKAT' & return.variance.explained, sqrt(1. / nullmod$total.var), 1) # Cov, Cor
		covariate <- cf * X
		if (test != 'famFLM') {
			P11 <- diag(n1) - covariate %*% solve(t(covariate) %*% covariate, t(covariate))
			P11CholInvCo <- cf * P11
		} else {
			CholInvCo <- diag(n1)
		}
	}
}

############ ANALYSIS

t00 <- proc.time()

lgt <- 4 # n.columns in output
if (test == 'famSKAT' & return.variance.explained) {
	if (reml) lgt <- 5 else lgt <- 8
} else if (test == 'famFLM') {
	lgt <- 5
} else if (test == 'famBT') { lgt <- 6 }

environment(analyze.region) <- environment()
environment(genotypes) <- environment()
pval.region <- as.function(get(paste('pval', test, sep='.')))
environment(pval.region) <- environment()
if (test == 'famSKAT' & rho) environment(pval.famSKATO) <- environment()
if (test == 'famFLM') environment(pval.MLR) <- environment()

#### SEQUENTIAL MODE

if (ncores == 1) {

	out <- as.data.frame(matrix(NA, nrow = nreg, ncol = lgt))

	for (i in 1:nreg) {
		r <- as.character(l[i])
		if (!is.null(regions)) {
			if (rtype == 1) { reg <- regions == r
			} else if (rtype == 2) {
				reg <- snvnames %in% regions[regions[, 2] == r, 1]
#				if (sum(reg) == 0) {
#					out[i, ] <- c(r, 0, 0, rep(NA, lgt - 3))
#					next
#				}
			} else if (rtype == 3) reg <- r
		} else { reg <- p1[i]:p2[i] }
		out[i, ] <- analyze.region()

	}

	#out[, 1] <- l

#### PARALLEL MODE

} else {

	nj <- floor(nreg / ncores)
	vj <- rep(nj, ncores)
	nadd <- nreg %% ncores
	if (nadd > 0) vj[1:nadd] <- nj + 1

	cl <- makeCluster(ncores)
	doParallel::registerDoParallel(cl)

	clusterExport(cl, varlist = ls(), envir = environment())
	out <- foreach::'%dopar%'(foreach::foreach(j = 1:ncores, .combine = rbind, .inorder = F), {

		i0 <- (sum(vj[1:(j - 1)]) * (j > 1) + 1)
		out <- c()

		for (i in i0:(i0 + vj[j] - 1)) {
			r <- as.character(l[i])
			if (!is.null(regions)) {
				if (rtype == 1) { reg <- regions == r 
				} else if (rtype == 2) {
					reg <- snvnames %in% regions[regions[, 2] == r, 1]
#					reg <- as.character(regions[regions[, 2] == r, 1])
#					if (sum(reg) == 0) {
#						out[i, ] <- c(r, 0, 0, rep(NA, lgt - 3))
#						next
#					}
				} else if (rtype == 3) reg <- r
			} else { reg <- p1[i]:p2[i] }
			out <- rbind(out, analyze.region())
#out <- rbind(out, environment())
		}
		out
#	}, env)
	})

	stopCluster(cl)
#browser()
	out <- as.data.frame(out)
	out <- out[match(l, out[, 1]), ]

}

#### LAST REGION IN SLIDING WINDOW

if (is.null(regions)) {
	if (p2[nreg] < k & (nreg * sliding.window[2] + 1) < k) {
		i <- nreg + 1
		p1 <- (i - 1) * sliding.window[2] + 1
		p2 <- k
		r <- paste(p1, p2, sep = "_")
		reg <- p1:p2
		out <- rbind(out, analyze.region())
	}
}

#### FINAL MAKE UP

out[, 2:4] <- sapply(2:4, function(x) as.numeric(as.character(out[, x])))
if (test != 'famFLM' & lgt > 4) out[, 5:lgt] <- sapply(5:lgt, function(x) as.numeric(as.character(out[, x])))

colnames(out) [1:4] <- c('region', 'markers', 'cleaned.markers', 'pvalue')
if (test == 'famSKAT' & return.variance.explained) {
	colnames(out) [5] <- 'prop.variance'
	if (!reml) {
		colnames(out) [6:8] <- c('h2', 'total.var', 'logLH')
	}
} else if (test == 'famFLM') {
	colnames(out) [5] <- 'model'
	out[, 5] <- gsub("fourier", "F", out[, 5])
	out[, 5] <- gsub("bspline", "B", out[, 5])
	out[, 2:4] <- sapply(2:4, function(x) as.numeric(as.character(out[, x])))
} else if (test == 'famBT') colnames(out) [5:6] <- c('beta', 'se.beta')

ttt <- (proc.time() - t00)

############ FINISH

out0 <- out
out <- c()
out$results <- out0
out$sample.size <- n1
out$nullmod <- nullmod
if (return.time) {
out$time$null <- ttt0
out$time$regions <- ttt
out$time$total <- (proc.time() - t0)
}

out

}