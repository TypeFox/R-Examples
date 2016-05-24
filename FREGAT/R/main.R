# FREGAT (c) 2016 Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

'FFBSKAT' <- function (formula, phenodata, genodata, kin = NULL, nullmod, regions = NULL,
sliding.window = c(20, 10), mode = 'add', ncores = 1, return.time = FALSE,
kernel = 'linear.weighted', beta.par = c(1, 25), weights = NULL, method = 'kuonen',
acc = 1e-8, lim = 1e+6, return.variance.explained = FALSE, reml = TRUE,
flip.genotypes = FALSE, impute.method = 'mean', rho = FALSE, ...) {

	do.call(do.regional, c(make.call(match.call()), test = 'famSKAT'), envir = parent.frame())

}

'famBT' <- function (formula, phenodata, genodata, kin = NULL, nullmod, regions = NULL,
sliding.window = c(20, 10), mode = 'add', ncores = 1, return.time = FALSE,
beta.par = c(1, 25), weights = NULL, flip.genotypes = FALSE, impute.method = 'mean', ...) {

	do.call(do.regional, c(make.call(match.call()), test = 'famBT'), envir = parent.frame())

}

'MLR' <- function (formula, phenodata, genodata, kin = NULL, nullmod, regions = NULL,
sliding.window = c(20, 10), mode = 'add', ncores = 1, return.time = FALSE,
stat = 'F', impute.method = 'mean', ...) {

	do.call(do.regional, c(make.call(match.call()), test = 'MLR'), envir = parent.frame())

}

'famFLM' <- function (formula, phenodata, genodata, kin = NULL, nullmod, regions = NULL,
sliding.window = c(20, 10), mode = 'add', ncores = 1, return.time = FALSE,
positions = NULL, GVF = FALSE, BSF = 'fourier', kg = 30, kb = 25, order = 4,
stat = 'F', flip.genotypes = FALSE, impute.method = 'mean', ...) {

	do.call(do.regional, c(make.call(match.call()), test = 'famFLM'), envir = parent.frame())

}

make.call <- function (cl) {
	cl <- cl0 <- as.list(cl)[-1L]
	names(cl) <- match.arg(names(cl0), names(formals(do.regional)), several.ok = TRUE)
	names(cl)[is.na(names(cl))] <- names(cl0)[is.na(names(cl))]
	cl
}

if (getRversion() >= "2.15.1") utils::globalVariables(c('r', 'lgt', 'test', 'rho', 'pval.region', 'genodata',
'n1', 'pheno', 'P11CholInvCo', 'nullmod', 'idnames', 'snvnames', 'bed', 'k', 'beta.par', 'positions', 'kb', 'kg',
'n', 'geneFile', 'annoType', 'rtype', 'nreg', 'l', 'j', 'gtype', 'reg', 'measured.ids', 'flip.genotypes', 'fweights',
'H2est', 'model0', 'g', 'kg0', 'genobasis0', 'J0', 'GVF', 'order0', 'kb0', 'betabasis0', 'BSF', 'stat', 'SIG_res',
'SPS', 'CholSigmaiP11', 'method', 'acc', 'lim', 'return.variance.explained', 'reml', 'X0', 'kin', 'y', 'SIGMA',
'rhos', 'sliding.window', 'P11CholInvCor', 'omit.linear.dependent', 'CholInvCo', 'impute.method', 'colInvOmega'))
