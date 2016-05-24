# FREGAT (c) 2016 Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

'null' <- function (formula, phenodata, kin = NULL, opt.method = 'optimize', ih2 = 0.3, eps = 1.e-04) {

############## CHECKS

ch <- check.ini(substitute(formula), phenodata, kin)
for(i in 1:length(ch)) assign(names(ch)[i], ch[[i]])

X <- check.covariates(n, formula, phenodata)

measured.ids <- as.logical(X$measured.ids)
X <- as.matrix(X$X)

############# NULL MODEL

suppressWarnings(NullMixedModel(X[, 1], X[, -1], kin * 2, opt.method = opt.method, ih2 = ih2, eps = 1.e-04, H2est = H2est))

}