`tTable` <-
function (model, ...) 	{
	.Deprecated("coefTable")
	coefTable(model, ...)
}

`coefTable` <-
function (model, ...) UseMethod("coefTable")

.makeCoefTable <- 
function(x, se, df = NA_real_, coefNames = names(x)) {
	if(n <- length(x)) {
		xdefined <- !is.na(x)
		ndef <- sum(xdefined)
		if(ndef < n) {
			if(length(se) == ndef) {
				y <- rep(NA_real_, n); y[xdefined] <- se; se <- y
			}
			if(length(df) == ndef) {
				y <- rep(NA_real_, n); y[xdefined] <- df; df <- y
			}
		}
	}
	if(n && n != length(se)) stop("length(x) is not equal to length(se)")
	ret <- matrix(NA_real_, ncol = 3L, nrow = length(x),
		dimnames = list(coefNames, c("Estimate", "Std. Error", "df")))
	if(n) ret[, ] <- cbind(x, se, rep(if(is.null(df)) NA_real_ else df,
		length.out = n), deparse.level = 0L)
	class(ret) <- c("coefTable", "matrix")
	ret
}

`coefTable.default` <-
function(model, ...) {
	dfs <- tryCatch(df.residual(model), error = function(e) NA_real_)
	cf <- summary(model, ...)$coefficients
	.makeCoefTable(cf[, 1L], cf[, 2L], dfs, coefNames = rownames(cf))
}

`coefTable.lm` <-
function(model, ...)
	.makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))), model$df.residual)


`coefTable.survreg` <- 
function(model, ...) {
.makeCoefTable(
	coeffs(model), 
	sqrt(diag(vcov(model, ...))), 
	NA,
	)
}

`coefTable.coxph` <-
function(model, ...) {
	.makeCoefTable(coef(model), if(all(is.na(model$var))) 
		rep(NA_real_, length(coef(model))) else sqrt(diag(model$var)), 
		model$df.residual)
}
	
`coefTable.glmmML` <- function(model, ...)
	.makeCoefTable(model$coefficients, model$coef.sd)

`coefTable.gls` <-
function (model, ...)
	.makeCoefTable(coef(model), sqrt(diag(as.matrix(model$varBeta))),
		model$dims$N - model$dims$p)

`coefTable.lme` <-
function(model, adjustSigma = TRUE, ...) {
	se <- sqrt(diag(as.matrix(model$varFix)))
	if (adjustSigma && model$method == "ML")
		se <- se * sqrt(model$dims$N / (model$dims$N - length(se)))
	.makeCoefTable(nlme::fixef(model), se, model$fixDF[["X"]])
}

`coefTable.multinom` <- 
function (model, ...) {
	.makeCoefTable(coeffs(model), sqrt(diag(vcov(model, ...))))
}

`coefTable.sarlm` <-
`coefTable.spautolm` <-
function(model, ...) {
	x <- coef(model)
	.makeCoefTable(x, sqrt(diag(summary(model, ...)$resvar))[names(x)])
}

`coefTable.coxme` <-
`coefTable.lmekin` <-
function(model, ...)  {
	# code from coxme:::print.coxme
	beta <- model$coefficients # for class coxme:
	if(is.list(beta) && !is.null(beta$fixed))
		beta <- beta$fixed # for class lmekin and older coxme
	nvar <- length(beta)
	if(nvar) {
		nfrail <- nrow(model$var) - nvar
		se <- sqrt(get("diag", getNamespace("Matrix"))(model$var)[nfrail + 1L:nvar])
	} else se <- NULL
	.makeCoefTable(beta, se)
}

`coefTable.rq` <- 
function(model, ...)
	.makeCoefTable(model$coefficients, rep(NA_real_, length(model$coefficients)))
	
`coefTable.zeroinfl` <-
function(model, ...)
	.makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))))

`coefTable.hurdle` <- 
function(model, ...) {
	cts <- summary(model)$coefficients
	ct <- do.call("rbind", unname(cts))
	cfnames <- paste0(rep(names(cts), vapply(cts, nrow, 1L)), "_", rownames(ct))
	.makeCoefTable(ct[, 1L], ct[, 2L], coefNames = cfnames)
	#.makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))))
}

`coefTable.aodql` <-
`coefTable.betareg` <- 
`coefTable.glimML` <-
`coefTable.unmarkedFit` <- 
function(model, ...)
	.makeCoefTable(coef(model), sqrt(diag(vcov(model, ...))))

`coefTable.gee` <-
`coefTable.geeglm` <-
function(model, ..., type = c("naive", "robust")) {
	cf <- summary(model, ...)$coefficients
	j <- if(match.arg(type) == "naive") 2L else 4L
	.makeCoefTable(cf[, 1L], cf[, j], coefNames = rownames(cf))
}

`coefTable.geem` <-
function(model, ..., type = c("naive", "robust")) {
	smr <- summary(model)
	.makeCoefTable(smr$beta, smr[[if(match.arg(type) == "naive")
								  "se.model" else "se.robust"]],
	               coefNames = smr$coefnames)
}

`coefTable.geese` <-
function(model, ..., type = c("naive", "robust")) {
	cf <- summary(model, ...)$mean
	type <- match.arg(type)
	j <- if(type == "naive") 2L else 4L
	.makeCoefTable(cf[, 1L], cf[, j], coefNames = rownames(cf))
}

`coefTable.yagsResult` <-
function(model, ..., type = c("naive", "robust")) {
	type <- match.arg(type)
	vcv <- slot(model, if(type == "naive") "naive.parmvar" else "robust.parmvar")
	.makeCoefTable(model@coefficients, sqrt(diag(vcv)), coefNames = model@varnames)
}

`coefTable.splm` <- 
function (model, ...) {
	cf <- sapply(c("coefficients", "arcoef", "errcomp"), function(i)
		if(is.matrix(model[[i]])) model[[i]][, 1L] else model[[i]],
		simplify = FALSE)
	
	ncf <- sapply(cf, length)
	vcovlab <- c(coefficients = "vcov", arcoef = "vcov.arcoef", errcomp = "vcov.errcomp")
	se <- sqrt(unlist(lapply(names(vcovlab), function(i) {
		vcv2 <- diag(model[[vcovlab[i]]])
		c(vcv2, rep(NA_real_, ncf[[i]] - length(vcv2)))
	})))
	
	.makeCoefTable(unlist(cf, use.names = FALSE), se,
		coefNames = unlist(lapply(cf, names), use.names = FALSE))
}

`coefTable.MCMCglmm` <-
function (model, ...) {
	cf <- coeffs(model)
	.makeCoefTable(cf, se = rep(NA_real_, length.out = length(cf)))
}

`coefTable.gamm` <-
function (model, ...) coefTable.lm(model$gam, ...)

`coefTable.mark` <- 
function (model, orig.names = FALSE, ...) {
    dfs <- model$results[['n']] - model$results[['npar']]
    beta <- model$results[['beta']]
    .makeCoefTable(beta[, 1L], beta[, 2L], dfs,
		coefNames = if(orig.names) rownames(beta) else
			gsub("^([a-zA-Z]+):(.*)$", "\\1(\\2)", rownames(beta), perl = TRUE))
}

`coefTable.logistf` <-
function (model, ...)
.makeCoefTable(model$coefficients, sqrt(diag(model$var)))

`coefTable.aodml` <-
function (model, ...) {
	.makeCoefTable(coeffs(model), sqrt(diag(vcov(model))))
	#.makeCoefTable(coeffs(model), sqrt(diag(model$varparam)))
}

## XXX: fixed effects coefficients only
`coefTable.asreml` <- 
function (model, ...)  {
	.makeCoefTable(
		x = model$coefficients$fixed, 
		se = sqrt(model$vcoeff$fixed * model$sigma2) ## ?
		) 
}

`coefTable.cplm` <-
function (model, ...) 
.makeCoefTable(coef(model), sqrt(diag(vcov(model))),
			   model@df.residual)

`coefTable.cpglmm` <-
function (model, ...) 
.makeCoefTable(coeffs(model), sqrt(diag(vcov(model))))

`coefTable.maxlikeFit` <-
function (model, ...)
.makeCoefTable(model$Est[, 1L], model$Est[, 2L])


coefTable.bic.glm <-
function (model, ...) {
	.makeCoefTable(model$condpostmean, model$condpostsd, NA_integer_,
		dimnames(model$mle)[[2L]])
}

# coefTable methods:

`print.coefTable` <-
function (x, ...)
stats::printCoefmat(x[, if(all(is.na(x[, 3L]))) -3L else TRUE, drop = FALSE],
	has.Pvalue = FALSE)

summary.coefTable <-
function (object, ...) {
	tvalue <- object[, 1L] / object[, 2L]
	if (all(is.na(object[, 3L]))) {
		pvalue <- 2 * pnorm(-abs(tvalue))
		rval <- cbind(object, tvalue, pvalue)
		cn <- c("z value", "Pr(>|z|)")
	} else if (any(is.finite(tvalue))) {
		pvalue <- 2 * pt(-abs(tvalue), object[, 3L])
		cn <- c("t value", "Pr(>|t|)")
    } else {
        pvalue <- tvalue <- NaN
		cn <- c("t value", "Pr(>|t|)")
    }
	rval <- cbind(object, tvalue, pvalue)
	colnames(rval)[4L:5L] <- cn
	class(rval) <- c("summary.coefTable", class(object))
	rval
}

`print.summary.coefTable` <-
function (x, signif.stars = getOption("show.signif.stars"), ...) {
	j <- if(all(is.na(x[, 3L]))) -3L else TRUE
	stats::printCoefmat(x[, j],
		has.Pvalue = any(is.finite(x[, 5L])),
		signif.stars = signif.stars,
		...)
}

plot.mcoefTable <-
function (x, y, labAsExpr = FALSE, n = 101, w = 5, ...) {
	lab_as_expr <- function(x) {
		x <- gsub(":", "%*%", x, perl = TRUE)
		x <- gsub("\\B_?(\\d+)(?![\\w\\._])", "[\\1]", x, perl = TRUE)
		parse(text = x)
	}
	
	xd <- function(z, n, w, x = NULL) {
		rval <- matrix(NA_real_, ncol = 3L, nrow = n)
		rval[, 1L] <- if(is.null(x)) seq(z[1L] - (w * z[2L]), z[1L] + (w * z[2L]),
			length.out = n) else x
		if(!is.na(z[3L]))
			rval[, 2L] <- dt((rval[, 1L] - z[1L]) / z[2L], z[3L]) / z[2L]
		rval[, 3L] <- dnorm(rval[, 1L], z[1L], z[2L])
		rval
	}
	
	m <- nrow(x)
	lab <- if(labAsExpr) lab_as_expr(rownames(x)) else rownames(x)
	nmodels <- dim(x)[3L]
	
	col <- 1L:nmodels

	par(mfrow = n2mfrow(m), ...)
	for(i in 1L:m) {
		#cat("--", dimnames(x)[[1]][i], "--\n")

		
		mat <- matrix(0, n, 2L * nmodels + 1L)
		for(k in 1L:nmodels) {
			j <- seq.int(length.out = 2 + (k == 1), from = 1 + (k - 1)* 2 + (k != 1))
			xlim <- c(min(x[i, 1L, ]) - w * min(x[i, 2L, ]), max(x[i, 1L,]) + w * max(x[i,2L,]))
			vx <- seq(xlim[1L], xlim[2L], length.out = n)
			v <- xd(x[i, , k], n = n, w = w, x = vx)
			mat[, j] <- v[, (2 - (k == 1)):3]
		}
		
		plot.new()
		plot.window(xlim = range(mat[, 1L]), ylim = range(mat[,-1L], na.rm = TRUE))
		j <- seq.int(2, length.out = nmodels, by = 2)
		if(any(!is.na(mat[, j])))
			matplot(mat[, 1], mat[, j], type = "l", lty = 2, add = TRUE, col = col)
		matplot(mat[, 1L], mat[, j + 1L], type = "l", lty = 1, add = TRUE, col = col)
		abline(v = x[i, 1L, ], lty = 1L, col = col)
		abline(v = 0, lty =3L, col = 8)
		axis(1L)
		axis(2L)
		box()
		
		title(lab[[i]])
		
	}
	invisible()
}

plot.coefTable <-
function (x, y, labAsExpr = FALSE, n = 101, w = 5,...) {
	lab_as_expr <- function(x) {
		x <- gsub(":", "%*%", x, perl = TRUE)
		x <- gsub("\\B_?(\\d+)(?![\\w\\._])", "[\\1]", x, perl = TRUE)
		parse(text = x)
	}
	xd <- function(z, n, w) {
		rval <- matrix(NA_real_, ncol = 3L, nrow = n)
		rval[, 1L] <- seq(z[1L] - w * z[2L], z[1L] + w * z[2L],
			length.out = n)
		if(!is.na(z[3L]))
			rval[, 2L] <- dt((rval[, 1L] - z[1L]) / z[2L], z[3L]) / z[2L]
		rval[, 3L] <- dnorm(rval[, 1L], z[1L], z[2L])
		rval
	}
	
	m <- nrow(x)
	lab <- if(labAsExpr) lab_as_expr(rownames(x)) else rownames(x)
	i <- 1
	nmodels <- dim(x)[3L]

	par(mfrow = n2mfrow(m))
	for(i in 1L:m) {	
		v <- xd(x[i, ], n = n, w = w)
		plot.new()
		plot.window(xlim = range(v[, 1L]), ylim = range(v[,-1L], na.rm = TRUE))
		lines(v[, 1L], v[, 2L])
		lines(v[, 1L], v[, 3L], col = "red")
		abline(v = c(x[i, 1L], 0), lty = c(1L, 3L))
		axis(1L)
		axis(2L)
		box()
		title(subst(expression(A == B %+-% C), A = lab[[i]], B = round(x[i, 1L], 1L), C = round(x[i, 2L], 2L)))
	}
	invisible()
}

