glmfit_aicloglik <-
function(object, y, x, wt) {
    fam <- object$family
    nobs <- NROW(x)
    n <- if (NCOL(y) == 1) 
        rep.int(1, nobs) else rowSums(y)
    mu <- fam$linkinv((x %*% object$coefficients)[, 1L])
    dev <- sum(fam$dev.resids(y, mu, wt))
    aic <- fam$aic(y, n, mu, wt, dev) + 2 * object$rank
    p <- object$rank
    if (fam$family %in% c("gaussian", "Gamma", "inverse.gaussian")) 
        p <- p + 1
    ll <- p - aic/2
    # c(aic = aic, loglik = ll, nobs = nobs, df = p)
    c(aic, ll, p)
} 


arm.glm <-
function(object, R = 250, weight.by = c("aic", "loglik"), trace = FALSE) {
	if(!inherits(object, c("glm", "lm")))
		stop("'object' must be a \"glm\" or \"lm\" object")
		
	maxtrials <- 10L
	weight.by <- switch(match.arg(weight.by), aic = 1L, loglik = 2L)
	allterms <- getAllTerms(object)
	ordtrm <- attr(allterms, "order")
	deps <- attr(allterms, "deps")[ordtrm, ordtrm]
	nterms <- length(allterms)
	mm <- model.matrix(object)
	n1 <- ceiling((nall <- nrow(mm))/2)
	n2 <- nall - n1 

	combinTerms <- lapply(0L:(2^nterms - 1L), function(j)
		as.logical(intToBits(j)[1L:nterms]))
	combinTerms <- combinTerms[vapply(combinTerms,
		formula_margin_check, FALSE, deps)]
	## combinPredictors
	combin <- vapply(combinTerms, function(x, assign)
		assign %in% c(0, which(x)), logical(ncol(mm)),
		assign = attr(mm, "assign"))
	combin <- combin[, i <- colSums(combin) < min(n1, n2) - 1L, drop = FALSE]
	combinTerms <- combinTerms[i]

	nModels <- ncol(combin)
	fam <- family(object)
	y <- object$y
	if(is.null(y)) y <- get.response(object)

	prwt <- weights(object, "prior")
	if(is.null(prwt)) prwt <- rep(1, nall) 
	
	if(fam$family == 'binomial') {
		y <- prwt * cbind(y, 1 - y, deparse.level = 0L)
		yvectorize <- function(y) y[, 1L] / rowSums(y)
		prwt <- rep(1, NROW(y))
		jresp <- c(1L, 2L)
	} else {
		yvectorize <- function(y) y[, 1L] ## XXX === drop  ?
		jresp <- 1L
	}
	Z <- cbind(y, mm, prwt, deparse.level = 0L)
	rownames(Z) <- NULL
	jprwts <- ncol(Z)
	jterms <- (1L:(ncol(Z) - 1L))[-jresp]
	
	traceinfo <- if(trace) { function(...) cat(..., sep = "") } else function(...) {}
	
	## begin iteration
	wts <- matrix(NA_real_, ncol = nModels, nrow = R)
	for(iter in 1L:R) {
		traceinfo("iteration=", iter, "\n")
		for (u in 1L:maxtrials) {
			traceinfo("  trial=", u, "\n    ")
			i <- sample(nall, n1)
			y1 <- Z[i, jresp, drop = FALSE]
			y2 <- Z[-i, jresp, drop = FALSE]
			vy2 <- yvectorize(y2)
			x1 <- Z[i, jterms, drop = FALSE]
			x2 <- Z[-i, jterms, drop = FALSE]
			prwts1 <- Z[i, jprwts]
			prwts2 <- Z[-i, jprwts]
			ic <- numeric(nModels)
			for (k in seq.int(nModels)) {
				traceinfo("k=", k, " ")
				fit1 <- glm.fit(x1[, combin[, k], drop = FALSE], y1, family = fam, weights = prwts1)
				if (problem <- (any(is.na(fit1$coefficients)) || !fit1$converged)) {
					traceinfo("<!>")
					break
				}
				ic[k] <- glmfit_aicloglik(fit1, vy2, x2[, combin[, k], drop = FALSE], 
					prwts2)[1L]
			}
			traceinfo("\n")
			if (!problem) 
				break
		}
		d <- exp(-ic/2)
		wts[iter, ] <- d/sum(d)
	} 

	wts[!is.finite(wts)] <- NA_real_
	wts <- wts[round(rowSums(wts, na.rm = TRUE)) == 1, ]
	wts <- wts/rowSums(wts, na.rm = TRUE)
	wtsmean <- colMeans(wts)
	
	cfnames <- colnames(Z[, jterms])
	x <- Z[, jterms, drop = FALSE]
	y <- Z[, jresp, drop = FALSE]
	msTable <- matrix(NA_real_, nrow = nModels, ncol = 6L, dimnames = list(1L:nModels, 
		c("df", "logLik", "AIC", "delta", "weight", "ARM weight")))
	## coefArray(model, c(coef, se, df), coefficients)
	coefArray <- array(NA_real_, dim = c(nModels, 3L, length(cfnames)), dimnames = list(1L:nModels, 
		c("Estimate", "Std. Error", "df"), cfnames))
	for (k in seq.int(nModels)) {
		fit1 <- glm.fit(x[, combin[, k], drop = FALSE], y, family = fam)
		coefArray[k, , combin[, k]] <- rbind(t(summary.glm(fit1)$coefficients[, 1L:2L, 
			drop = FALSE]), fit1$df.residual)
		msTable[k, c(3L, 2L, 1L)] <- glmfit_aicloglik(fit1, fit1$y, x[, combin[, k], 
			drop = FALSE], fit1$prior.weights)
	}
	msTable[, 4L] <- msTable[, 3L] - min(msTable[, 3L])
	msTable[, 5L] <- Weights(msTable[, 3L])
	msTable[, 6L] <- wtsmean
	
	cfmat <- coefArray[, 1L, ]
	cfmat[is.na(cfmat)]<- 0
	coefMat <- array(dim = c(2L, ncol(cfmat)),
		dimnames = list(c("full", "subset"), colnames(cfmat)))
	coefMat[1L, ] <- drop(wtsmean %*% cfmat)
	
	#debug <- list(wtsmean = wtsmean, cfmat = cfmat, coefArray = coefArray)
	
	ass <- attr(mm, "assign")
	bp <- !is.na(coefArray[, 1L, ass != 0L & !duplicated(ass)])
	tenm <- allterms[ordtrm]
	allmodelnames <- .modelNames(allTerms = apply(bp, 1L, function(z) tenm[z]),
						uqTerms = tenm)
	rownames(msTable) <- c(allmodelnames)
	
	ordmod <- order(msTable[,4L], decreasing = FALSE)
		
	rval <- list(
		msTable = structure(as.data.frame(msTable[ordmod, ]),
			term.codes = attr(allmodelnames, "variables")),
		coefficients = coefMat,
		coefArray = coefArray[ordmod, , ],
		importance = {
			structure(wtsmean %*% bp,
				n.models = structure(colSums(bp), names = tenm),
				names = tenm, class = "importance")
		 },
		formula = object$formula,
		call = match.call() 
	    #, debug = debug
	)
		
	attr(rval, "rank") <- .getRank(AIC)  ## TODO
	attr(rval, "nobs") <- nrow(x)
	attr(rval, "beta") <- "none"
	attr(rval, "revised.var") <- TRUE
	attr(rval, "arm") <- TRUE

	class(rval) <- "averaging"
	rval 
}