
`r.squaredGLMM` <-
function(x)
	UseMethod("r.squaredGLMM")

`r.squaredGLMM.default` <-
function(x) .NotYetImplemented()
	
`r.squaredGLMM.lme` <-
function(x) {
	VarFx <- var(fitted(x, level = 0L))
	mmRE <- model.matrix(x$modelStruct$reStruct,
						 data = x$data[rownames(x$fitted), , drop = FALSE])
	n <- nrow(mmRE)

	sigma2 <- x$sigma^2
	reStruct <- x$modelStruct$reStruct
	if((m <- length(reStruct)) > 1L) {
		nams <- names(reStruct)
		for(i in seq.int(m))
			attr(reStruct[[i]], "Dimnames")[[2L]] <-
				paste(nams[[i]], attr(reStruct[[i]], "Dimnames")[[2L]], sep = ".")
	}
	
	varRe <- sum(sapply(reStruct, function(z) {
		sig <- nlme::pdMatrix(z) * sigma2
		mm1 <-  mmRE[, rownames(sig), drop = FALSE]
		#sum(diag(mm1 %*% sig %*% t(mm1))) / n
		sum(matmultdiag(mm1 %*% sig, ty = mm1)) / n
	}))

	varTot <- sum(VarFx, varRe)
	res <- c(VarFx, varTot) / (varTot + sigma2)
	names(res) <- c("R2m", "R2c")
	res
}

## extracts random effect formula. e.g:
## ~ ... + (a | ...) + (b + c | ...) --> ~ a + b + c
ranform <- function (form) {
	ans <- update.formula(reformulate(vapply(lapply(.findbars(form),
		"[[", 2L), deparse, "")), ~ . + 1)
	environment(ans) <- environment(form)
	ans
}

`r.squaredGLMM.merMod` <-
function(x) {
	fam <- family(x)
	useObsLevVar <- (fam$family == "poisson" && fam$link == "log") || fam$family == "binomial"
	## for poisson(log) and binomial(*), update 'x' to include individual-level
	## variance (1 | 1:nobs(x)):
    if (useObsLevVar && !any(sapply(x@flist, nlevels) == nobs(x))) {
		cl <- get_call(x)
        frm <- formula(x)
		nRowData <- eval(call("eval", as.expression(call("NROW", cl$formula[[2L]])),
							  envir = cl$data), envir = environment(frm),
						 enclos = parent.frame())
		fl <- length(frm)
		frx <- . ~ . + 1
		frx[[3L]][[3L]] <- call("(", call("|", 1, call("gl", nRowData, 1)))
		cl$formula <- update.formula(frm, frx)		
		x <- tryCatch(eval(cl, envir = environment(frm), enclos = parent.frame()),
			error = function(e) {
				cry(conditionCall(e), conditionMessage(e), warn = TRUE)
				cry(cl, "fitting model with the observation-level random effect term failed. Add the term manually")
		})
		message("The result is correct only if all data used by the model ",
				"has not changed since model was fitted.", domain = "R-MuMIn")
    }


	mmAll <- model.matrix(ranform(formula(x)), data = model.frame(x))
		##Note: Argument 'contrasts' can only be specified for fixed effects
		##contrasts.arg = eval(cl$contrasts, envir = environment(formula(x))))	
	
	vc <- lme4::VarCorr(x)

	n <- nrow(mmAll)
	fx <- lme4::fixef(x) # fixed effect estimates
	fxpred <- as.vector(model.matrix(x) %*% fx)
	
	if (useObsLevVar) {
        vname <- names(x@flist)[sapply(x@flist, nlevels) == n][1L]
        varResid <- vc[[vname]][1L]
        beta0 <- mean(fxpred)
        vc <- vc[names(vc) != vname]
    } else {
        varResid <- attr(vc, "sc")^2
        beta0 <- NULL
    }
	
	if(!all(unlist(sapply(vc, rownames), use.names = FALSE) %in% colnames(mmAll)))
		stop("random term names do not match those in model matrix. \n",
			 "Have 'options(contrasts)' changed since the model was fitted?")
	
	varRe <- if(length(vc) == 0L) 0L else
		sum(sapply(vc, function(sig) {
			mm1 <-  mmAll[, rownames(sig), drop = FALSE]
			# sum(matmult(mm1 %*% sig, t(mm1), diag.only = TRUE)) / n
			sum(matmultdiag(mm1 %*% sig, ty = mm1)) / n
			#sum(diag(mm1 %*% sig %*% t(mm1))) / n
		}))
	
	#varRe0 <- if(length(vc) == 0L) 0L else
	#          sum(sapply(vc, function(sig) sig[[1]]))
		
	.rsqGLMM(fam = family(x), varFx = var(fxpred), varRe = varRe,
			 varResid = varResid, beta0 = beta0)
}

`r.squaredGLMM.glmmML` <-
function(x) {
	if(is.null(x$x))
		stop("glmmML must be fitted with 'x = TRUE'")
		
	fam <- family(x)
	useObsLevVar <- (fam$family == "poisson" && fam$link == "log") || fam$family == "binomial"
		if(useObsLevVar) {
			cry(, "cannot calculate 'unit variance' in glmmML")
	} 
	fxpred <- as.vector(x$x %*% coef(x))
	.rsqGLMM(family(x), varFx = var(fxpred), varRe = x$sigma^2, varResid = NULL,
			 beta0 = mean(fxpred))
}

`r.squaredGLMM.lm` <-
function(x) {
	fam <- family(x)
	.rsqGLMM(fam,
		 varFx = var(as.vector(model.matrix(x) %*% coef(x))),
		 #varFx = var(fitted(x)),
		 varRe = 0,
		 varResid = sum(if(is.null(x$weights)) resid(x)^2 else
					   resid(x)^2 * x$weights) / df.residual(x),
		 beta0 = if(fam$family == "poisson" && fam$link == "log") 
			log(mean(model.response(model.frame(x)))) else 
			NULL
		 )
}

`.rsqGLMM` <-
function (fam, varFx, varRe, varResid, beta0) {
    varDistr <- switch(paste(fam$family, fam$link, sep = "."), 
        gaussian.identity = 0,
		binomial.logit = 3.28986813369645, 
        binomial.probit = 1,
		poisson.log = {
            expBeta0 <- exp(beta0)
            if (expBeta0 < 6) cry(sys.call(-1L), "exp(beta0) of %0.1f is too close to zero, estimate may be unreliable \n", 
                expBeta0, warn = TRUE)
            log1p(1 / expBeta0)
        },
		poisson.sqrt = 0.25,
		cry(sys.call(-1L), "do not know how to calculate variance for this family/link combination")
		)
	
	#print(c(Sf = varFx, Sl = varRe, Se = varResid, Sd = varDistr))
	#  total.var <- [Sf + Sl] + [Se + Sd]
    varTot <- sum(varFx, varRe)
    res <- c(varFx, varTot) / (varTot + varDistr + varResid)
    names(res) <- c("R2m", "R2c")
    res
}
