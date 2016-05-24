#---------------------------------------------------------------------------------------
# Revision history:
#   2009-01-16: replaced unlist(options("foo")) with getOption("foo")
#   2009-09-16: optionally allow models with aliased coefficients. J. Fox
#   2009-12-10: modification by A. Zeileis to allow wider range of coef. names.
#   2009-12-22: small changes to linearHypothesis.mlm() to handle user-specified
#               within-subjects designs in Anova()
#   2010-05-21: linearHypothesis.default() and .lm() changed so that differences
#               in df, etc. will be postive.
#   2010-06-12: linearHypothesis.mlm() changed to allow observation weights
#	2010-06-22: fixed bug in linearHypothesis.lm caused by 2010-05-21 revision
#   2010-01-21: added methods for mixed models; added matchCoefs() and methods. J. Fox
#   2011-05-03: fixed bug in displaying numbers starting with "-1" or "+1" in printed representation. J. Fox
#   2011-06-09: added matchCoefs.mlm(). J. Fox
#   2011-11-27: added linearHypothesis.svyglm(). John
#   2011-12-27: fixed printing bug in linearHypothesis(). John
#   2012-02-28: added F-test to linearHypothesis.mer(). John
#   2012-03-07: singular.ok argument added to linearHypothesis.mlm(). J. Fox
#   2012-08-20: Fixed p-value bug for chisq test in .mer method. John
#   2012-09-17: updated linearHypothesis.mer for pkrtest 0.3-2. John
#   2012-11-21: test for NULL rhs to avoid warning in R 2.16.0. John
#   2013-01-28: hypotheses can now contain newlines and tabs
#   2013-02-14: fixed bug in printing constants of the form 1.x*. John
#   2013-06-20: added .merMod() method. John
#   2013-06-22: tweaks for lme4. John
#   2013-06-22: test argument uniformly uses "Chisq" rather than "chisq". J. Fox
#   2013-08-19: removed calls to unexported functions in stats. J. Fox
#   2014-08-17: added call to requireNamespace() and :: as needed (doesn't work for pbkrtest). J. Fox
#   2014-08-18: fixed bug in linearHypothesis.survreg(). J. Fox
#   2014-09-23: added linearHypothesis.rlm. J. Fox
#   2014-12-18: check that residual df nonzero in Anova.lm() and Anova.default
#               and residual SS nonzero in Anova.lm(). John
#   2015-01-27: KRmodcomp() and methods now imported from pbkrtest. John
#   2015-02-03: Check for NULL df before 0 df in default method. John
#---------------------------------------------------------------------------------------

vcov.default <- function(object, ...){
	stop(paste("there is no vcov() method for models of class",
					paste(class(object), collapse=", ")))
}

has.intercept.matrix <- function (model, ...) {
	"(Intercept)" %in% colnames(model)
}


makeHypothesis <- function(cnames, hypothesis, rhs = NULL){
	parseTerms <- function(terms){
		component <- gsub("^[-\\ 0-9\\.]+", "", terms)
		component <- gsub(" ", "", component, fixed=TRUE)
		component
	}
	stripchars <- function(x) {
	  x <- gsub("\\n", " ", x)
	  x <- gsub("\\t", " ", x)
		x <- gsub(" ", "", x, fixed = TRUE)
		x <- gsub("*", "", x, fixed = TRUE)
		x <- gsub("-", "+-", x, fixed = TRUE)
		x <- strsplit(x, "+", fixed = TRUE)[[1]]
		x <- x[x!=""]
		x
	}
	char2num <- function(x) {
		x[x == ""] <- "1"
		x[x == "-"] <- "-1"
		as.numeric(x)
	}
	constants <- function(x, y) { 
		with.coef <- unique(unlist(sapply(y,
								function(z) which(z == parseTerms(x)))))
		if (length(with.coef) > 0) x <- x[-with.coef]
		x <- if (is.null(x)) 0 else sum(as.numeric(x))
		if (any(is.na(x)))
			stop('The hypothesis "', hypothesis,
					'" is not well formed: contains bad coefficient/variable names.')
		x
	}
	coefvector <- function(x, y) {
		rv <- gsub(" ", "", x, fixed=TRUE) ==
				parseTerms(y)
		if (!any(rv)) return(0)
		if (sum(rv) > 1) stop('The hypothesis "', hypothesis,
					'" is not well formed.')
		rv <- sum(char2num(unlist(strsplit(y[rv], x, fixed=TRUE))))
		if (is.na(rv))
			stop('The hypothesis "', hypothesis,
					'" is not well formed: contains non-numeric coefficients.')
		rv
	}
	
	if (!is.null(rhs)) rhs <- rep(rhs, length.out = length(hypothesis))
	if (length(hypothesis) > 1)
		return(rbind(Recall(cnames, hypothesis[1], rhs[1]),
						Recall(cnames, hypothesis[-1], rhs[-1])))
	
	cnames_symb <- sapply(c("@", "#", "~"), function(x) length(grep(x, cnames)) < 1)
	
	if(any(cnames_symb)) {
		cnames_symb <- head(c("@", "#", "~")[cnames_symb], 1)
		cnames_symb <- paste(cnames_symb, seq_along(cnames), cnames_symb, sep = "")
		hypothesis_symb <- hypothesis
		for(i in order(nchar(cnames), decreasing = TRUE))
			hypothesis_symb <- gsub(cnames[i], cnames_symb[i], hypothesis_symb, fixed = TRUE)
	} else {
		stop('The hypothesis "', hypothesis,
				'" is not well formed: contains non-standard coefficient names.')
	}
	
	lhs <- strsplit(hypothesis_symb, "=", fixed=TRUE)[[1]] 
	if (is.null(rhs)) {
		if (length(lhs) < 2) rhs <- "0"
		else if (length(lhs) == 2) {
			rhs <- lhs[2]
			lhs <- lhs[1]
		}
		else stop('The hypothesis "', hypothesis,
					'" is not well formed: contains more than one = sign.')
	}
	else {
		if (length(lhs) < 2) as.character(rhs)
		else stop('The hypothesis "', hypothesis,
					'" is not well formed: contains a = sign although rhs was specified.')
	}
	lhs <- stripchars(lhs)
	rhs <- stripchars(rhs)
	rval <- sapply(cnames_symb, coefvector, y = lhs) - sapply(cnames_symb, coefvector, y = rhs) 
	rval <- c(rval, constants(rhs, cnames_symb) - constants(lhs, cnames_symb)) 
	names(rval) <- c(cnames, "*rhs*")
	rval
}

printHypothesis <- function(L, rhs, cnames){
	hyp <- rep("", nrow(L))
	for (i in 1:nrow(L)){
		sel <- L[i,] != 0
		h <- L[i, sel]
		h <- ifelse(h < 0, as.character(h), paste("+", h, sep=""))
		nms <- cnames[sel]
		h <- paste(h, nms)
		h <- gsub("-", " - ", h)
		h <- gsub("+", "  + ", h, fixed=TRUE)
		h <- paste(h, collapse="")
		h <- gsub("  ", " ", h, fixed=TRUE)
		h <- sub("^\\ \\+", "", h)
		h <- sub("^\\ ", "", h)
		h <- sub("^-\\ ", "-", h)	
		h <- paste(" ", h, sep="")
		h <- paste(h, "=", rhs[i])		
		h <- gsub(" 1([^[:alnum:]_.]+)[ *]*", "", 
				gsub("-1([^[:alnum:]_.]+)[ *]*", "-", 
						gsub("- +1 +", "-1 ", h)))
		h <- sub("Intercept)", "(Intercept)", h)		
		h <- gsub("-", " - ", h)
		h <- gsub("+", "  + ", h, fixed=TRUE)
		h <- gsub("  ", " ", h, fixed=TRUE)
		h <- sub("^ *", "", h)
		hyp[i] <- h
	}
	hyp
}

linearHypothesis <- function (model, ...)
	UseMethod("linearHypothesis")

lht <- function (model, ...)
	UseMethod("linearHypothesis")
	
linearHypothesis.nlsList <- function(model,  ..., vcov., coef.){
   vcov.nlsList <- function(object, ...) {
       vlist <- lapply(object, vcov)
       ng <- length(vlist)
       nv <- dim(vlist[[1]])[1]
       v <- matrix(0, nrow=ng*nv, ncol=ng*nv)
       for (j in 1:ng){
          cells <- ((j-1)*nv + 1):(j*nv)
          v[cells, cells] <- vlist[[j]]
        }
      v
      }
   linearHypothesis.default(model, vcov.=vcov.nlsList(model), 
       coef.=unlist(lapply(model, coef)), ...)}


linearHypothesis.default <- function(model, hypothesis.matrix, rhs=NULL,
		test=c("Chisq", "F"), vcov.=NULL, singular.ok=FALSE, verbose=FALSE, 
    coef. = coef(model), ...){
	df <- df.residual(model)
	if (is.null(df)) df <- Inf ## if no residual df available
    if (df == 0) stop("residual df = 0")
	V <- if (is.null(vcov.)) vcov(model)
			else if (is.function(vcov.)) vcov.(model) else vcov.
	b <- coef.
	if (any(aliased <- is.na(b)) && !singular.ok)
		stop("there are aliased coefficients in the model")
	b <- b[!aliased]
	if (is.null(b)) stop(paste("there is no coef() method for models of class",
						paste(class(model), collapse=", ")))
	if (is.character(hypothesis.matrix)) {
		L <- makeHypothesis(names(b), hypothesis.matrix, rhs)
		if (is.null(dim(L))) L <- t(L)
		rhs <- L[, NCOL(L)]
		L <- L[, -NCOL(L), drop = FALSE]
		rownames(L) <- hypothesis.matrix
	}
	else {
		L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix)
				else hypothesis.matrix
		if (is.null(rhs)) rhs <- rep(0, nrow(L))
	}
	q <- NROW(L)
	if (verbose){
		cat("\nHypothesis matrix:\n")
		print(L)
		cat("\nRight-hand-side vector:\n")
		print(rhs)
		cat("\nEstimated linear function (hypothesis.matrix %*% coef - rhs)\n")
		print(drop(L %*% b - rhs))
		cat("\n")
	}
	SSH <- as.vector(t(L %*% b - rhs) %*% solve(L %*% V %*% t(L)) %*% (L %*% b - rhs))
	test <- match.arg(test)
	if (!(is.finite(df) && df > 0)) test <- "Chisq"
	name <- try(formula(model), silent = TRUE)
	if (inherits(name, "try-error")) name <- substitute(model)
	title <- "Linear hypothesis test\n\nHypothesis:"
	topnote <- paste("Model 1: restricted model","\n", "Model 2: ", 
			paste(deparse(name), collapse = "\n"), sep = "")
	note <- if (is.null(vcov.)) ""
			else "\nNote: Coefficient covariance matrix supplied.\n"
	rval <- matrix(rep(NA, 8), ncol = 4)
	colnames(rval) <- c("Res.Df", "Df", test, paste("Pr(>", test, ")", sep = ""))
	rownames(rval) <- 1:2
	rval[,1] <- c(df+q, df)
	if (test == "F") {
		f <- SSH/q
		p <- pf(f, q, df, lower.tail = FALSE)
		rval[2, 2:4] <- c(q, f, p)
	}
	else {
		p <- pchisq(SSH, q, lower.tail = FALSE)
		rval[2, 2:4] <- c(q, SSH, p)
	}
	if (!(is.finite(df) && df > 0)) rval <- rval[,-1]
	structure(as.data.frame(rval),
			heading = c(title, printHypothesis(L, rhs, names(b)), "", topnote, note),
			class = c("anova", "data.frame"))
}

linearHypothesis.glm <- function(model, ...)
	linearHypothesis.default(model, ...)

linearHypothesis.lm <- function(model, hypothesis.matrix, rhs=NULL,
		test=c("F", "Chisq"), vcov.=NULL,
		white.adjust=c(FALSE, TRUE, "hc3", "hc0", "hc1", "hc2", "hc4"),
		singular.ok=FALSE, ...){
    if (df.residual(model) == 0) stop("residual df = 0")
    if (deviance(model) < sqrt(.Machine$double.eps)) stop("residual sum of squares is 0 (within rounding error)")
	if (!singular.ok && is.aliased(model))
		stop("there are aliased coefficients in the model.")
	test <- match.arg(test)
	white.adjust <- as.character(white.adjust)
	white.adjust <- match.arg(white.adjust)
	if (white.adjust != "FALSE"){
		if (white.adjust == "TRUE") white.adjust <- "hc3"
		vcov. <- hccm(model, type=white.adjust)
	}
	rval <- linearHypothesis.default(model, hypothesis.matrix, rhs = rhs,
			test = test, vcov. = vcov., singular.ok=singular.ok, ...)
	if (is.null(vcov.)) {
		rval2 <- matrix(rep(NA, 4), ncol = 2)
		colnames(rval2) <- c("RSS", "Sum of Sq")
		SSH <- rval[2,test]
		if (test == "F") SSH <- SSH * abs(rval[2, "Df"])
		df <- rval[2, "Res.Df"]
		error.SS <- deviance(model)
		rval2[,1] <- c(error.SS + SSH * error.SS/df, error.SS)
		rval2[2,2] <- abs(diff(rval2[,1]))
		rval2 <- cbind(rval, rval2)[,c(1, 5, 2, 6, 3, 4)]
		class(rval2) <- c("anova", "data.frame")
		attr(rval2, "heading") <- attr(rval, "heading")
		rval <- rval2
	}
	rval
}


check.imatrix <- function(X, terms){ 
# check block orthogonality of within-subjects model matrix
	XX <- crossprod(X)
	if (missing(terms)) terms <- attr(X, "assign")
	for (term in unique(terms)){
		subs <- term == terms
		XX[subs, subs] <- 0
	}
	if (any(abs(XX) > sqrt(.Machine$double.eps)))
		stop("Terms in the intra-subject model matrix are not orthogonal.")
}

linearHypothesis.mlm <- function(model, hypothesis.matrix, rhs=NULL, SSPE, V,
		test, idata, icontrasts=c("contr.sum", "contr.poly"), idesign, iterms,
		check.imatrix=TRUE, P=NULL, title="", singular.ok=FALSE, verbose=FALSE, ...){
	if (missing(test)) test <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
	test <- match.arg(test, c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
			several.ok=TRUE)
	df.residual <- df.residual(model)
	wts <- if (!is.null(model$weights)) model$weights else rep(1,nrow(model.matrix(model)))
	# V = (X'WX)^{-1}
	if (missing (V)) V <- solve(wcrossprod(model.matrix(model), w=wts))
	B <- coef(model)
	if (is.character(hypothesis.matrix)) {
		L <- makeHypothesis(rownames(B), hypothesis.matrix, rhs)
		if (is.null(dim(L))) L <- t(L)
		L <- L[, -NCOL(L), drop = FALSE]
		rownames(L) <- hypothesis.matrix
	}
	else {
		L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix)
				else hypothesis.matrix
	}
	# SSPE = E'WE
	if (missing(SSPE)) SSPE <- wcrossprod(residuals(model),w=wts)
	if (missing(idata)) idata <- NULL
	if (missing(idesign)) idesign <- NULL
	if (!is.null(idata)){
		for (i in 1:length(idata)){
			if (is.null(attr(idata[,i], "contrasts"))){
				contrasts(idata[,i]) <- if (is.ordered(idata[,i])) icontrasts[2]
						else icontrasts[1]
			}
		}
		if (is.null(idesign)) stop("idesign (intra-subject design) missing.")
		X.design <- model.matrix(idesign, data=idata)
		if (check.imatrix) check.imatrix(X.design)
		intercept <- has.intercept(X.design)
		term.names <- term.names(idesign)
		if (intercept) term.names <- c("(Intercept)", term.names)
		which.terms <- match(iterms, term.names)
		if (any(nas <- is.na(which.terms))){
			if (sum(nas) == 1)
				stop('The term "', iterms[nas],'" is not in the intrasubject design.')
			else stop("The following terms are not in the intrasubject design: ",
						paste(iterms[nas], collapse=", "), ".")
		}
		select <- apply(outer(which.terms, attr(X.design, "assign") + intercept, "=="),
				2, any)
		P <- X.design[, select, drop=FALSE]
	}
	if (!is.null(P)){
		rownames(P) <- colnames(B)
		SSPE <- t(P) %*% SSPE %*% P
		B <- B %*% P
	}
	rank <- sum(eigen(SSPE, only.values=TRUE)$values >= sqrt(.Machine$double.eps))
	if (!singular.ok && rank < ncol(SSPE))
		stop("The error SSP matrix is apparently of deficient rank = ",
				rank, " < ", ncol(SSPE))
	r <- ncol(B)
	if (is.null(rhs)) rhs <- matrix(0, nrow(L), r)
	rownames(rhs) <- rownames(L)
	colnames(rhs) <- colnames(B)
	q <- NROW(L)
	if (verbose){
		cat("\nHypothesis matrix:\n")
		print(L)
		cat("\nRight-hand-side matrix:\n")
		print(rhs)
		cat("\nEstimated linear function (hypothesis.matrix %*% coef - rhs):\n")
		print(drop(L %*% B - rhs))
		cat("\n")
	}
	SSPH <- t(L %*% B - rhs) %*% solve(L %*% V %*% t(L)) %*% (L %*% B - rhs)
	rval <- list(SSPH=SSPH, SSPE=SSPE, df=q, r=r, df.residual=df.residual, P=P,
			title=title, test=test, singular=rank < ncol(SSPE))
	class(rval) <- "linearHypothesis.mlm"
	rval
}

#linearHypothesis.mlm <- function(model, hypothesis.matrix, rhs=NULL, SSPE, V,
#   test, idata, icontrasts=c("contr.sum", "contr.poly"), idesign, iterms,
#   check.imatrix=TRUE, P=NULL, title="", verbose=FALSE, ...){
#   if (missing(test)) test <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
#   test <- match.arg(test, c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
#       several.ok=TRUE)
#   df.residual <- df.residual(model)
#   if (missing (V)) V <- solve(crossprod(model.matrix(model)))
#   B <- coef(model)
#   if (is.character(hypothesis.matrix)) {
#       L <- makeHypothesis(rownames(B), hypothesis.matrix, rhs)
#       if (is.null(dim(L))) L <- t(L)
#       L <- L[, -NCOL(L), drop = FALSE]
#       rownames(L) <- hypothesis.matrix
#   }
#   else {
#       L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix)
#           else hypothesis.matrix
#   }
#   if (missing(SSPE)) SSPE <- crossprod(residuals(model))
#   if (missing(idata)) idata <- NULL
#   if (missing(idesign)) idesign <- NULL
#   if (!is.null(idata)){
#       for (i in 1:length(idata)){
#           if (is.null(attr(idata[,i], "contrasts"))){
#               contrasts(idata[,i]) <- if (is.ordered(idata[,i])) icontrasts[2]
#                   else icontrasts[1]
#           }
#       }
#       if (is.null(idesign)) stop("idesign (intra-subject design) missing.")
#       X.design <- model.matrix(idesign, data=idata)
#       if (check.imatrix) check.imatrix(X.design)
#       intercept <- has.intercept(X.design)
#       term.names <- term.names(idesign)
#       if (intercept) term.names <- c("(Intercept)", term.names)
#       which.terms <- match(iterms, term.names)
#       if (any(nas <- is.na(which.terms))){
#           if (sum(nas) == 1)
#               stop('The term "', iterms[nas],'" is not in the intrasubject design.')
#           else stop("The following terms are not in the intrasubject design: ",
#                   paste(iterms[nas], collapse=", "), ".")
#       }
#       select <- apply(outer(which.terms, attr(X.design, "assign") + intercept, "=="),
#           2, any)
#       P <- X.design[, select, drop=FALSE]
#   }
#   if (!is.null(P)){
#       rownames(P) <- colnames(B)
#       SSPE <- t(P) %*% SSPE %*% P
#       B <- B %*% P
#   }
#   rank <- sum(eigen(SSPE, only.values=TRUE)$values >= sqrt(.Machine$double.eps))
#   if (rank < ncol(SSPE))
#       stop("The error SSP matrix is apparently of deficient rank = ",
#           rank, " < ", ncol(SSPE))
#   r <- ncol(B)
#   if (is.null(rhs)) rhs <- matrix(0, nrow(L), r)
#   rownames(rhs) <- rownames(L)
#   colnames(rhs) <- colnames(B)
#   q <- NROW(L)
#   if (verbose){
#       cat("\nHypothesis matrix:\n")
#       print(L)
#       cat("\nRight-hand-side matrix:\n")
#       print(rhs)
#       cat("\nEstimated linear function (hypothesis.matrix %*% coef - rhs):\n")
#       print(drop(L %*% B - rhs))
#       cat("\n")
#   }
#   SSPH <- t(L %*% B - rhs) %*% solve(L %*% V %*% t(L)) %*% (L %*% B - rhs)
#   rval <- list(SSPH=SSPH, SSPE=SSPE, df=q, r=r, df.residual=df.residual, P=P,
#       title=title, test=test)
#   class(rval) <- "linearHypothesis.mlm"
#   rval
#}

print.linearHypothesis.mlm <- function(x, SSP=TRUE, SSPE=SSP,
		digits=getOption("digits"), ...){
	test <- x$test
	if (!is.null(x$P) && SSP){
		P <- x$P
		cat("\n Response transformation matrix:\n")
		attr(P, "assign") <- NULL
		attr(P, "contrasts") <- NULL
		print(P, digits=digits)
	}
	if (SSP){
		cat("\nSum of squares and products for the hypothesis:\n")
		print(x$SSPH, digits=digits)
	}
	if (SSPE){
		cat("\nSum of squares and products for error:\n")
		print(x$SSPE, digits=digits)
	}
	if ((!is.null(x$singular)) && x$singular){
		warning("the error SSP matrix is singular; multivariate tests are unavailable")
		return(invisible(x))
	}
	SSPE.qr <- qr(x$SSPE)
	# the following code is adapted from summary.manova
	eigs <- Re(eigen(qr.coef(SSPE.qr, x$SSPH), symmetric = FALSE)$values)
	tests <- matrix(NA, 4, 4)
	rownames(tests) <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
	if ("Pillai" %in% test)
		tests[1, 1:4] <- Pillai(eigs, x$df, x$df.residual)
	if ("Wilks" %in% test)
		tests[2, 1:4] <- Wilks(eigs, x$df, x$df.residual)
	if ("Hotelling-Lawley" %in% test)
		tests[3, 1:4] <- HL(eigs, x$df, x$df.residual)
	if ("Roy" %in% test)
		tests[4, 1:4] <- Roy(eigs, x$df, x$df.residual)
	tests <- na.omit(tests)
	ok <- tests[, 2] >= 0 & tests[, 3] > 0 & tests[, 4] > 0
	ok <- !is.na(ok) & ok
	tests <- cbind(x$df, tests, pf(tests[ok, 2], tests[ok, 3], tests[ok, 4],
					lower.tail = FALSE))
	colnames(tests) <- c("Df", "test stat", "approx F", "num Df", "den Df", "Pr(>F)")
	tests <- structure(as.data.frame(tests),
			heading = paste("\nMultivariate Test",
					if (nrow(tests) > 1) "s", ": ", x$title, sep=""),
			class = c("anova", "data.frame"))
	print(tests, digits=digits)
	invisible(x)
}

linearHypothesis.survreg <- function(model, hypothesis.matrix, rhs=NULL,
		test=c("Chisq", "F"), vcov., verbose=FALSE, ...){
	if (missing(vcov.)) {
		vcov. <- vcov(model)
		p <- which(rownames(vcov.) == "Log(scale)")
		if (length(p) > 0) vcov. <- vcov.[-p, -p]
	}
	linearHypothesis.default(model, hypothesis.matrix, rhs, test, vcov., verbose=verbose, ...)
}

linearHypothesis.polr <- function (model, hypothesis.matrix, rhs=NULL, vcov., verbose=FALSE, ...){
	k <- length(coef(model))
	V <- vcov(model)[1:k, 1:k]
	linearHypothesis.default(model, hypothesis.matrix, rhs, vcov.=V, verbose=verbose, ...)
}

coef.multinom <- function(object, ...){
    # the following local function is copied from nnet:::coef.multinom
    coef.m <- function (object, ...) {
            r <- length(object$vcoefnames)
            if (length(object$lev) == 2L) {
                coef <- object$wts[1L + (1L:r)]
                names(coef) <- object$vcoefnames
            }
            else {
                coef <- matrix(object$wts, nrow = object$n[3L], byrow = TRUE)[, 
                                                                              1L + (1L:r), drop = FALSE]
                if (length(object$lev)) 
                    dimnames(coef) <- list(object$lev, object$vcoefnames)
                if (length(object$lab)) 
                    dimnames(coef) <- list(object$lab, object$vcoefnames)
                coef <- coef[-1L, , drop = FALSE]
            }
            coef
        }
    
	b <- coef.m(object, ...)
	cn <- colnames(b)
	rn <- rownames(b)
	b <- as.vector(t(b))
	names(b) <- as.vector(outer(cn, rn, function(c, r) paste(r, c, sep=":")))
	b
}

## functions for mixed models

linearHypothesis.merMod <- function(model, hypothesis.matrix, rhs=NULL,
                                 vcov.=NULL, test=c("Chisq", "F"), 
                                 singular.ok=FALSE, verbose=FALSE, ...){
    linearHypothesis.mer(model=model, hypothesis.matrix=hypothesis.matrix,
                         vcov.=vcov., test=test, singular.ok=singular.ok,
                         verbose=verbose, ...)
}

linearHypothesis.mer <- function(model, hypothesis.matrix, rhs=NULL,
                                 vcov.=NULL, test=c("Chisq", "F"), singular.ok=FALSE, verbose=FALSE, ...){
    test <- match.arg(test)
    V <- as.matrix(if (is.null(vcov.))vcov(model)
                   else if (is.function(vcov.)) vcov.(model) else vcov.)
    b <- fixef(model)
    if (any(aliased <- is.na(b)) && !singular.ok)
        stop("there are aliased coefficients in the model")
    b <- b[!aliased]
    if (is.character(hypothesis.matrix)) {
        L <- makeHypothesis(names(b), hypothesis.matrix, rhs)
        if (is.null(dim(L))) L <- t(L)
        rhs <- L[, NCOL(L)]
        L <- L[, -NCOL(L), drop = FALSE]
        rownames(L) <- hypothesis.matrix
    }
    else {
        L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix)
        else hypothesis.matrix
        if (is.null(rhs)) rhs <- rep(0, nrow(L))
    }
    q <- NROW(L)
    if (verbose){
        cat("\nHypothesis matrix:\n")
        print(L)
        cat("\nRight-hand-side vector:\n")
        print(rhs)
        cat("\nEstimated linear function (hypothesis.matrix %*% coef - rhs)\n")
        print(drop(L %*% b - rhs))
        cat("\n")
    }
    if (test == "Chisq"){
        df <- Inf
        SSH <- as.vector(t(L %*% b - rhs) %*% solve(L %*% V %*% t(L)) %*% (L %*% b - rhs))
    }
    else {
        if (!requireNamespace("lme4")) stop("lme4 package is missing")
#        if (!require("pbkrtest") || packageVersion("pbkrtest") < "0.3.2") stop("pbkrtest package version >= 0.3.2 required for F-test on linear mixed model")
        if (!lme4::isREML(model)) 
            stop("F test available only for linear mixed model fit by REML")
#        res <- pbkrtest::KRmodcomp(model, L)$test
        res <- KRmodcomp(model, L)$test
        df <- res["Ftest", "ddf"]
        F <- res["Ftest", "stat"]
        p <- res["Ftest", "p.value"]
    }
    name <- try(formula(model), silent = TRUE)
    if (inherits(name, "try-error")) name <- substitute(model)
    title <- "Linear hypothesis test\n\nHypothesis:"
    topnote <- paste("Model 1: restricted model","\n", "Model 2: ", 
                     paste(deparse(name), collapse = "\n"), sep = "")
    note <- if (is.null(vcov.)) ""
    else "\nNote: Coefficient covariance matrix supplied.\n"
    rval <- matrix(rep(NA, 8), ncol = 4)
    if (test == "Chisq"){
        colnames(rval) <- c("Res.Df", "Df", "Chisq",  paste("Pr(>Chisq)", sep = ""))
        rownames(rval) <- 1:2
        rval[,1] <- c(df+q, df)
        p <- pchisq(SSH, q, lower.tail = FALSE)
        rval[2, 2:4] <- c(q, SSH, p)
        rval <- rval[,-1]
    }
    else{
        colnames(rval) <- c("Res.Df", "Df", "F",  paste("Pr(>F)", sep = ""))
        rownames(rval) <- 1:2
        rval[,1] <- c(df+q, df)
        rval[2, 2:4] <- c(q, F, p)
    }
    structure(as.data.frame(rval),
              heading = c(title, printHypothesis(L, rhs, names(b)), "", topnote, note),
              class = c("anova", "data.frame"))
}

linearHypothesis.lme <- function(model, hypothesis.matrix, rhs=NULL,
		vcov.=NULL, singular.ok=FALSE, verbose=FALSE, ...){
	V <- as.matrix(if (is.null(vcov.))vcov(model)
					else if (is.function(vcov.)) vcov.(model) else vcov.)
	b <- fixef(model)
	if (any(aliased <- is.na(b)) && !singular.ok)
		stop("there are aliased coefficients in the model")
	b <- b[!aliased]
	if (is.character(hypothesis.matrix)) {
		L <- makeHypothesis(names(b), hypothesis.matrix, rhs)
		if (is.null(dim(L))) L <- t(L)
		rhs <- L[, NCOL(L)]
		L <- L[, -NCOL(L), drop = FALSE]
		rownames(L) <- hypothesis.matrix
	}
	else {
		L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix)
				else hypothesis.matrix
		if (is.null(rhs)) rhs <- rep(0, nrow(L))
	}
	q <- NROW(L)
	if (verbose){
		cat("\nHypothesis matrix:\n")
		print(L)
		cat("\nRight-hand-side vector:\n")
		print(rhs)
		cat("\nEstimated linear function (hypothesis.matrix %*% coef - rhs)\n")
		print(drop(L %*% b - rhs))
		cat("\n")
	}
	df <- Inf
	SSH <- as.vector(t(L %*% b - rhs) %*% solve(L %*% V %*% t(L)) %*% (L %*% b - rhs))
	name <- try(formula(model), silent = TRUE)
	if (inherits(name, "try-error")) name <- substitute(model)
	title <- "Linear hypothesis test\n\nHypothesis:"
	topnote <- paste("Model 1: restricted model","\n", "Model 2: ", 
			paste(deparse(name), collapse = "\n"), sep = "")
	note <- if (is.null(vcov.)) ""
			else "\nNote: Coefficient covariance matrix supplied.\n"
	rval <- matrix(rep(NA, 8), ncol = 4)
	colnames(rval) <- c("Res.Df", "Df", "Chisq",  paste("Pr(>Chisq)", sep = ""))
	rownames(rval) <- 1:2
	rval[,1] <- c(df+q, df)
	p <- pchisq(SSH, q, lower.tail = FALSE)
	rval[2, 2:4] <- c(q, SSH, p)
	rval <- rval[,-1]
	structure(as.data.frame(rval),
			heading = c(title, printHypothesis(L, rhs, names(b)), "", topnote, note),
			class = c("anova", "data.frame"))
}

## for svyglm

linearHypothesis.svyglm <- function(model, ...) linearHypothesis.default(model, ...)

## for rlm

df.residual.rlm <- function(object, ...){
  p <- length(coef(object))
  wt.method <- object$call$wt.method
  if (!is.null(wt.method) && wt.method == "case") {
    sum(object$weights) - p
  }
  else length(object$wresid) - p
}

linearHypothesis.rlm <- function(model, ...) linearHypothesis.default(model, test="F", ...)


## matchCoefs

matchCoefs <- function(model, pattern, ...) UseMethod("matchCoefs")

matchCoefs.default <- function(model, pattern, coef.=coef, ...){
	names <- names(coef.(model))
	grep(pattern, names, value=TRUE)
}

matchCoefs.mer <- function(model, pattern, ...) NextMethod(coef.=fixef)

matchCoefs.merMod <- function(model, pattern, ...) NextMethod(coef.=fixef)

matchCoefs.lme <- function(model, pattern, ...) NextMethod(coef.=fixef)

matchCoefs.mlm <- function(model, pattern, ...){
	names <- rownames(coef(model))
	grep(pattern, names, value=TRUE)
}


