terms.averaging <-
function (x, ...) {
	terms(formula(x))
}

## helper function
toarray <- function(x)
	array(unlist(x), dim = c(dim(x[[1L]]), length(x)),
	 dimnames = c(dimnames(x[[1L]]),  list(names(x))))

mergeContrasts <-
function(models) {
	ctr <- lapply(models, get.contrasts)
	ctrnm <- unique(unlist(lapply(ctr, names)))
	ret <- structure(vector("list", length = length(ctrnm)), names = ctrnm)
	for(x in ctr) {
		for(i in names(x)) {
			if(is.null(ret[[i]])) {
				ret[[i]] <- x[[i]]
			} else 
				if(!identical(ret[[i]], x[[i]]))
					stop(gettextf("inconsistent contrasts in '%s'", i) )
		}
	}
	ret
}

mergeMF <-
function(models) {
	mf <- model.frame(models[[1L]])
	mfNames <- colnames(mf)
	lhs <- asChar(getResponseFormula(mf))

	f <- attr(fixTermsObject(terms(mf)), "term.labels")
	#m <- models[[2]]
	for(m in models[-1L]) {
		mf1 <- model.frame(m)
		if(!identical(lhs, lhs1 <- asChar(getResponseFormula(terms(mf1)))))
			stop("response differs between models: ", sQuote(c(asChar(lhs), lhs1)))
		mf <- cbind(mf, mf1[, !(colnames(mf1) %in% mfNames), drop = FALSE])
		tt1 <- fixTermsObject(terms(mf1))
		f <- c(f, attr(tt1, "term.labels"))
		if(!is.null(attr(tt1,"offset")))
			f <- c(f, sapply(as.list(attr(tt1,"variables")[attr(tt1,"offset") + 1L]), asChar))
	}
	
	f <- reformulate(f[!duplicated(f)])
	f <- as.formula(call("~", parse(text = lhs, keep.source = FALSE)[[1L]], f[[length(f)]]))
	environment(f) <- environment(formula(models[[1L]]))
	tt <- fixTermsObject(terms(f))
	mf <- mf[, rownames(attr(tt, "factors"))]
	attr(tt, "dataClasses") <- vapply(mf, .MFclass, "")
	attr(mf, "terms") <- tt
	mf
}

offsetTermNames <- function(x)
vapply(as.list(attr(x, "variables")[attr(x,"offset") + 1L]), deparse, "", control = NULL)

offsetWeights <-
function(wts, Terms, models) {
	if(is.null(off <- attr(Terms, "offset")))
		return(NULL)
	offnames <- rownames(attr(Terms, "factors"))[off]
	n <- length(offnames)
	v <-  matrix(vapply(models, function(x) {
		offnames %in% offsetTermNames(terms(x))
		}, logical(n)), nrow = n)
	offwts <- vapply(1L:n, function(i) sum(wts[v[i, ]]), 0)
	names(offwts) <- offnames
	offwts
}

## orders terms alphabetically in 'terms' object
fixTermsObject <-
function(x, peel = TRUE) {
	
	peelfun <- function(z) if(is.call(z))
		paste0(as.character(z[-1L]), collapse = " ") else as.character(z)
	
	factors <- attr(x, "factors")
    if (length(factors) != 0L) {
		
		z <- rep(1L, nrow(factors))
		z[attr(x, "response")] <- 0L
		if(hasOff <- !is.null(attr(x, "offset")))
			z[attr(x, "offset")] <- 2L
		charvar <- if (peel) sapply(as.list(attr(x, "variables")[-1L]), peelfun) else
			rownames(factors)
		ov <- order(z, charvar)
		factors <- factors[ov, , drop = FALSE]
		charvar <- charvar[ov]
		v <- rownames(factors)
		
		lab <- lab_ord <- character(ncol(factors))
		lfac <- factors != 0L
		for(i in 1L:ncol(factors)) {
			j <- lfac[, i]
			lab[i] <- paste0(v[j], collapse = ":")
			lab_ord[i] <- paste0(charvar[j], collapse = ":")
		}
		o <- order(attr(x, "order"), lab_ord)
		lab <- lab[o]
		ans <- reformulate(c(lab, offsetTermNames(x)), intercept = attr(x, "intercept"))
		
		if(attr(x,"response") != 0) ans <-
			as.formula(call("~", attr(x, "variables")[[attr(x, "response") + 1L]], ans[[2L]]))
		attributes(ans) <- attributes(x)
		attr(ans, "factors") <- factors[, o, drop = FALSE]
		colnames(attr(ans,"factors")) <- attr(ans, "term.labels") <- lab
		
		if(hasOff) attr(ans, "offset") <- which(z[ov] == 2L)
		for(j in c("variables", "predvars"))
			if(!is.null(attr(ans, j))) attr(ans, j) <- attr(ans, j)[c(1L, ov + 1L)]
		if(!is.null(attr(ans, "dataClasses")))
			attr(ans, "dataClasses") <- attr(ans, "dataClasses")[ov]
		ans		
	} else x
}

get.contrasts <-
function(x) UseMethod("get.contrasts")

get.contrasts.lm <-
function(x) x$contrasts

get.contrasts.averaging <-
function(x) {
	mergeContrasts(getModelList(x))
}


getModelList <-
function(object, error = TRUE) {
	if(is.null(models <- attr(object, "modelList")))
		if(error) stop("component models not included in this 'averaging' object")
	invisible(models)
}

