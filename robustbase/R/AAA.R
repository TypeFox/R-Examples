
## if(getRversion() < "2.13") {
##     nobs <- function (object, ...) UseMethod("nobs")
##     ## also used for mlm fits  *and* lmrob :
##     nobs.lm <- function(object, ...)
## 	if(!is.null(w <- object$weights)) sum(w != 0) else NROW(object$residuals)
##     ## for glmrob :
##     nobs.glm <- function(object, ...) sum(!is.na(object$residuals))
## }

## Here and in NAMESPACE:
if(getRversion() < "3.1.0") {

## cut'n'paste from R's source src/library/stats/R/confint.R
format.perc <- function(probs, digits)
    ## Not yet exported, maybe useful in other contexts:
    ## quantile.default() sometimes uses a version of it
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
	  "%")

confint.lm <- function(object, parm, level = 0.95, ...)
{
    cf <- coef(object)
    pnames <- names(cf)
    if(missing(parm)) parm <- pnames
    else if(is.numeric(parm)) parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    fac <- qt(a, object$df.residual) # difference from default method
    pct <- format.perc(a, 3)
    ci <- array(NA, dim = c(length(parm), 2L),
		dimnames = list(parm, pct))
    ses <- sqrt(diag(vcov(object)))[parm] # gives NA for aliased parms
    ci[] <- cf[parm] + ses %o% fac
    ci
}

## cut'n'paste from R's source src/library/stats/R/dummy.coef.R
dummy.coef.lm <- function(object, use.na=FALSE, ...)
{
    Terms <- terms(object)
    tl <- attr(Terms, "term.labels")
    int <- attr(Terms, "intercept")
    facs <- attr(Terms, "factors")[-1, , drop=FALSE]
    Terms <- delete.response(Terms)
    vars <- all.vars(Terms)
    xl <- object$xlevels
    if(!length(xl)) {			# no factors in model
	return(as.list(coef(object)))
    }
    nxl <- setNames(rep.int(1, length(vars)), vars)
    tmp <- unlist(lapply(xl, length)) ## ?? vapply(xl, length, 1L)
    nxl[names(tmp)] <- tmp
    lterms <- apply(facs, 2L, function(x) prod(nxl[x > 0]))
    nl <- sum(lterms)
    args <- setNames(vector("list", length(vars)), vars)
    for(i in vars)
	args[[i]] <- if(nxl[[i]] == 1) rep.int(1, nl)
	else factor(rep.int(xl[[i]][1L], nl), levels = xl[[i]])
    dummy <- do.call("data.frame", args)
    pos <- 0
    rn <- rep.int(tl, lterms)
    rnn <- rep.int("", nl)
    for(j in tl) {
	i <- vars[facs[, j] > 0]
	ifac <- i[nxl[i] > 1]
	if(length(ifac) == 0L) {        # quantitative factor
	    rnn[pos+1] <- j
	} else if(length(ifac) == 1L) {	# main effect
	    dummy[ pos+1L:lterms[j], ifac ] <- xl[[ifac]]
	    rnn[ pos+1L:lterms[j] ] <- as.character(xl[[ifac]])
	} else {			# interaction
	    tmp <- expand.grid(xl[ifac])
	    dummy[ pos+1L:lterms[j], ifac ] <- tmp
	    rnn[ pos+1L:lterms[j] ] <-
		apply(as.matrix(tmp), 1L, function(x) paste(x, collapse=":"))
	}
	pos <- pos + lterms[j]
    }
    ## some terms like poly(x,1) will give problems here, so allow
    ## NaNs and set to NA afterwards.
    mf <- model.frame(Terms, dummy, na.action=function(x)x, xlev=xl)
    mm <- model.matrix(Terms, mf, object$contrasts, xl)
    if(any(is.na(mm))) {
        warning("some terms will have NAs due to the limits of the method")
        mm[is.na(mm)] <- NA
    }
    coef <- object$coefficients
    if(!use.na) coef[is.na(coef)] <- 0
    asgn <- attr(mm,"assign")
    res <- setNames(vector("list", length(tl)), tl)
    for(j in seq_along(tl)) {
	keep <- asgn == j
	ij <- rn == tl[j]
	res[[j]] <-
	    setNames(drop(mm[ij, keep, drop=FALSE] %*% coef[keep]), rnn[ij])
    }
    if(int > 0) {
	res <- c(list("(Intercept)" = coef[int]), res)
    }
    class(res) <- "dummy_coef"
    res
}

}# if R <= 3.1.0

## Not exported, and only used because CRAN checks must be faster
doExtras <- function() {
    interactive() || nzchar(Sys.getenv("R_robustbase_check_extra")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}

if(getRversion() < "3.3") {
    sigma <- function(object, ...) UseMethod("sigma")

    ## For completeness, and when comparing with nlrob() results:
    sigma.nls <- function(object, ...)
	## sqrt (  sum( R_i ^ 2) / (n - p) ) :
	sqrt( deviance(object) / (nobs(object) - length(coef(object))) )
}


## shortcut -- used often in print() etc:
pasteK <- function(...) paste(..., collapse = ", ")

## stopifnot(..) helper :
is.1num <- function(x) is.numeric(x) && length(x) == 1L
