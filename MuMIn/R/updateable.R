`updateable` <-
function (FUN, eval.args = NULL, Class) {
    FUN <- match.fun(FUN)
	env <- environment(FUN)
	if(!isNamespace(env) && all(c("Class", "eval.args", "FUN", "FUNV") %in% names(env))) {
		warning("it looks that 'FUN' is already a result of 'updateable'. Using the original function instead")
		FUN <- env[['FUN']]
	}
	rm(env, inherits = FALSE)

    FUNV <- function() {
		rval <- do.call(FUN, as.list(match.call())[-1L])
		cl <- match.call()
		parentframe <- parent.frame()
        for (i in eval.args) cl[[i]] <- eval(cl[[i]], parentframe)
        if (!isS4(rval) && is.list(rval)) 
            rval$call <- cl
			else attr(rval, "call") <- cl
        class(rval) <- Class
        rval
    }
    body(FUNV)[c(if (missing(eval.args)) c(4L, 5L), if (missing(Class)) 7L)] <- NULL
    formals(FUNV) <- formals(FUN)
    FUNV
}

`updateable2` <-
function (FUN, Class) .Defunct("updateable")


`get_call` <- function(x) {
	rval <-
	if(isS4(x)) {
		if(any(i <- (sln <- c("call", "CALL", "Call")) %in% methods::slotNames(x)))
			slot(x, sln[i][1L]) else
			if(!is.null(attr(x, "call")))
				attr(x, "call") else NULL
	} else {
		if(!is.atomic(x) && (i <- match("call", names(x), nomatch = 0L)) != 0L) {
			x[[i]]
		} else if(!is.null(attr(x, "call"))) {
			attr(x, "call")
		} else
			NULL
	}
	if(is.null(rval)) stats::getCall(x) else rval
}

##==============================================================================

`uGamm` <-
 function(formula, random = NULL, ..., lme4 = inherits(random, "formula")) {
	pkg <- if(lme4) "gamm4" else "mgcv"
	if (!require(pkg, character.only = TRUE)) stop("cannot load package ", sQuote(pkg))
	funcl <- call("get", if(lme4) "gamm4" else "gamm", ns <- asNamespace(pkg))
	clx <- cl <- match.call()
	clx$lme4 <- NULL	
	clx <- match.call(clx, definition = eval(funcl, envir = ns))
	clx[[1L]] <- funcl
	res <- eval.parent(clx)
	res$call <- cl
	class(res) <- c(if(lme4) "gamm4", "gamm")
	res
}


`gamm` <- 
function(...) .Deprecated("uGamm", old = "MuMIn::gamm")
