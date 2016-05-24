as.function.formula <- function(x, ...) {
	vars <- setdiff(all.vars(x[[2]]), c("letters", "LETTERS", "pi"))
	dotdot <- grepl("^[.][.][1-9.]$", vars)
	if (any(dotdot)) vars <- c(setdiff(vars[!dotdot], "..."), "...")
	if ("&" %in% vars) vars <- c("&", setdiff(vars, c("...", "&")), "...")
	if (length(vars) == 0) { 
		f0 <- function() {}
		body(f0) <- x[[length(x)]]
		environment(f0) <- environment(x)
		f0
	} else {
		f <- function(x) {}
		formals(f) <- rep(as.list(formals(f)), length(vars))
		names(formals(f)) <- vars
		body(f) <- x[[length(x)]]
		environment(f) <- environment(x)
		f
	}
}

match.funfn <- function(FUN, descend = TRUE) {
    if (is.function(FUN)) 
        return(FUN)
	if (inherits(FUN, "formula"))
		return(as.function.formula(FUN))
    if (!(is.character(FUN) && length(FUN) == 1 || is.symbol(FUN))) {
        FUN <- eval.parent(substitute(substitute(FUN)))
        if (!is.symbol(FUN)) 
            stop(gettextf("'%s' is not a function, character or symbol", 
                deparse(FUN)), domain = NA)
    }
    envir <- parent.frame(2)
    if (descend) 
        FUN <- get(as.character(FUN), mode = "function", envir = envir)
    else {
        FUN <- get(as.character(FUN), mode = "any", envir = envir)
        if (!is.function(FUN)) 
            stop(gettextf("found non-function '%s'", FUN), domain = NA)
    }
    return(FUN)
}

fn <- structure(NA, class = "fn")
"$.fn" <- function(x, FUN) {
	env <- parent.frame()
	mf <- match.fun(FUN)

	function(...) {
		args <- list(...)
		mc <- if (is.primitive(mf)) match.call()
		else match.call(mf)
		mc1 <- mc[-1]
		nm <- names(mc1)
		if (is.null(nm)) nm <- rep("", length(args))

		mcList <- as.list(mc1)
		p <- parent.frame()
		mcListE <- lapply(mcList, eval, p)

		# if simplify found set it and remove it from lists
		simplify <- NULL
		idx <- match("simplify", tolower(nm), nomatch = 0)
		if (idx > 0) {
			if (!is.logical(mcListE[[idx]])) {
				simplify <- mcListE[[idx]]
				mcListE <- mcListE[-idx]
				mcList <- mcList[-idx]
				nm <- nm[-idx]
			}
		}

		# is.fo is a logical vector indicating whether
		#    each list element is or is not a formula
		# is.fo2 is a logical vector indicating whether each
		#    list element has or does not have a ~~ (double ~)

		is.fo <- sapply(mcListE, function(x) inherits(x, "formula"))
		any.fo <- any(is.fo)
		is.fo2 <- sapply(mcListE, function(x) inherits(x, "formula") &&
			length(x[[length(x)]]) > 1 &&
			identical(x[[length(x)]][[1]], as.name("~")))
		# change ~~ to ~
		any.fo2 <- any(is.fo2)
		if (any.fo2)
		   for(i in seq(along = mcListE))
			if (is.fo2[i]) {
			   len <- length(mcListE[[i]])
			   mcListE[[i]][[len]] <- mcListE[[i]][[len]][[2]]
			   mcListE[[i]] <- as.function(mcListE[[i]])
			} 
					
		is.char <- sapply(mcListE, is.character)
		any.char <- any(is.char)
		is.chara <- sapply(mcListE, function(x) 
			is.character(x) && substring(x, 1, 1) == "\1")
		# remove leading \1 on character strings
		any.chara <- any(is.chara)
		if (any.chara)
		   for(i in seq(along = mcListE))
		      if (is.chara[i])
			mcListE[[i]] <- gsubfn(x = substring(mcListE[[i]], 2), env = p)

		# if no ~~ formulas and no \1 strings use default strategy
		# of converting all formulas to functions and if no formulas
		# performing perl-style interpolation on all strings
		if (!any.fo2 && !any.chara) {
		   if (any.fo) {
		      for(i in seq(along = mcListE))
		         if (is.fo[i])
			    mcListE[[i]] <- as.function(mcListE[[i]])
		   } else {
		      if (any.char)
		         for(i in seq(along = mcListE))
		            if (is.char[i])
			       mcListE[[i]] <- gsubfn(x = mcListE[[i]], env = p)
		   }
		}
			
		# out <- do.call(FUN, args)
		# out <- withVisible(FUN, mcListE, env=p))
		out <- withVisible(do.call(FUN, mcListE, envir=p))
		vis <- out$visible
		out <- out $value
		if (!is.null(simplify)) {
			if(!is.list(out)) out <- list(out) 
			out <- withVisible(do.call(simplify, out))
			vis <- out$visible
			out <- out$value
		}
		if (vis) out else invisible(out)
	}
}


matrixfn <- function (data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL) {
	if (!is.function(data) & !inherits(data, "formula")) return(matrix(data = data, nrow = nrow, ncol = ncol, byrow = byrow, dimnames = dimnames))

	data. <- match.funfn(data)

	out <- matrix(nrow = nrow, ncol = ncol, dimnames = dimnames)
	for(i in seq_len(nrow)) for(j in seq_len(ncol)) out[i, j] <- data.(i, j)
	out
}

# test
# fn$list(x ~ 2*x)
# fn$mapply(~ x + y, 1:10, 21:30)

cat0 <- function(..., sep = "") cat(..., sep = sep)
# paste0 <- function(..., sep = "") paste(..., sep = sep)




