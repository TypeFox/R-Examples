
anova.dyn <- function(object, ...) {

	series <- function(x) attr(model.frame(x), "series")
	args <- list(object, ...)

	# eliminate non-tsreg objects
	L <- lapply(args, series)
	L <- L[!sapply(L, is.null)]

	# subset= is not supported currently
	stopifnot(all(sapply(L, function(x) is.null(x$call$subset))))

	L <- lapply(L, function(x) na.omit(do.call("merge.zoo", x)))
	common.times <- time(na.omit(do.call("merge.zoo", L)))

	# for each tsreg arg, set subset= in its call to common times
	for(i in seq(along = args))
	   if (!is.null(series(args[[i]]))) {
	      times. <- time(do.call("merge.zoo", series(args[[i]])))
	      subset. <- !is.na(MATCH(times., common.times))
	      args[[i]] <- eval(substitute(update(args[[i]], subset = subset.), 
				list(subset. = subset.)))
	    }

	do.call(anova, args)
}

