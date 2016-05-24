# Created by David Bickel on 18 July 2007.

qvaluen <- function(object, ...){printGeneric("qvaluen")}#browser
#removeMethods("qvaluen")
setMethod("qvaluen", signature(object = "numeric"), function(object, ...)
{
	if(any(object < 0 | object > 1, na.rm = TRUE))
		stop("bad argument in qvaluen")
	qvaluen(Numeric(object), ...)
})
setMethod("qvaluen", signature(object = "Numeric"), function(object, threshold, verbose, ...)
{
	if(missing(verbose))
		verbose <- FALSE
	if(missing(threshold))
		threshold <- numeric(0)
	stopifnot(is.numeric(threshold) && (length(threshold) %in% c(0, 1) || all(threshold >= 0 & threshold <= 1)))
	stopifnot(all(object <= 1, na.rm = TRUE))
	boo <- is.finite(object)
	if(all(boo))
	{
		ra <- rank(object, ...)
		reject <- function(q) # BH 1995
		{
			object <= ra * q / length(object)
		}
		qval <- rep(1, length(object))
		if(length(threshold) %in% c(0, 1))
		{
			get.threshold <- function(length.out)
			{
				if(length(length.out) == 0)
				{
					length.out <- 1000
					message("qvaluen ", class(object), " threshold = ", length.out)
					get.threshold(length.out = length.out)
				}
				else if(length.out >= 2 && floor(length.out) == length.out)
				{
					z <- zvalue(object)
					z <- z[is.finite(z)]
					z.threshold <- seq(from = min(z), to = max(z), length.out = length.out)
					pnorm(z.threshold)
				}
				else
					stop("bad threshold argument")
			}
			threshold <- get.threshold(length.out = threshold)
		}
		if(any(threshold < 0 | threshold > 1))
			stop("bad threshold")
		uq <- unique(sort(threshold, decreasing = TRUE))
		if(verbose)
			cat("  performing BH 1995 for q = ")
		for(q in uq)
		{
			if(verbose)
				cat(q, " ")
			qval <- ifelse(reject(q = q), q, qval) # lowest FDR at which rejected
		}
		if(verbose)
			cat("\n")
		stopifnot(length(qval) == length(object) && all(is.finite(qval)))
		stopifnot(all(qval >= object))
		names(qval) <- names(object)
		Numeric(qval)
	}	
	else if(any(boo))
	{
		object[boo] <- qvaluen(object = object[boo], threshold = threshold, verbose = verbose, ...)
		object
	}
	else
	{
		warning("cannot compute q-value")
		object
	}
})
setMethod("qvaluen", signature(object = "biasEstimate"), function(object, alternative, correct, ...)
{
	if(missing(alternative))
	{
		alternative <- "two.sided"
		message("qvaluen ", class(object), " alternative = ", alternative)
	}
	if(missing(correct))
	{
		correct <- FALSE
		message("qvaluen ", class(object), " correct = ", correct)
	}
	pval <- pvalue(object = object, alternative = alternative, correct = correct, ...)
	qvaluen(object = pval)
})

change.pvalue <- function(object, old.alternative, new.alternative) # added 7 May 2008; see also "pvalue", signature(object = "eFdr"); low-level
{
	stopifnot(is(object, "numeric") && all(object >= 0 & object <= 1, na.rm = TRUE))
	if(old.alternative == "less")
	{
		if(new.alternative == "less")
			object
		else if(new.alternative == "greater")
		{
			warning('consider using change.zvalue based on "pvalue", signature(object = "eFdr") for less p-value roundoff error')
			1 - object
		}
		else if(new.alternative == "two.sided")
			2 * pmin(object, 1 - object)
		else
			stop("new.alternative not recognized")
	}
	else
		stop(paste("change.pvalue not implemented for old.alternative ", old.alternative, sep = ""))
}
change.zvalue <- function(object, old.alternative, new.alternative) # added 9 May 2008; see also "pvalue", signature(object = "eFdr")
{
	stopifnot(is(object, "numeric"))
	if(old.alternative == "less")
	{
		if(new.alternative == "less")
			object
		else if(new.alternative == "greater")
			-object
		else if(new.alternative == "two.sided")
			pmin(object, -object)
		else
			stop("new.alternative not recognized")
	}
	else
		stop(paste("change.zvalue not implemented for old.alternative ", old.alternative, sep = ""))
}
pvalue.via.zvalue <- function(object, alternative, ...) # low-level
{
	if(missing(alternative))
		alternative <- "less"
	q <- zvalue(object = object, alternative = alternative, ...) #type = "pvalue", ...)
	p <- pnorm(q = q)
	if(alternative == "two.sided")
		2 * p
	else
		p
}

pvalue <- function(object, ...){printGeneric("pvalue")}#browser # in genetics.r prior to 12 February 2008
#removeMethods("pvalue")
setMethod("pvalue", signature(object = "numeric"), function(object)
{
	if(is(object, "Numeric"))
		object
	else
		Numeric(object)
})
setMethod("pvalue", signature(object = "list"), function(object)
{
	pvalue(object$p.value)
})
setMethod("pvalue", signature(object = "Fdr"), function(object, alternative)
{
	old.alternative <- "less"
	if(missing(alternative))
		alternative <- old.alternative
	change.pvalue(object = object@p.value, old.alternative = old.alternative, new.alternative = alternative)
})
setMethod("pvalue", signature(object = "eFdr"), function(object, alternative, rounded)
{
	if(missing(alternative))
		alternative <- "less"
	if(missing(rounded))
		rounded <- FALSE
	p <- if(rounded) # sometimes rounds off
	{
		pvalue(as(object, "basicFdr"), alternative = alternative) # prior to 25 March 2008
	}
	else # best for volcano plots
	{
		pvalue.via.zvalue(object = object, alternative = alternative, type = "pvalue")
	}
	pvalue(p)
})
setMethod("pvalue", signature(object = "biasEstimate"), function(object, ...)
{
	pvalue.via.zvalue(object = object, ...)
})

zvalue <- function(object, type, ...){printGeneric("zvalue")}#browser # new 24 April 2008
#removeMethods("zvalue")
setMethod("zvalue", signature(object = "Fdr", type = "missing"), function(object, type, ...)
{
	type <- "pvalue" # [sic]
	message("zvalue type = ", type)
	zvalue(object = object, type = type, ...)
})
setMethod("zvalue", signature(object = "eFdr", type = "character"), function(object, type, alternative, ...)
{
	stopifnot(length(type) == 1)
	if(type == "pvalue") # moved from "pvalue", signature(object = "eFdr")
	{
		old.alternative <- "less"
		if(missing(alternative))
			alternative <- old.alternative
		change.zvalue(object = object@zz, old.alternative = old.alternative, new.alternative = alternative)
	}
	else if(missing(alternative))
		zvalue(object = as(object, "basicFdr"), type = type)
	else
		stop(paste("specification of alternative incompatible with type = ", type, sep = ""))
})
setMethod("zvalue", signature(object = "Numeric", type = "missing"), function(object, type)
{
	qnorm(as(object, "numeric"))
})
setMethod("zvalue", signature(object = "Fdr", type = "character"), function(object, type, ...)
{
	stopifnot(length(type) == 1)
	Num <- do.call(type, list(object, ...))
	p <- if(type == "fdr")
	{
		pval <- pvalue(object) # assumed to come from a one-sided test
		lfdr <- Num
		stopifnot(any(pval != lfdr, na.rm = TRUE))
		lower <- lfdr / 2
		nonsig <- 1 / 2 # one-sided p-value or lfdr.p always considered nonsignificant
		higher <- 1 - lfdr / 2 
		lfdr.p <- ifelse(pval < nonsig, lower, ifelse(pval > nonsig, higher, nonsig))
		names(lfdr.p) <- names(Num)
		Numeric(lfdr.p)
	}
	else
		Num
	stopifnot(is(p, "Numeric"))
	tolerance <- 0
	if(any(p < 0 - tolerance | p > 1 + tolerance, na.rm = TRUE))
	{ message("p for zvalue out of range")}#; browser()
	zvalue(object = p) # [sic]
})
setMethod("zvalue", signature(object = "biasEstimate", type = "missing"), function(object, type, ...)
{
	uncorrected <- object@uncorrected
	class.ok <- is(uncorrected, "unsortableEstimate")
	if(!class.ok)
		stop(paste("qvaluen ", class(object), " not implemented for uncorrected of class ", class(uncorrected), sep = ""))
	type <- uncorrected@estimator@type
	zvalue(object = object, type = type, ...)
})
setMethod("zvalue", signature(object = "biasEstimate", type = "character"), function(object, type, correct, alternative, ...)
{
	stopifnot(is.logical(correct))
	uncorrected <- object@uncorrected
	z <- if(type %in% c("pvalue.z", "fdr.z"))
	{
		if(correct)
			uncorrected - object # error fixed 9 May 2008 4 pm
		else
			uncorrected
	}
	else
		stop(paste("qvaluen ", class(object), " not implemented for type ", type, sep = ""))
	old.alternative <- "less"
	if(missing(alternative))
	{
		alternative <- old.alternative
		message("zvalue ", class(object), " alternative = ", alternative)
	}
#	p.less <- pnorm(z)
#	change.pvalue(object = p.less, old.alternative = old.alternative, new.alternative = alternative)
	change.zvalue(object = z, old.alternative = old.alternative, new.alternative = alternative)
})

prior0 <- function(object, ...){printGeneric("prior0")}#browser
#removeMethods("prior0")
setMethod("prior0", signature(object = "numeric"), function(object, ...)
{
	priorFDR(p.values = object, ...)
})
setMethod("prior0", signature(object = "eFdr"), function(object, ...)
{
	p0 <- as(object, "list")[[1]]$fp0[1, "p0"]
	stopifnot(is.numeric(p0) && length(p0) == 1)
	p0
})
setMethod("prior0", signature(object = "bFdr"), function(object, ...)
{
	object@prior.fdr
})

fdr <- function(object, name) # extracts numeric fdr or prob0 from object
{
	if(missing(name)) name <- "fdr"
	f <- Slot(object = object, name = name)
	stopifnot(is.numeric(f))
	f
}
#removeMethods("fdr")
setMethod("fdr", signature(object = "Prob0", name = "missing"), function(object, name)
{
	as(object, "Numeric")
})

new.basicFdr <- function(fdr, p.value)
{
	fnam <- names(fdr)
	pnam <- names(p.value)
	if(is.character(fnam) && is.character(pnam))
		stopifnot(all(fnam == pnam))
	else if(is.null(fnam))
		names(fdr) <- names(p.value) # new 24 April 2008
	new("basicFdr", fdr = Numeric(fdr), p.value = Numeric(p.value))
}

eFdr <- function(object, ...){message("generic eFdr")}#browser # e = Efron
#removeMethods("eFdr")
setMethod("eFdr", signature(object = "numeric"), function(object, max.abs.z, verbose, save.time, trivial, ...)
{
	if(missing(save.time))
		save.time <- TRUE
	if(missing(trivial))
		trivial <- FALSE
	time0 <- Sys.time()
	if(!trivial)
		message("Beginning main eFdr computations on ", date(), ".")
	if(missing(verbose))
	  verbose <- FALSE
	if(missing(max.abs.z))
	{
		max.abs.z <- max(abs(object[is.finite(object)]), na.rm = TRUE)
		if(verbose) message("max.abs.z = ", max.abs.z)
	}
	stopifnot(length(max.abs.z) == 1 && max.abs.z >= 0)
	nam <- names(object)
	stopifnot(is.character(names(object)))
	object <- ifelse(abs(object) > max.abs.z, sign(object) * max.abs.z, object)
	names(object) <- nam
	if(all(object >= 0, na.rm = TRUE) || all(object <= 0, na.rm = TRUE))
	{ message("no element of zz has different direction")}#; browser()
	boo <- !is.na(object)
	zz <- object[boo]
	get.blank.s3 <- function(){list()}
	s3 <- if(trivial)
		get.blank.s3()
	else
		try(mylocfdr(zz = zz, ...))
	add.NA <- function(x)
	{
		x2 <- sapply(nam, function(na)
		{
			if(na %in% names(x))
				x[na]
			else
				as.numeric(NA)
		})
		if(length(x2) == length(nam))
			names(x2) <- nam
		else
			stop("bad x2")
		x2
	}
	addNAs <- function(x)
	{
		x2 <- as.numeric(rep(NA, length(object)))
		stopifnot(sum(boo) == length(x))
		x2[boo] <- x
		if(length(x2) == length(nam))
			names(x2) <- nam
		else
			stop("bad x2")
		x2 # error corrected 27 Feb. 2008
	}
	add.na <- function(x){if(save.time) addNAs(x) else add.NA(x)}
	get.fdr <- function(vec)
	{
		names(vec) <- names(zz)
		vec2 <- add.na(vec)
		if(verbose) message("New fdr slot has ", length(vec2), " values, ", length(vec), " of which were returned by locfdr.")
		Numeric(vec2) # Slot(from, "fdr")
	}
	fdrNum <- if(is(s3, "try-error") || trivial) #	{ message("bad s3 in eFdr")}#browser
	{
		if(!trivial)
			warning("bad s3 in eFdr")
		s3 <- get.blank.s3()
		get.fdr(vec = rep(1, length(zz)))
	}
	else
	{
		get.fdr(vec = s3$fdr)
	}
	if(verbose)
		message(length(fdrNum), " local false discovery rates computed in ", (Sys.time() - time0) / 60, " minutes, ending on ", date(), ".")
	if(length(fdrNum) != length(object))
	{ message("fdrNum of bad length")}#; browser()
	new("eFdr", list(s3 = s3), fdr = fdrNum, zz = add.na(zz))
})
setMethod("eFdr", signature(object = "Numeric"), function(object, direction, ...)
{
  if(!all(object <= 1, na.rm = TRUE))
  { message("not all p-values are sufficiently small")} #; browser()
  if(is.null(names(object)))
  { message("bad names(object) in eFdr")} #; browser()
  if(!missing(direction))
 	{
 		stopifnot(is.numeric(direction))
	  stopifnot(all(names(object) == names(direction)))
 		p1 <- object / 2
 		nam <- names(object)
 		object <- ifelse(is.na(direction), direction, ifelse(direction > 0, p1, 1 - p1)) # one-sided p-value
 		if(!(is.null(nam) || length(object) == length(nam)))
 		{ message("direction killed length")} # browser()
 		names(object) <- nam
 	}
 	zz <- zvalue(object = object) # qnorm(as(object, "numeric"))
 	names(zz) <- names(object)
 	present.zz <- zz[!is.na(zz)]
 	if(is(zz, "Numeric") || all(present.zz[1] == present.zz))
 	{
 		message("bad zz")
 		#browser()
 	}
	if(all(zz >= 0, na.rm = TRUE) || all(zz <= 0, na.rm = TRUE))
	{
		message("no element of zz has different direction; press c to continue anyway")
		#browser()
	}
	eFdr(object = zz, ...)
})
setAs(from = "eFdr", to = "Numeric", function(from)
{
	from@fdr
})
setAs(from = "eFdr", to = "numeric", function(from)
{
	as(as(from, "Numeric"), "numeric")
})
eFdr0 <- function(object, ...)
{
	eFdr(object = object, nulltype = 0, ...)
}
trivialFdr <- function(object, ...) # for p-value computation without local FDR
{
	eFdr(object = object, trivial = TRUE, ...)
}
is.eFdr.FUN <- function(object) # new 24 April 2008
{
	stopifnot(is.function(object))
	funs <- list(eFdr, eFdr0, trivialFdr)
	stopifnot(all(sapply(funs, is.function)))
	any(sapply(funs, function(fun){identical(fun, object)}))
}

# begin added 15 Feb. 2008
bFdr <- function(object, ...){message("generic bFdr")}#browser # b = Bickel
#removeMethods("bFdr")
setMethod("bFdr",  signature(object = "numeric"), function(object, prior.fdr, verbose, ...)
{
	if(missing(verbose))
	  verbose <- FALSE
	if(verbose)
		message("\n\nBeginning main bFdr computations on ", date(), ".\n")
	p.value <- Numeric(object)
	if(missing(prior.fdr))
	{
		if(verbose)
			message("Calling priorFDR on ", date(), ".")
		prior.fdr <- priorFDR(p.values = p.value, verbose = verbose)
		message("prior.fdr = ", prior.fdr)
	}
	fdrNum <- Numeric(localFDR(p.values = p.value, prior.fdr = prior.fdr, verbose = verbose, ...))
	if(verbose)
		message("\n", length(fdrNum), " local false discovery rates ending on ", date(), ".")
	new("bFdr", fdr = fdrNum, p.value = p.value, prior.fdr = prior.fdr)
})
setMethod("bFdr", signature(object = "Numeric"), function(object, direction, ...)
{
  if(!all(object <= 1, na.rm = TRUE))
  { message("not all p-values are sufficiently small")}#; browser()
  if(is.null(names(object)))
  { message("bad names(object) in bFdr")} #; browser()
  if(!missing(direction))
 	{
 		stopifnot(is.numeric(direction))
	  stopifnot(all(names(object) == names(direction)))
 		p1 <- object / 2
 		object <- ifelse(direction > 0, p1, 1 - p1) # one-sided p-value
 	}
	bFdr(object = as(object, "numeric"), ...)
})
setAs(from = "bFdr", to = "Numeric", function(from)
{
	from@fdr
})
setAs(from = "bFdr", to = "numeric", function(from)
{
	as(as(from, "Numeric"), "numeric")
})
# end added 15 Feb. 2008

pseudoFdr.FUN.default <- function(x)
{
	dexp(exp(abs(x)) - 1)
}
get.pseudoFdr.FUN <- function(foldchange.threshold)
{
	stopifnot(is.numeric(foldchange.threshold))
	if(length(foldchange.threshold) == 0)
		pseudoFdr.FUN.default
	else
	{
		stopifnot(length(foldchange.threshold) == 1)
		function(x)
		{
			foldchange <- exp(abs(x))
			ifelse(foldchange >= foldchange.threshold, 0, 1)
		}
	}
}
new.pseudoFdr <- function(estimate, FUN, foldchange.threshold)
{
	if(missing(foldchange.threshold))
	{
		foldchange.threshold <- numeric(0)
		if(missing(FUN))
			message("new.pseudoFdr foldchange.threshold = ", foldchange.threshold)
	}
	if(missing(FUN))
	{
		FUN <- get.pseudoFdr.FUN(foldchange.threshold = foldchange.threshold)
	}
	else
		stopifnot(length(foldchange.threshold) == 0)
	stopifnot(is.function(FUN))
	fdr <- FUN(estimate)
	if(any(fdr < 0 | fdr > 1, na.rm = TRUE))
	{	message("bad pseudo fdr")} #; browser()
	if(!is(fdr, "Numeric"))
		fdr <- Numeric(fdr)
	new("pseudoFdr", fdr = fdr, estimate = estimate, FUN = FUN)
}
pseudoFdr <- function(FUN, estimator, foldchange.threshold, ...)
{
	if(missing(foldchange.threshold))
	{
		foldchange.threshold <- numeric(0)
		message("pseudoFdr foldchange.threshold = ", foldchange.threshold)
	}
	if(missing(FUN))
	{
		FUN <- get.pseudoFdr.FUN(foldchange.threshold = foldchange.threshold)
	}
	else
		stopifnot(length(foldchange.threshold) == 0)
	if(missing(estimator))
		estimator <- meanHat
	stopifnot(is(estimator, "Estimator"))
	estimate <- unsortableEstimate(estimator = estimator, ...)
	new.pseudoFdr(estimate = estimate, FUN = FUN, foldchange.threshold = numeric(0))
}
setMethod("p.value", signature(object = "pseudoFdr"), function(object, ...)
{
	vec <- rep(1, length(object))
	names(vec) <- names(object)
	Numeric(vec)
})

colnames.design.default <- "Grp2vs1"
sFdr <- function(x, y, ...){message("generic sFdr")}#browser # s = Smyth
#removeMethods("sFdr")
setMethod("sFdr",  signature(x = "Numeric", y = "missing"), function(x, y, verbose, ...)
{
	stop("functon not sufficiently modified from eFdr")
	if(missing(verbose))
	  verbose <- FALSE
	nam <- names(object)
	stopifnot(is.character(nam))
	p.value <- object[!is.na(object)]
	s3 <- try(mylocfdr(p.value = p.value, ...))#XXX:locfdr does not have p.value as input!
	add.NA <- function(x)
	{
		x2 <- sapply(nam, function(na)
		{
			if(na %in% names(x))
				x[na]
			else
				as.numeric(NA)
		})
		if(length(x2) == length(nam))
			names(x2) <- nam
		else
			stop("bad x2")
		x2
	}
	get.fdr <- function(vec)
	{
		names(vec) <- names(p.value)
		vec2 <- add.NA(vec)
		if(verbose) message("New fdr slot has ", length(vec2), " values, ", length(vec), " of which were returned by locfdr.")
		Numeric(vec2) # Slot(from, "fdr")
	}
	fdrNum <- if(is(s3, "try-error")) #	{ message("bad s3 in sFdr")}#browser
	{
		warning("bad s3 in sFdr")
		s3 <- list()
		get.fdr(vec = rep(1, length(p.value)))
	}
	else
	{
		get.fdr(vec = s3$fdr)
	}
#	new("sFdr", list(s3 = s3), fdr = fdrNum, p.value = add.NA(p.value))
})
setMethod("sFdr", signature(x = "numeric", y = "missing"), function(x, y, ...)
{
	sFdr(x = Numeric(x), ...)
})
setMethod("sFdr", signature(x = "xprnSet", y = "missing"), function(x, y, design, ...)
{
	if(missing(design))
		design <- NULL
	get.sfdr <- function(es)
	{
		stopifnot(is(es, "xprnSet"))
		fit <- try(lmFit(object = exprs(es), design = design))
		if(is(fit, "try-error"))
		{ message("error in lmFit")} #; browser()
		warn <- getOption("warn")
		options(warn = 2)
		sfdr <- try(sFdr(x = fit, ...))
		options(warn = warn)
		if(is(sfdr, "try-error"))
			message('error in "sFdr", signature(x = "MArrayLM", y = "missing")')			
		sfdr
	}
	small.es <- removeMissing(x)
	ann <- annotation(x)
	if(!is.character(ann) || length(ann) != 1 || ann == '')
		ann <- class(x)
	nfeatures <- length(featureNames(x))
	nremoved <- nfeatures - length(featureNames(small.es))
	if(nremoved > 0)
		message("    Removed ", nremoved, " features from ", ann, ".")
	sfdr <- get.sfdr(es = small.es)
	if(is(sfdr, "try-error"))
	{
		nsamples <- ncol(exprs(x))
		vec <- as.numeric(rep(NA, nfeatures))
		names(vec) <- featureNames(x)
		Num <- Numeric(vec)
		sfdr <- new("sFdr", list(s3 = list()), fdr = Num, p.value = Num, t = Num, proportion = Scalar(NA))
		stopifnot(validObject(sfdr))
		mess <- paste("    get.sfdr generated an error for ", ann, " of ", nfeatures, " features and ", nsamples, " samples; returning blank sFdr object of length ", length(sfdr), ".", sep = "")
		warning(mess)
		message(mess)
	}
#	{ message("error in get.sfdr; to debug, try something like class(get.sfdr(x[1:50, ]))")}#browser
	nmissing <- nfeatures - length(sfdr)
	if(nmissing == 0)
		sfdr
	else
	{
		message("    Adding ", nmissing, " features back to sfdr.")
		sfdr <- sFdr(x = sfdr, y = x)
		stopifnot(nfeatures == length(sfdr))
		sfdr
	}
})
setMethod("sFdr", signature(x = "xprnSet", y = "xprnSet"), function(x, y, ...)
{
	if(is(x, "XprnSet") || is(y, "XprnSet")) stop("x and y were transformed differently")
	sFdr(x = xprnSetPair(x = x, y = y), ...)
})
setMethod("sFdr", signature(x = "XprnSet", y = "missing"), function(x, y, ...)
{
	sFdr(x = logb(x), ...)
})
setMethod("sFdr", signature(x = "XprnSet", y = "XprnSet"), function(x, y, ...)
{
	sFdr(x = logb(x), y = logb(y), ...)
})
setMethod("sFdr", signature(x = "xprnSetPair", y = "missing"), function(x, y, design, ...)
{
	if(missing(design))
	{
		half.col <- function(object, binary)
		{
			rep(binary, ncol(exprs(object)))
		}
		design <- cbind(Grp1 = 1, Grp2vs1 = c(half.col(x@x, 0), half.col(x@y, 1)))
		colnames(design) <- c("Grp1", colnames.design.default)
	}
	big.es <- as(x, "xprnSet")
	small.es <- removeMissing(big.es)
	nmissing <- length(featureNames(big.es)) - length(featureNames(small.es))
	if(nmissing > 0)
		message(" Subtracted ", nmissing, " missing features from ", class(x), ".")
	small.sfdr <- sFdr(x = small.es, design = design, ...)
	big.sfdr <- if(nmissing > 0)
	{
		message(" Adding ", nmissing, " blank features back to ", class(small.sfdr), ".")
		sFdr(x = small.sfdr, y = big.es)
	}
	else
		small.sfdr
	if(length(big.sfdr) != length(featureNames(big.es)))
	{ message("big.sfdr is too small")} #; browser()
	big.sfdr
})
setMethod("sFdr", signature(x = "sFdr", y = "xprnSet"), function(x, y, ...)
{
	nam <- featureNames(y)
	sfdr <- sFdr(x = x, y = nam, ...)
	if(length(sfdr) != length(nam))
	{ message('error in "sFdr", signature(x = "sFdr", y = "xprnSet")')} #; browser()
	sfdr
})
setMethod("sFdr", signature(x = "sFdr", y = "character"), function(x, y)
{
	stopifnot(is.character(names(x)) && all(names(x) %in% y))
	nam.x <- y[y %in% names(x)]
	if(length(nam.x) != length(x))
	{
		message("problem with nam.x")
		#browser()
	}
	x <- x[nam.x]
	sfdr <- x
	add.na <- function(vec)
	{
		if(is.null(names(vec)))
			names(vec) <- names(x)
		get.present <- function(object){object[!is.na(object)]}
		present.vec <- get.present(vec)
		new.vec <- as.numeric(rep(NA, length(y)))
		names(new.vec) <- y
		new.vec[nam.x] <- vec
		present.new.vec <- get.present(new.vec)
		len.ok <- length(new.vec) == length(y)
		other.ok <- all(names(new.vec) == y) && all(present.new.vec == present.vec) && all(names(present.new.vec) == names(present.vec))
		ok <- len.ok && other.ok
		if(!ok)
		{ message("add.na failed")} #; browser()
		new.vec
	}
	sfdr@fdr <- Numeric(add.na(x@fdr))
	sfdr@p.value <- Numeric(add.na(x@p.value))
	sfdr@t <- add.na(x@t)
	if(!validObject(sfdr))
	{ message("lengthed ", class(sfdr), " is invalid")} #; browser()
	sfdr
})
setMethod("sFdr", signature(x = "MArrayLM", y = "missing"), function(x, y, max.niter, stdev.coef.lim, proportion, verbose, j, tolerance, ...)
{
	if(missing(verbose))
	{
		verbose <- TRUE
		if(verbose)
			message("sFdr verbose = ", verbose)
	}
	if(missing(max.niter))
	{
		max.niter <- if(missing(proportion)) 1000 else 1
		if(verbose)
			message("max.niter = ", max.niter)
	}
	if(missing(tolerance))
	{
		tolerance <- 0.001
		if(verbose)
			message("tolerance = ", tolerance)
	}
	if(missing(stdev.coef.lim)) stdev.coef.lim <- c(0.1,4)
	if(missing(proportion))
	{
		proportion <- 0.1 # 1 - convest(p = p.value, max.niter = max.niter)
		if(verbose)
			message("proportion = ", proportion)
	}
	stopifnot(max.niter >= 1)
	old.proportion.change <- 0
	if(missing(j))
		j <- NULL
	for(i in 1:max.niter)
	{
		eb <- eBayes(fit = x, proportion = proportion, stdev.coef.lim = stdev.coef.lim) # try()
		if(is(eb, "try-error"))
		{ message("error in eBayes")} #; browser()
		get.vec <- function(object, nam)
		{
			if(is.matrix(object))
			{
				if(is.null(j))
				{
					if(colnames.design.default %in% colnames(object))
						j <- colnames.design.default
					else
					{
						j <- 2
						warning(paste(colnames.design.default, "is not in design, so j =", j))
					}
				}
				object <- object[, j]
			}
			if(!is.numeric(object))
			{ message(nam, " is non-numeric")}#; browser()
			object
		}
		p.value <- Numeric(get.vec(eb$p.value, nam = "p.value"))
		stopifnot(is.character(names(p.value)))
		stat <- get.vec(eb$t, nam = "stat")
		names(stat) <- names(p.value)
		odds <- exp(eb$lods)
		odds <- get.vec(odds, nam = "odds")
		names(odds) <- names(p.value)
		fdrNum <- Numeric(1 / (1 + odds))
		old.proportion <- proportion
		proportion <- 1 - mean(fdrNum, na.rm = TRUE)
		proportion.change <- proportion - old.proportion
		sign.changed <- old.proportion.change * proportion.change < 0 # change in sign
		small.change <- is.numeric(tolerance) && length(tolerance) == 1 && tolerance > 0 && abs(proportion.change) <= tolerance
		converged <- sign.changed || small.change
		if(converged)
		{
			if(verbose)
				message("  Converged; direction ", if(sign.changed) "changed" else "did not change", "; proportion.change: ", proportion.change, ".")
			break
		}
		else
			old.proportion.change <- proportion.change
	}
	if(verbose)
		message("  Estimate ", i, ": ", round(proportion, 3) * 100, "% of features differ.")
	if(!converged && max.niter > 1)
	{
		mess <- paste("  portion of differing features failed to converge after", max.niter, "iterations")
		warning(mess)
		if(verbose) message(mess)
	}
	if(length(p.value) != length(stat))
	{ message("length(p.value) != length(stat)")}#; browser()
	new("sFdr", list(s3 = x), fdr = fdrNum, p.value = p.value, t = stat, proportion = Scalar(proportion))
})
setAs(from = "MArrayLM", to = "sFdr", function(from)
{
	sFdr(x = from)
})
setAs(from = "sFdr", to = "Numeric", function(from)
{
	from@fdr
})
setAs(from = "sFdr", to = "numeric", function(from)
{
	as(as(from, "Numeric"), "numeric")
})
setAs(from = "sFdr", to = "MArrayLM", function(from)
{
	malm <- from$s3
	stopifnot(is(malm, "MArrayLM"))
	malm
})

statistic <- function(object, ...){printGeneric("statistic")}#browser
#removeMethods("statistic")
setMethod("statistic", signature(object = "sFdr"), function(object, ...)
{
#	statistic(object = as(object, "MArrayLM"), ...)
	object@t
})
setMethod("statistic", signature(object = "MArrayLM"), function(object, ...)
{
#	stop(paste("statistic method not yet implemented for", class(object)))
	stat <- object$t
	if(is.numeric(stat))
		stat
	else if(is.matrix(stat))
	{
		message(class(object), " statistic not yet implemented for matrix")
		#browser()
	}
	else
	{
		message("bad moderated statistic: class(stat) is ", class(stat), ".")
		#browser()
	}
	stat
})

Prob0 <- function(x, y, ...){printGeneric("Prob0")}#browser
#removeMethods("Prob0")
setMethod("Prob0", signature(x = "xprnSet", y = "missing"), function(x, y, probabilityNull, PValue, verbose, na.rm, mle, return.class, ...)
{
	if(missing(return.class))
	{
		return.class <- "Prob0"
		message("Prob0 (x) return.class = ", return.class)
	}
	if(missing(mle))
	{
		mle <- FALSE
		if(return.class == "Prob0")
			message("Prob0 (x) mle = ", mle)
	}
	if(missing(na.rm)) na.rm <- TRUE
	if(missing(verbose)) verbose <- FALSE
	if(missing(probabilityNull) || is.null(probabilityNull))
	{
		pn.name <- if(missing(probabilityNull))
			"eFdr"
		else
			"trivialFdr"
		probabilityNull <- eval(parse(text = pn.name)) # probabilityNullFUN(PValue = PValue)
		message("probabilityNull = ", pn.name)
	}
	nam <- featureNames(x)
	mat <- exprs(x)
	if(ncol(mat) == 1) # added 5 October 2007 for leave-1-out cross-validation
	{
		PValue <- Numeric(numeric(0))
		probabilityNull <- Numeric(rep(1, nrow(mat)))
	}
	else
	{
		if(identical(probabilityNull, sFdr) || identical(probabilityNull, pseudoFdr))
		{
			if(!missing(PValue) && !is.null(PValue))
				stop("PValue specified for incompatible probabilityNull (x)")
			probabilityNull <- probabilityNull(x = x, ...)
			PValue <- p.value(probabilityNull) # probabilityNull@p.value
		}
		else if(is.function(probabilityNull)) # then get probabilityNull from PValue
		{
			if(missing(PValue))
				PValue <- PValueFUN() # function(...){PValueFUN(FUN = t.test, ...)}
			if(identical(PValue, sFdr))
			{
				if(verbose) message("computing p-values based on moderated t statistics for non-Smyth local FDR (x)")
				sfdr <- PValue(x = x, max.niter = 1, ...)
				stopifnot(is(sfdr, "sFdr"))
				two.sided.p <- sfdr@p.value
				side.sign <- sign(statistic(sfdr)) # sign(statistic(as(sfdr, "MArrayLM")))
				PValue <- ifelse(side.sign < 0, two.sided.p / 2, 1 - two.sided.p / 2)
			} # do not forget "else" and "(x, y)"
			else if(is.function(PValue))
				PValue <- apply(mat, 1, PValue)
			else if(is.numeric(PValue))
			{
				if(length(PValue) != nrow(mat) || is.null(names(PValue)) || any(names(PValue) != nam))
				{ message("bad PValue argument")} #; browser()
			}
			if(length(PValue) != nrow(mat))
			{ message("bad PValue length")}#; browser()
			names(PValue) <- nam
			if(!is(PValue, "Numeric"))
				PValue <- Numeric(PValue)
			if(any(PValue > 1, na.rm = TRUE) || !validObject(PValue) || !is(PValue, "Numeric"))
			{ message("PValue has element too large or other problem")}#browser
			if(verbose)
			{ 
				message("PValue:")
				print(class(PValue))
				print(summary(PValue))
			}
			probabilityNull <- probabilityNull(PValue, ...)
		}
		else if(missing(PValue))
			PValue <- Numeric(numeric(0))
		else
			stop("PValue was specified even though probabilityNull is not a function")
	}
	to.return <- if(is(probabilityNull, return.class))
		probabilityNull
	else if(return.class %in% c("Prob0", "Fdr"))
	{
		if(!is.numeric(probabilityNull))
			probabilityNull <- fdr(probabilityNull)
		pn.ok <- is.numeric(probabilityNull) && nrow(mat) == length(probabilityNull) && all(names(probabilityNull) == nam)
		if(!pn.ok)
		{
			message("probabilityNull is not a valid vector of probabilities that the nulls are true, a function returning such, or an object such that the fdr method returns such.")
			#browser()
		}
		get.vec <- function(FUN, ...)
		{
			ok <- try(is.function(FUN))
			if(is(ok, "try-error") || !ok)
			{ message("bad FUN in get.vec")}#browser
			vec <- apply(mat, 1, FUN, ...)
			names(vec) <- nam
			vec
		}
		Scale0 <- Numeric(get.vec(Sd, na.rm = na.rm, mu = 0, mle = mle))
		Scale1 <- Numeric(get.vec(Sd, na.rm = na.rm, mle = mle))
		Location1 <- get.vec(mean, na.rm = na.rm)
		xSize <- Numeric(get.vec(function(object, ...){sum(!is.na(object), ...)}, na.rm = na.rm))
		ySize <- Numeric(rep(0, length(xSize)))
		suppressWarnings(names(ySize) <- names(xSize))
		if(!is(PValue, "Numeric"))
		{
			PValue <- if(is.numeric(PValue))
				Numeric(PValue)
			else
			{
				warning(class(PValue), "PValue coerced to Numeric of length 0")
				Numeric(numeric(0))
			}
		}
		new("Prob0", probabilityNull, Scale0 = Scale0, Scale1 = Scale1, Location1 = Location1, xSize = xSize, ySize = ySize, PValue = PValue, mle = mle, annotation = annotation(x))
	}
	else
	{
		message("cannot return ", return.class)
		#browser()
	}
	if(is(to.return, "Prob0") && return.class == "Fdr")
		to.return <- as(to.return, "Fdr")
	stopifnot(is(to.return, return.class))
	to.return
})
setMethod("Prob0", signature(x = "xprnSet", y = "xprnSet"), function(x, y, probabilityNull, PValue, verbose, na.rm, mle, normalization, return.class, ...)
{
	if(missing(return.class))
	{
		return.class <- "Prob0"
		message("(x, y) return.class = ", return.class)
	}
	if(missing(normalization))
	{
		normalization.name <- "normalize"
		normalization <- eval(parse(text = normalization.name))
		message("Prob0 normalization = ", normalization.name)
	}
	if(missing(verbose)) verbose <- FALSE
	if(is.function(normalization))
	{
		if(verbose) message("normalizing arrays")
		x <- normalization(x)
		y <- normalization(y)
	}
	else if(verbose)
		message("not normalizing arrays since normalization = ", normalization)
	nam <- featureNames(x)
	stopifnot(length(nam) >= 2 && all(nam == featureNames(y)))
	if(missing(mle))
	{
		mle <- FALSE
		if(return.class == "Prob0")
			message("Prob0 (x, y) mle = ", mle)
	}
	if(missing(na.rm)) na.rm <- TRUE
	if(missing(probabilityNull) || is.null(probabilityNull))
	{
		pn.name <- if(missing(probabilityNull))
			"eFdr"
		else
			"trivialFdr"
		probabilityNull <- eval(parse(text = pn.name))
		message("probabilityNull = ", pn.name)
	}
	x.mat <- exprs(x)
	y.mat <- exprs(y)
	if(ncol(x.mat) == 1 && ncol(y.mat) == 1) # added 4 October 2007 for leave-1-out cross-validation
	{
		PValue <- Numeric(numeric(0))
		probabilityNull <- Numeric(rep(1, nrow(y.mat)))
	}
	else
	{
		if(identical(probabilityNull, sFdr) || identical(probabilityNull, pseudoFdr))
		{
			if(!missing(PValue) && !is.null(PValue))
				stop("PValue specified for incompatible probabilityNull (x, y)")
			probabilityNull <- probabilityNull(x = x, y = y, ...)
			PValue <- p.value(probabilityNull) # probabilityNull@p.value
			if(length(PValue) != nrow(x.mat))
			{ message("bad (x, y) PValue from ", class(probabilityNull))}#browser
		}
		else if(is.function(probabilityNull)) # then get probabilityNull from PValue
		{
			if(missing(PValue))
				PValue <- PValueFUN()
			if(identical(PValue, sFdr))
			{
				if(verbose) message("computing p-values based on moderated t statistics for non-Smyth local FDR (x, y)")
				sfdr <- PValue(x = x, y = y, max.niter = 1, ...)
				stopifnot(is(sfdr, "sFdr"))
				two.sided.p <- sfdr@p.value
				side.sign <- sign(statistic(sfdr)) # sign(statistic(as(sfdr, "MArrayLM")))
				PValue <- ifelse(side.sign < 0, two.sided.p / 2, 1 - two.sided.p / 2)
			}
			else if(is.function(PValue))
				PValue <- sapply(1:nrow(x.mat), function(i)
				{
					PValue(x.mat[i, ], y.mat[i, ])
				}) # apply(mat, 1, PValue)
			if(length(PValue) != nrow(y.mat))
			{ message("bad PValue length (2-sample)")}#browser
			names(PValue) <- nam
			if(!is(PValue, "Numeric"))
				PValue <- Numeric(PValue)
			if(any(PValue > 1, na.rm = TRUE) || !validObject(PValue) || !is(PValue, "Numeric"))
			{ message("PValue has element too large or other problem (2-sample)")}#browser
			if(verbose)
			{ 
				message("2-sample PValue:")
				print(class(PValue))
				print(summary(PValue))
			}
			present.PValue <- PValue[!is.na(PValue)]
			if(all(present.PValue[1] == present.PValue))
			{ message("all p-values are equal!")}#browser
			probabilityNull <- probabilityNull(PValue, ...)
		}
		else if(missing(PValue))
			PValue <- Numeric(numeric(0))
		else
			stop("PValue was specified even though probabilityNull is not a function (2-sample)")
	}
	to.return <- if(is(probabilityNull, return.class))
		probabilityNull
	else if(return.class %in% c("Prob0", "Fdr")) # from before 25 April 2008
	{
		if(!is.numeric(probabilityNull))
			probabilityNull <- fdr(probabilityNull)
		pn.ok <- is.numeric(probabilityNull) && nrow(x.mat) == length(probabilityNull) && all(names(probabilityNull) == nam)
		if(!pn.ok)
		{
			message("probabilityNull is not a valid vector of probabilities that the nulls are true, a function returning such, or an object such that the fdr method returns such (x, y).")
			#browser()
		}
		get.vec <- function(mat, FUN, ...)
		{
			ok <- try(is.function(FUN))
			if(is(ok, "try-error") || !ok)
			{ message("bad FUN in get.vec")}#browser
			vec <- apply(mat, 1, FUN, ...)
			names(vec) <- nam
			vec
		}
		combinedScale <- function(x.Scale, y.Scale)
		{
			stopifnot(all(names(x.Scale) == names(y.Scale)))
			Numeric(sqrt(x.Scale ^ 2 + y.Scale ^ 2))
		}
		subScale <- function(mat, ...)
		{
			Numeric(get.vec(mat = mat, FUN = Sd, na.rm = na.rm, mle = mle, ...))
		}
		superScale <- function(...){combinedScale(x.Scale = subScale(mat = x.mat, ...), y.Scale = subScale(mat = y.mat, ...))}
		Scale0 <- superScale(mu = 0)
		Scale1 <- superScale()
		Location1 <- get.vec(mat = x.mat, FUN = mean, na.rm = na.rm) - get.vec(mat = y.mat, FUN = mean, na.rm = na.rm)
		subSize <- function(mat)
		{
			Numeric(get.vec(mat = mat, FUN = function(object, ...){sum(!is.na(object), ...)}, na.rm = na.rm))
		}
		xSize <- subSize(mat = x.mat)
		ySize <- subSize(mat = y.mat)
		stopifnot(all(names(ySize) == names(xSize)))
		if(!is(PValue, "Numeric"))
		{
			PValue <- if(is.numeric(PValue))
				Numeric(PValue)
			else
			{
				warning(class(PValue), "PValue coerced to Numeric of length 0")
				Numeric(numeric(0))
			}
		}
		new("Prob0", probabilityNull, Scale0 = Scale0, Scale1 = Scale1, Location1 = Location1, xSize = xSize, ySize = ySize, PValue = PValue, mle = mle, annotation = paste(annotation(x), annotation(y), sep = " / "))
	}
	else
	{
		message("(x, y) cannot return ", return.class)
		#browser()
	}
	if(is(to.return, "Prob0") && return.class == "Fdr")
		to.return <- as(to.return, "Fdr")
	stopifnot(is(to.return, return.class))
	to.return
})
setMethod("Prob0", signature(x = "XprnSet", y = "missing"), function(x, y, ...)
{
	Prob0(x = logb(x), ...)
})
setMethod("Prob0", signature(x = "XprnSet", y = "XprnSet"), function(x, y, ...)
{
	Prob0(x = logb(x), y = logb(y), ...)
})
setMethod("Prob0", signature(x = "xprnSetPair", y = "missing"), function(x, y, ...) # added 21 April 2008
{
	Prob0(x = x@x, y = x@y, ...)
})
setMethod("Scale", signature(object = "Prob0"), function(object, null.hypothesis)
{
	Num <- if(missing(null.hypothesis) || length(null.hypothesis) == 0)
	{
		vec <- object * Scale(object = object, null.hypothesis = TRUE) + (1 - object) * Scale(object = object, null.hypothesis = FALSE)
		ifelse(is.nan(vec), Inf, vec)
	}
	else if(null.hypothesis)
		object@Scale0
	else
		object@Scale1
	Num <- as(Num, "Numeric")
	names(Num) <- names(object)
	Num
})
setMethod("Location", signature(object = "Prob0"), function(object, null.hypothesis, na.keep)
{
	if(missing(na.keep)) na.keep <- FALSE
	num <- if(missing(null.hypothesis) || length(null.hypothesis) == 0)
	{
		loc0 <- Location(object = object, null.hypothesis = TRUE)
		loc1 <- Location(object = object, null.hypothesis = FALSE)
		ifelse(object == 0, loc1, ifelse(object == 1, loc0, object * loc0 + (1 - object) * loc1)) # ifelse needed due to NAs
	}
	else if(null.hypothesis)
		rep(0, length(object@Location1))
	else
		object@Location1
	num <- as(num, "numeric")
	if(na.keep)
		num <- ifelse(is.na(object), as.numeric(NA), num)
	names(num) <- names(object)
	num
})
Order.Prob0 <- function(object, decreasing, by.PValue, by.Location1)
{
	p <- object@PValue
	if(missing(by.Location1))
	{
		by.Location1 <- FALSE
		message("Order.Prob0 by.Location1 = ", by.Location1)
	}
	if(missing(decreasing)) decreasing <- by.Location1 # FALSE
	boo <- !is.na(p)
	can.by.PValue <- sum(boo) >= 2 && any(p[boo] != p[boo][1])
	if(missing(by.PValue) || (by.PValue && !can.by.PValue))
	{
		by.PValue <- can.by.PValue && !by.Location1
		message("by.PValue = ", by.PValue)
	}
	if(by.PValue && by.Location1)
		stop("by.PValue and by.Location1 both TRUE")
	vec <- if(by.PValue) # order by p-value
		p
	else if(by.Location1)
		abs(object@Location1)
	else # order by local FDR
		as(object, "numeric")
	ord <- order(vec, decreasing = decreasing)
	stopifnot(length(ord) == length(object))
	ord
}
setMethod("Order", signature(object = "Prob0"), Order.Prob0)

arbitrarize <- function(object, ...){printGeneric("arbitrarize")}#browser
#removeMethods("arbitrarize")
setMethod("arbitrarize", signature(object = "Prob0"), function(object, ...)
{
	estimate <- object@Location1
	nam <- names(object)
	object@.Data <- new.pseudoFdr(estimate = estimate, ...)@fdr
	names(object) <- nam
	object
})

Probs0.biasEstimates <- function(P0, bias)
{
	new("Probs0.biasEstimates", P0 = P0, bias = bias)
}

loss <- function(x, y, ...){printGeneric("loss")}#browser
#removeMethods("loss")
setMethod("loss", signature(x = "xprnSet", y = "xprnSet"), function(x, y, normalization, ...)
{
	if(missing(normalization))
	{
		normalization.name <- "normalize"
		normalization <- eval(parse(text = normalization.name))
		message("loss normalization = ", normalization.name)
	}
	loss(x = xprnSetPair(x = x, y = y), normalization = normalization, ...)
})
setMethod("loss", signature(x = "xprnSetObject", y = "missing"), function(x, y, FUN, mle, hypothesis.test, PValue, probabilityNull, ...)
{
	message("\n\nBeginning computation of loss of ", annotation(x), " on ", date(), ".")
	if(identical(probabilityNull, sFdr) || identical(probabilityNull, pseudoFdr))
	{
		if(!missing(PValue) && !is.null(PValue))
			stop("PValue specified with incompatible probabilityNull")
		PValue <- NULL
	}
  else
  {
		if(missing(PValue))
		{
			if(missing(hypothesis.test))
				hypothesis.test <- t.test
			stopifnot(is.function(hypothesis.test))
			PValue <- PValueFUN(FUN = hypothesis.test)
		}
		else if(!missing(hypothesis.test))
			stop("hypothesis.test and PValue were both specified")
	}
	if(missing(mle))
	{
		mle <- FALSE
		message("loss mle = ", mle)
	}
	if(missing(probabilityNull))
	{
		pn.name <- "eFdr"
		probabilityNull <- eval(parse(text = pn.name)) # probabilityNullFUN(PValue = PValue)
		message("loss probabilityNull = ", pn.name)
	}
	if(missing(FUN))
		FUN <- squaredError
	stopifnot(is.function(FUN))
	olo <- oneLeftOut(x, FUN = Prob0, mle = mle, PValue = PValue, probabilityNull = probabilityNull, ...)
	P0 <- noneLeftOut(olo)
	stopifnot(is(P0, "Prob0"))
	lossBayes <- loss0 <- loss1 <- numeric(length(P0))
	get.loc <- function(x, null.hypothesis) {Location(object = x, null.hypothesis = null.hypothesis, na.keep = TRUE)}
	get.loss <- function(training.null.hypothesis, gene)
	{
		stopifnot(is.logical(training.null.hypothesis))
		mat <- sapply(1:length(olo@test), function(i)
		{
			predicted <- get.loc(olo@training[[i]], null.hypothesis = training.null.hypothesis)
			observed <- get.loc(olo@test[[i]], null.hypothesis = FALSE)
			FUN(predicted, observed)
		})
		if(!missing(gene))
		{
			ro <- mat[gene, ]
			names(ro) <- as.character(1:ncol(mat))
			print(ro)
		}
		vec <- as.numeric(apply(mat, 1, mean, na.rm = TRUE))
		stopifnot(length(vec) == length(P0))
		names(vec) <- names(P0)
		vec
	}
	loss0 <- get.loss(training.null.hypothesis = TRUE)
	loss1 <- get.loss(training.null.hypothesis = FALSE)
	lossBayes <- get.loss(training.null.hypothesis = logical(0))
#	bad.boo <- (P0 == 1 & (loss0 != lossBayes)) | (P0 == 0 & (loss1 != lossBayes))
#	if(any(bad.boo, na.rm = TRUE))
#	{ message(sum(bad.boo, na.rm = TRUE), " genes have unexpected losses")}#browser
	new("loss", lossBayes = lossBayes, loss0 = loss0, loss1 = loss1, FUN = FUN, P0 = P0)
})
setMethod("accumulate", signature(object = "loss"), function(object, indices, nfeatures, relative, ...)
{
	if(missing(relative))
	{
		relative <- FALSE
		message("accumulate loss relative = ", relative)
	}
	len <- length(object@P0)
	if(missing(indices))
	{
		if(missing(nfeatures))
		{
			nfeatures <- len
			message(class(object), " nfeatures = ", nfeatures)
		}
		indices <- 1:nfeatures
	}
	stopifnot(all(indices >= 1 & indices <= len))
	get.vec <- function(x, rel)
	{
		if(missing(rel)) rel <- relative
		if(rel)
			x <- x / Location(object@P0, null.hypothesis = FALSE) ^ 2
		vec <- accumulate(object = x, indices = indices, ...)
#		if(rel)
#			vec / get.vec(x = Location(object@P0, null.hypothesis = FALSE) ^ 2, rel = FALSE)
#		else
#			vec
		vec
	}
	vo1 <- try(validObject(object))
	if(is(vo1, "try-error") || !vo1)
	{ message("bad vo1")}#browser
	object0 <- object
	object@lossBayes <- get.vec(object@lossBayes)
	object@loss0 <- get.vec(object@loss0)
	object@loss1 <- get.vec(object@loss1)
	object@P0 <- object@P0[indices]
	vo2 <- try(validObject(object))
	if(is(vo2, "try-error") || !vo2)
	{ message("bad vo2")}#browser
	object
})
setMethod("Order", signature(object = "loss"), function(object, decreasing, ...)
{
	if(missing(decreasing)) decreasing <- FALSE
	Order(object@P0, decreasing = decreasing, ...)
})
setMethod("arbitrarize", signature(object = "loss"), function(object, foldchange.threshold, ...)
{
	if(missing(foldchange.threshold))
	{
		foldchange.threshold <- numeric(0)
		message(class(object), " foldchange.threshold = ", foldchange.threshold)
	}
	object@P0 <- arbitrarize(object = object@P0, foldchange.threshold = foldchange.threshold, ...)
	object@lossBayes <- if(length(foldchange.threshold) == 0)
		object@loss1
	else
	{
		lb <- -object@P0 # dummy
		boo <- !is.na(lb)
		p0 <- object@P0[boo]
		if(!all(p0 == 0 | p0 == 1))
		{ message("foldchange.threshold failed")}#browser
		lb[boo] <- ifelse(p0 == 1, object@loss0[boo], object@loss1[boo])
		if(sum(!is.na(lb)) > sum(boo) || !all(is.na(lb[!boo])))
		{ message("bad lb")}#browser
		names(lb) <- names(object@lossBayes)
		lb
	}
	object
})

losses <- function(x, ...){printGeneric("losses")}#browser
#removeMethods("losses")
setMethod("losses", signature(x = "ANY"), function(x, test.names, ...)
{
	if(missing(test.names))
	{
		test.names <- c("t.test", "wilcox.test")
		cat("test.names = ")
		print(test.names)
	}
	stopifnot(is.character(test.names) && length(test.names) >= 1)
	arglis <- list(...)
	ann <- if("y" %in% names(arglis))
	{
		paste(annotation(x), annotation(arglis$y), sep = " & ")
	}
	else
		annotation(x)
	lis <- lapply(test.names, function(test.name)
	{
		hypothesis.test <- eval(parse(text = test.name))
		stopifnot(is.function(hypothesis.test))
		message("\n\nComputing loss of ", ann, " from ", test.name, " on ", date(), ".")
		loss(x = x, hypothesis.test = hypothesis.test, ...)
	})
	names(lis) <- test.names
	losses(lis, test.names = test.names)
})
setMethod("losses", signature(x = "list"), function(x, test.names)
{
	if(missing(test.names))
		as(x, "losses")
	else
	{
		if(length(test.names) == 1)
			test.names <- rep(test.names, length(x))
		new("losses", x, test.names = test.names)
	}
})
losses0 <- function(...)
{
	losses(..., probabilityNull = eFdr0)
}


ProbabilityEstimate <- function(object, ...){printGeneric("ProbabilityEstimate")}#browser
#removeMethods("ProbabilityEstimate") # cf. probabilityEstimate of "reprod.s"
setMethod("ProbabilityEstimate", signature(object = "Prob0"), function(object, qRatio, lower.tail, null.hypothesis)
{
	if(missing(null.hypothesis))
		null.hypothesis <- logical(0)
	if(missing(qRatio))
	{
		qRatio <- 1
		message("ProbabilityEstimate qRatio = ", qRatio)
	}
	if(missing(lower.tail))
	{
		lower.tail <- logical(0)
		message("lower.tail not specified, so both tails are used")
	}
	stopifnot(is.logical(lower.tail) && length(lower.tail) <= 2)
	stopifnot(is.numeric(qRatio) && length(qRatio) == 1)
	qRatio <- Scalar(qRatio)
	get.p <- function(qRatio, lower.tail)
	{
		stopifnot(qRatio >= 1)
		as(ProbabilityEstimate(object = object, qRatio = if(lower.tail) 1 / qRatio else qRatio, lower.tail = lower.tail, null.hypothesis = null.hypothesis), "numeric")
	}
	PE <- if(length(lower.tail) == 0)
	{
		Numeric(get.p(qRatio = qRatio, lower.tail = TRUE) + get.p(qRatio = qRatio, lower.tail = FALSE))
	}
	else if(length(lower.tail) == 2)
	{
		stopifnot(qRatio >= 1)
		p1 <- get.p(qRatio = qRatio, lower.tail = lower.tail[1])
		p2 <- get.p(qRatio = qRatio, lower.tail = lower.tail[2])
		p1.p2 <- pmax(p1, p2, na.rm = TRUE)
		stopifnot(all(names(p1) == names(p2)))
		names(p1.p2) <- names(p1)
		Numeric(p1.p2)
	}
	else if(length(null.hypothesis) == 0) # added 070828 17:13
	{
		conditionalPE <- function(null.hypothesis)
		{
			ProbabilityEstimate(object = object, qRatio = qRatio, lower.tail = lower.tail, null.hypothesis = null.hypothesis)
		}
		PE0 <- conditionalPE(null.hypothesis = TRUE)
		PE1 <- conditionalPE(null.hypothesis = FALSE)
		vec <- object * PE0 + (1 - object) * PE1
		ifelse(is.nan(vec), Inf, vec)
		names(vec) <- names(object)
		if(any(is.nan(vec)))
		{ message("some vec elements NaN 1")}#browser
		Numeric(vec)
	}
	else
	{
		mean <- Location(object, null.hypothesis = null.hypothesis)
		sd <- Scale(object, null.hypothesis = null.hypothesis)
		vec <- pnorm(q = logb(qRatio), mean = mean, sd = sd, lower.tail = lower.tail, log.p = FALSE)
		names(vec) <- names(object)
		if(any(is.nan(vec)))
		{ message("some vec elements NaN 2")}#browser
		Numeric(vec)
	}
	if(any(PE > 1.0001, na.rm = TRUE))
	{ message("probability estimate exceeds 1")}#browser
	pe <- new("ProbabilityEstimate", PE, qRatio = qRatio, lower.tail = lower.tail, null.hypothesis = null.hypothesis)
	names(pe) <- names(PE)
	pe
})

