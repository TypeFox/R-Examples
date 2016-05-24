# estimate.s created by David Bickel on 30 August 2007.
# new.unsortableEstimate modified 7 April 2008, after Zahra Montazeri changed log to logb.

estimate <- function(x, y, ...){printGeneric("estimate"); browser()}
#removeMethods("estimate")
setMethod("estimate", signature(x = "matrix", y = "missing"), function(x, y, verbose, na.rm, ann)
{
  if(missing(na.rm))
  	na.rm <- FALSE
	if(missing(verbose))
	{
		verbose <- na.rm
		if(verbose)
			message("(x) verbose = ", verbose)
	}
  Location <- unsortableEstimate(x = x, estimator = meanHat, ann = ann) # get.vec(FUN = mean, na.rm = minor.na.rm)
  Scale <- unsortableEstimate(x = x, estimator = SeMeanHat, ann = ann) # Numeric(get.vec(FUN = Sd, na.rm = minor.na.rm, se.of.mean = TRUE))
  if(na.rm)
  {
  	warning("na.rm = TRUE not tested after implementing unsortableEstimate '['; try na.rm = FALSE if needed")
  	boo <- is.finite(Location) #& is.finite(log(Scale))
  	if(!all(boo))
  	{
  		if(verbose)
  			message("Only the ", sum(boo), " features that have finite locations are retained out of the original ", length(boo), " features.")
  		Location <- Location[boo]
  		Scale <- Scale[boo]
  	}
  }
  if(!is.character(ann))
  { message("bad (x) ann"); browser()}
  new("estimate", Location = Location, Scale = Scale, annotation = ann)
})
setMethod("estimate", signature(x = "xprnSet", y = "missing"), function(x, y, ...)
{
  estimate(x = exprs(x), ann = annotation(x), ...)
})
setMethod("estimate", signature(x = "XprnSet", y = "missing"), function(x, y, ...)
{
  estimate(x = logb(x), ...)
})
setMethod("estimate", signature(x = "ANY", y = "ANY"), function(x, y, verbose, na.rm, ...)
{
  if(missing(na.rm))
  	na.rm <- FALSE
	if(missing(verbose))
	{
		verbose <- na.rm
		if(verbose)
			message("(ANY, ANY) verbose = ", verbose)
	}
	get.est <- function(x){estimate(x = x, verbose = verbose, na.rm = na.rm, ...)}
  x.est <- get.est(x = x)
  y.est <- get.est(x = y)
  estimate(x = x.est, y = y.est, verbose = verbose, ...)
})
setMethod("estimate", signature(x = "estimate", y = "estimate"), function(x, y, verbose, ...)
{
	if(missing(verbose))
	{
		verbose <- FALSE
		if(verbose)
			message("(estimate, estimate) verbose = ", verbose)
	}
	x.nam <- names(x)
	y.nam <- names(y)
	nam <- intersect(x.nam, y.nam)
	stopifnot(length(nam) >= 1)
	if(verbose && length(nam) < max(length(x.nam), length(y.nam)))
		message(length(nam), " features in common among groups of ", length(x.nam), " and ", length(y.nam), " features.")
	x <- x[nam]
	y <- y[nam]
	Location <- unsortableEstimate(x = Location(x), y = Location(y)) # Location(x) - Location(y)
	get.Scale <- function(object){Scale(object, shrinkage = 0, ...)}
	Scale <- unsortableEstimate(x = get.Scale(x), y = get.Scale(y)) # Numeric(sqrt(get.Scale(x) ^ 2 + get.Scale(y) ^ 2))
  ann <- paste(annotation(x), annotation(y), sep = " / ")
  if(!is.character(ann))
  { message("bad (x, y) ann"); browser()}
  new("estimate", Location = Location, Scale = Scale, annotation = ann)
})
setMethod("Location", signature(object = "estimate"), function(object)
{
	Slot(object = object, name = "Location")
})
setMethod("Scale", signature(object = "estimate"), function(object, shrinkage)
{
	if(missing(shrinkage))
	{
		shrinkage <- 0.5
		message("shrinkage = ", shrinkage)
	}
	stopifnot(all(shrinkage >= 0 & shrinkage <= 1))
	unshrunkScale <- Slot(object = object, name = "Scale")
	if(length(shrinkage) == 1)
		shrinkage <- rep(shrinkage, length(unshrunkScale))
	stopifnot(length(unshrunkScale) == length(shrinkage))
	if(all(shrinkage == 0))
		unshrunkScale
	else
	{
		shrunkScale <- median(unshrunkScale, na.rm = TRUE)
		vec <- sqrt(shrinkage * shrunkScale ^ 2 + (1 - shrinkage) * unshrunkScale ^ 2)
		sca.ok <- length(vec) == length(unshrunkScale) && all(vec >= min(unshrunkScale, na.rm = TRUE) & vec <= max(unshrunkScale, na.rm = TRUE), na.rm = TRUE)
		if(!sca.ok)
		{ message("Scale estimate error"); browser()}
		names(vec) <- names(unshrunkScale)
		vec <- Numeric(vec)
		if(is.null(names(vec)))
		{ message("vec nameless"); browser()}
		new.unsortableEstimate(vec, annotation = annotation(object), estimator = unshrunkScale@estimator)
	}
})
setMethod("Order", signature(object = "estimate"), function(object, na.last, decreasing, ...)
{
	if(missing(decreasing))
		decreasing <- FALSE
	if(missing(na.last)) na.last <- logical(0)
	mod.t <- moderatedT(object = object, ...)
	if(length(na.last) == 0)
	{
		mod.t <- ifelse(is.finite(mod.t), mod.t, 0)
		if(any(is.na(mod.t))) stop("Order estimate error")
		na.last <- NA
	}
	else
		warning(paste("na.last = ", na.last, " is not recommended for an ", class(object), sep = ""))
	order(mod.t, na.last = na.last, decreasing = decreasing)
})
setMethod("moderatedT", signature(object = "estimate"), function(object, ...)
{
	Location(object = object) / Scale(object = object, ...)
})

new.unsortableEstimate <- function(vec, annotation, estimator)
{
	nam <- names(vec)
	stopifnot(is.character(nam))
	if(!is.character(annotation))
	{ message("scope problem"); browser()}
	est <- new("unsortableEstimate", vec, annotation = annotation, estimator = estimator)
	names(est) <- nam
#	names(est@.Data) <- nam
#	stopifnot(is.character(names(est@.Data)))
	est
}
unsortableEstimate <- function(x, y, estimator, ...){printGeneric("unsortableEstimate"); browser()}
#removeMethods("unsortableEstimate")
setMethod("unsortableEstimate", signature(x = "matrix", y = "missing", estimator = "Estimator"), function(x, y, estimator, verbose, na.rm, ann)
{
	stopifnot(1 %in% estimator@nsamples)
  if(missing(na.rm))
  	na.rm <- FALSE
	if(missing(verbose))
	{
		verbose <- na.rm
		if(verbose)
			message("(x) verbose = ", verbose)
	}
  get.vec <- function(FUN, ...)
  {
  	vec <- as.numeric(apply(x, 1, FUN = FUN, ...))
  	names(vec) <- rownames(x)
  	vec
  }
  estimate.vec <- if(is(estimator, "Estimator"))
  	get.vec(FUN = estimator)
  else
  {
	  minor.na.rm <- TRUE
  	get.vec(FUN = estimator, na.rm = minor.na.rm)
  }
  if(na.rm)
  {
  	boo <- is.finite(estimate.vec) #& is.finite(log(Scale))
  	if(!all(boo))
  	{
  		if(verbose)
  			message("Only the ", sum(boo), " features that have finite locations are retained out of the original ", length(boo), " features (unsortableEstimate).")
  		estimate.vec <- estimate.vec[boo]
  	}
  }
  if(!is.character(ann))
  { message("bad (x) ann"); browser()}
  if(is.null(names(estimate.vec)))
  { message("estimate.vec nameless"); browser()}
  new.unsortableEstimate(estimate.vec, annotation = ann, estimator = estimator)
})
setMethod("unsortableEstimate", signature(x = "xprnSet", y = "missing", estimator = "Estimator"), function(x, y, estimator, ann, ...)
{
	if(missing(ann))
		ann <- annotation(x)
  unsortableEstimate(x = exprs(x), ann = ann, estimator = estimator, ...)
})
setMethod("unsortableEstimate", signature(x = "XprnSet", y = "missing", estimator = "Estimator"), function(x, y, estimator, ...)
{
  unsortableEstimate(x = logb(x), estimator = estimator, ...)
})
setMethod("unsortableEstimate", signature(x = "matrix", y = "matrix", estimator = "Estimator"), function(x, y, estimator, verbose, na.rm, ann)
{
	stopifnot(2 %in% estimator@nsamples && nrow(x) == nrow(y) && all(rownames(x) == rownames(y)))
	ngenes <- nrow(x)
  if(missing(na.rm))
  	na.rm <- FALSE
	if(missing(verbose))
	{
		verbose <- na.rm
		if(verbose)
			message("unsortableEstimate(x, y) verbose = ", verbose)
	}
  get.vec <- function(FUN, ...)
  {
  	vec <- sapply(1:ngenes, function(i)
  	{
  		get.row <- function(object){object[i, ]}
  		FUN(x = get.row(x), y = get.row(y), ...)
  	}) # as.numeric(apply(x, 1, FUN = FUN, ...))
  	names(vec) <- rownames(x)
  	vec
  }
  estimate.vec <- if(is(estimator, "Estimator"))
  	get.vec(FUN = estimator)
  else
  {
	  minor.na.rm <- TRUE
  	get.vec(FUN = estimator, na.rm = minor.na.rm)
  }
  if(na.rm)
  {
  	boo <- is.finite(estimate.vec) #& is.finite(log(Scale))
  	if(!all(boo))
  	{
  		if(verbose)
  			message("Only the ", sum(boo), " features that have finite locations are retained out of the original ", length(boo), " features (unsortableEstimate).")
  		estimate.vec <- estimate.vec[boo]
  	}
  }
  if(!is.character(ann))
  { message("bad (x, y) ann"); browser()}
  if(is.null(names(estimate.vec)))
  { message("estimate.vec nameless"); browser()}
  new.unsortableEstimate(estimate.vec, annotation = ann, estimator = estimator)
})
setMethod("unsortableEstimate", signature(x = "xprnSet", y = "xprnSet", estimator = "Estimator"), function(x, y, estimator, ...)
{
	unsortableEstimate(x = exprs(x), y = exprs(y), ann = paste(annotation(x), annotation(y), sep = " vs. "), estimator = estimator, ...)
})
setMethod("unsortableEstimate", signature(x = "XprnSet", y = "XprnSet", estimator = "Estimator"), function(x, y, estimator, ...)
{
	unsortableEstimate(x = logb(x), y = logb(y), estimator = estimator, ...)
})
setMethod("unsortableEstimate", signature(x = "ANY", y = "ANY", estimator = "Estimator"), function(x, y, estimator, ...)
{
	get.est <- function(object){unsortableEstimate(x = object, estimator = estimator, ...)}
	if(2 %in% estimator@nsamples)
	{
		stop("code me and two-sample 'estimate'!")
	}
	else if(1 %in% estimator@nsamples)
	{
		message("estimator might be improved; press c to keep it anyway")
		unsortableEstimate(x = get.est(x), y = get.est(y))
	}
	else
		stop("bad estimator!")
})
# setMethod("unsortableEstimate", signature(x = "ANY", y = "ANY", estimator = "Estimator"), function(x, y, estimator, ...)
# {
# 	unsortableEstimate(x = x, y = y, estimator = estimator, estimatorType = estimator@type, ...)
# })
setMethod("unsortableEstimate", signature(x = "unsortableEstimate", y = "unsortableEstimate", estimator = "missing"), function(x, y, estimator, ...)
{
	stopifnot(length(x) == length(y))
	stopifnot(identical(x@estimator, y@estimator))
	combine.location.estimates <- function(x, y)
	{
		x - y
	}
	combine.variance.estimates <- function(x, y)
	{
		Numeric(x + y)
	}
	combine.scale.estimates <- function(x, y)
	{
		Numeric(sqrt(combine.variance.estimates(x = x ^ 2, y = y ^ 2)))
	}
	estimatorType <- x@estimator@type
	stopifnot(estimatorType == y@estimator@type)
	if(estimatorType == "location")
	{
		FUN <- combine.location.estimates
		sep <- " - "
	}
	else if(estimatorType == "variance")
	{
		FUN <- combine.variance.estimates
		sep <- " + "
	}
	else if(estimatorType == "scale")
	{
		FUN <- combine.scale.estimates
		sep <- " '+' "
	}
	else
		stop(paste("Cannot combine estimates of type", estimatorType))
	if(!is.function(FUN))
		stop("FUN not specified as a function that combines an estimate from x with an estimate from y")
	estimate.vec <- FUN(x = x, y = y)
	ann <- if(length(annotation(x)) == 1 && annotation(x) == annotation(y))
		annotation(x)
	else
		paste(annotation(x), annotation(y), sep = sep)
  if(is.null(names(estimate.vec)))
  { message("estimate.vec nameless 2"); browser()}
  new.unsortableEstimate(estimate.vec, annotation = ann, estimator = x@estimator)
})

new.sortableEstimate <- function(estimate.vec, LocationScale)
{
	nam <- names(estimate.vec)
	nam.ok <- is.character(nam) #&& is.character(names(estimate.vec@.Data))
	if(!nam.ok)
	{ message("error naming sortableEstimate"); browser()}
	est <- new("sortableEstimate", estimate.vec, LocationScale = LocationScale)
	names(est) <- nam
	est
}
sortableEstimate <- function(x, y, estimator, ...){printGeneric("sortableEstimate"); browser()}
#removeMethods("sortableEstimate")
setMethod("sortableEstimate", signature(x = "ANY", y = "missing", estimator = "Estimator"), function(x, y, estimator, ...)
{
	estimate.vec <- unsortableEstimate(x = x, estimator = estimator, ...)
	LocationScale <- estimate(x)
	new.sortableEstimate(estimate.vec, LocationScale = LocationScale)
})
setMethod("sortableEstimate", signature(x = "ANY", y = "ANY", estimator = "Estimator"), function(x, y, estimator, ...)
{
	estimate.vec <- unsortableEstimate(x = x, y = y, estimator = estimator, ...)
	LocationScale <- estimate(x, y)
	new.sortableEstimate(estimate.vec, LocationScale = LocationScale)
})
setMethod("sortableEstimate", signature(x = "ANY", y = "ANY", estimator = "ANY"), function(x, y, estimator, ...)
{
	stop("No estimator argument of class 'Estimator' was specified:")
})
setMethod("Order", signature(object = "sortableEstimate"), function(object, na.last, decreasing, ...)
{
	if(missing(decreasing))
		decreasing <- FALSE
	if(missing(na.last)) na.last <- logical(0)
	Order(object@LocationScale, na.last = na.last, decreasing = decreasing, ...)
})
setMethod("moderatedT", signature(object = "sortableEstimate"), function(object, ...)
{
	moderatedT(object@LocationScale, ...)
})

EstimateLeftOut <- function(object, estimator, shrinkage, ...){printGeneric("EstimateLeftOut"); browser()}
#removeMethods("EstimateLeftOut")
setMethod("EstimateLeftOut", signature(object = "ANY", estimator = "missing", shrinkage = "numeric"), function(object, estimator, shrinkage, ...)
{
	estimator.name <- "relative.frequencyHat"
	estimator <- eval(parse(text = estimator.name))
	message("estimator = ", estimator.name)
	EstimateLeftOut(object = object, estimator = estimator, shrinkage = shrinkage, ...)
})
setMethod("EstimateLeftOut", signature(object = "ANY", estimator = "Estimator", shrinkage = "numeric"), function(object, estimator, shrinkage, save.memory, ...)
{
	shrinkage <- Scalar(shrinkage)
	unsorted.elo <- EstimateLeftOut(object = object, estimator = estimator, call.sort = FALSE, ...)
#	e <- oneLeftOut(olo, estimate
#	olo <- sort(olo, 
#	order <- oneLeftOut(object = object, FUN = 
	is.sorted <- function(x) {sorted(x)}
	if(is.sorted(unsorted.elo))
		stop("unsorted.elo was sorted")
	if(is.na(shrinkage))
		unsorted.elo
	else
	{
		if(missing(save.memory))
		{
			save.memory <- TRUE
			message("save.memory = ", save.memory)
		}
		stopifnot(shrinkage <= 1)
		sorted.elo <- sort(unsorted.elo, shrinkage = shrinkage)
		if(save.memory)
		{
			message("Even more memory might be saved if oneLeftOut were only called once instead of the three times of the present implementation. ", date())
			sorted.olo <- oneLeftOut(sorted.elo, FUN = as, Class = "unsortableEstimate")
			if(is(sorted.olo, "EstimateLeftOut"))
				stop("sorted.olo is wrong class")
			sorted.elo <- new("EstimateLeftOut", sorted.olo, shrinkage = shrinkage)
			stopifnot(is(noneLeftOut(sorted.elo), "unsortableEstimate"))
		}
		else
			stopifnot(is(noneLeftOut(sorted.elo), "sortableEstimate"))
		if(!is.sorted(sorted.elo))
			stop("sorted.elo was not sorted")
		sorted.elo
	}
})
setMethod("EstimateLeftOut", signature(object = "ANY", estimator = "Estimator", shrinkage = "missing"), function(object, estimator, shrinkage, call.sort, ...)
{
	if(missing(call.sort) || call.sort)
		stop("Either specify shrinkage argument (recommended) or specify call.sort = FALSE.")
#	sortable.estimate <- function()
#	{
#		sortableEstimate(x = object, estimator = estimator)
#	}
	olo <- oneLeftOut(object, FUN = sortableEstimate, estimator = estimator, ...) # stop("missing shrinkage not yet implemented!!!!")
	new("EstimateLeftOut", olo, shrinkage = Scalar(as.numeric(NA)))
})
setMethod("sorted", signature(object = "EstimateLeftOut"), function(object)
{
	!is.na(object@shrinkage)
})

setMethod("biasEstimate", signature(object = "xprnSetObject", uncorrected = "missing"), function(object, uncorrected, estimator, shrinkage, verbose, factor.name, ...)
# object consists of ordered estimates
{
	if(missing(verbose)) verbose <- FALSE
	get.elo <- function(...)
	{
		EstimateLeftOut(object = object, estimator = estimator, shrinkage = shrinkage, verbose = verbose, ...)
	}
	elo <- if(missing(factor.name))
		get.elo()
	else
		get.elo(factor.name = factor.name)
	biasEstimate(object = elo, ...)
})
setMethod("biasEstimate", signature(object = "EstimateLeftOut", uncorrected = "missing"), function(object, uncorrected, jackknife, ...)
# object consists of ordered estimates
{
	if(missing(jackknife))
	{
		estimator <- noneLeftOut(object)@estimator
		jackknife <- FALSE # estimator@functional
		message("jackknife argument is missing, so jackknife = ", jackknife)
		# message("jackknife argument is missing and ", annotation(estimator), if(estimator@functional) " is " else " is not ", "a functional statistic, so jackknife = ", jackknife, "; cf. Efron (1982), Theorem 2.1 and pp. 10, 33-34, 44-45.")
	}
	biasEstimate(as(object, "oneLeftOut"), sorted = sorted(object), jackknife = jackknife, ...)
})

biasEstimates <- function(object, shrinkage, ...) {printGeneric("biasEstimates"); browser()}
#removeMethods("biasEstimates")
setMethod("biasEstimates", signature(object = "ANY", shrinkage = "numeric"), function(object, shrinkage, verbose, ...)
{
	stopifnot(length(shrinkage) >= 1)
	if(missing(verbose))
	{
		verbose <- TRUE
		if(verbose)
			message("biasEstimates verbose = ", verbose)
	}
	if(verbose)
	{
		ann <- try(annotation(object))
		if(is(ann, "try-error") || length(ann) == 0)
			ann <- "unannotated"
		message("\n\n\nBeginning biasEstimates iterations for ", ann, " on ", date())
	}
	lis <- lapply(shrinkage, function(shr)
	{
		if(verbose)
		{
			message("\n\n  Calling biasEstimate for ", ann, " with ", round(shr, 3), " shrinkage on ", date(), ".")
		}
		be <- biasEstimate(object = object, shrinkage = shr, verbose = verbose, ...)
		if(is(be, "try-error"))
		{ message("be error"); browser()}
		be
	})
	if(verbose)
		message("Ending biasEstimates iterations for ", ann, " on ", date(), "\n\n\n")
	new("biasEstimates", lis, shrinkage = shrinkage)
})
setMethod("biasEstimates", signature(object = "biasEstimates", shrinkage = "missing"), function(object, ...)
{
	lis <- lapply(object, biasEstimate, ...)
	new("biasEstimates", lis, shrinkage = object@shrinkage)
})

setMethod("Rank", signature(object = "biasEstimates"), function(object, ...)
{
	lis <- lapply(object, Rank, ...)
	ra <- lis[[1]]
  names(ra) <- NULL
  if(!all(sapply(lis, function(elem){all(ra == elem)})))
  	stop("different ranks")
  ra
})

setMethod("smoothData", signature(x = "biasEstimates", y = "missing"), function(x, y, class2, ...)
{
	if(missing(class2))
	{
		class2 <- "Lowess"
		message(class(x), "'s class2 = ", class2)
	}
	stopifnot(class2 %in% smoothData.classes)
	lis <- lapply(x, smoothData, class2 = class2, ...)
	new("biasEstimates", lis, shrinkage = x@shrinkage)
})

new.sim.biasEstimates <- function(gene.rank, uncorrected.estimate, rough.bias.estimate, smooth.bias.estimate, sample.size, means, sds, rxprnSet, estimator, shrinkage, true.rank, ...)
{
	new("sim.biasEstimates", gene.rank = Numeric(gene.rank), uncorrected.estimate = uncorrected.estimate, rough.bias.estimate = rough.bias.estimate, smooth.bias.estimate = smooth.bias.estimate, sample.size = Scalar(sample.size), means = means, sds = Numeric(sds), rxprnSet = rxprnSet, estimator = estimator, shrinkage = Scalar(shrinkage), true.rank = Scalar(true.rank), ...)
}

normal.rxprnSet <- function(sample.size, means, sds)
{
	stopifnot(length(means) == length(sds) && length(sample.size) == 1)
	mat <- t(sapply(means, function(mean){rnorm(mean = mean, sd = 1, n = sample.size)}))
	datf <- data.frame(random = factor(rep("a", sample.size)))
	colnames(mat) <- rownames(datf) <- sample.names <- as.character(paste("sample", 1:sample.size, sep = ""))
	# rownames(mat) <- gene.names <- paste("m", means, "s", sds, sep = "")
	rownames(mat) <- names(means)
	stopifnot(is.character(rownames(mat)) && !any(duplicated(rownames(mat))))
	xprnSet(phenoData = datf, exprs = mat, annotation = "normal.rxprnSet")
}
Bickel2004.means <- function(ngenes, max.mean, positive, exponent) # based on 2004 Bioinformatics paper
{
	if(missing(max.mean))
	{
		max.mean <- 5
		message("Bickel2004.means max.mean = ", max.mean)
	}
	if(missing(exponent))
	{
		exponent <- 8
		message("exponent = ", exponent)
	}
	recall <- function(ngenes, positive)
	{
		Bickel2004.means(ngenes = ngenes, max.mean = max.mean, positive = positive, exponent = exponent)
	}
	means <- if(missing(positive))
	{
		recurse <- function(positive)
		{
			half.ngenes <- ngenes / 2
			stopifnot(half.ngenes == floor(ngenes / 2))
			recall(ngenes = half.ngenes, positive = positive)
		}
		vec <- c(recurse(FALSE), recurse(TRUE))
		stopifnot(length(vec) == ngenes)
		vec
	}
	else if(positive)
	{
		max.mean * ((1:ngenes) / ngenes) ^ exponent
	}
	else
		-rev(recall(ngenes = ngenes, positive = TRUE))
	names(means) <- paste("m", as.character(means), sep = "")
	means
}

# normal.auc <- function(object, ...){printGeneric("auc"); browser()}

sim.biasEstimates <- function(nsamples, sample.size, ...){printGeneric("sim.biasEstimates"); browser()}
#removeMethods("sim.biasEstimates")
setMethod("sim.biasEstimates", signature(nsamples = "numeric", sample.size = "numeric"), function(nsamples, sample.size, ngenes, max.mean, means, sds, estimator, shrinkage, call.browser, f, true.rank, ...)
{
	if(missing(f))
	{
		f <- 1/100
		message("f = ", f)
	}
	if(missing(call.browser))
		call.browser <- FALSE
	if(missing(means))
	{
		if(missing(max.mean))
		{
			max.mean <- 1
			message("sim.biasEstimates max.mean = ", max.mean)
		}
		means <- Bickel2004.means(ngenes = ngenes, max.mean = max.mean)
	}
	if(missing(sds)) sds <- 1
	if(length(sds) == 1)
		sds <- rep(sds, length(means))
	if(missing(true.rank))
	{
		true.rank <- length(means)
		message("true.rank = ", true.rank)
	}
	gene.name <- names(means)[rank(means / sds) == true.rank]
	stopifnot(length(gene.name) == 1)
	if(missing(estimator))
		estimator <- relative.frequencyHat
	if(missing(shrinkage))
	{
		shrinkage <- 0
		message("shrinkage = ", shrinkage)
	}
	rxprnSet <- function()
	{
		normal.rxprnSet(sample.size = sample.size, means = means, sds = sds)
	}
	nam <- c("gene.rank", "uncorrected.estimate", "rough.bias.estimate", "smooth.bias.estimate")
	lis <- lapply(1:nsamples, function(i)
	{
		message("\nBeginning sample ", i, " of ", nsamples, " on ", date())
		get.bias.estimate <- function(object)
		{
			stopifnot(is(object, "biasEstimates") && length(object) == 1)
			as(object[[1]], "numeric")[gene.name]
		}
		es <- rxprnSet()
		x <- biasEstimates(es, estimator = estimator, shrinkage = shrinkage)
		gene.rank <- Rank(x[[1]][gene.name])
		uncorrected.estimate <- as(x[[1]][gene.name]@uncorrected, "numeric")
		rough.bias.estimate <- get.bias.estimate(object = x)
		get.smooth.bias.estimate <- function(x, smooth, f) # based on "plot", signature(x = "biasEstimates", y = "missing")
		{
			smooth.x <- function(smooth, x)
			{
				if(is.logical(smooth) && !smooth)
					x # cannot smooth x
				else
				{
					if(!is.function(smooth))
						smooth <- smoothDataFUN(f = f)
					smooth(x)
				}
			}
			graph.object <- smooth.x(smooth = smooth, x = x)
			if(FALSE) #call.browser)
			{
				message("after calling smooth.x; ", date())
				browser()
			}
#			nfeatures <- ceiling(length(means) / 2)
#			biasEstimates(object = graph.object, nfeatures = nfeatures, use.lower.ranks = true.rank <= nfeatures)
			get.bias.estimate(object = graph.object)
		}
		smooth.bias.estimate <- get.smooth.bias.estimate(x = x, smooth = TRUE, f = f) # replace me
		if(call.browser)
		{
			message("after calling get.smooth.bias.estimate; ", date())
			browser()
		}
		vec4 <- c(gene.rank = unname(gene.rank), uncorrected.estimate = unname(uncorrected.estimate), rough.bias.estimate = unname(rough.bias.estimate), smooth.bias.estimate = unname(smooth.bias.estimate))
		vec4.ok <- is.numeric(vec4) && length(vec4) == length(nam) && all(names(vec4) == nam)
		if(!vec4.ok)
		{ message("!vec4.ok"); browser()}
		vec4
	})
	get.vec <- function(name)
	{
		sapply(lis, function(elem){elem[name][[1]]})
	}
	gene.rank <- get.vec("gene.rank")
	uncorrected.estimate <- get.vec("uncorrected.estimate")
	rough.bias.estimate <- get.vec("rough.bias.estimate")
	smooth.bias.estimate <- get.vec("smooth.bias.estimate")
	new.sim.biasEstimates(gene.rank = gene.rank, uncorrected.estimate = uncorrected.estimate, rough.bias.estimate = rough.bias.estimate, smooth.bias.estimate = smooth.bias.estimate, sample.size = sample.size, means = means, sds = sds, rxprnSet = rxprnSet, estimator = estimator, shrinkage = shrinkage, true.rank = true.rank)
})

setMethod("Combine", signature(x = "sim.biasEstimates", y = "sim.biasEstimates"), function(x, y, ...)
{
	join <- function(name)
	{
		get.vec <- function(object){slot(object, name = name)}
		xvec <- get.vec(x)
		yvec <- get.vec(y)
		vec <- Combine(x = xvec, y = yvec)
		if(!is(vec, class(xvec)))
		{ message("bad Combine"); browser()}
		vec
	}
	x@gene.rank <- join("gene.rank")
	x@uncorrected.estimate <- join("uncorrected.estimate")
	x@rough.bias.estimate <- join("rough.bias.estimate")
	x@smooth.bias.estimate <- join("smooth.bias.estimate")
	x
})

true.bias <- function(object, ...)
{
	stopifnot(is(object, "sim.biasEstimates"))
	boo <- object@gene.rank == object@true.rank
	if(sum(boo, na.rm = TRUE) == 0)
	{
		false.bias <- 0
		message("The rank was never correctly estimated, so returning ", false.bias, " bias; press c to continue.")
		browser()
		false.bias
	}
	else
	{
		est <- object@uncorrected.estimate
		mean(est[boo], na.rm = TRUE) - mean(est, na.rm = TRUE)
	}
}

estimand.value <- function(object, gene, ...)
{
	stopifnot(is(object, "sim.biasEstimates"))
	if(object@estimator@annotation == "relative frequency")
	{
		if(missing(gene))
		{
			gene <- names(object@means)[object@true.rank == rank(object@means / object@sds)]
			message("gene = ", gene)
		}
		mean <- object@means[gene]
		sds <- object@sds
		if(is.null(names(sds)))
			names(sds) <- names(object@means)
		sd <- sds[gene]
#		stopifnot(length(mean) == 1 && length(sd) == 1)
		pnorm(q = 0, mean = mean, sd = sd, lower.tail = FALSE)
	}
	else
	{
		warning("Approximating estimand.value by estimates.")
		mean(object@uncorrected.estimate, na.rm = TRUE)
	}
}

# used by extended.r and/or extended.s (in estimate.s prior to 090716):
setMethod("mean", signature(x = "xprnSet"), function(x, ...) # simpler version: stat of data.s
{
	unsortableEstimate(x = x, estimator = meanHat, ...) # ann = annotation(x), 
})
setMethod("mean", signature(x = "XprnSet"), function(x, ...)
{
	mean(x = logb(x), ...)
})
setMethod("Density", signature(x = "xprnSet"), function(x, ...)
{
	Density(x = mean(x), ...)
})
