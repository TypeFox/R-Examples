

EstimatorTypes <- c("scale", "relative.scale", "inverse.relative.scale", "variance", "relative.variance", "location", "association", "distribution") # modified 24 April 2008 and 20 November 2008


#
setClass("Estimator", representation("function", functional = "logical", type = "character", annotation = "character", nsamples = "Numeric"))
setValidity("Estimator", function(object)
{
	sam.ok <- length(object@nsamples) >= 1
	len.ok <- length(object@functional) == 1 && length(object@type) == 1 && length(object@annotation) == 1
	type.ok <- len.ok && (object@type %in% EstimatorTypes)
	oks <- c(len.ok, type.ok, sam.ok)
	ok <- all(oks==T)
	if(!ok)
	{ printInvalid(object)}
	ok
})

setClass("unsortableEstimate", representation("numeric", annotation = "character", estimator = "Estimator"))
setValidity("unsortableEstimate", function(object)
{
	len.ok <- length(object@annotation) <= 1
	nam <- names(object)
	nam.ok <- is.character(nam)
	if(!nam.ok) warning(class(object), " object has no names, possibly due to creation of sortableEstimate")
	ok <- len.ok
	if(!ok)
	{ printInvalid(object)}
	ok
})
setClass("estimate", representation(Location = "unsortableEstimate", Scale = "unsortableEstimate", annotation = "character")) #, Size = "Numeric"))
setValidity("estimate", function(object)
{
	vec.names <- setdiff(slotNames(object), "annotation")
	len.ok <- all(length(object@Location) == sapply(vec.names, function(name){length(slot(object, name))})) && length(object@annotation) <= 1
	nam <- names(object)
	nam.ok <- len.ok && is.character(nam) && all(sapply(vec.names, function(name){all(nam == names(slot(object, name)))}))
	ok <- nam.ok
	if(!ok)
	{ printInvalid(object)}
	ok
})







setClass("sortableEstimate", representation("unsortableEstimate", LocationScale = "estimate"))
setValidity("sortableEstimate", function(object)
{
	len.ok <- length(object) == length(object@LocationScale@Scale)
	nam.ok <- len.ok && all(names(object) == names(object@LocationScale))
	ok <- nam.ok
	if(!ok)
	{ printInvalid(object)}
	ok
})

#
setClass("basicFdr", representation(fdr = "Numeric", p.value = "Numeric"))
setValidity("basicFdr", function(object)
{
	nam <- names(object) # as(object, "Numeric"))
	nam.ok <- !is.null(nam) && sameNames(object@fdr, object@p.value)
	range.ok <- all(c(object@fdr, object@p.value) <= 1.001, na.rm = TRUE)
	ok <- nam.ok && range.ok
	if(!ok)
	{ message("invalid basicFdr")}
	ok
})

setClass("eFdr", representation("list", fdr = "Numeric", zz = "numeric")) # e is for Efron
setValidity("eFdr", function(object)
{
	nam <- names(object) # as(object, "Numeric"))
	nam.ok <- !is.null(nam) && sameNames(object@fdr, object@zz)
#	if(nam.ok)
#		message("Some names: ", paste(nam[1:10], collapse = "; "))
	ok <- nam.ok
	if(!ok)
	{ message("invalid eFdr")}
	ok
})


setClass("sFdr", representation("list", fdr = "Numeric", p.value = "Numeric", t = "numeric", proportion = "Scalar")) # s is for Smyth
setValidity("sFdr", function(object)
{
	nam <- names(object) # as(object, "Numeric"))
	nam.ok <- !is.null(nam) && sameNames(object@fdr, object@p.value, object@t)
#	if(nam.ok)
#		message("Some names: ", paste(nam[1:10], collapse = "; "))
	ok <- nam.ok
	if(!ok)
	{ message("invalid sFdr")}
	ok
})
setClass("bFdr", representation(fdr = "Numeric", p.value = "Numeric", prior_fdr = "Scalar")) # b is for Bickel
setValidity("bFdr", function(object)
{
	nam <- names(object) # as(object, "Numeric"))
	nam.ok <- !is.null(nam) && sameNames(object@fdr, object@p.value)
	range.ok <- all(c(object@fdr, object@p.value, object@prior_fdr) <= 1.001, na.rm = TRUE)
	ok <- nam.ok && range.ok
	if(!ok)
	{ message("invalid bFdr")}
	ok
})

setClass("pseudoFdr", representation(fdr = "Numeric", estimate = "numeric", FUN = "function"))
setValidity("pseudoFdr", function(object)
{
	nam <- names(object)
	nam.ok <- !is.null(nam) && sameNames(object@fdr, object@estimate)
	ok <- nam.ok
	if(!ok)
	{ message("invalid pseudoFdr")}
	ok
})


#
setClass("Prob0", representation("Numeric", Scale0 = "Numeric", Scale1 = "Numeric", Location1 = "numeric", xSize = "Numeric", ySize = "Numeric", PValue = "Numeric", mle = "logical", annotation = "character"))
setValidity("Prob0", function(object)
{
	sn <- setdiff(slotNames(object), c("mle", ".Data", "annotation"))
	if(length(object@PValue) == 0)
		sn <- setdiff(sn, c("PValue"))
	len.ok <- all(length(object) == sapply(sn, function(name){length(slot(object = object, name = name))})) && length(object@mle) == 1 && length(object@annotation) <= 1
	nam.ok <- len.ok && all(sapply(sn, function(name){all(names(object) == names(slot(object = object, name = name)))}))
	val.ok <- all(object <= 1.001, na.rm = TRUE)
	ok <- nam.ok && val.ok
	if(!ok){printInvalid(object)}
	ok
})
setClass("Probs0.biasEstimates", representation(P0 = "list", bias = "list"))
setValidity("Probs0.biasEstimates", function(object)
{
	len.ok <- length(object@P0) == length(object@bias)
	P0.nam <- names(object@P0)
	bias.nam <- names(object@bias)
	nam.ok <- len.ok && ((is.null(P0.nam) && is.null(bias.nam)) || all(P0.nam == bias.nam))
	cla.ok <- len.ok && nam.ok
	ok <- nam.ok && cla.ok
	if(!ok)
	{	printInvalid(object)}
	ok
})
#
setClass("loss", representation(lossBayes = "numeric", loss0 = "numeric", loss1 = "numeric", P0 = "Prob0", FUN = "function"))
setValidity("loss", function(object)
{
	is_name.ok <- function(x)
	{
		nam <- names(x)
		is.character(nam) && !any(is.na(nam))
	}
	loss.ok <- function(vec){is_name.ok(vec) && length(vec) == length(object@P0) && all(names(vec) == names(object@P0))}
	nam.ok <- is_name.ok(object@P0)
	l.ok <- loss.ok(object@lossBayes) && loss.ok(object@loss0) && loss.ok(object@loss1)
	nl.ok <- nam.ok && l.ok # && length(object@annotation) == 1
	if(!identical(nl.ok, TRUE))
	{ printInvalid(object)}
	nl.ok
})


setClass("losses", representation("list", test.names = "character"))
setValidity("losses", function(object)
{
	len.ok <- length(object) == length(object@test.names)
	cla.ok <- all(sapply(object, is, class2 = "loss"))
	ok <- len.ok && cla.ok
	if(!ok) printInvalid(object)
	ok
})

#

setClass("oneLeftOut", representation(training = "list", test = "list", noneLeftOut = "list"))
setValidity("oneLeftOut", function(object)
{
  len.ok <- length(object@training) >= 2 && length(object@training) == length(object@test) && length(object@noneLeftOut) == 1
  lens.ok <- all(length(object@noneLeftOut) == sapply(object@test, length)) && all(length(object@noneLeftOut) == sapply(object@training, length))
  cla.ok <- all(class(object@noneLeftOut) == sapply(object@training, class) && class(object@noneLeftOut) == sapply(object@test, class))
  ok <- len.ok && lens.ok && cla.ok
  if(!ok){message("invalid oneLeftOut")}
  ok
})

#

setClassUnion("Fdr", c("basicFdr"))
setAs(from = "Fdr", to = "basicFdr", function(from)
{
	new_basicFdr(fdr = fdr(from), p.value = pvalue(from))
})
setClassUnion("numericEstimate", c("numeric", "unsortableEstimate", "sortableEstimate")) # used in slot of biasEstimate class
setClassUnion("Fdr", c("sFdr", "eFdr", "bFdr", "pseudoFdr", "basicFdr"))
#
#
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
setAs(from = "eFdr", to = "Numeric", function(from)
{
	from@fdr
})
setAs(from = "eFdr", to = "numeric", function(from)
{
	as(as(from, "Numeric"), "numeric")
})
setAs(from = "bFdr", to = "Numeric", function(from)
{
	from@fdr
})
setAs(from = "bFdr", to = "numeric", function(from)
{
	as(as(from, "Numeric"), "numeric")
})
#
Probs0.biasEstimates <- function(P0, bias){new("Probs0.biasEstimates", P0 = P0, bias = bias)}
p.value <- function(object, ...){printGeneric("p.value")}#removeMethods("p.value")
Scale <- function(object, ...){printGeneric("Scale")}#removeMethods("Scale")
Location <- function(object, ...){printGeneric("Location")}#removeMethods("Location")
Order <- function(object, ...){printGeneric("Order");message("generic Order")}
#

setMethod("p.value", signature(object = "list"), function(object, ...)
{
	object$p.value#Slot(object = object, name = "p.value", ...)
})
setMethod("p.value", signature(object = "Fdr"), function(object, ...)
{
	if(is(object,"list")){object$p.value}
        else{Slot(object = object, name = "p.value", ...)}
})

#
setMethod("nannotation", signature(object = "Estimator"), function(object){object@annotation})
#
setAs(from = "unsortableEstimate", to = "numeric", function(from)
{
	vec <- as.numeric(from)
	names(vec) <- names(from)
	vec
})
setMethod("[", signature(x = "unsortableEstimate", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	vec <- as(x, "numeric")[i]
#	names(vec) <- names(x)[i]
	if(is.null(names(vec)))
	{ message("'[' names error")}
	new_unsortableEstimate(vec, annotation = x@annotation, estimator = x@estimator)
})


setAs(from = "sortableEstimate", to = "unsortableEstimate", function(from)
{
	vec <- as.numeric(from)
	if(is.null(names(vec)))
	{
		names(vec) <- if(is.null(names(from))) # this occurs with new("unsortableEstimate", ...)
			names(from@LocationScale)
		else
			names(from)
	}
	if(is.null(names(vec)))
	{ warning("from = 'sortableEstimate' vec lacks names")}
	new_unsortableEstimate(vec, annotation = from@annotation, estimator = from@estimator)
})
setAs(from = "sortableEstimate", to = "numeric", function(from)
{
	as(as(from, "unsortableEstimate"), "numeric")
})
#
#
Extract.sortableEstimate <- function(x, i, j, drop){
	stopifnot(is(x, "sortableEstimate"))
	if(is.null(names(x)))
	{ message("'[' sortableEstimate names error")}
	vec <- as(x, "unsortableEstimate")[i]
	LocationScale <- x@LocationScale[i]
	new_sortableEstimate(vec, LocationScale = LocationScale)
}

setMethod("[", signature(x = "sortableEstimate", i = "ANY", j = "missing"), Extract.sortableEstimate)
setMethod("sort", signature(x = "sortableEstimate"), function(x, decreasing = logical(0), ...)
{#XXX|:removed decreasing = "ANY"
	if(missing(decreasing) || length(decreasing) == 0)
		decreasing <- FALSE
	if(is.null(names(x)))
	{ message("sort sortableEstimate names error")}
	x[Order(object = x, decreasing = decreasing, ...)]
})

#
#
SeMeanHat <- new("Estimator", function(...){
    Sd(..., se.of.mean = TRUE)}, functional = FALSE, type = "scale", annotation = "SE of mean", nsamples = nNumeric(1))
meanHat <- new("Estimator", function(x, y, ...){
	sample.mean <- function(object) {mean(object, na.rm = TRUE, ...)}
	if(missing(y))
		sample.mean(x)
	else
		sample.mean(x) - sample.mean(y)
}, functional = TRUE, type = "location", annotation = "sample mean", nsamples = nNumeric(c(1, 2)))
meanHat0 <- new("Estimator", function(x, y, ...){
	0
}, functional = TRUE, type = "location", annotation = "zero mean", nsamples = nNumeric(c(1, 2)))


#
setMethod("names", signature(x = "estimate"), function(x)
{
	names(Location(x))
})
setMethod("[", signature(x = "estimate", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	x@Location <- x@Location[i]
	x@Scale<- x@Scale[i]
	x
})
setMethod("sort", signature(x = "estimate"), function(x, decreasing = logical(0), ...)
{#XXX|:removed decreasing = "ANY"
	if(missing(decreasing) || length(decreasing) == 0)
		decreasing <- FALSE
	x[Order(object = x, decreasing = decreasing, ...)]
})



#---
 ProbabilityEstimate <- function(object, ...){printGeneric("ProbabilityEstimate")}#browser#removeMethods("ProbabilityEstimate") # cf. probabilityEstimate of "reprod.s"
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
	qRatio <- nScalar(qRatio)
	get.p <- function(qRatio, lower.tail)
	{
		stopifnot(qRatio >= 1)
		as(ProbabilityEstimate(object = object, qRatio = if(lower.tail) 1 / qRatio else qRatio, lower.tail = lower.tail, null.hypothesis = null.hypothesis), "numeric")
	}
	PE <- if(length(lower.tail) == 0)
	{
		nNumeric(get.p(qRatio = qRatio, lower.tail = TRUE) + get.p(qRatio = qRatio, lower.tail = FALSE))
	}
	else if(length(lower.tail) == 2)
	{
		stopifnot(qRatio >= 1)
		p1 <- get.p(qRatio = qRatio, lower.tail = lower.tail[1])
		p2 <- get.p(qRatio = qRatio, lower.tail = lower.tail[2])
		p1.p2 <- pmax(p1, p2, na.rm = TRUE)
		stopifnot(all(names(p1) == names(p2)))
		names(p1.p2) <- names(p1)
		nNumeric(p1.p2)
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
		ifelse(is.na(vec), Inf, vec)
		names(vec) <- names(object)
		if(any(is.na(vec)))
		{ message("some vec elements NaN 1")}#browser
		nNumeric(vec)
	}
	else
	{
		mean <- Location(object, null.hypothesis = null.hypothesis)
		sd <- Scale(object, null.hypothesis = null.hypothesis)
		vec <- pnorm(q = logb(qRatio), mean = mean, sd = sd,  lower.tail = lower.tail, log.p = FALSE) 
		names(vec) <- names(object)
		if(any(is.na(vec)))
		{ message("some vec elements NaN 2")}#browser
		nNumeric(vec)
	}
	if(any(PE > 1.0001, na.rm = TRUE))
	{ message("probability estimate exceeds 1")}#browser
	pe <- new("ProbabilityEstimate", PE, qRatio = qRatio, lower.tail = lower.tail, null.hypothesis = null.hypothesis)
	names(pe) <- names(PE)
	pe
})

#
squaredError <- function(predicted, observed, object, ...){printGeneric("squaredError")}#removeMethods("squaredError")
setMethod("squaredError", signature(predicted = "numeric", observed = "numeric", object = "missing"), function(predicted, observed, object, relative, decay.scale)
{
	if(missing(relative)) relative <- !missing(decay.scale)
	stopifnot(is(predicted, "numeric") && is(observed, "numeric"))
	if(length(predicted) == length(observed))
		stopifnot(all(names(predicted) == names(observed)))
	else
		stopifnot(any(1 == c(length(predicted), length(observed))))
	err <- (predicted - observed) ^ 2
	if(relative)
	{
		stopifnot(length(observed) %in% c(1, length(err)))
		err <- err / observed ^ 2
	}
	if(!missing(decay.scale))
	{
		stopifnot(is.numeric(decay.scale) && length(decay.scale) %in% c(1, length(err)))
		err <- exp(- err / decay.scale)
		if(any(is.na(err)) || any(err < 0 | err > 1))
		{ message("bad err; decay.scale: ", decay.scale)}
	}
	err
})


#

noneLeftOut <- function(object, ...){message("generic noneLeftOut")}#removeMethods("noneLeftOut")
setMethod("noneLeftOut", signature(object = "oneLeftOut"), function(object, ...)
{
  (object@noneLeftOut)#[[1]]
})

oneLeftOut <- function(object, ...){message("generic oneLeftOut")}#removeMethods("oneLeftOut")
setMethod("oneLeftOut", signature(object = "numeric"), function(object, FUN, ...)
{	#XXX|:multiple local function definitions for "fun" with different formal arguments
	#XXX|:removed , FUN = "function", 
	assert.is(FUN,"function")
        arglis_FUN<-list(...)
#fun <- function(x){FUN(x,  ...)}XXX|:changed	
  fun <- function(x){do.call(FUN,c(list(x),arglis_FUN))}#XXX|:changed
  noneLeftOut <- list(fun(object))
  js <- 1:length(object)
  training <- lapply(js, function(j)
  {
    boo <- js != j
    stopifnot(sum(!boo) == 1 && length(boo) == length(object))
    fun(object[boo])
  })
  test <- lapply(js, function(j)
  {
    fun(object[j])
  })
  new("oneLeftOut", training = training, test = test, noneLeftOut = noneLeftOut)
})
setMethod("oneLeftOut", signature(object = "nxprnSet"), function(object, FUN, verbose, factor_name, ...)
{#XXX|:removed , FUN = "function"
    
	assert.is(FUN,"function")
        arglis_FUN<-list(...)
	if(missing(verbose)) verbose <- FALSE
	if(!missing(factor_name) && is.character(factor_name) && length(factor_name) == 1) # two-sample
	{
		oneLeftOut(object = new_nxprnSetPair(x = object, factor_name = factor_name), FUN = FUN, verbose = verbose, ...)
	}
	else if(missing(factor_name)) # single-sample; stopping condition
	{	fun <- function(x){do.call(FUN,c(list(x),arglis_FUN))}# , XXX|:changed
		#fun <- function(x){FUN(x, verbose = verbose, ...)} verbose passed 31 October 2007
		noneLeftOut <- list(fun(object))
		js <- 1:ncol(nexprs(object))
		training <- lapply(js, function(j)
		{
			boo <- js != j
			stopifnot(sum(!boo) == 1 && length(boo) == ncol(nexprs(object)))
			fun(object[, boo])
		})
		test <- lapply(js, function(j)
		{
			fun(object[, j])
		})
		if(verbose)
		{
			message("oneLeftOut mapped ", class(object), " to ", 2 * length(training) + 1, " objects of class ", class(noneLeftOut[[1]]))
			print(FUN)
		}
		new("oneLeftOut", training = training, test = test, noneLeftOut = noneLeftOut)
	}
	else
		stop(paste("bad oneLeftOut factor_name:", factor_name))
})
setMethod("oneLeftOut", signature(object = "nxprnSetPair"), function(object, FUN, verbose,  ...)
{#XXX|:removed , FUN = "function", changed "..." by arglis_FUN=list()
	assert.is(FUN,"function")
        arglis_FUN<-list(...)
	if(missing(verbose)) verbose <- FALSE
	x <- object@x
	y <- object@y
	fun <- function(x, y)
	{
		if(missing(y))
			stop("y is missing in the two-sample case")
	  #FUN(x = x, y = y, ...)
	  do.call(FUN,c(list(x=x,y=y),arglis_FUN))# , XXX|:changed
	

	}
	noneLeftOut <- list(fun(x = x, y = y))
	get.js <- function(object) {1:ncol(nexprs(object))}
	x.js <- get.js(x)
	y.js <- get.js(y)
	training <- test <- list()
	get.es <- function(object, j, negate.j)
	{
		if(length(j) != 1) stop("bad j, get.es")
		js <- get.js(object)
		boo <- if(negate.j)
			js != j
		else
			js == j
		size <- if(negate.j)
			ncol(nexprs(object)) - 1
		else
			1
		stopifnot(sum(boo) == size && length(boo) == ncol(nexprs(object)))
		object[, boo]
	}
	for(x.j in x.js)
	{
		for(y.j in y.js)
		{
			if(verbose)
			{
				paste.iteration <- function(i, j){paste("(", i, ", ", j, ")", sep = "")}
				message("      objects ", paste.iteration(x.j, y.j), " of ", paste.iteration(max(x.js), max(y.js)), " began on ", date())
			}
			get.lis <- function(negate)
			{
				elem <- fun(x = get.es(object = x, j = x.j, negate = negate), y = get.es(object = y, j = y.j, negate = negate))
				ok <- is(elem, class(noneLeftOut[[1]])) && length(elem) == length(noneLeftOut[[1]])
				if(!ok){message("get.lis_elem error")}
				list(elem)
			}
			training <- c(training, get.lis(negate = TRUE))
			test <- c(test, get.lis(negate = FALSE))
			if(length(training) != length(test))
			{ message("training and test error")}
		}
	}
	if(verbose)
	{
		message("two-sample oneLeftOut mapped ", class(object), " to ", 2 * length(training) + 1, " objects of class ", class(noneLeftOut[[1]]))
		print(FUN)
	}
	olo <- try(new("oneLeftOut", training = training, test = test, noneLeftOut = noneLeftOut))
	if(is(olo, "try-error") || !validObject(olo))
	{ message('"oneLeftOut", signature(object = "nxprnSetPair", FUN = "function") error')}
	olo
})
setMethod("oneLeftOut", signature(object = "oneLeftOut"), function(object, FUN, ...)
{#XXX|:multiple local function definitions for "fun" with different formal arguments
	##XXX|:removed , FUN = "function",
	assert.is(FUN,"function")
        
	arglis <- arglis_FUN<-list(...)
	alo.nam <- "anotherLeftOut"
	if(alo.nam %in% names(arglis))
	{
		alo <- arglis[alo.nam][[1]]
		if(!is(alo, "oneLeftOut") || length(arglis) != 1)
		{ message('cannot use "oneLeftOut", signature(object = "oneLeftOut", FUN = "function") this way')}
		# arglis <- arglis[names(arglis) != alo]
		#funx <- function(x, y){FUN(x, y)}
                funx <- function(x,y){do.call(FUN,c(list(x,y),arglis))}#XXX|:changed
		change.lis <- function(name)
		{
			ok <- name %in% slotNames(object) && name %in% slotNames(alo)
			if(!ok){ message("change.lis error")}
			xlis <- slot(object = object, name = name)
			ylis <- slot(object = alo, name = name)
			compat <- length(xlis) == length(ylis) && sameNames(xlis, ylis)
			if(!compat)
			{ message("inner problem of oneLeftOut")}
			new_lis <- lapply(1:length(xlis), function(i)
			{
				funx(x = xlis[[i]], y = ylis[[i]])
			})
			names(new_lis) <- names(xlis)
			new_lis
		}
		training <- change.lis("training")
		test <- change.lis("test")
		noneLeftOut <- change.lis("noneLeftOut")
	}
	else
	{
		#fun <- function(x){FUN(x, ...)}
		fun <- function(x){do.call(FUN,c(list(x),arglis))}#XXX|:changed
		training <- lapply(object@training, fun)
		test <- lapply(object@test, fun)
		noneLeftOut <- lapply(object@noneLeftOut, fun)
	}
	len.ok <- length(training) == length(object@training) && length(test) == length(object@test) && length(noneLeftOut) == length(object@noneLeftOut)
	if(!len.ok)
	{ message("length problem in oneLeftOut")}
	new("oneLeftOut", training = training, test = test, noneLeftOut = noneLeftOut)
})

#--
estimate <- function(x, y, ...){printGeneric("estimate")}#removeMethods("estimate")
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
  Location <- unsortableEstimate(x = x, estimator = meanHat, ann = ann) # get.vec(FUN = mean, na.rm = minor_na.rm)
  Scale <- unsortableEstimate(x = x, estimator = SeMeanHat, ann = ann) # nNumeric(get.vec(FUN = Sd, na.rm = minor_na.rm, se.of.mean = TRUE))
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
  { message("bad (x) ann")}
  new("estimate", Location = Location, Scale = Scale, annotation = ann)
})
setMethod("estimate", signature(x = "nxprnSet", y = "missing"), function(x, y, ...)
{
  estimate(x = nexprs(x), ann = nannotation(x), ...)
})
setMethod("estimate", signature(x = "nXprnSet", y = "missing"), function(x, y, ...)
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
	Scale <- unsortableEstimate(x = get.Scale(x), y = get.Scale(y)) # nNumeric(sqrt(get.Scale(x) ^ 2 + get.Scale(y) ^ 2))
  ann <- paste(nannotation(x), nannotation(y), sep = " / ")
  if(!is.character(ann))
  { message("bad (x, y) ann")}
  new("estimate", Location = Location, Scale = Scale, annotation = ann)
})
#
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
		{ message("Scale estimate error")}
		names(vec) <- names(unshrunkScale)
		vec <- nNumeric(vec)
		if(is.null(names(vec)))
		{ message("vec nameless")}
		new_unsortableEstimate(vec, annotation = nannotation(object), estimator = unshrunkScale@estimator)
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
#
moderatedT <- function(object, ...){printGeneric("moderatedT")}
setMethod("moderatedT", signature(object = "estimate"), function(object, ...)
{
	Location(object = object) / Scale(object = object, ...)
})

#
new_unsortableEstimate <- function(vec, annotation, estimator){
	nam <- names(vec)
	stopifnot(is.character(nam))
	if(!is.character(annotation))
	{ message("scope problem")}
	est <- new("unsortableEstimate", vec, annotation = annotation, estimator = estimator)
	names(est) <- nam
#	names(est@.Data) <- nam
#	stopifnot(is.character(names(est@.Data)))
	est
}
unsortableEstimate <- function(x, y, estimator, ...){printGeneric("unsortableEstimate")}#removeMethods("unsortableEstimate")
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
	  minor_na.rm <- TRUE
  	get.vec(FUN = estimator, na.rm = minor_na.rm)
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
  { message("bad (x) ann")}
  if(is.null(names(estimate.vec)))
  { message("estimate.vec nameless")}
  new_unsortableEstimate(estimate.vec, annotation = ann, estimator = estimator)
})
setMethod("unsortableEstimate", signature(x = "nxprnSet", y = "missing", estimator = "Estimator"), function(x, y, estimator, ann, ...)
{
	if(missing(ann))
		ann <- nannotation(x)
  unsortableEstimate(x = nexprs(x), ann = ann, estimator = estimator, ...)
})
setMethod("unsortableEstimate", signature(x = "nXprnSet", y = "missing", estimator = "Estimator"), function(x, y, estimator, ...)
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
	  minor_na.rm <- TRUE
  	get.vec(FUN = estimator, na.rm = minor_na.rm)
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
  { message("bad (x, y) ann")}
  if(is.null(names(estimate.vec)))
  { message("estimate.vec nameless")}
  new_unsortableEstimate(estimate.vec, annotation = ann, estimator = estimator)
})
setMethod("unsortableEstimate", signature(x = "nxprnSet", y = "nxprnSet", estimator = "Estimator"), function(x, y, estimator, ...)
{
	unsortableEstimate(x = nexprs(x), y = nexprs(y), ann = paste(nannotation(x), nannotation(y), sep = " vs. "), estimator = estimator, ...)
})
setMethod("unsortableEstimate", signature(x = "nXprnSet", y = "nXprnSet", estimator = "Estimator"), function(x, y, estimator, ...)
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
		nNumeric(x + y)
	}
	combine.scale.estimates <- function(x, y)
	{
		nNumeric(sqrt(combine.variance.estimates(x = x ^ 2, y = y ^ 2)))
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
	ann <- if(length(nannotation(x)) == 1 && nannotation(x) == nannotation(y))
		nannotation(x)
	else
		paste(nannotation(x), nannotation(y), sep = sep)
  if(is.null(names(estimate.vec)))
  { message("estimate.vec nameless 2")}
  new_unsortableEstimate(estimate.vec, annotation = ann, estimator = x@estimator)
})

#
new_sortableEstimate <- function(estimate.vec, LocationScale)
{
	nam <- names(estimate.vec)
	nam.ok <- is.character(nam) #&& is.character(names(estimate.vec@.Data))
	if(!nam.ok)
	{ message("error naming sortableEstimate")}
	est <- new("sortableEstimate", estimate.vec, LocationScale = LocationScale)
	names(est) <- nam
	est
}
sortableEstimate <- function(x, y, estimator, ...){printGeneric("sortableEstimate")}#removeMethods("sortableEstimate")
setMethod("sortableEstimate", signature(x = "ANY", y = "missing", estimator = "Estimator"), function(x, y, estimator, ...)
{
	estimate.vec <- unsortableEstimate(x = x, estimator = estimator, ...)
	LocationScale <- estimate(x)
	new_sortableEstimate(estimate.vec, LocationScale = LocationScale)
})
setMethod("sortableEstimate", signature(x = "ANY", y = "ANY", estimator = "Estimator"), function(x, y, estimator, ...)
{
	estimate.vec <- unsortableEstimate(x = x, y = y, estimator = estimator, ...)
	LocationScale <- estimate(x, y)
	new_sortableEstimate(estimate.vec, LocationScale = LocationScale)
})
setMethod("sortableEstimate", signature(x = "ANY", y = "ANY", estimator = "ANY"), function(x, y, estimator, ...)
{
	stop("No estimator argument of class 'Estimator' was specified:")
})
#
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


#----



setAs(from = "Fdr", to = "basicFdr", function(from)
{
	new_basicFdr(fdr = fdr(from), p.value = pvalue(from))
})
setMethod("length", signature(x = "basicFdr"), function(x)
{
	length(names(x))
})
setMethod("names", signature(x = "basicFdr"), function(x)
{
	names(x@p.value)
})
setMethod("[", signature(x = "basicFdr", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@fdr <- (x@fdr)[i]
  x@p.value <- (x@p.value)[i]
  x
})
setMethod("nannotation", signature(object = "basicFdr"), function(object)
{
	"basic"
})

#

setMethod("length", signature(x = "eFdr"), function(x)
{
	length(names(x))
})
setMethod("names", signature(x = "eFdr"), function(x)
{
	names(x@zz) # names(as(x, "Numeric"))
})
setMethod("[", signature(x = "eFdr", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@fdr <- (x@fdr)[i]
  x@zz <- (x@zz)[i] # error fixed 27 October 2007
  x
})
setMethod("nannotation", signature(object = "eFdr"), function(object)
{
	"locfdr"
})

# bFdr added 15 Feb. 2008:

setMethod("length", signature(x = "bFdr"), function(x)
{
	length(names(x))
})
setMethod("names", signature(x = "bFdr"), function(x)
{
	names(x@p.value) # names(as(x, "Numeric"))
})
setMethod("[", signature(x = "bFdr", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@fdr <- (x@fdr)[i]
  x@p.value <- (x@p.value)[i]
  x
})
setMethod("nannotation", signature(object = "bFdr"), function(object)
{
	"empiricalBayes"
})

#
setMethod("length", signature(x = "pseudoFdr"), function(x)
{
	length(names(x))
})
setMethod("names", signature(x = "pseudoFdr"), function(x)
{
	names(x@estimate)
})
setMethod("[", signature(x = "pseudoFdr", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@fdr <- (x@fdr)[i]
  x@estimate <- (x@estimate)[i]
  x
})
setMethod("nannotation", signature(object = "pseudoFdr"), function(object)
{
	"fold-change"
})

#
setMethod("length", signature(x = "sFdr"), function(x)
{
	length(names(x))
})
setMethod("names", signature(x = "sFdr"), function(x)
{
	names(x@p.value) # names(as(x, "Numeric"))
})
setMethod("[", signature(x = "sFdr", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@fdr <- (x@fdr)[i]
  x@p.value <- (x@p.value)[i]
  x@t <- (x@t)[i]
  x
})
setMethod("nannotation", signature(object = "sFdr"), function(object)
{
	"limma"
})

#
setMethod("[", signature(x = "Prob0", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@.Data <- as(x, "Numeric")[i]
  x@Scale0 <- x@Scale0[i]
  x@Scale1 <- x@Scale1[i]
  x@Location1 <- x@Location1[i]
  x@xSize <- x@xSize[i]
  x@ySize <- x@ySize[i]
  x@PValue <- x@PValue[i]
#	stop('"[", signature(x = "Prob0", i = "ANY", j = "missing") not yet implemented')
  x
})
setMethod("names", signature(x = "Prob0"), function(x)
{
	names(x@Location1)
})
setAs(from = "Prob0", to = "numeric", function(from)
{
	vec <- from@.Data
	names(vec) <- names(from)
	vec
})
setAs(from = "Prob0", to = "Numeric", function(from)
{
	vec <- as(from, "numeric")
	nNumeric(vec)
})
setAs(from = "Prob0", to = "Fdr", function(from) # new 24 April 2008
{
	Num <- from@.Data
	if(is(Num, "Fdr"))
		Num
	else
	{
		warning("Prob0 may have lost Fdr info. such as eFdr's zz slot")
		new_basicFdr(fdr = Num, p.value = from@PValue)
	}
})
setMethod("nannotation", signature(object = "Prob0"), function(object){object@annotation})
#
#setAs(from = "loss", to = "predictionError", function(from)
#{
#	mat <- as.matrix(data.frame(from@loss0, from@lossBayes, from@loss1))
#	colnames(mat) <- c("null hypotheses", "empirical Bayes", "unbiased estimation")
# 	pE <- try(new("predictionError", mat, parameter = colnames(mat), parameter_lab = character(0)))
# 	if(is(pE, "try-error"))
# 	{ message("pE err"); browser()}
# 	pE
#})
setMethod("nannotation", signature(object = "loss"), function(object){nannotation(object@P0)})
setMethod("[", signature(x = "loss", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@lossBayes <- x@lossBayes[i]
  x@loss0 <- x@loss0[i]
  x@loss1 <- x@loss1[i]
  x@P0 <- x@P0[i]
  x
})
setAs(from = "losses", to = "character", function(from)
{
	get.main <- function(x)
	{
		if(x == "t.test")
			"Student t test"
		else if(x == "wilcox.test")
			"Wilcoxon test"
		else if(x == "moderated.t.test")
			"moderated t test; linear model" # Smyth (2004)
		else if(x == "moderated.t.test.eFdr0")
			"moderated t test; theoretical null" # Efron (2004)
		else
			x
	}
	sapply(from@test.names, get.main)
})

#
pvalue <- function(object, ...){printGeneric("pvalue")}#browser # in genetics.r prior to 12 February 2008
#removeMethods("pvalue")
setMethod("pvalue", signature(object = "numeric"), function(object)
{
	if(is(object, "Numeric"))
		object
	else
		nNumeric(object)
})
setMethod("pvalue", signature(object = "list"), function(object)
{
	pvalue(object$p.value)
})
#setMethod("pvalue", signature(object = "biasEstimate"), function(object, ...)
#{
#	pvalue.via.zvalue(object = object, ...)
#})
#
zvalue <- function(object, type, ...){printGeneric("zvalue")}#browser # new 24 April 2008
#removeMethods("zvalue")
setMethod("zvalue", signature(object = "Numeric", type = "missing"), function(object, type)
{
	qnorm(as(object, "numeric"))
})
#
change.pvalue <- function(object, old.alternative, new_alternative) # added 7 May 2008; see also "pvalue", signature(object = "eFdr"); low-level
{
	stopifnot(is(object, "numeric") && all(object >= 0 & object <= 1, na.rm = TRUE))
	if(old.alternative == "less")
	{
		if(new_alternative == "less")
			object
		else if(new_alternative == "greater")
		{
			warning('consider using change.zvalue based on "pvalue", signature(object = "eFdr") for less p-value roundoff error')
			1 - object
		}
		else if(new_alternative == "two.sided")
			2 * pmin(object, 1 - object)
		else
			stop("new_alternative not recognized")
	}
	else
		stop(paste("change.pvalue not implemented for old.alternative ", old.alternative, sep = ""))
}
change.zvalue <- function(object, old.alternative, new_alternative) # added 9 May 2008; see also "pvalue", signature(object = "eFdr")
{
	stopifnot(is(object, "numeric"))
	if(old.alternative == "less")
	{
		if(new_alternative == "less")
			object
		else if(new_alternative == "greater")
			-object
		else if(new_alternative == "two.sided")
			pmin(object, -object)
		else
			stop("new_alternative not recognized")
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


#setMethod("zvalue", signature(object = "biasEstimate", type = "missing"), function(object, type, ...)
#{
#	uncorrected <- object@uncorrected
#	class.ok <- is(uncorrected, "unsortableEstimate")
#	if(!class.ok)
#		stop(paste("qvaluen ", class(object), " not implemented for uncorrected of class ", class(uncorrected), sep = ""))
#	type <- uncorrected@estimator@type
#	zvalue(object = object, type = type, ...)
#})
#setMethod("zvalue", signature(object = "biasEstimate", type = "character"), function(object, type, correct, alternative, ...)
#{
#	stopifnot(is.logical(correct))
#	uncorrected <- object@uncorrected
#	z <- if(type %in% c("pvalue.z", "fdr_z"))
#	{
#		if(correct)
#			uncorrected - object # error fixed 9 May 2008 4 pm
#		else
#			uncorrected
#	}
#	else
#		stop(paste("qvaluen ", class(object), " not implemented for type ", type, sep = ""))
#	old.alternative <- "less"
#	if(missing(alternative))
#	{
#		alternative <- old.alternative
#		message("zvalue ", class(object), " alternative = ", alternative)
#	}
##	p.less <- pnorm(z)
##	change.pvalue(object = p.less, old.alternative = old.alternative, new_alternative = alternative)
#	change.zvalue(object = z, old.alternative = old.alternative, new_alternative = alternative)
#})

#=======================
