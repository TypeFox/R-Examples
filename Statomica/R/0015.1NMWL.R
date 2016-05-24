# NMWL.r created by David Bickel on 11 September 2010 by moving material from MDL.r. 

#Source("MDL.r")

setClass("DensityFunUNML", representation("function", weight = "Weight"))
setClass("Complexity", representation("Scalar", weight = "Weight"))

setClass("DensityFunNML", representation("function", family = "Family", COMP = "Complexity"))
setValidity("DensityFunNML", function(object)
{
	len.ok <- length(object@COMP) == 1
	ok <- len.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("Plot", signature(x = "DensityFunNML", y = "missing"), function(x, y, xlim, main, ...)
{
	fam <- x@family
	if(missing(xlim))
		xlim <- default(c(max(fam@lower.unknown.param, 1e-6), fam@upper.unknown.param), "xlim")
	if(missing(main))
	{
		main <- if(length(fam@known.param) == 1)
		{
			val <- fam@known.param
			nam <- names(fam@known.param)
			if(nam == "df")
			{
				nam <- if(val == 1) "degree of freedom" else "degrees of freedom"
			}
			default(paste(val, nam), "main")
		}
		else
			"optimal predictive density"
	}
	plot(vectorize(x), xlab = "statistic (reduced data)", ylab = "predictive density", xlim = xlim, main = main, ...)
})
setMethod("plot", signature(x = "DensityFunNML", y = "missing"), function(x, y, param0 = stop("param0 was not specified"), x.converges = FALSE, xlim = stop("xlim was not specified"), base = 2, xlab, ...)
{
	message("consider calling Plot instead")
	lf <- likFun(object = x, param0 = param0, base = base) # not really a pseudo-likelihood function, but this is used for plotting purposes
	if(missing(xlab))
	{
		xlab <- if(x.converges)
		{
			if(identical(x@family@Distr.fun, Td))
				"sample coefficient of variation"
			else
				"consistent estimate" # "sample coefficient of variation"
		}
		else # if(stat.name == "sample.cv")
			"statistic"
	}
#	else
#		stat.name
	test <- try(c(lf(x = xlim[1], x.converges = x.converges), lf(x = xlim[2], x.converges = x.converges)))
	if(is.err(test))
	{ message("cannot plot ", class(x), ); browser()}
	print(test)
	plot(lf, xlim = xlim, xlab = xlab, x.converges = x.converges, ...)
	if(x.converges && identical(x@family@Distr.fun, Td))
		abline(v = c(-1, 1), col = "gray")
})
setClassUnion("familyContaining", c("Family", "DensityFunNML"))
setMethod("family", signature(object = "DensityFunNML"), function(object)
{
	family(object@family)
})

setClass("DensityFunNMLs", representation("list", statistic.fun = "function", reduced.x = "statistic", weight = "Weight", family = "Family"))
setValidity("DensityFunNMLs", function(object)
{
	len.ok <- length(object) >= 1 && length(object@reduced.x) == length(object)
	nam.ok <- len.ok && sameNames(object@reduced.x, object)
	cla.ok <- are(object, "DensityFunNML")
	weight.ok <- all(object@weight@incidental.weight == object@weight@incidental.weight[1]) # same weights for all incidental features or comparisons
	ok <- cla.ok && len.ok && nam.ok && weight.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("family", signature(object = "DensityFunNML"), function(object)
{
	family(object@family)
})
setAs(from = "DensityFunNMLs", to = "lik.ratio", function(from)
{
	lik.ratio(object1 = from, object2 = 0)
})


setAs(from = "xprnSet", to = "Weight", function(from)
{
	Weight(x = from)
})
setAs(from = "xprnSetPair", to = "Weight", function(from)
{
	Weight(x = from)
})
setAs(from = "DensityFunUNML", to = "Weight", function(from)
{
	from@weight
})
setAs(from = "Complexity", to = "Weight", function(from)
{
	from@weight
})
setAs(from = "DensityFunNML", to = "Weight", function(from)
{
	as(from@COMP, "Weight")
})


# near end of file:

#Source(file = "NMWL.s") # functions

