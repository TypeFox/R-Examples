# Created by David Bickel on 17 August 2007.
# "Rank", signature(object = "numeric") changed and "sort", signature(x = "biasEstimate", decreasing = "ANY") added 8 April 2008.
#  More changes 16 April 2008, 24 April 2008, and later.



print_stats <- function(object, name, ...){message(name, ":"); print(stats(object, ...))}
stats <- function(object, ...)
{
	if(is(object, "numeric"))
		object <- as(object, "numeric")
	summ <- try(summary(object, ...))
	if(is.err(summ))
		summ <- c(mean = mean(object), median = median(object), n.na = sum(is.na(object)))
	if(is(object, "numeric"))
	{
		nam <- c("n.finite", "size", "mode", names(summ))
		vec <- c(n.finite = sum(is.finite(object)), size = length(object), mode = hsm(object), as.numeric(summ))
		names(vec) <- nam
		vec
	}
	else
		summ

}

load.file <- function(file, ...)
{
	refresh()#XXX|:no visible global function definition for ÔrefreshÕ
	paste(file, "RData", sep = ".")
}
saveo <- function(..., file, use.elem1 = FALSE)
{
	refresh()
	ext <- ".RData"
	if(missing(file))
	{
		warning("file argument must be explicit to be used")
		arglis <- list(...)
		cla <- sapply(arglis, class)
		boo <- cla == "character" & sapply(arglis, length) == 1
		if(FALSE) # any(boo))
		{
			whi <- which(boo)[1]
			lis <- try(c(arglis, list(file = arglis[[whi]]))) # [-whi]
			if(is(lis, "try-error"))
			{ message("bad lis"); browser()}
			dc <- try(do.call("saveo", lis))
			if(is(dc, "try-error"))
			{ message("bad dc"); browser()}
			return()
		}
		else
			file <- paste(Sys.Date(), " ", round(runif(1) * 1000), ext, sep = "")
	}
	sav <- function(...)
	{
		save(...)
	}
	if(is.character(file) && length(file) > 1)
	{
		warning("using only the first element of file")
		saveo(file = file[1], ...)
	}
	else if(is.character(file) && length(file) == 1)
	{
		nc <- nchar(file)
		exten <- substr(file, start = nc - 5, stop = nc)
		if(exten != ".RData")
			file <- paste(file, ext, sep = "")
		message("saving ", file, " on ", date())
		arglis <- list(...)
		if(use.elem1 && length(arglis) != 1)
			stop("bad num args")
		if(!use.elem1)
			sa <- try(sav(file = file, ...), silent = TRUE)
		if(use.elem1 || is(sa, "try-error"))
		{
			if(length(arglis) == 1)
			{
				elem1 <- arglis[[1]]
				s <- try(save(file = file, elem1))
				if(is(s, "try-error"))
				{ message("cannot save first element of arglis"); browser()}
			}
			else
			{
				message("cannot save to ", file)
				browser()
			}
		}
	}
	else
	{
		warning(paste(file, "is not a file name"))
		sav(file, ...)
	}
}
saveh <- function()
{
	refresh()
	file <- paste(Sys.Date(), "Rhistory", sep = ".")
	message("saving ", file, " in ", pwd(), " on ", date())
	savehistory(file = file)
}
savei <- function(file.base = "", ...)
{
	refresh()
	file <- paste(file.base, "RData", sep = ".")
	message("saving ", file, " in ", pwd(), " on ", date())
	save.image(file = file, ...)
}
savej <- function(...)
{
	saveh()
	savei(...)
	getwd()
}
pwd <- function(..., n = 20)
{
	path <- getwd(...)
	nc <- nchar(path)
	substr(path, max(1, nc - n + 1), nchar(path))
}

assign.text <- function(x.name, value.name)
{
	paste(x.name, value.name, sep = " <- ")
}
get.expr <- function(text, verbose = TRUE)
{
	text <- as.character(text)
	if(verbose)
		message("\n > ", text, "\n    (", date(), ")")
	parse(text = text)
}
assign.expr <- function(..., verbose = TRUE){get.expr(text = assign.text(...), verbose = verbose)} #parse(text = assign.text(...))}
Assign <- function(x.name, value, verbose = FALSE)
{
	assert.is(x.name, "character")
	if(verbose)
		message("\n", x.name, " <- object of class ", class(value))
	sa <- try(saveo(value, file = x.name, use.elem1 = TRUE))
	if(is(sa, "try-error"))
	{ message("cannot save ", class(value), " to ", x.name); browser()}
	invisible(sa)
}
Get <- function(x.name, verbose = FALSE)
{
	assert.is(x.name, "character")
#	try(rm(elem1), silent = TRUE)
	elem1<-load(file = load.file(x.name))#XXX|elem1<-load(file = load.file(x.name))  ????
	#before April2014:load(file = load.file(x.name)) #XXX|no visible binding for global variable Ôelem1Õ
	if(verbose)
		message("\n <- ", x.name, " # object of class ", class(elem1))
	
	elem1# assumes Assign was used to save file
}

objectSize <- function(x, ...){printGeneric("objectSize"); browser()}
#removeMethods("objectSize")
setMethod("objectSize", signature(x = "character"), function(x, i, verbose)
{
	if(missing(verbose))
	{
		verbose <- TRUE
		message("objectSize verbose = ", verbose)
	}
	nam <- x
	if(missing(i))
	{
		text <- paste("i = ", length(nam) - 9, ":", length(nam), sep = "")
		message(text)
		eval(parse(text = text))
	}
	sizes <- sapply(nam, function(nam)
	{
		eval(parse(text = paste("object.size(", nam, ")", sep = "")))
	}) / 1e6
	names(sizes) <- nam
	ord <- order(sizes)
	if(verbose)
	{
#		args.with.commas <- if(length(i) == 0)
#			character(0)
#		else
#			paste(nam[ord][i], ", ", sep = "")
		args.with.commas <- paste(nam[ord][i], collapse = ", ")
		message("For copy and paste: rm(", args.with.commas, ")")
		message(sum(sizes), " MBs of ", length(nam), " objects")
	}
	datf <- data.frame(object = nam, MB = sizes, index = 1:length(sizes))[ord, ]
	if(verbose)
		print(datf[i, ])
	invisible(datf)
})
global.environment <- environment()
setMethod("objectSize", signature(x = "missing"), function(x, name = global.environment, ...)
{
	nam <- objects(name = name)
	objectSize(x = nam, ...)
})
setClassUnion("index", c("logical", "numeric", "character"))
# end

plot_pch <- function()
{
	Pex <- 3 ## good for both .Device=="postscript" and "nx11"
	ipch <- 0:35; np <- length(ipch); k <- floor(sqrt(np)); dd <- c(-1,1)/2
	rx <- dd + range(ix <- ipch %/% k)
	ry <- dd + range(iy <- 3 + (k-1)- ipch %% k)
	pch <- as.list(ipch)
	pch[26+ 1:10] <- as.list(c("*",".", "o","O","0","+","-","|","%","#"))
	plot(rx, ry, type="n", axes = FALSE, xlab = "", ylab = "",
			 main = paste("plot symbols :  points (...  pch = *, cex =", Pex,")"))
	abline(v = ix, h = iy, col = "lightgray", lty = "dotted")
	for(i in 1:np) {
		pc <- pch[[i]]
		points(ix[i], iy[i], pch = pc, col = "red", bg = "yellow", cex = Pex)
		## red symbols with a yellow interior (where available)
		text(ix[i] - .3, iy[i], pc, col = "brown", cex = 1.2)
	}
}

Stripchart <- function(x, ...)
{
	stripchart(x = x, ...)
}
#removeMethods("Stripchart")
setMethod("Stripchart", signature(x = "matrix"), function(x, group.names, pch, ylab, circle, triangle, log, text.x, ...)
{
	if(missing(log))
	{
		log <- ""
		message("log = ", log)
	}
	if(missing(ylab)) ylab <- ""
	stopifnot(length(pch) == nrow(x))
	if(missing(circle) || length(circle) == 0)
	{
		circle <- numeric(0)
		message("matrix circle = ", circle)
	}
	if(is.character(circle) && is.character(pch))
	{
		stopifnot(circle %in% pch)
		circle <- which(circle == pch)
	}
	if(missing(triangle) || length(triangle) == 0)
	{
		triangle <- circle
		message("matrix triangle = ", triangle)
	}
	if(is.character(triangle) && is.character(pch))
	{
		stopifnot(triangle %in% pch)
		triangle <- which(triangle == pch)
	}
	shape.ok <- function(shape){length(shape) == 0 || shape %in% 1:length(pch)}
	stopifnot(shape.ok(circle) && shape.ok(triangle))
	y <- rep(ncol(x):1, nrow(x))
	xvec <- as.numeric(x)
	plot(x = xvec, y = y, type = "n", ylim = c(1, max(y) + 0.7), yaxt = "n", ylab = ylab, log = log, ...)
	stopifnot(length(group.names) == ncol(x))
	for(j in 1:ncol(x))
	{
		xv <- x[, j]
		yv <- rep(y[j], length(xv))
		points(x = xv, y = yv, pch = pch, ...)
		shape.points <- function(shape, pch)
		{
			if(length(shape) != 0)
				points(x = xv[shape], y = yv[shape], pch = pch, cex = 3)
		}
		shape.points(circle, pch = 1)
		shape.points(triangle, pch = 0)
#		points(x = xv[triangle], y = yv[triangle], pch = 0, cex = 3)
#		ti <- ceiling(length(xvec) / 2)
		midpoint <- function(vec)
		{
			ra <- range(vec, na.rm = TRUE)
			if(log != "")
				ra <- log(ra)
			mp <- mean(ra)
			if(log != "")
				mp <- exp(mp)
			mp
		}
		if(missing(text.x) || length(text.x) == 0)
		{
			text.x <- midpoint(vec = xvec)
			message(class(x), " text.x = ", text.x)
		}
		text(x = text.x, y = yv[1] + 0.3, labels = group.names[j])
	}
})
setMethod("Stripchart", signature(x = "list"), function(x, pch, ...)
{
	if(missing(pch) || length(pch) == 1 || any(length(pch) != sapply(x, length)))
	{
		message("Stripchart calling stripchart; press c to continue")
		browser()
		stripchart(x = x, pch = pch, ...)
	}
	else
	{
		mat <- try(sapply(x, as.numeric))
		if(is.matrix(mat))
		{
			message("Stripchart (list) calling Stripchart (matrix)")
			tr <- try(Stripchart(x = mat, pch = pch, ...))
			if(is(tr, "try-error"))
			{ message("Stripchart (matrix) failed"); browser()}
			tr
		}
		else
		{
			message("mat is not a matrix")
			browser()
		}
	}
})

scalar <- function(object, ...){message("scalar generic"); browser()}
#removeMethods("scalar")
setMethod("scalar", signature(object = "scalar"), function(object)
{
	as(object, "scalar")
})
setMethod("scalar", signature(object = "ANY"), function(object)
{
	scalar(as.numeric(object))
})
setMethod("scalar", signature(object = "numeric"), function(object) # changed 16 April 2008
{
#  if(length(object) == 0)
#    stop("cannot convert vector of length 0 to scalar")
#  else
  if(length(object) > 1)
  {
#    message("extra vector elements will be lost in coersion of object to scalar")
    print(stats(object))
    stop("extra vector elements prevents coersion of object to scalar")
  	object <- object[1]
  }
  new("scalar", object)
})

Scalar <- function(object, ...){message("Scalar generic"); browser()}
#removeMethods("Scalar")
setMethod("Scalar", signature(object = "ANY"), function(object) # changed 16 April 2008
{
#  if(length(object) == 0)
#    stop("cannot convert vector of length 0 to Scalar")
#  else if(length(object) > 1)
#    warning("extra vector elements lost in coersion to Scalar")
  Sca <- scalar(object)
  if(length(Sca) == 1 && !is.na(Sca) && Sca < 0)
    stop("cannot convert negative number to Scalar")
  new("Scalar", Sca)
})

pNumeric <- function(object, ...){message('pNumeric generic'); browser()}
#removeMethods("pNumeric")
setMethod("pNumeric", signature(object = "numeric"), function(object)
{
  x <- new("pNumeric", as.numeric(object))
  names(x) <- names(x@.Data) <- names(object)
  x
})
Numeric <- function(object, ...){message('Numeric generic'); browser()}
#removeMethods("Numeric")
setMethod("Numeric", signature(object = "numeric"), function(object)
{
  x <- new("Numeric", as.numeric(object))
  names(x) <- names(x@.Data) <- names(object)
  x
})
Matrix <- function(object, ...){message('Matrix generic'); browser()}
#removeMethods("Matrix")
setMethod("Matrix", signature(object = "matrix"), function(object)
{
  new("Matrix", object)
})
setMethod("Matrix", signature(object = "numeric"), function(object)
{
  Numeric(object) # [sic]
})

Aggregate <- function(object, ...){printGeneric("Aggregate"); browser()}
#removeMethods("Aggregate")
setMethod("Aggregate", signature(object = "data.frame"), function(object, by, FUN, ...)
{
	paste.dim <- function(x)
	{
		main.str <- paste(class(x), "of", nrow(x), "rows and", ncol(x), "columns")
		numeric.str <- paste("(", sum(sapply(x, is.numeric)), " of which are numeric)", sep = "")
		if(is(x, "data.frame"))
			paste(main.str, numeric.str)
		else
			main.str
	}
	dim0 <- paste.dim(object)
	get.fac.boo <- function(x)
	{
		boo <- sapply(x, is.factor)
		n.fac <- sum(boo)
		ok <- n.fac < ncol(x)
		if(!ok)
		{ message("fac.boo error"); browser()}
		boo
	}
	if(missing(by))
	{
		fac.boo <- get.fac.boo(x = object)
		if(sum(fac.boo) == 0)
		{
			message("aggregate not called because 'by' not specified and could not be extracted from object")
			return(object)
		}
		by <- object[, fac.boo]
		by <- if(is.data.frame(by))
			as.list(by)
		else if(is.factor(by))
			list(by)
		else
			stop("bad by")
		names(by) <- names(object)[fac.boo]
		object <- object[, !fac.boo]
	}
	if(missing(FUN)) FUN <- function(x){mean(x, na.rm = TRUE)}
	datf <- aggregate(x = object, by = by, FUN = FUN, ...)
	datf.fac.boo <- get.fac.boo(datf)
	if(sum(datf.fac.boo) >= 1)
	{
		get.rownames <- function()
		{
			paste.names <- function(...){paste(..., sep = ".")}
			nam.lis <- lapply(datf[, datf.fac.boo, drop = FALSE], as.character)
			ok <- all(sapply(nam.lis, is.character))
			if(!ok)
			{ message("get.rownames error"); browser()}
			nam <- do.call("paste.names", nam.lis)
			if(length(nam) != nrow(datf))
			{ message("nam of incorrect length"); browser()}
			nam
		}
		rownames(datf) <- get.rownames()
		datf <- datf[, !datf.fac.boo]
	}
	if(all(sapply(datf, is.numeric)))
		datf <- as.matrix(datf)
	message(dim0, " converted to ", paste.dim(datf))
	datf
})

xprnSubset <- function(object, ...){printGeneric("xprnSubset"); browser()}
#removeMethods("xprnSubset")
setMethod("xprnSubset", signature(object = "ExpressionSet"), function(object, level, factor.name, ...)
{
	if(missing(factor.name))
	{
		fac.boo <- sapply(pData(object), is.factor)
		if(!any(fac.boo))
			stop("there are no factors in pData(object)")
		factor.name <- names(pData(object))[fac.boo][1]
		message("factor.name = ", factor.name)
	}
	fac <- pData(object)[, factor.name] # [,  == "MyoT"]
	stopifnot(is.factor(fac))
	if(!all(level %in% as.character(fac)))
	{ message("cannot complete xprnSubset"); browser()}
	boo <- fac %in% level
	stopifnot(sum(boo) >= 1 && length(boo) == ncol(exprs(object)))
	es <- object[, boo]
	ann <- paste(level, collapse = "&")
	annotation(es) <- if(length(annotation(es)) == 1)
		paste(ann, annotation(es), sep = " of ")
	else
		ann
	es
})
setMethod("xprnSubset", signature(object = "xprnSet"), function(object, ...)
{
	object@es <- xprnSubset(object = object@es, ...)
	object
})

unpair <- function(object, ...){printGeneric("unpair"); browser()}
#removeMethods("unpair")
setMethod("unpair", signature(object = "xprnSet"), function(object, factor.name, pairing.name, change.sign, ...)
{
	if(missing(change.sign))
	{
		change.sign <- FALSE
		message("change.sign = ", change.sign)
	}
	datf <- pData(object)
	fac <- datf[, factor.name]
	stopifnot(is.factor(fac))
	levs <- levels(fac)
	stopifnot(length(levs) == 2)
	if(change.sign)
		levs <- rev(levs)
	get.pairing <- function(es)
	{
		pData(es)[, pairing.name]
	}
	get.sub <- function(level)
	{
		es <- xprnSubset(object = object, level = level, factor.name = factor.name)
		pairing <- get.pairing(es = es)
		if(any(duplicated(pairing)))
		{ message("bad pairing"); browser()}
		ord <- order(pairing)
		stopifnot(length(ord) == ncol(exprs(es)))
		es[, ord]
	}
	x <- get.sub(level = levs[1])
	y <- get.sub(level = levs[2])
	stopifnot(all(get.pairing(x) == get.pairing(y)))
	es <- x
	exprs(es@es) <- exprs(x) - exprs(y)
	pData(es@es)[, factor.name] <- factor(paste(levs, collapse = ".minus."))
	colnames(exprs(es@es)) <- rownames(pData(es@es)) <- paste(rownames(pData(x)), rownames(pData(y)), sep = ".")
	es
})
setMethod("unpair", signature(object = "XprnSet"), function(object, ...)
{
	unpair(logb(object), ...)
})

removeMissing <- function(object, ...){printGeneric("removeMissing"); browser()}
#removeMethods("removeMissing")

xprnSet <- function(phenoData, exprs, featureData,...){message("xprnSet generic"); browser()}
#removeMethods("xprnSet")
setMethod("xprnSet", signature(phenoData = "AnnotatedDataFrame", exprs = "matrix",featureData = "missing"), function(phenoData, exprs, ...)
{
  ed <- new("ExpressionSet", phenoData = phenoData, exprs = exprs, ...)
  new("xprnSet", es = ed)
})
setMethod("xprnSet", signature(phenoData = "AnnotatedDataFrame", exprs = "matrix",featureData = "AnnotatedDataFrame"), function(phenoData, exprs,featureData, ...)
{
  ed <- new("ExpressionSet", phenoData = phenoData, exprs = exprs, featureData = featureData, ...)
  new("xprnSet", es = ed)
})
setMethod("xprnSet", signature(phenoData = "data.frame", exprs = "matrix",featureData = "missing"), function(phenoData, exprs, ...)
{
	stopifnot(ncol(exprs) == nrow(phenoData) && all(colnames(exprs) == rownames(phenoData)))
  xprnSet(phenoData = as(phenoData, "AnnotatedDataFrame"), exprs = exprs, ...)
})
setMethod("xprnSet", signature(phenoData = "missing", exprs = "matrix",featureData = "missing"), function(phenoData, exprs, ...)
{
	if(is.null(rownames(exprs)))
		rownames(exprs) <- make.names(1:nrow(exprs))
	if(is.null(colnames(exprs)))
		colnames(exprs) <- make.names(1:ncol(exprs))
	phenoData <- data.frame(dummy = factor(rep(0, ncol(exprs))))
	rownames(phenoData) <- colnames(exprs)
	xprnSet(phenoData = phenoData, exprs = exprs, ...)
})
setMethod("xprnSet", signature(phenoData = "missing", exprs = "numeric",featureData = "missing"), function(phenoData, exprs, nrow, ...)
{
	if(missing(nrow))
		nrow <- 1
	xprnSet(exprs = matrix(exprs, nrow = nrow), ...)
})
setMethod("xprnSet", signature(phenoData = "missing", exprs = "xprnSet",featureData = "missing"), function(phenoData, exprs, ...)
{
	mat <- exprs(exprs)
	feature.okay = sapply(1:nrow(mat), function(i)
	{
		boo <- try(!any(is.na(mat[i, ])))
		if(is.logical(boo) && length(boo) == 1)
			boo
		else
		{ message("boo error"); browser()}
	})
	stopifnot(any(feature.okay))
	exprs[feature.okay, ]
})
#setMethod("xprnSet", signature(phenoData = "ANY", exprs = "missing",featureData = "missing"), function(phenoData, exprs, ...)
#{
#	xprnSet(exprs = phenoData, ...)
#})

as.xprnSet <- function(object, ...){printGeneric("as.xprnSet"); browser()}
#removeMethods("as.xprnSet")
setMethod("as.xprnSet", signature(object = "data.frame"), function(object, samples1.j, samples2.j, translation = 0, ...)
{
	exprs <- delta.ln.from.datf(x = object, samples1.j = samples1.j, samples2.j = samples2.j, translation = translation)
	repnames <- paste(as.character(samples1.j), as.character(samples2.j), sep = ".")
	stopifnot(all(length(repnames) == c(length(samples1.j), length(samples2.j))))
	phenoData <- data.frame(as.factor(repnames))
	stopifnot(ncol(phenoData) == 1)
	names(phenoData) <- "replicate"
	stopifnot(length(repnames) == nrow(phenoData))
	rownames(phenoData) <- repnames
	stopifnot(length(repnames) == ncol(exprs))
	colnames(exprs) <- repnames
	xprnSet(phenoData = phenoData, exprs = exprs)
})

setMethod("removeMissing", signature(object = "xprnSet"), function(object)
{
  missing.gene <- as.logical(apply(exprs(object), 1, function(ro){all(is.na(ro))}))
  stopifnot(length(missing.gene) == nrow(exprs(object)))
  if(all(missing.gene))
  {
  	message("all genes are missing in ", class(object))
  	browser()
  }
  if(any(missing.gene))
  	object[!missing.gene]
  else
	  object
})

RepNames <- function(object, times, unique = TRUE, ...)
{
	nam0 <- if(is.character(names(object)))
		names(object)
	else if(is(object, "matrix") || is(object, "data.frame"))
	{
		rownam0 <- rownames(object)
		if(is.character(rownam0))
			rownam0
		else
			1:nrow(object)
	}
	else if(length(object) >= 1)
		1:length(object)
	else
		stop("bad RepNames object")
#	if(is(object, "numeric") && length(object) == 1)
#		1:object
	assert.is(nam0, c("numeric", "character"))
	assert.is(times, "numeric")
	stopifnot(length(times) == 1 && times >= 1 && floor(times) == times)
	nam <- make.names(rep(nam0, times), unique = unique, ...)
	stopifnot(length(nam) == times * length(nam0))
	nam
}

Rep <- function(object, times, ...){stop(paste("bad class for Rep (", class(object), ") or times (", class(times), ")", sep = ""))}
setMethod("Rep", signature(object = "Vector"), function(object, times, ...)
{
	nam <- RepNames(object = object, times = times, ...)
	vec <- rep(object, times)
	stopifnot(length(vec) == length(nam))
	names(vec) <- nam
	vec
})
setMethod("Rep", signature(object = "matrix"), function(object, times, ...)
{
	nam <- RepNames(object = object, times = times, ...)
	mat <- matrix(rep(as.numeric(object), times), ncol = ncol(object), byrow = TRUE)
	stopifnot(nrow(mat) == times * nrow(object))
	stopifnot(nrow(mat) == length(nam))
	rownames(mat) <- nam
	colnames(mat) <- colnames(object)
	mat
})
setMethod("Rep", signature(object = "xprnSet"), function(object, times, ...)
{
	datf <- pData(object)
	stopifnot(is(datf, "data.frame"))
	mat <- Rep(object = exprs(object), times = times, ...)
	stopifnot(is.character(rownames(mat)) && !any(duplicated(rownames(mat))))
	es <- xprnSet(phenoData = datf, exprs = mat, annotation = "normal.rxprnSet")
	nam <- names(es)
	stopifnot(length(nam) == times * length(names(object)) && all(nam == names(es)))
	es
})
setMethod("Rep", signature(object = "Bernoulli.normal.rxprnSet"), function(object, times, ...)
{
	get.slot <- function(slot0){Rep(slot0, times = times, ...)}
	es <- get.slot(as(object, "xprnSet"))
	alternative <- get.slot(object@alternative)
	param <- get.slot(object@param)
	stopifnot(sameNames(es, alternative, param))
	new("Bernoulli.normal.rxprnSet", es, alternative = alternative, param = param)
})

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
new.Bernoulli.normal.rxprnSet <- function(object, alternative, param)
{
	add.nam <- function(x)
	{
		if(is.null(names(x)))
			names(x) <- names(object)
		x
	}
	new("Bernoulli.normal.rxprnSet", object, alternative = add.nam(alternative), param = add.nam(param))
}
Bernoulli.normal.rxprnSet <- function(size, nfeature, P0, mean0 = 0, mean1 = 2, sd0 = 1, sd1 = 1.5 * sd0, nfeatureIID = nfeature)
{
#	nfeatureActual <- nfeature
#	nfeature
#	if(nfeature != nfeatureIID)
	assert.are(list(size, nfeature, mean1, sd0, sd1, P0, nfeatureIID), "numeric")
	times <- nfeature / nfeatureIID
	stopifnot(times > 0 && floor(times) == times)
	rm(nfeature)
	stopifnot(all(lapply(list(size, nfeatureIID, P0, mean0, mean1, sd0, sd1), length) == 1))
	stopifnot(is.prob(P0))
	alternative <- rbinom(nfeatureIID, 1, 1 - P0) == 1
	if(mean0 != 0)
		stop("param1 implemented only for mean0 = 0")
	param1 <- mean1 * ifelse(rbinom(nfeatureIID, 1, 0.5) == 1, 1, -1)
	param <- ifelse(alternative, param1, mean0)
	names(param) <- make.names(paste("feature", 1:length(param), sep = "."))
	sds <- ifelse(alternative, sd1, sd0)
	es <- normal.rxprnSet(sample.size = size, means = param, sds = sds)
	stopifnot(Length(es) == nfeatureIID)
	res <- new.Bernoulli.normal.rxprnSet(es, alternative = alternative, param = param)
	Rep(object = res, times = times)
}

leave.out <- function(object, not.j, ...){printGeneric("leave.one.out"); browser()}
#removeMethods("leave.out")
setMethod("leave.out", signature(object = "xprnSet", not.j = "Scalar"), function(object, not.j, ...)
{
	nsample <- ncol(exprs(object))
	stopifnot(nsample > 1)
	stopifnot(not.j >= 1 && not.j <= nsample)
	j <- setdiff(1:nsample, not.j)
	es <- object[, j]
	prefix <- annotation(object)
	suffix <- paste("without sample", not.j)
	annotation(es@es) <- if(is.character(prefix))
		paste(prefix, suffix)
	else
		suffix
	es
})
setMethod("leave.out", signature(object = "xprnSet", not.j = "missing"), function(object, not.j, random, ...)
{
	if(missing(random))
	{
		random <- FALSE
		message('"leave.out", signature(object = "ExpressionSet", not.j = "missing") random = ', random)
	}
	nsample <- ncol(exprs(object))
	lo <- function(not.j){leave.out(object = object, not.j = Scalar(not.j), ...)}
	if(random)
		lo(not.j = sample(1:nsample, 1))
	else
	{
		lapply(1:nsample, function(not.j){lo(not.j = not.j)})
	}
})

XprnSet <- function(phenoData, exprs,featureData, ...){message("XprnSet generic"); browser()}
#removeMethods("XprnSet")
#setMethod("XprnSet", signature(phenoData = "ANY", exprs = "Matrix"), function(phenoData, exprs)
#{
#})
#setMethod("XprnSet", signature(phenoData = "ANY", exprs = "matrix",featureData = "missing"), function(phenoData, exprs, ...)
#{
#  if(!all(as.numeric(exprs) > 0, na.rm = TRUE))
#  { message("XprnSet requires all positive; consider xprnSet"); browser()}
##  XprnSet(phenoData = phenoData, exprs = Matrix(exprs))
#  new("XprnSet", xprnSet(phenoData = phenoData, exprs = exprs, ...))
#})


# D. R. Bickel, "Microarray gene expression analysis: Data transformation and multiple-comparison bootstrapping," Computing Science and Statistics 34, 383-400, Interface Foundation of North America (Proceedings of the 34th Symposium on the Interface, Montréal, Québec, Canada, April 17-20, 2002).
Transform.term.default <- function(object)
{
	if(is(object, "numeric"))
		median(object, na.rm = TRUE)
	else
		Transform.term.default(as(object, "numeric"))
}
Transform.number <- function(object, term = default(Transform.term.default(object), "Transform.number term"), ...)
{
	stopifnot(is(object, "numeric") || is(object, "matrix"))
	total <- object + term
	stopifnot(all(total > 0, na.rm = TRUE))
	logb(total)
}
Transform.xprnSet <- function(object, ...)
{
	assert.is(object, "xprnSet")
	if(FALSE) # is(object, "XprnSet"))
	{
		warning("taking the log before applying Transform, appropriate for converting ratios to generalized ratios")
		object <- logb(object)
	}
	es <- as(object, "ExpressionSet")
	exprs(es) <- Transform(exprs(es), ...)
	prefix <- "Transformed"
	annotation(es) <- if(length(annotation(es) == 0))
		prefix
	else
		paste(prefix, annotation(es))
	as(es, "xprnSet")
}
Transform <- function(object, ...)
{
	class.err <- function(){stop(paste("Transform not yet implemented for class", class(object)))}
	fun <- if(is(object, "xprnSet"))
		Transform.xprnSet
	else if(is(object, "numeric") || is(object, "matrix"))
		Transform.number
	else
		class.err()
	fun(object = object, ...)
}

xprnSetPair <- function(x, y, factor.name, ...){printGeneric("xprnSetPair"); browser()}
#removeMethods("xprnSetPair")
setMethod("xprnSetPair", signature(x = "xprnSet", y = "missing", factor.name = "character"), function(x, y, factor.name, level, verbose)
{
	if(missing(verbose))
		verbose <- FALSE
	fac <- pData(x)[, factor.name]
	if(!is.factor(fac)) stop("!is.factor(fac)")
	if(missing(level))
		level <- levels(unique((fac)))#as.character
	if(length(level) != 2) stop("length(level) != 2; try specifying different level argument.")
	get.subset <- function(level)
	{
		if(!all(level %in% as.character(fac)))
		{ message("cannot complete xprnSetPair"); browser()}
		if(verbose)
			message("calling xprnSubset with level = ", level, ", factor.name = ", factor.name)
		su <- xprnSubset(object = x, level = level, factor.name = factor.name)
		if(verbose)
			message("finished xprnSubset with level = ", level, ", factor.name = ", factor.name)
		su
	}
	x.sub <- get.subset(level = level[1])
	y.sub <- get.subset(level = level[2])
	xprnSetPair(x = x.sub, y = y.sub)
})
setMethod("xprnSetPair", signature(x = "xprnSet", y = "xprnSet", factor.name = "missing"), function(x, y, factor.name)
{
	if(is(x, "XprnSet") || is(y, "XprnSet"))
		stop("x, y inconsistency")
	else
		new("xprnSetPair", x = x, y = y)
})
setMethod("xprnSetPair", signature(x = "XprnSet", y = "XprnSet", factor.name = "missing"), function(x, y, factor.name)
{
	xprnSetPair(x = logb(x), y = logb(y))
})
setMethod("removeMissing", signature(object = "xprnSetPair"), function(object)
{
  x <- removeMissing(object@x)
  y <- removeMissing(object@y)
  nam <- intersect(featureNames(object@x), featureNames(y))
  if(!is.character(nam) || length(nam) == 0)
  { message("x and y have no non-missing genes in common"); browser()}
  object@x <- x[nam]
  object@y <- y[nam]
  if(!validObject(object))
  { message(class(object), " is no longer valid"); browser()}
  object
})

new.xprnSetObjectPair <- function(training, test){new("xprnSetObjectPair", training = training, test = test)} # new 24 April 2008
xprnSetObjectPair <- function(x, y, ...){printGeneric("xprnSetObjectPair"); browser()}
#removeMethods("xprnSetObjectPair")
setMethod("xprnSetObjectPair", signature(x = "xprnSetObject", y = "missing"), function(x, y, parametric, ...)
{
	if(missing(parametric))
	{
		parametric <- TRUE
		message("parametric = ", parametric)
	}
	get.xprnSetObject <- if(parametric)
	{
		function()
		{
			parametricBootstrap(object = x, ...)
		}
	}
	else
		function()
		{
			Sample(object = x, replace = TRUE, ...) # see nonparametricBootstrap
		}
	pair <- new.xprnSetObjectPair(training = get.xprnSetObject(), test = get.xprnSetObject())
	stopifnot(is(pair, "xprnSetObjectPair"))
	pair
})

new.xprnSetObjects <- function(object) # new 28 April 2008
{
	stopifnot(is.list(object))
	if(is.null(names(object)))
		names(object) <- sapply(object, annotation)
	x <- new("xprnSetObjects", object)
	names(x) <- names(object)
	x
}

normalize <- function(object, ...){printGeneric("normalize"); browser()}
#removeMethods("normalize")
setMethod("normalize", signature(object = "numeric"), function(object,  na.rm=T, location.fun, scale.fun,  ...)
{
	if(missing(location.fun))
		location.fun <- default(median, "location.fun {numeric}")
	if(missing(scale.fun))
		scale.fun <- default(mad, "scale.fun {numeric}")
	if(missing(na.rm)) na.rm <- TRUE
	me <- if(is.null(location.fun))
		0
	else if(is.function(location.fun))
		location.fun(object, na.rm = na.rm)
	else
		stop("bad location")
	stopifnot(is.finite(me))
	ma <- if(is.null(scale.fun))
		1
	else if(is.function(scale.fun))
		scale.fun(object, na.rm = na.rm)
	else
		stop("bad scale")
	stopifnot(ma > 0)
	(object - me) / ma
})
setMethod("normalize", signature(object = "matrix"), function(object, na.rm=T, location.fun, scale.fun, ...)
{
	if(missing(location.fun))
		location.fun <- default(median, "location.fun {matrix}")
	if(missing(scale.fun))
		scale.fun <- default(mad, "scale.fun {matrix}")
	mat <- apply(object, 2, normalize, location.fun = location.fun, scale.fun = scale.fun,na.rm=na.rm, ...)
	stopifnot(all(rownames(mat) == rownames(object)) && all(colnames(mat) == colnames(object)))
	mat
})
setMethod("normalize", signature(object = "ExpressionSet"), function(object, ...)
{
	exprs(object) <- normalize(exprs(object), ...)
	object
})
setMethod("normalize", signature(object = "xprnSet"), function(object, ...)
{
	object@es <- normalize(object = object@es, ...)
	object
})
setMethod("normalize", signature(object = "XprnSet"), function(object, ...)
{
	norm <- normalize(logb(object), ...)
	warning(paste("normalize converted ", class(object), " to ", class(norm), sep = ""))
#	stop(paste("normalize not yet implemeted for", class(object)))
	norm
})
setMethod("normalize", signature(object = "xprnSetPair"), function(object, ...)
{
	object@x <- normalize(object = object@x, ...)
	object@y <- normalize(object = object@y, ...)
	object
})

standardize <- function(..., location.fun = mean, scale.fun = sd)
{
	normalize(..., location.fun = location.fun, scale.fun = scale.fun)
}
recenter <- function(..., location.fun = mean)
{
	normalize(..., location.fun = location.fun, scale.fun = NULL)
}

identity <- function(x){x}



Order <- function(object, ...){message("generic Order"); browser()}
#removeMethods("Order")

nondecreasing <- function(object, ...)
{
	stopifnot(is.numeric(object))
	for(j in (length(object) - 1):1)
	{
		subsequent <- object[j + 1]
		current <- object[j]
		object[j] <- if(is.na(subsequent) && is.na(current))
			as.numeric(NA)
		else
			min(subsequent, current, na.rm = TRUE)
	}
	object
}

smoothSpline <- function(x, y, ...){message("generic smoothSpline"); browser()}
#removeMethods("smoothSpline")
setMethod("smoothSpline", signature(x = "numeric", y = "numeric"), function(x, y, df, nondecreasing, FUN, na.rm, ...)
{
	if(missing(na.rm))
		na.rm <- TRUE
	if(missing(nondecreasing))
	{
		nondecreasing <- FALSE
		message("nondecreasing = ", nondecreasing)
	}
	if(missing(df) && missing(FUN))
	{
		df <- 7
		message("df = ", df)
	}
	boo <- !is.na(x) & !is.na(y)
	if(!missing(df))
		stopifnot(sum(boo) >= df)
#	if(!na.rm)
#	{
#		all.x <- x
#		all.y <- y
#	}
	x <- x[boo]
	y <- y[boo]
	if(missing(FUN))
		FUN <- function(x, y, ...)
		{
			smooth.spline(x = x, y = y, df = df, ...)
		}
	ss <- FUN(x, y, ...)
	if(nondecreasing)
		ss$y <- nondecreasing(object = ss$y)
	if(!na.rm)
	{
		add.miss <- function(small)
		{
			tall <- as.numeric(rep(NA, length(boo)))
			if(length(small) != sum(boo))
				stop("add.miss error")
			tall[boo] <- small
			tall
		}
		ss$x <- add.miss(small = ss$x)
		ss$y <- add.miss(small = ss$y)
	}
	sS <- list(s3 = ss)
	new("smoothSpline", sS)
})
setMethod("smoothSpline", signature(x = "Density", y = "missing"), function(x, y, ...)
{
	smoothSpline(x = xi(x), y = eta(x), ...)
})

Lowess <- function(x, y, ...)#{printGeneric("Lowess"); browser()}
{
	if(missing(x)) stop("Lowess error")
	if(!is.numeric(x)) {message("Lowess error: x is not numeric"); browser()}
	FUN <- function(x, y, ...)
	{
		lowess(x = x, y = y, ...)
	}
	arglis <- list(x = x, FUN = FUN, ...)
	if(!missing(y)) arglis <- c(arglis, list(y = y))
	smoo <- do.call("smoothSpline", arglis)
	if(is(smoo, "smoothSpline"))
		as(smoo, "Lowess")
	else
		smoo
}
setMethod("Lowess", signature(x = "smoothSpline", y = "missing"), function(x, y, ...)
{
	arglis <- list(...)
	if(length(arglis) > 0)
	{ message("unused Lowess arguments:"); print(names(arglis)); browser()}
	as(x, "Lowess")
})
#setMethod("Lowess", signature(x = "ANY", y = "missing"), function(x, y, ...)
#{
#	Lowess(x = smoothSpline(x = x, ...))
#})
#setMethod("Lowess", signature(x = "ANY", y = "ANY"), function(x, y, ...)

movingAverage <- function(x, y, ...)#{printGeneric("movingAverage"); browser()}
{
	if(missing(x)) stop("movingAverage error")
	FUN <- function(x, y, ...)
	{
		s4 <- movingLocation(x = x, y = y, ...)
		list(x = xi(s4), y = eta(s4))
	}
	arglis <- list(x = x, FUN = FUN, ...)
	if(!missing(y)) arglis <- c(arglis, list(y = y))
	smoo <- do.call("smoothSpline", arglis)
	if(is(smoo, "smoothSpline"))
		movingAverage(x = smoo)
	else
		smoo
}
setMethod("movingAverage", signature(x = "smoothSpline", y = "missing"), function(x, y, ...)
{
	arglis <- list(...)
	if(length(arglis) > 0)
	{ message("unused movingAverage arguments:"); print(names(arglis)); browser()}
	#before April2014: new("movingAverage", from)
	as(x,'movingAverage')
})

smoothData <- function(x, y, ...){printGeneric("smoothData"); browser()}
#removeMethods("smoothData")
setMethod("smoothData", signature(x = "biasEstimate", y = "missing"), function(x, y, class2, ...)
{
	if(missing(class2))
	{
		class2 <- "Lowess"
		message("x of ", class(x), "; class2 = ", class2)
	}
	stopifnot(class2 %in% smoothData.classes)
	fu <- eval(parse(text = class2))
	if(!is.function(fu)) stop(paste("bad fu", class2))
	xvec <- Rank(x)
	yvec <- as(x, "numeric")
	if(!is.numeric(xvec) || length(xvec) != length(yvec))
	{ message("smoothData error"); browser()}
	sD <- fu(x = xvec, y = yvec, na.rm = FALSE, ...)
	if(!is(sD, class2))
	{ message("sD error"); browser()}
	if(length(xi(sD)) != length(Rank(x)) || !all(xi(sD) == Rank(x), na.rm = TRUE))
	{ message("smoothData x problem"); browser()}
	bias <- eta(sD)
	if(length(bias) != length(Rank(x)))
	{ message("smoothData y problem"); browser()}
	names(bias) <- names(x)
	if(length(annotation(x@uncorrected)) != 1)
	{ message(class(x), " annotation problem"); browser()}
	bE <- new.biasEstimate(bias = bias, uncorrected = x@uncorrected, sorted = x@sorted, jackknife = x@jackknife, Rank = x@Rank)
	smooth.ann <- try(annotation(bE), silent = TRUE)
	if(length(smooth.ann) == 0 || annotation(x) != smooth.ann)
	{
		message("smoothData bE has no annotation")
		browser()
	}
	bE
})
setMethod("smoothData", signature(x = "predictionError", y = "missing"), function(x, y, class2, ...)
{
	if(missing(class2))
	{
		class2 <- "Lowess"
		message("x of ", class(x), "; class2 = ", class2)
	}
	stopifnot(class2 %in% smoothData.classes)
	fu <- eval(parse(text = class2))
	if(!is.function(fu)) stop(paste("bad fu", class2))
	xvec <- 1:nrow(x)
	smooth.vec <- function(yvec)
	{
		if(!is.numeric(xvec) || length(xvec) != length(yvec))
		{ message(class(x), " smoothData error"); browser()}
		sD <- fu(x = xvec, y = yvec, na.rm = FALSE, ...)
		if(!is(sD, class2))
		{ message("sD error"); browser()}
		eta(sD)
	}
	x@.Data <- apply(x@.Data, 2, smooth.vec)
	stopifnot(validObject(x))
	x
})
smoothDataFUN <- function(...)
{
	function(x)
	{
		smoothData(x = x, ...)
	}
}

first.difference <- function(object, ...)
{
	if(!is.numeric(object))
		object <- as.numeric(object)
	object[2:length(object)] - object[1:(length(object) - 1)]
}

is.smooth <- function(object, maxDiff, verbose = FALSE)
{
	stopifnot(is.numeric(object))
	maxDiff <- Scalar(maxDiff)
	object.forward <- c(object[2:length(object)], object[length(object)])
	Diff <- abs(object.forward - object)
	if(verbose)
	{ message("maxDiff==", maxDiff, "; Diff:"); print(summary(Diff))}
	Diff <= maxDiff
}

movingLocation <- function(x, y, width, ...)
{
  message("generic movingLocation")
  browser()
}
#removeMethods("movingLocation")
setMethod("movingLocation", signature(x = "ANY", y = "ANY", width = "missing"), function(x, y, width, ...)
{
	stop("width missing")
})

# methods moved to reprod.s 4 September 2007

s3 <- function(object, name){message("generic s3"); browser()}
#removeMethods("s3")
s3.recursive <- function(object, name, first.call)
{
	if(missing(first.call))
		first.call <- TRUE
  x <- if(is.list(object) && "s3" %in% names(object))
    object$s3
  else if("s3" %in% slotNames(object))
    object@s3
  else if(first.call)
	  stop("s3 error")
  else
  	object
	if(first.call)
		s3.recursive(x, first.call = FALSE)
	else
		x
}
setMethod("s3", signature(object = "ANY", name = "missing"), function(object, name)
{
	s3.recursive(object = object, name = name)
})
setMethod("s3", signature(object = "ANY", name = "character"), function(object, name)
{
  s3(object)[name][[1]]
})


Slot <- function(object, name)
{
  if(name %in% slotNames(object))
    slot(object = object, name = name)
	else if((is.list(object) && "s3" %in% names(object)) || "s3" %in% slotNames(object))
		s3(object = object, name = name)
	else if(is.list(object) && name %in% names(object))
		object[name][[1]]
	else
	  stop(paste("Slot cannot extract", name, "from object of class", class(object)))
}

p.value <- function(object, ...){printGeneric("p.value"); browser()}
#removeMethods("p.value")
setMethod("p.value", signature(object = "list"), function(object, ...)
{
	Slot(object = object, name = "p.value", ...)
})

xi <- function(object, ...)
{
  Slot(object = object, name = "x", ...)
}
eta <- function(object, ...)
{
  Slot(object = object, name = "y", ...)
}

Lines <- function(x, y, ...){message("generic Lines"); browser()}
#removeMethods("Lines")
setMethod("Lines", signature(x = "ANY", y = "missing"), function(x, y, ...)
{
  lines(x = xi(x), y = eta(x), ...)
})
Plot <- function(x, y, ...){message("generic Plot; calling plot"); plot(x = x, y = y, ...)}
#removeMethods("Plot")
#setMethod("Plot", signature(x = "numeric", y = "function"), function(x, y, ...)
#{
#  plot(x = x, y = y, ...)
#})
setMethod("Plot", signature(x = "ANY", y = "missing"), function(x, y, ...)
{
  plot(x = xi(x), y = eta(x), ...)
})
setMethod("Plot", signature(x = "numeric", y = "numeric"), function(x, y, ...)
{
	nam <- intersect(names(x), names(y))
	stopifnot(length(nam) > 0)
  plot(x = x[nam], y = y[nam], ...)
})
Plot.XprnSet <- function(alpha = NULL, i, significance, es, x = 1, y = 1, call.par = TRUE, main = NULL, estimate, ...)
{
	if(missing(i))
	{
		significance <- significance[is.finite(significance)]
		Diff <- abs(significance - alpha)
		i <- names(Diff)[which(Diff == min(Diff))]
	}
	print(i)
	get.main <- function(lg)
	{
		if(is.null(main))
		{
			if(is.null(alpha))
				""
			else
			{
		#		pr <- function(tex, x){paste("P(", tex, ") = ", percent(x), "%", sep = "")}
		#		paste(pr("ratio > 1", 1 - alpha), " or ", pr("ratio = 1", 1), " ?")
		#		co <- function(tex, x){paste("confidence(", tex, ") = ", percent(x), "%", sep = "")}
				co <- function(tex, x){paste(percent(x), "% certain that ", if(lg) "log " else "", "ratio > ", if(lg) 0 else 1, sep = "")}
				paste(co("ratio > 1", 1 - alpha))
			}
		}
		else
			main
	}
	main <- get.main(FALSE)
	stopifnot(length(i) == 1)
	assert.is(es, "XprnSet")
	Xp <- as.numeric(exprs(es)[i, ])
	Xp <- Xp[is.finite(Xp)]
	print(Xp)
	xp <- log10(Xp)
	lab <- "observed expression ratio"
	lim <- c(1 / max(Xp), max(Xp))
	print(lim)
	if(call.par)
	{
		Mfrow()
		plot(Xp, Xp, log = "y", xlab = lab, ylab = lab, xlim = lim, ylim = lim, main = main)
		abline(v = 1, h = 1, col = "gray")
	#	text(x = 1, y = 1, labels = main)
	#	stripchart(Xp, xlab = lab, xlim = lim)
	#	hist(Xp, xlab = lab, ylab = "number of ratios", main = main)
		barplot(Xp, ylab = lab, main = main, log = "y", ...)
		abline(a = 1, b = 0)
		barplot(Xp, ylab = lab, main = main, ...)
		abline(a = 1, b = 0)
	}
	xp2 <- log2(Xp)
	sub = if(missing(estimate))
		""
	else
	{
		assert.is(estimate, "numeric")
		if(length(estimate) > 1 && is.character(i))
			estimate <- estimate[i]
		if(length(estimate) == 1)
			paste("local false discovery rate estimate:", percent(estimate))
		else
			stop("problematic estimate")
	}
	barplot(xp2, ylab = paste(lab, "(log2)"), main = get.main(TRUE), ylim = c(-max(abs(xp2)), max(abs(xp2))), sub = sub, ...)
	message(paste(round(sort(Xp), 1), collapse = ", "))
	t.test(xp, alternative = "greater")
}
setMethod("Plot", signature(x = "XprnSet", y = "missing"), function(x, ...)
{
	Plot.XprnSet(es = x, ...)
})



setMethod("annotation", signature(object = "ANY"), function(object)
{
	if(!isS4(object))
	{
		warning("getting annotation of non-S4 object")
		"S3 object"
	}
	else
	{
		if("annotation" %in% slotNames(object))
			object@annotation
		else if("ann" %in% slotNames(object))
			object@ann
		else
			stop("bad annotation call")
	}
})
setMethod("annotation", signature(object = "xprnSet"), function(object)
{
	annotation(object@es)
})

sameSet <- function(x, y){length(x) == length(y) && all(x %in% y) && all(y %in% x)}
disjoint <- function(x, y){length(intersect(x, y)) == 0}

sameNames <- function(object, ...){message("generic sameNames"); browser()}
#removeMethods("sameNames")
setMethod("sameNames", signature(object = "ANY"), function(object, ...)
{
	arglis <- list(...) # c(list(object), list(...))
	if("order.sensitive" %in% names(arglis))
	{
		order.sensitive <- arglis$order.sensitive
		arglis <- arglis[names(arglis) != "order.sensitive"]
		if("order.sensitive" %in% names(arglis))
		{ message("bad arglis in sameNames"); browser()}
	}
	else
		order.sensitive <- TRUE
	na.rm <- FALSE
	verbose <- FALSE
	nam0 <- names(object) # arglis[[1]])
	nam <- try(if(order.sensitive) nam0 else sort(nam0))
	if(is(nam, "try-error"))
	{
		message("nam err")
		browser()
	}
	if(verbose)
		print(nam)
	all(sapply(arglis, function(elem)
	{
	  elem.nam <- names(elem)
	  if(verbose)
	  	print(elem.nam)
	  if(!order.sensitive)
	  	elem.nam <- sort(elem.nam)
		identical(elem.nam, nam)
	}), na.rm = na.rm)
})

sameLengths <- function(...)
{
	lis <- list(...)
	all(length(lis[[1]]) == sapply(lis, length))
}
compatible <- function(...){sameLengths(...) && sameNames(...)}

printInvalid <- function(object, ...)
{
	mess <- paste("Invalid ", class(object), ..., " on ", date(), ".", sep = "")
	message(mess)
	warning(mess)
}

printGeneric <- function(object, ...)
{
	message("Generic ", if(is.character(object)) object else class(object), ..., " on ", date(), ".")
}

sorted <- function(object, ...)
{
	Slot(object = object, name = "sorted", ...)
}
#removeMethods("sorted")

sample.size <- function(object, ...){printGeneric("sample.size"); browser()}
#removeMethods("sample.size")
setMethod("sample.size", signature(object = "numeric"), function(object)
{
	sum(!is.na(object))
})
setMethod("sample.size", signature(object = "xprnSet"), function(object)
{
	mat <- exprs(object)
	size <- sapply(1:nrow(mat), function(i){sample.size(as.numeric(mat[i, ]))})
	names(size) <- featureNames(object)
	size
})

wilkinson.test <- function(x, mu = 0, alternative = "greater")
{
	assert.is(x, "numeric")
	p.value <- if(alternative == "greater")
	{
		prob <- function(q, lower.tail)
		{
			pr <- pt(q = q, df = length(x) - 1, lower.tail = lower.tail)
			stopifnot(length(pr) == 1)
			pr
		}
		q <- abs(mean(x)) / Sd(x, se.of.mean = TRUE)
		stopifnot(length(q) == 1 && q >= 0)
		sum(c(prob(q = q, lower.tail = FALSE), prob(q = -q, lower.tail = TRUE))) # pchisq(q = sum(x ^ 2), df = length(x), lower.tail = FALSE)
	}
	else if(alternative == "less")
		1 - wilkinson.test(x = x, mu = mu, alternative = "greater")$p.value
	else
		stop(paste("alternative '", alternative, "' not recognized", sep = ""))
	list(p.value = p.value)
}
PValueFUN <- function(FUN, alternative, ...)
{
	arglis <- list(...)
	null.arg.name <- "mu"
	mu <- if(null.arg.name %in% names(arglis))
		arglis[null.arg.name][[1]]
	else
		0
	stopifnot(is.numeric(mu) && length(mu) == 1)
  if(missing(FUN))
    FUN <- t.test
  if(!is.function(FUN))
  { message("bad FUN"); browser()}
  if(missing(alternative))
  	alternative <- character(0)
  function(x, y)
  {
  	if(missing(y))
  		y <- NULL
  	insignificant.p <- if(length(alternative) == 0)
     	0.5 # 1 before 20 August 2007
    else
     	1
  	get.p <- function(alternative)
  	{
  		extract.p <- function(...){FUN(alternative = alternative, ...)$p.value}
  		if(length(alternative) == 0 || alternative == "Greater") # one-sided, but conservative toward 0.5 rather than toward 1
  		{
  			less <- get.p(alternative = "less")
  			greater <- get.p(alternative = "greater")
  			if(length(alternative) == 0)
  			{
					if(less < greater)
						less
					else if(greater < less)
						1 - greater
					else
						insignificant.p
				}
				else if(alternative == "Greater") # added 090716; corrects some conservatism
					mean(c(1 - less, greater))
				else
					stop("bad alternative")
  		}
  		else if(is.null(y))
  			extract.p(x = x)
  		else
  			extract.p(x = x, y = y)
  	}
  	effective.sample.size <- function()
  	{
  		ss <- sample.size # function(vec){sum(!is.na(vec))}; generalized 12 May 2008
  		if(is.null(y))
  			ss(x)
  		else if(identical(FUN, wilcox.test) && min(ss(x), ss(y)) >= 1)
  			max(ss(x), ss(y))
  		else
  			min(ss(x), ss(y))
  	}
    p <- if(effective.sample.size() >= 2)
		{
			present.part <- function(vec){vec[!is.na(vec)]}
			x <- present.part(x)
			if(!is.null(y)) y <- present.part(y)
			is.essentially.constant <- function(vec){all(vec[1] == vec)}
			t.test.like <- identical(FUN, t.test) || identical(FUN, wilkinson.test)
			if(t.test.like && (is.essentially.constant(x) && (is.null(y) || is.essentially.constant(y))))
			{
				sample.mean <- if(is.null(y))
					x[1]
				else
					x[1] - y[1]
				if(sample.mean == mu)
					insignificant.p
				else
				{
					significant.p <- if(sample.mean > mu && length(alternative) == 0) 1 else 0 # 0 before 20 August 2007
					warning(paste("returning p-value of ", significant.p, " for ", paste(x, collapse = " | ")), " since t.test fails.", sep = "")
					significant.p
				}
			}
			else
				get.p(alternative = alternative) # typical computation of the p-value
		}
		else
			as.numeric(NA) # insignificant.p before 070828 15:57
    Scalar(p)
  } # end function(x, y)
}
RatioFUN <- function(FUN, na.rm, ...)
{
  if(missing(FUN))
    FUN <- function(x, na.rm){exp(mean(x, na.rm = na.rm))}
  if(!is.function(FUN))
  { message("bad FUN"); browser()}
  if(missing(na.rm))
    na.rm <- TRUE
  function(x)
  {
    Rat <- FUN(x, na.rm = na.rm, ...)
    Scalar(Rat)
  }
}
PValue <- function(object, get.PValue, alternative = default("Greater", "alternative"), verbose = TRUE, ...)
{
	if(missing(get.PValue))
		get.PValue <- PValueFUN(alternative = alternative, ...)
	recall <- function(object){PValue(object = object, get.PValue = get.PValue, alternative = alternative, ...)}
	if(is(object, "numeric"))
	{
		pval <- get.PValue(object)
		names(pval) <- names(object)
	}
	else if(is(object, "matrix"))
	{
		if(verbose)
			message("\n  began p-value computations ", date())
		pval <- apply(object, 1, get.PValue)
		if(verbose)
			message("\  completed p-value computations ", date(), "\n")
		stopifnot(nrow(object) == length(pval))
		names(pval) <- rownames(object)
	}
	else if(is(object, "xprnSetPair"))
	{
		get.mat <- function(es)
		{
			if(is(es, "XprnSet"))
				es <- logb(es)
			exprs(es)
		}
		xmat <- exprs(object@x)
		ymat <- exprs(object@y)
		stopifnot(all(rownames(xmat) == rownames(ymat)))
		nam <- rownames(xmat)
		if(verbose)
			message("\n  began p-value computations (2-sample test) ", date())
		pval <- sapply(1:length(nam), function(i)
		{
			get.vec <- function(mat)
			{
				vec <- try(as.numeric(mat[i, ]))
				if(is.err(vec))
				{ message("vec problem"); browser()}
				vec
			}
			get.PValue(x = get.vec(xmat), y = get.vec(ymat))
		})
		if(verbose)
			message("\  completed p-value computations ", date(), "\n")
		stopifnot(length(nam) == length(pval))
		names(pval) <- nam
	}
	else if(is(object, "XprnSet"))
		pval <- recall(object = logb(object))
	else if(is(object, "xprnSet"))
	{
		pval <- recall(object = exprs(object))
		nam <- featureNames(object)
		stopifnot(length(pval) == length(nam))
		names(pval) <- nam
	}
	else
		stop("bad class for PValue")
	if(!is.na(pval[1]) && all(pval == pval[1], na.rm = TRUE))
	{ message("all p-values same"); browser()}
	Numeric(pval)
}

exp_root <- function(vec)
{
	root <- sqrt(vec)
	exp(root)
}

Mean <- function(x, ...){printGeneric("Mean"); browser()}
#removeMethods("Mean")
setMethod("Mean", signature(x = "numeric"), function(x, na.rm, geometric = FALSE, ...)
{
	if(missing(na.rm)) na.rm <- TRUE
	fun <- if(geometric)
		geomean
	else
		mean
	fun(x, na.rm = na.rm, ...)
})
setMethod("Mean", signature(x = "matrix"), function(x, ...)
{
	vec <- sapply(1:nrow(x), function(i){Mean(x[i, ], ...)})
	stopifnot(nrow(x) == length(vec) && is.numeric(vec))
	names(vec) <- rownames(x)
	vec
})

stat <- function(x, y = numeric(0), FUN = mean, na.rm = TRUE, ...) # modified 090716 to non-numeric classes; this replaces "mean", signature(x = "xprnSet") of estimate.s
{
	assert.is(FUN, "function")
	arglis <- list(...)
	if(length(arglis) > 0)
	{
		fun0 <- FUN
		FUN <- function(...){do.call(fun0, c(list(...), arglis))}
	}
	recall <- function(X){stat(x = X, FUN = FUN)}
	recall.xy <- function(X, Y){stat(x = X, y = Y, FUN = FUN)}
	add.featureNames <- function(vec)
	{
		stopifnot(is(vec, "numeric"))
		nam <- featureNames(x)
		stopifnot(is.numeric(vec) && length(vec) == length(nam))
		names(vec) <- featureNames(x)
		vec
	}
	if(is(x, "xprnSetPair") && is.nothing(y))
	{
		vec <- recall.xy(X = exprs(x@x), Y = exprs(x@y))
		add.featureNames(vec = vec)
	}
	else if(is(x, "matrix") && !is.nothing(y) && is(y, "matrix"))
	{
		nam <- rownames(x)
		stopifnot(all(nam == rownames(y)))
		vec <- sapply(1:nrow(x), function(i)#april 2013 corrected by marta: nrow(x) <->length(nam)
		{
			get.num <- function(mat){mat[i, ]}
			recall.xy(X = get.num(x), Y = get.num(y))
		})
		stopifnot(is.numeric(vec) && length(vec) == nrow(y))
		names(vec) <- nam
		vec
	}
	else if(is(x, "matrix") && is.nothing(y))
	{
		vec <- apply(x, 1, recall)
		stopifnot(is.numeric(vec) && length(vec) == nrow(x))
		names(vec) <- rownames(x)
		vec
	}
	else if(is(x, "XprnSet") && is.nothing(y))
		recall(logb(x))
	else if(is(x, "xprnSet") && is.nothing(y))
	{
		vec <- recall(exprs(x))
		add.featureNames(vec = vec)
	}
	else if(is(x, "numeric"))
	{
		get.num <- function(num)
		{
			na.boo <- is.na(num)
			if(all(na.boo))
			{
				warning("there is nothing for stat to compute; returning NA")
				return(as.numeric(NA))
			}
			else if(any(na.boo))
				warning(paste(sum(na.boo), "of", length(num), "missing values removed to compute stat"))
			num
		}
		x <- get.num(num = x)
		if(is.nothing(y))
			FUN(x, na.rm = na.rm)
		else if(is(y, "numeric"))
		{
			y <- get.num(num = y)
			FUN(x = x, y = y, na.rm = na.rm)
		}
		else
			stop("cannot compute stat from numeric arguments given")
	}
	else
		stop("cannot compute stat from arguments given")
} # end stat

Stats <- function(x)
{
	Q <- function(...){stat(x, FUN = function(y, na.rm){quantile(y, na.rm = na.rm, ...)})}
	Q1 <- Q(probs = 0.25)
	average <- Q2 <- Q(probs = 0.5)
	Q3 <- Q(probs = 0.75)
	mean.abs <- stat(x, FUN = function(y, na.rm){mean(abs(y), na.rm = na.rm)})
	vec <- c(mean.abs = mean.abs, Q1 = Q1, Q2 = Q2, Q3 = Q3)
	names(vec) <- c("mean.abs", "Q1", "Q2", "Q3")
	if(any(is.na(vec)))
	{ warning("missing value in Stats")}
	vec
}

Length <- function(object)
{
	if(is(object, "xprnSetObject"))
	{
		Length(names(object))
	}
	else
	{
		length(object)
	}
}
Size <- function(object)
{
	nfinite <- function(x){sum(is.finite(x))}
	if(is(object, "xprnSetObject"))
	{
		sizes <- sapply(1:Length(object), function(i)
		{
			get.size <- function(es)
			{
				assert.is(es, "xprnSet")
				Size(as(exprs(es)[i, ], "numeric"))
			}
			if(is(object, "xprnSet"))
				get.size(object)
			else if(is(object, "xprnSetPair"))
				get.size(object@x) + get.size(object@y)
			else
				stop("unrecognized expression object 2")
		})
		assert.is(sizes, "numeric")
		stopifnot(length(sizes) == Length(object))
		names(sizes) <- names(object)
		sizes
	}
	else if(is(object, "numeric"))
		nfinite(object)
	else
		stop("unrecognized object 3")
}


geoMean <- function(x, ...)
{
	printGeneric("geoMean"); browser()
}
#try(removeMethods("geoMean"), silent = TRUE)
setMethod("geoMean", signature(x = "numeric"), function(x, ...) # na.rm corrected 18 March 2010
{
	Mean(x = x, geometric = TRUE, ...)
})

Sd <- function(object, ...){printGeneric("Sd"); browser()}
#removeMethods("Sd")
setMethod("Sd", signature(object = "numeric"), function(object, na.rm, mle, mu, se.of.mean)
{
	if(missing(na.rm)) na.rm <- TRUE
	if(missing(mle)) mle <- !missing(mu)
	if(missing(se.of.mean)) se.of.mean <- FALSE
	sample.mean <- mean(object, na.rm = na.rm)
	ss <- sum(!is.na(object))
	Sca <- if(ss == 0)
		as.numeric(NA)
	else if(ss == 1 && !mle && (na.rm || !any(is.na(object))))
		Inf
	else
	{
		va <- if(missing(mu) || length(mu) == 0)
			var(x = object, na.rm = na.rm)
		else
			sum((object - mu) ^ 2, na.rm = na.rm) / (ss - 1)
		if(mle)
			va <- va * (ss - 1) / ss
		sqrt(va)
	}
	if(se.of.mean)
		Sca <- Sca / sqrt(ss)
	Scalar(Sca)
})
se.mean <- function(...)
{
	Sd(..., se.of.mean = TRUE)
}
Var <- function(...){Sd(...) ^ 2}

cv <- function(object = numeric(0), y = numeric(0), x = numeric(0), na.rm, mle, mu, ...)
{
	if(missing(na.rm))
		na.rm <- TRUE
	assert.is(na.rm, "logical")
	if(missing(mle))
		mle <- !missing(mu)
	if(missing(mu))
		mu <- numeric(0)
	if(is.nothing(x) && !is.nothing(y))
	{
		x <- object
		object <- numeric(0)
	}
	else if(is.nothing(y) && is.nothing(object))
	{
		object <- x
		x <- numeric(0)
	}
	get.mu.hat <- function(vec)
	{
		if(length(mu) == 0)
			mean(x = vec, na.rm = na.rm)
		else
			mu
	}
	if(is.nothing(x) && is.nothing(y))
	{
		mu.hat <- get.mu.hat(vec = object)
		sigma.hat <- Sd(object = object, na.rm = na.rm, mle = mle, mu = mu, ...)
	}
	else if(is.nothing(object) && is.nothing(mu) && !mle)
	{
		mu.hat <- get.mu.hat(vec = x) - get.mu.hat(vec = y)
		get.size <- function(vec){sum(is.finite(vec))}
		size <- get.size(c(x, y))
#		sigma.hat <- Sd(object = c(x, y), na.rm = na.rm, mle = mle, ...) * sqrt((size - 1) / (size - 2)) # fails
		get.ss <- function(vec){(get.size(vec) - 1) * var(vec, na.rm = na.rm)}
		sigma.hat <- sqrt((get.ss(x) + get.ss(y)) / (size - 2))
	}
	else
		stop("bad cv arguments")
	scalar(as(sigma.hat, "numeric") / mu.hat)
}
t_size <- function(nx, ny){1 / (1 / nx + 1 / ny)}
t_statistic <- function(x, y, na.rm = TRUE)
{
	get.size <- function(vec){sum(is.finite(vec))}
	size <- get.size(c(x, y))
	nx <- get.size(x)
	ny <- get.size(y)
	stopifnot(size == nx + ny)
	get.ss <- function(vec){(get.size(vec) - 1) * var(vec, na.rm = na.rm)}
	(mean(x) - mean(y)) / sqrt(1 / t_size(nx = nx, ny = ny) * (get.ss(x) + get.ss(y)) / (size - 2))
}
icv <- function(...){1 / cv(...)}
test.t_statistic <- function(x, y, ...) # test.t_statistic(x = 1:3, y = 4:7); test.t_statistic(x = 1:3, y = 4:7, alternative = "greater")
{
	get.t <- function(var.equal){t.test(x = x, y = y, var.equal = var.equal, ...)$stat}
	c(t_statistic = t_statistic(x = x, y = y), var.eq = get.t(var.equal = TRUE), var.ne = get.t(var.equal = FALSE))
}
test.icv <- function(x, y) # test.icv(x = 1:3, y = 4:7)
{
	get.size <- function(vec){sum(is.finite(vec))}
	size.x <- get.size(x)
	size.y <- get.size(y)
	fac <- sqrt(1 / size.x + 1 / size.y)
	t.test.icv <- fac * t.test(x = x, y = y, var.equal = TRUE)$stat
	t_statistic.icv <- fac * t_statistic(x = x, y = y)
	c(icv = icv(x = x, y = y), t.test.icv = t.test.icv, t_statistic.icv = t_statistic.icv)
}
functional.cv <- function(...)
{
	cv(..., mle = TRUE)
}
relativeVar <- function(...){Scalar(cv(...) ^ 2)}
functional.relativeVar <- function(...)
{
	relativeVar(..., mle = TRUE)
}
t.stat <- function(...)
{
	1 / cv(..., se.of.mean = TRUE)
}
relative.frequency <- function(x, threshold, offset.p = default(0.5, "offset.p")) # offset.p added 090401 (previous default was effectively NULL)
{
	if(missing(threshold))
		threshold <- 0
	x <- x[!is.na(x)]
	if(length(x) == 0)
		as.numeric(NA)
	else
	{
		offset.size <- if(is.null(offset.p))
		{
			offset.p <- 0
			0
		}
		else
			1
		(sum(x > threshold) + sum(x == threshold) / 2 + offset.p) / (length(x) + offset.size)
	}
}
auc <- function(x, y, na.rm, ...)
{
	if(missing(na.rm))
		na.rm <- TRUE
	get.vec <- function(object){object[!is.na(object)]}
	x <- get.vec(x)
	y <- get.vec(y)
	if(length(x) == 0 || length(y) == 0)
		as.numeric(NA)
	else
	{
		compare <- length(x) * length(y)
		check.w <- function(w){if(w < 0 || w > 1) stop("bad w statistic")}
		get.w <- function(x, y)
		{
			options(warn = -1)
			w <- wilcox.test(x = x, y = y, na.rm = na.rm, alternative = "two.sided", ...)$statistic / compare
			options(warn = 0)
			check.w(w)
			w
		}
		ww <- get.w(x, y) # mean(get.w(x, y), 1 - get.w(y, x))
		check.w(ww)
		ww
	}
}

hsm <- function(x, na.rm) # half-sample mode of a vector [D. R. Bickel and R. Frühwirth (contributed equally), "On a Fast, Robust Estimator of the Mode: Comparisons to Other Robust Estimators with Applications," Computational Statistics and Data Analysis 50, 3500-3530 (2006)]
{
	stopifnot(is.numeric(x))
	if(missing(na.rm))
		na.rm <- FALSE
	if(na.rm)
		x <- x[!is.na(x)]
  y <- sort(x);
  while (length(y)>=4)
  {
    m <- ceiling(length(y)/2);
    w.min <- y[length(y)]-y[1];
    for(i in 1:(length(y)-m+1))
    {
      w <- y[i+m-1]-y[i];
      if(w<=w.min)
      {
        w.min <- w;
        j <- i
      }
    }
    if(w==0)
      y <- y[j]
    else
      y <- y[j:(j+m-1)]
  }
  if(length(y) == 3)
  {
    z <- 2*y[2]-y[1]-y[3];
    if(!is.finite(z))
    {
      print('ERROR: z is not finite; x, y, and z follow:');
      print(x);
      print(y);
      print(z);
    }
    if(z < 0)
      mean(y[1:2])
    else if(z > 0)
      mean(y[2:3])
    else
      y[2]
  }
  else
    mean(y)
}

assymetry <- function(x, ...)
{
	if(!is.numeric(x))
		x <- as.numeric(x)
	Mean <- mean(x, ...)
	Median <- median(x, ...)
	2 * (Mean - Median) / abs(Mean + Median)
}

Scale <- function(object, ...){printGeneric("Scale")}
#removeMethods("Scale")
Location <- function(object, ...){printGeneric("Location")}
#removeMethods("Location")

oneLeftOut <- function(object, ...){message("generic oneLeftOut"); browser()}
#removeMethods("oneLeftOut")
setMethod("oneLeftOut", signature(object = "numeric"), function(object, FUN, arglis.FUN=list())
{	#XXX|:multiple local function definitions for ÔfunÕ with different formal arguments
	#XXX|:removed , FUN = "function", changed "..." by arglis.FUN=list()
	assert.is(FUN,"function")
#fun <- function(x){FUN(x,  ...)}XXX|:changed	
  fun <- function(x){do.call(FUN,c(list(x),arglis.FUN))}#XXX|:changed
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
setMethod("oneLeftOut", signature(object = "xprnSet"), function(object, FUN, verbose, factor.name,  arglis.FUN=list())
{#XXX|:removed , FUN = "function", changed "..." by arglis.FUN=list()
	assert.is(FUN,"function")
	if(missing(verbose)) verbose <- FALSE
	if(!missing(factor.name) && is.character(factor.name) && length(factor.name) == 1) # two-sample
	{
		oneLeftOut(object = xprnSetPair(x = object, factor.name = factor.name), FUN = FUN, verbose = verbose, ...)
	}
	else if(missing(factor.name)) # single-sample; stopping condition
	{	fun <- function(x){do.call(FUN,c(list(x),arglis.FUN))}# , XXX|:changed
		#fun <- function(x){FUN(x, verbose = verbose, ...)} verbose passed 31 October 2007
		noneLeftOut <- list(fun(object))
		js <- 1:ncol(exprs(object))
		training <- lapply(js, function(j)
		{
			boo <- js != j
			stopifnot(sum(!boo) == 1 && length(boo) == ncol(exprs(object)))
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
		stop(paste("bad oneLeftOut factor.name:", factor.name))
})
setMethod("oneLeftOut", signature(object = "xprnSetPair"), function(object, FUN, verbose,  arglis.FUN=list())
{#XXX|:removed , FUN = "function", changed "..." by arglis.FUN=list()
	assert.is(FUN,"function")
	if(missing(verbose)) verbose <- FALSE
	x <- object@x
	y <- object@y
	fun <- function(x, y)
	{
		if(missing(y))
			stop("y is missing in the two-sample case")
	  #FUN(x = x, y = y, ...)
	  do.call(FUN,c(list(x=x,y=y),arglis.FUN))# , XXX|:changed
	

	}
	noneLeftOut <- list(fun(x = x, y = y))
	get.js <- function(object) {1:ncol(exprs(object))}
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
			ncol(exprs(object)) - 1
		else
			1
		stopifnot(sum(boo) == size && length(boo) == ncol(exprs(object)))
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
				if(!ok){message("get.lis.elem error"); browser()}
				list(elem)
			}
			training <- c(training, get.lis(negate = TRUE))
			test <- c(test, get.lis(negate = FALSE))
			if(length(training) != length(test))
			{ message("training and test error"); browser()}
		}
	}
	if(verbose)
	{
		message("two-sample oneLeftOut mapped ", class(object), " to ", 2 * length(training) + 1, " objects of class ", class(noneLeftOut[[1]]))
		print(FUN)
	}
	olo <- try(new("oneLeftOut", training = training, test = test, noneLeftOut = noneLeftOut))
	if(is(olo, "try-error") || !validObject(olo))
	{ message('"oneLeftOut", signature(object = "xprnSetPair", FUN = "function") error'); browser()}
	olo
})
setMethod("oneLeftOut", signature(object = "oneLeftOut"), function(object, FUN,  arglis.FUN=list())
{#XXX|:multiple local function definitions for ÔfunÕ with different formal arguments
	##XXX|:removed , FUN = "function", changed "..." by arglis.FUN=list()
	assert.is(FUN,"function")
	arglis <- arglis.FUN
	alo.nam <- "anotherLeftOut"
	if(alo.nam %in% names(arglis))
	{
		alo <- arglis[alo.nam][[1]]
		if(!is(alo, "oneLeftOut") || length(arglis) != 1)
		{ message('cannot use "oneLeftOut", signature(object = "oneLeftOut", FUN = "function") this way'); browser()}
		# arglis <- arglis[names(arglis) != alo]
		fun <- function(x, y)
		{
			FUN(x, y)
		}
		change.lis <- function(name)
		{
			ok <- name %in% slotNames(object) && name %in% slotNames(alo)
			if(!ok){ message("change.lis error"); browser()}
			xlis <- slot(object = object, name = name)
			ylis <- slot(object = alo, name = name)
			compat <- length(xlis) == length(ylis) && sameNames(xlis, ylis)
			if(!compat)
			{ message("inner problem of oneLeftOut"); browser()}
			new.lis <- lapply(1:length(xlis), function(i)
			{
				fun(x = xlis[[i]], y = ylis[[i]])
			})
			names(new.lis) <- names(xlis)
			new.lis
		}
		training <- change.lis("training")
		test <- change.lis("test")
		noneLeftOut <- change.lis("noneLeftOut")
	}
	else
	{
		#fun <- function(x){FUN(x, ...)}
		fun <- function(x){do.call(FUN,c(list(x),arglis.FUN))}#XXX|:changed
		training <- lapply(object@training, fun)
		test <- lapply(object@test, fun)
		noneLeftOut <- lapply(object@noneLeftOut, fun)
	}
	len.ok <- length(training) == length(object@training) && length(test) == length(object@test) && length(noneLeftOut) == length(object@noneLeftOut)
	if(!len.ok)
	{ message("length problem in oneLeftOut"); browser()}
	new("oneLeftOut", training = training, test = test, noneLeftOut = noneLeftOut)
})

noneLeftOut <- function(object, ...){message("generic noneLeftOut"); browser()}
#removeMethods("noneLeftOut")
setMethod("noneLeftOut", signature(object = "oneLeftOut"), function(object, ...)
{
  (object@noneLeftOut)[[1]]
})

sign_abs_error <- function(predicted, observed, sign.err.factor)
{
	err <- abs(predicted - observed)
	if(missing(sign.err.factor))
	{
		message("no penalty for sign errors")
		err
	}
	else
		err * ifelse(sign(predicted) * sign(observed) >= 0, 1, sign.err.factor)
}
sign_abs_error.fun <- function(sign.err.factor)
{
	function(predicted, observed){sign_abs_error(predicted = predicted, observed = observed, sign.err.factor = sign.err.factor)}
}

squaredError <- function(predicted, observed, object, ...){printGeneric("squaredError"); browser()}
#removeMethods("squaredError")
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
		{ message("bad err; decay.scale: ", decay.scale); browser()}
	}
	err
})

mean.squaredError <- function(...)
{
	error <- squaredError(...)
	FUN <- function(x){mean(x, na.rm = TRUE)}
	mse <- if(is.matrix(error))
		apply(error, 2, FUN)
	else if(is.numeric(error))
		FUN(error)
	if(!is.numeric(mse))
	{ message("mse is not numeric"); browser()}
	mse
}

coherence <- function(..., decay.scale) # a wrapper of squaredError
{
	if(missing(decay.scale))
	{
		decay.scale <- 1
		message("coherence decay.scale = ", decay.scale)
	}
	stopifnot(is.numeric(decay.scale) && length(decay.scale) == 1 && decay.scale > 0)
	co <- squaredError(..., relative = TRUE, decay.scale = decay.scale) # exp(-squaredError(..., relative = TRUE) / decay.scale)
	if(any(is.na(co)) || any(co < 0 | co > 1))
	{ message("coherence bad co"); browser()}
	co
}

mean.coherence <- function(..., decay.scale) # a wrapper of mean.squaredError
{
	if(missing(decay.scale))
	{
		decay.scale <- 1
		message("mean.coherence decay.scale = ", decay.scale)
	}
	stopifnot(is.numeric(decay.scale) && length(decay.scale) == 1 && decay.scale > 0)
	co <- mean.squaredError(..., relative = TRUE, decay.scale = decay.scale) # exp(-squaredError(..., relative = TRUE) / decay.scale)
	if(any(is.na(co)) || any(co < 0 | co > 1))
	{ message("mean.coherence bad co"); browser()}
	co
}

new.log <- function(old.log, base, old.base, new.base)
{
	if(missing(new.base))
		new.base <- base
	else
		stopifnot(missing(base))
	if(missing(old.base))
		old.log / logb(new.base)
	else if(all(old.base == new.base))
		old.log
	else if(length(old.base) == 0)
		logb(old.log, base = new.base)
	else if(length(new.base) == 0)
		old.base ^ old.log
	else
		old.log * logb(old.base, base = new.base)
}
new.log2 <- function(...){new.log(..., base = 2)}

sampleSize <- function(object, ...){printGeneric("sampleSize"); browser()}
#removeMethods("sampleSize")

predictionError <- function(object, ...){message("generic predictionError"); browser()}
#removeMethods("predictionError")
setMethod("predictionError", signature(object = "oneLeftOut"), function(object, error.fun, prediction.fun, call.Order, shrinkage, is.relative, ...)
{
	if(missing(error.fun))
	{
		error.fun.name <- "squaredError"
		error.fun <- eval(parse(text = error.fun.name))
		message("error.fun = ", error.fun.name)
	}
	stopifnot(is.function(error.fun))
	if(missing(is.relative))
	{
		is.relative <- identical(error.fun, squaredError)
		message("is.relative = ", is.relative)
	}
	nlo <- noneLeftOut(object)
	prediction.fun.stop <- function(){stop("prediction.fun should be specified as a function of a ", class(nlo), " returning a numeric.")}
	if(missing(prediction.fun))
		prediction.fun <- if(is(nlo, "estimate")) # estimate is in estimate.r
			Location
		else
			NULL
  if(is.null(prediction.fun)) # object consists of predictions
  {
  	if(!is(nlo, "numeric"))
  		prediction.fun.stop()
  	errors <- sapply(1:length(object@test), function(i)
  	{
  		error.fun(predicted = object@training[[i]], observed = object@test[[i]])
  	})
  	ef.ok <- is.matrix(errors) && ncol(errors) == length(object@test)
  	if(!ef.ok)
  	{ message("error applying error.fun"); browser()}
  	vec <- as.numeric(apply(errors, 1, mean, na.rm = TRUE))
  	stopifnot(length(vec) == length(nlo))
  	if(is.relative)
  	{
  		vec <- vec / nlo ^ 2
  		names(vec) <- names(nlo)
  	}
  	vec
  }
  else if(is.function(prediction.fun))
  {
		if(missing(call.Order))
		{
			call.Order <- TRUE
			message("call.Order = ", call.Order)
		}
		pe.vec <- function(...)
		{
			obj <- if(call.Order)
				sort(x = object, ...)
			else
				object
			predictionError(object = oneLeftOut(object = obj, FUN = prediction.fun), error.fun = error.fun, prediction.fun = NULL, call.Order = FALSE, is.relative = is.relative)
		}
		if(missing(shrinkage) && !is(nlo, "estimate"))
	  	pe.vec(...)
	  else if(!missing(shrinkage) && is(shrinkage, "numeric") && length(shrinkage) >= 1 && is(nlo, "estimate"))
	  {
	  	mat <- sapply(shrinkage, function(shri){pe.vec(shrinkage = shri, ...)})
	  	if(!is.matrix(mat))
	  		mat <- matrix(mat, ncol = 1)
	  	if(length(shrinkage) != ncol(mat))
	  	{ message("mat fails to match shrinkage"); browser()}
	  	new("predictionError", mat, parameter = shrinkage, parameter.lab = "amount of shrinkage")
	  }
	  else
	  	stop("incompatible predictionError arguments")
  }
  else
  	prediction.fun.stop()
})

accumulate <- function(object, ...){printGeneric("accumulate"); browser()}
#removeMethods("accumulate")
setMethod("accumulate", signature(object = "numeric"), function(object, indices, cumulative, nfeatures, two.sided, ...)
{
	if(missing(two.sided)) two.sided <- FALSE
	if(missing(cumulative))
		cumulative <- TRUE
	if(missing(indices))
	{
		if(missing(nfeatures))
		{
			nfeatures <- length(object)
			message(class(object), " nfeatures = ", nfeatures)
		}
		indices <- 1:nfeatures
	}
	stopifnot(all(indices >= 1 & indices <= length(object)))
	vec <- if(cumulative)
		sapply(indices, function(i)
		{
			sub.indices <- if(!two.sided || i < length(object) / 2)
				1:i
			else
				i:length(object)
			mean(object[sub.indices], na.rm = TRUE)
		})
	else
		object[indices]
	if(is.character(names(object)))
		names(vec) <- names(object)[indices]
	vec
})
setMethod("accumulate", signature(object = "matrix"), function(object, ...)
{
	mat <- sapply(1:ncol(object), function(j){accumulate(object = object[, j], ...)})
#	dimnames(mat) <- dimnames(object)
	mat
})
setMethod("accumulate", signature(object = "predictionError"), function(object, ...)
{
	dn <- dimnames(object)
	mat <- accumulate(object = as(object, "matrix"), ...)
	object@.Data <- mat
#	dimnames(object) <- dn
	vo <- try(validObject(object))
	if(is(vo, "try-error") || !vo)
	{ message("vo problem"); browser()}
	object
})

moderatedT <- function(object, ...){printGeneric("moderatedT"); browser()}
#removeMethods("moderatedT")

annotated <- function(object)
{
	tr <- try(annotation(object), silent = TRUE)
	tr != "try-error"
}

new.biasEstimate <- function(bias, uncorrected, sorted, jackknife, Rank, Weight) # changed 16 April 2008
{
	if(!annotated(uncorrected)){message("bad uncorrected for new.biasEstimate"); browser()}
	ann <- try(annotation(uncorrected))
	if(is(ann, "try-error")){stop("bad annotation(uncorrected)")}#; browser()}
	if(missing(Rank))
		Rank <- Numeric(numeric(0))
	stopifnot(length(names(bias)) == length(names(uncorrected)) && all(names(bias) == names(uncorrected)))
	if(missing(Weight)) Weight <- numeric(0)
	if(!is(Weight, "Scalar"))
		Weight <- Scalar(Weight)
	bE <- new("biasEstimate", bias, uncorrected = uncorrected, sorted = sorted, jackknife = jackknife, Rank = Rank, Weight = Weight)
	smooth.ann <- try(annotation(bE), silent = TRUE)
	if(length(smooth.ann) == 0)
	{
		warning("new.biasEstimate bE has no annotation")
#		browser()
	}
	bE
}
biasEstimate <- function(object, uncorrected, ...){message("generic biasEstimate"); browser()}
#removeMethods("biasEstimate")
setMethod("biasEstimate", signature(object = "numeric", uncorrected = "missing"), function(object, uncorrected) # object consists of observations
{
	stop("not yet implemented for object that consists of observations rather than leave-1-out estimates")
})
setMethod("biasEstimate", signature(object = "numeric", uncorrected = "scalar"), function(object, uncorrected, jackknife)
{
	ok <- is.logical(jackknife) && length(jackknife) == 1
	if(!ok)
	{
		mess <- paste('jackknife argument must be logical of length 1 in method for class "', class(object), '"; ', 'this method is better called via the method "biasEstimate", signature(object = "estimateLeftOut") of file estimate.s', sep = '')
		stop(mess)
	}
	leave.1.out <- object[is.finite(object)]
	size <- length(leave.1.out)
	if(size == 0 || !is.finite(uncorrected))
		as.numeric(NA)
	else
	{
		sample.mean <- mean(leave.1.out, na.rm = FALSE)
		stopifnot(length(sample.mean) == 1)
		if(jackknife) # object consists of leave-1-out training-set estimates
		{
			jack.est <- (size - 1) * (sample.mean - uncorrected) # Theorem 2.1 of Efron: good for functional
			jack.est # jackknife
		}
		else # object consists of leave-1-out difference training-set estimates and test-set estimates
		{
			sample.mean # uncorrected is not used here
		}
			# jack.est * (size - 1) / size # bootstrap, Efron p. 10; cf. pp. 33-34, 44-45
	}
})
setMethod("biasEstimate", signature(object = "matrix", uncorrected = "numeric"), function(object, uncorrected, ...) # object consists of estimates
{
	stopifnot(nrow(object) == length(uncorrected))
	vec <- sapply(1:nrow(object), function(i)
	{
		biasEstimate(object = object[i, ], uncorrected = scalar(uncorrected[i]), ...)
	})
	names(vec) <- names(uncorrected)
	vec
})
setMethod("biasEstimate", signature(object = "oneLeftOut", uncorrected = "missing"), function(object, uncorrected, sorted, jackknife, ...) # object consists of ordered estimates
{
	ok <- is.logical(jackknife) && length(jackknife) == 1
	if(!ok)
	{
		mess <- paste('jackknife argument must be logical of length 1 in method for class "', class(object), '"; ', 'this method is better called via the method "biasEstimate", signature(object = "estimateLeftOut") of file estimate.s', sep = '')
		stop(mess)
	}
	if(missing(sorted))
		stop('sorted missing; this method is better called by the method "biasEstimate", signature(object = "estimateLeftOut") of file estimate.s')
	if(!sorted)
	{
		message("object of estimates is not sorted; press c to continue anyway")
		browser()
		warning("object of estimates is not sorted")
	}
#	if(missing(estimator) || !is.function(estimator))
#	{
#		stop("estimator is not the function used to generate object of ordered estimates")
#		estimator.name <- "cv"
#		estimator <- eval(parse(text = estimator.name))
#		message("estimator = ", estimator.name)
#	}
	nlo <- noneLeftOut(object)
	if(is(nlo, "numeric"))
		uncorrected <- nlo
	else
		stop("object does not consist of numeric estimates")
	get.mat <- function(name) {sapply(slot(object, name = name), identity)}
	training.mat <- get.mat(name = "training")
	mat <- if(jackknife)
		training.mat
	else
	{
		test.mat <- get.mat(name = "test")
		stopifnot(all(dim(training.mat) == dim(test.mat)))
		training.mat - test.mat
	}
	stopifnot(nrow(mat) == length(uncorrected))
	bias <- biasEstimate(object = mat, uncorrected = uncorrected, jackknife = jackknife, ...) # as.numeric(apply(mat, 1, biasEstimate, ...))
	names(bias) <- names(uncorrected)
	new.biasEstimate(bias = bias, uncorrected = uncorrected, sorted = sorted, jackknife = jackknife)
})
setMethod("biasEstimate", signature(object = "biasEstimate", uncorrected = "missing"), function(object, uncorrected, ranks, nfeatures, use.lower.ranks, verbose, ...)
{
	total.nfeatures <- length(object)
	if(missing(verbose)) verbose <- FALSE
	if(missing(ranks))
	{
		if(missing(nfeatures))
		{
			nfeatures <- min(200, total.nfeatures)
			message("nfeatures = ", nfeatures)
		}
		if(missing(use.lower.ranks))
		{
			use.lower.ranks <- logical(0)
			message("use.lower.ranks = ", use.lower.ranks)
		}
		nf <- min(nfeatures, total.nfeatures)
		if(length(use.lower.ranks) == 0)
			nf <- ceiling(nf / 2)
		lower <- 1:nf
		upper <- (total.nfeatures + 1 - nf):total.nfeatures
		ranks <- if(length(use.lower.ranks) == 0)
			union(lower, upper)
		else if(use.lower.ranks)
			lower
		else
			upper
		if(verbose)
		{ cat("ranks = "); print(summary(ranks))}
	}
	else if(!(missing(nfeatures) && missing(use.lower.ranks)))
		stop('"biasEstimate", signature(object = "biasEstimate", uncorrected = "missing") has incompatible arguments')
	if(length(ranks) == 0 || length(ranks) > length(object))
	{ message('"biasEstimate", signature(object = "biasEstimate", uncorrected = "missing") has bad ranks'); browser()}
	stopifnot(all(ranks >= 1 & ranks <= max(Rank(object))))
	boo <- Rank(object) %in% ranks
	if(length(boo) != length(object) || sum(boo) != length(ranks))
	{ message('"biasEstimate", signature(object = "biasEstimate", uncorrected = "missing") has bad boo'); browser()}
	object[boo]
})

corrected <- function(object, ...){printGeneric("corrected"); browser()}
#removeMethods("corrected")
setMethod("corrected", signature(object = "biasEstimate"), function(object)
{
	object@uncorrected - object
})

Rank <- function(object, ...){printGeneric("Rank"); browser()}
#removeMethods("Rank")
setMethod("Rank", signature(object = "numeric"), function(object, factor, ...)
{
	if(missing(factor))
	{
		factor <- "warn"
		message("Rank factor = ", factor)
	}
	if(factor == "warn")
	{
		recurse <- function(factor){Rank(object = object, factor = factor, ...)}
		ra <- recurse(factor = FALSE) # no jitter
		if(all(ra == floor(ra)))
			ra
		else
		{
			warning(paste("ties broken randomly on ", date()))
			recurse(factor = 1)
		}
	}
	else
	{
		if(!is.logical(factor) || factor)
			object <- jitter(object, factor = factor)
		ra <- rank(object, ...)
		stopifnot(all(ra >= 0, na.rm = TRUE))
		Numeric(ra)
	}
})
setMethod("Rank", signature(object = "biasEstimate"), function(object)
{
	if(!object@sorted)
		stop(paste(class(object), "object is not sorted"))
	if(length(object@Rank) == 0)
	{
		vec <- Numeric(1:length(object))
		names(vec) <- names(object)
		Numeric(vec)
	}
	else
		object@Rank
})
setMethod("sort", signature(x = "biasEstimate"), function(x, decreasing = F, call.abs, save.memory, ...)#XXX|:removed decreasing = "ANY"
{
	ok <- is.logical(decreasing) && !decreasing
	if(!ok)
		stop(paste(class(x), "sort not implemented for non-default decreasing"))
	message('"sort", signature(x = "biasEstimate", decreasing = "missing")')
	if(x@sorted)
	{
		warning(class(x), " is already sorted")
		y <- x
	}
	else
	{
		warning("sorting ", class(x), " without additional information may be futile")
		if(missing(call.abs))
		{
			call.abs <- FALSE
			message("sort call.abs = ", call.abs)
		}
		if(missing(save.memory))
		{
			save.memory <- length(x >= 1e5)
			message("sort ", class(x), " save.memory = ", save.memory)
		}
		object <- x@uncorrected
		if(call.abs)
			object <- abs(object)
		message("Ranking ", class(x), " by", if(call.abs) " abs of " else " ","its ", class(object), " on ", date())
		ra <- Rank(object = object, ...)
		x@sorted <- TRUE
		if(any(is.na(ra)))
		{ message("missing rank value"); browser()}
		ord <- order(ra, decreasing = decreasing)
		sorted.ra <- ra[ord]
		stopifnot(all(sorted.ra == (1:length(ra))))
		if(!save.memory)
			x@Rank <- sorted.ra
		sort.slot <- function(name)
		{
			value <- slot(object = x, name = name)
			stopifnot(length(value) == length(ord))
			value[ord]
		}
		nam <- names(x)
		if(is.character(nam))
			nam <- nam[ord]
		x@.Data <- sort.slot(name = ".Data")
		stopifnot(length(x) == length(nam))
		names(x) <- nam
		x@uncorrected <- sort.slot(name = "uncorrected")
		stopifnot(validObject(x))
		x
	}
	stopifnot(y@sorted)
	y
})

blank.plot <- function(legend, ...)
{
	plot(x = 0, type = "n", col.axis = "white", xlab = "", ylab = "", tcl = 0, xaxt = "n", yaxt = "n", bty = "n", ...)
	if(!missing(legend))
		legend(x = "center", legend = legend, col = "white", bty = "n")
}

# Sample moved from reprod.s & numeric method added on 24 October 2007:
Sample <- function(object, ...){printGeneric("Sample"); browser()}
#removeMethods("Sample")
setMethod("Sample", signature(object = "numeric"), function(object, fac, sample.FUN = NULL, ...)
{
	if(is.null(sample.FUN))
		sample.FUN <- sample
	assert.is(sample.FUN, "function")
	if(missing(fac) || is.null(fac))
		sample.FUN(object, ...) # sample(x = object, ...)
	else if(length(object) == length(fac))
	{
		levs <- unique(fac)
		lis <- lapply(levs, function(lev)
		{
			stopifnot(length(lev) == 1)
			Sample(object[fac == lev], ...)
		})
		vec <- do.call("c", lis)
		stopifnot(length(vec) == length(object))
		vec
	}
	else
		stop("fatal fac")
})
setMethod("Sample", signature(object = "ExpressionSet"), function(object, factor.name, parametric = FALSE, ...)
{
	fac <- if(missing(factor.name) || is.null(factor.name))
		NULL
	else if(is.character(factor.name))
		pData(object)[, factor.name]
	else
		stop("fatal factor.name!")
	get.sample <- function(x, sample.FUN)
	{
		Sample(object = x, fac = fac, sample.FUN = sample.FUN, ...)
	}
	if(parametric)
	{
#		if(is.null(fac))
#			parametricBootstrap()
#		else
#			stop("fac and parametric")
		mat <- t(sapply(1:nrow(exprs(object)), function(i)
		{
			get.sample(x = exprs(object)[i, ], sample.FUN = parametricBootstrap)
		}))
		if(!all(dim(mat) == dim(exprs(object))))
		{ message("bad dim(mat) == dim(exprs(object))"); print(list(dim(mat), dim(exprs(object)))); browser()}
		rownames(mat) <- rownames(exprs(object))
		colnames(mat) <- colnames(exprs(object))
		exprs(object) <- mat
		stopifnot(validObject(object))
		object
	}
	else
		object[, get.sample(x = 1:ncol(exprs(object)), sample.FUN = NULL)] # sample
})
setMethod("Sample", signature(object = "xprnSet"), function(object, ...)
{
	object@es <- Sample(object@es, ...)
	object
})
setMethod("Sample", signature(object = "XprnSet"), function(object, ...)
{
	Sample(logb(object), ...)
})
setMethod("Sample", signature(object = "xprnSetPair"), function(object, ...) # new 25 April 2008
{
	sample.es <- function(es){Sample(object = es, ...)}
	object@x <- sample.es(es = object@x)
	object@y <- sample.es(es = object@y)
	stopifnot(validObject(object))
	object
})

mean_and_sd <- function(object, i) # not for the user
{
	if(is(object, "XprnSet"))
		object <- logb(object)
	stopifnot(is(object, "xprnSet"))
	if(missing(i))
	{
		i <- 1:10
		cat("i = "); print(i)
	}
	object <- object[i, ]
	apply(exprs(object), 1, function(ro){c(mean = mean(ro, na.rm = TRUE), sd = sd(ro, na.rm = TRUE))})
}
plot_mean_and_sd <- function(object, i, ...) # not for the user
{
	stopifnot(is(object, "list"))
	if(missing(i))
	{
		i <- 1:100
		cat("i = "); print(i)
	}
	mats <- lapply(object, function(es){mean_and_sd(object = es, i = i)})
	histogram <- length(i) >= 20
	if(length(i) <= 10) print(mats)
	plot.moment <- function(moment.name)
	{
		moments.mat <- sapply(mats, function(mat){mat[moment.name, ]})
		colnames(moments.mat) <- names(mats)
		get.y <- function(j){moments.mat[, j]}
		cat("\n", moment.name, " diagnostics for ", length(i), " genes:\n", sep = "")
		print.diagnostics <- function(reference.name, mat.name)
		{
			mat.minus.reference <- if(mat.name == "sd")
			{
				get.y(j = mat.name) ^ 2 - get.y(j = reference.name) ^ 2
				mat.name <- "var"
			}
			else
				get.y(j = mat.name) - get.y(j = reference.name)
			diff.name <- paste(mat.name, " minus ", reference.name, sep = "")
			message(diff.name, ":")
			print(summary(mat.minus.reference))
			print.portion <- function(ineq)
			{
				boo <- if(ineq == "less")
					mat.minus.reference < 0
				else if(ineq == "greater")
					mat.minus.reference > 0
				boo <- boo[!is.na(boo)]
				message(round(sum(boo) * 100, 1) / length(boo), "% are ", ineq, " than 0")
			}
			print.portion("less")
			print.portion("greater")
			if(histogram) hist(mat.minus.reference, xlab = diff.name, main = moment.name)
		}
		get.name <- function(index){names(mats)[index]}
		reference.index <- 1
		ylim <- range(as.numeric(moments.mat), na.rm = TRUE)
		for(index in 1:length(object))
		{
			if(index != reference.index)
			{
				print.diagnostics(reference.name = get.name(index = reference.index), mat.name = get.name(index = index))
			}
			graph <- function(fun, ...)
			{
				x <- i
				y <- get.y(j = index) # mats[[index]][moment.name, ]
				if(length(x) != length(y))
				{ message("bad news"); browser()}
				fun(x = x, y = y, ...)
			}
			if(!histogram)
				graph(fun = if(index == 1) function(...){plot(main = "diagnostics", xlab = "index", ylab = moment.name, ylim = ylim, ...)} else points, pch = index, col = index, ...)
		}
		if(!histogram) legend(legend = names(object), x = "topleft", pch = i, col = i)
	}
	par(mfrow = c(2, 2))
	plot.moment(moment.name = "mean")
	plot.moment(moment.name = "sd")
}

parametricBootstrap <- function(object, ...){printGeneric("parametricBootstrap"); browser()}
#removeMethods("parametricBootstrap")
setMethod("parametricBootstrap", signature(object = "numeric"), function(object, ...)
{
	sample.mean <- mean(object, na.rm = TRUE)
	sample.sd <- Sd(object, mle = FALSE, na.rm = TRUE) # sd(object, na.rm = TRUE) changed 25 and 28 (before 6:00 pm) April 2008 (cf. pp. 7, 9 of Peter Hall)
	normality.failed <- !is.finite(sample.mean) || !is.finite(sample.sd) || sample.sd == 0
	get.vec <- function()
	{
		num <- object
		boo <- !is.na(num)
		bootstrap.sample <- if(normality.failed)
		{
			function(n)
			{
				warning("missing bootstrap values returned")
				as.numeric(rep(NA, n)) # this might be replaced by a nonparametric bootstrap, but see nonparametricBootstrap
			}
		}
		else
		{
			function(n)
			{
				rnorm(mean = sample.mean, sd = sample.sd, n = n) # added 28 April 2008
			}
		}
		num[boo] <- bootstrap.sample(n = sum(boo)) # added 28 April 2008
		num
	}
	vec <- get.vec()
	stopifnot(length(vec) == length(object))
	names(vec) <- names(object)
#	if(normality.failed)
#		stopifnot(all(is.na(vec)))
	if(!normality.failed && !all(is.na(vec) == is.na(object)))
	{
		message("bad vec")
		browser()
	}
	vec
})
setMethod("parametricBootstrap", signature(object = "matrix"), function(object, ...)
{
	mat <- t(sapply(1:nrow(object), function(i)
	{
		parametricBootstrap(object = object[i, ], ...)
	}))
	stopifnot(all(dim(mat) == dim(object)))
	rownames(mat) <- rownames(object)
	colnames(mat) <- colnames(object)
	mat
})
setMethod("parametricBootstrap", signature(object = "ExpressionSet"), function(object, ...)
{
	exprs(object) <- parametricBootstrap(exprs(object))
	object
})
setMethod("parametricBootstrap", signature(object = "xprnSet"), function(object, ...)
{
	object@es <- parametricBootstrap(object@es, ...)
	object
})
setMethod("parametricBootstrap", signature(object = "XprnSet"), function(object, ...)
{
	parametricBootstrap(logb(object), ...)
})
setMethod("parametricBootstrap", signature(object = "xprnSetPair"), function(object, ...)
{
	get.es <- function(es){parametricBootstrap(es)}
	new("xprnSetPair", x = get.es(es = object@x), y = get.es(es = object@y))
})

beep <- function(n, pause)
{
	if(missing(n))
		n <- 1
	if(missing(pause))
	{
		pause <- if(n == 1)
			0
		else
			1e7
	}
	for(i in 1:n)
	{
		if(i > 1)
		{
			for(j in 1:pause)
				"a"
		}
		alarm()
	}
}

CombinedNames <- function(object, xlab = "x", ylab = "y")
{
	if(is(object, "character"))
	{
		if(xlab == ylab)
			stop("CombinedNames xlab == ylab")
		dup0 <- duplicated(object)
		if(any(dup0))
		{ message("duplication in object:"); print(dup0); browser()}
		subnam <- function(lab)
		{
			stopifnot(length(lab) == 1)
			paste(lab, object, sep = ".")
		}
		nam <- c(subnam(xlab), subnam(ylab))
		stopifnot(length(nam) == 2 * length(object))
		dup <- duplicated(nam)
		if(any(dup))
		{ message("duplication in nam:"); print(dup); browser()}
		nam
	}
	else if(is.null(object))
	{
		warning("No names were found; returning the object")
		NULL
	}
	else # adds combined names to object
		stop("bad class")
#	{
#		names(object) <- CombinedNames(object = names(object), xlab = xlab, ylab = ylab)
#		object
#	}
}
Combine <- function(x, y, ...){printGeneric("Combine"); browser()}
#removeMethods("Combine")
setMethod("Combine", signature(x = "Vector", y = "Vector"), function(x, y, ...)
{
	combo <- c(x, y, ...)
	if(sameNames(x, y))
		names(combo) <- CombinedNames(names(x))
	if(is(names(combo), "character"))
		names(combo) <- make.names(names(combo), unique = TRUE)
	combo
})
setMethod("Combine", signature(x = "Numeric", y = "Numeric"), function(x, y, ...)
{
	num <- function(object){as(object, "numeric")}
	Numeric(Combine(x = num(x), y = num(y), ...))
})
setMethod("Combine", signature(x = "ttest", y = "ttest"), function(x, y, xlab = stop("no xlab"), ylab = stop("no ylab"), ...) # like x = "cvalue", y = "cvalue"
{
	assert.is(xlab, "character")
	assert.is(ylab, "character")
	stopifnot(length(xlab) == 1)
	stopifnot(length(ylab) == 1)
	if(xlab == ylab)
		stop("Combine xlab == ylab")
	same <- function(name)
	{
		identical(slot(x, name), slot(y, name))
	}
	stopifnot(same("alternative"))
	stopifnot(same("level2"))
	len <- length(x) + length(y)
	nam <- if(all(names(x) == names(y)))
		CombinedNames(names(x), xlab = xlab, ylab = ylab)
	else
		stop("x and y have different names")
	get.vec <- function(name)
	{
		get.subvec <- function(object, lab)
		{
			z <- try(slot(object, name = name))
			if(is(z, "try-error"))
			{ message("slot(object, name = name) error!"); browser()}
			names(z) <- paste(lab, names(z), sep = ".") # sapply(names(z), function(na){paste(lab, na, sep = ".")})
			z
		}
		p1 <- get.subvec(x, lab = xlab)
		p2 <- get.subvec(y, lab = ylab)
		zvalue <- c(p1, p2)
		dup <- duplicated(names(zvalue))
		if(any(dup))
		{ message("duplication in names(zvalue):"); print(dup); browser()}
		if(any(nam != names(zvalue)))
		{ message("failed CombinedNames"); browser()}
		zvalue
	}
	x@pvalue <- Numeric(get.vec(name = "pvalue"))
	x@stat <- get.vec(name = "stat")
	x@df <- Numeric(get.vec(name = "df"))
	x@level1 <- paste(x@level1, y@level1, sep = ".")
	stopifnot(validObject(x))
	if(length(x) != len)
	{ message("length(x) != len"); browser()}
	x
})



Smooth <- function(object, smooth, f) # in "plot", signature(x = "numeric", y = "biasEstimate") before 9 May 2008
{
	stopifnot(is(object, "biasEstimate"))
	if(missing(smooth) || is.null(smooth))
	{
		smooth <- TRUE
		message("smooth = ", smooth)
	}
	if(!is.logical(smooth) || smooth)
	{
		if(!is.function(smooth))
		{
			if(missing(f) || is.null(f))
			{
				f <- 1 / 100
				message("f = ", f)
			}
			smooth <- smoothDataFUN(f = f) # smoothData
			message("smooth = [default function]")
		}
		ann <- annotation(object)
		object <- smooth(object)
		if(!is(object, "biasEstimate")) stop("smooth object lost its class")
		smooth.ann <- try(annotation(object), silent = TRUE)
		if(length(smooth.ann) == 0 || ann != smooth.ann) stop("smooth object lost its annotation, ", ann)
	}
	object
}

Domain <- function(object, ...){printGeneric("Domain"); browser()}
#removeMethods("Domain")
setMethod("Domain", signature(object = "Density"), function(object, ...)
{
	den <- s3(object)
	x <- den$x[is.finite(den$x)]
	dom <- range(x, na.rm = TRUE)
	ok <- is.numeric(dom) && length(dom) == 2 && dom[1] < dom[2]
	if(!ok)
	{ message("bad Domain"); browser()}
	dom
})
Range <- function(object, ...){printGeneric("Range"); browser()}
#removeMethods("Range")
setMethod("Range", signature(object = "Density"), function(object, positive, ...)
{
	den <- s3(object)
	y <- den$y[is.finite(den$y)]
	if(positive)
		y <- y[y > 0]
	ra <- range(y, na.rm = TRUE)
	ok <- is.numeric(ra) && length(ra) == 2 && ra[1] < ra[2]
	if(!ok)
	{ message("bad Range"); browser()}
	ra
})

new.Density <- function(s3, ann, ...)
{
	if(length(s3) > 1)
		s3 <- list(s3 = s3)
	new("Density", s3 = s3, annotation = ann, ...)
}
Density <- function(x, ...){printGeneric("Density"); browser()}
#removeMethods("Density")
setMethod("Density", signature(x = "numeric"), function(x, ann, ...)
{
	s3 <- density(x = x[is.finite(x)], ...)
	if(missing(ann))
	{
		ann <- try(annotation(x))
	}
	if(is(ann, "try-error") || ann == "S3 object")
	{
		ann <- date()
		message("Density ann = ", ann)
	}
	new.Density(s3 = s3, ann = ann)
})
Density.fun <- function(object, reflection.point = 0, fdr = default(FALSE, "fdr"), ...)
{
	assert.is(object, "Density")
	fun <- as(object, "function")
	if(fdr)
		function(v)
		{
			lfdr <- Density.fun(object = object, reflection.point = reflection.point, fdr = FALSE)(v) / fun(v)
			ifelse(lfdr <= 1, lfdr, 1)
		}
	else if(!is.numeric(reflection.point))
		fun
	else
		function(v)
		{
			reflected <- ifelse(v > reflection.point, -v, v)
			fun(reflected)
		}
}

# used by extended.r and/or extended.s (in estimate.s prior to 100316):
setMethod("mean", signature(x = "xprnSet"), function(x, ...) # simpler version: stat of data.s
{
	#if(is(try(when.estimate.last.loaded), "try-error"))
	#	stop("estimate.r is needed; consider implementing this using stat of data.s")
	#else
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

new.functions <- function(object)
{
	if(is(object, "function"))
		functions(object)
	else if(is(object, "functions"))
		object
	else if(is.list(object))
	{
		funs <- new("functions", object)
		names(funs) <- names(object)
		funs
	}
	else
		stop("new.functions error")
}
functions <- function(...)
{
	lis <- list(...)
	if(length(lis) == 0 || is(lis[[1]], "function"))
		new.functions(lis)
	else if(is(lis[[1]], "list"))
		new.functions(...)
	else
		stop("functions error")
}

closest <- function(object, target, ...){printGeneric("closest"); browser()}
#removeMethods("closest")
setMethod("closest", signature(object = "numeric", target = "numeric"), function(object, target)
{
	if(length(object) == 1 && length(target) > 1)
		closest(target, object)
	else if(length(target) == 1)
	{
		dis <- abs(object - target)
		boo <- dis == min(dis, na.rm = TRUE)
		boo <- boo & !is.na(boo)
		if(is.null(names(object)))
			names(object) <- as.character(which(boo))
		object[boo]
	}
	else if(length(target) > 1)
	{
		sapply(target, closest, object = object)
	}
	else
		stop("bad closest")
})

test <- function(object, ...){printGeneric("test"); browser()}
#removeMethods("test")

setSeed <- function(seed, kind = "Knuth-TAOCP-2002", verbose = TRUE)
{
	if(missing(seed) || length(seed) == 0 || is.null(seed))
	{
		seed <- round(1000 * runif(1), 0) + 1
		verbose <- TRUE
	}
	if(verbose)
		message("setSeed seed = ", seed, " ", kind, " on ", date())
	set.seed(seed = seed, kind = kind)
}
seed.random <- function(...){message("calling setSeed"); setSeed(...)}

memory <- function(x, bytes, total.only, verbose)
{
	nam <- slotNames(x)
	if(missing(verbose))
	{
		verbose <- TRUE
		if(verbose)
			message("verbose = ", verbose)
	}
	if(missing(total.only))
	{
		total.only <- FALSE
		if(verbose && is.character(nam))
			message("total.only = ", total.only)
	}
	if(missing(bytes))
		bytes <- 1e6
	if(verbose && bytes == 1e6)
		message("object sizes given in MB")
	else if(verbose && bytes == 1e9)
		message("object sizes given in GB")
	lis <- list(total = round(object.size(x) / bytes, 1))
	if(!total.only && is.character(nam) && lis[[1]] >= 1)
	{
		lis2 <- lapply(nam, function(name)
		{
			memory(slot(x, name = name), bytes = bytes, total.only = total.only, verbose = FALSE)
		})
		stopifnot(length(lis2) == length(nam))
		names(lis2) <- nam
		lis <- c(lis, lis2)
	}
	lis
}

Verbose <- FALSE
Sum <- function(object, ...){printGeneric("Sum"); browser()}
#removeMethods("Sum")
setMethod("Sum", signature(object = "numeric"), function(object, na.rm)
{
	if(missing(na.rm))
		na.rm <- TRUE
	fun <- function(...){sum(..., na.rm = na.rm)}
	do.call("fun", as.list(object))
})
Sum.list <- function(object)
{
	if(Verbose)
		message("calling Sum.list for ", length(object), " elements each of class ", class(object[[1]]))
	tot <- if(length(object) == 1)
		object[[1]]
	else if(length(object) > 1)
		try(object[[1]] + Sum(object[2:length(object)]))
	else
		stop("Sum error 1")
	if(is(tot, "try-error"))
	{ message("Sum error 2"); browser()}
	tot
} # setMethod("Sum", signature(object = "list"), Sum.list)
setMethod("Sum", signature(object = "list"), function(object, na.rm)
{
	if(missing(na.rm))
		na.rm <- FALSE
	if(na.rm)
		object <- object(!is.na(object))
	tot <- object[[1]]
	for(i in 2:length(object))
		tot <- add(tot, object[[i]]) # tot + object[[i]]
	if(is(object[[1]], "function"))
		assert.is(tot, "function")
	tot
})

multiply <- function(e1, e2, ...){e1 * e2}
#removeMethods("multiply")
setMethod("multiply", signature(e1 = "numeric", e2 = "function"), function(e1, e2)
{
	function(param)
	{
		e1 * e2(param)
	}
})
setMethod("multiply", signature(e1 = "function", e2 = "function"), function(e1, e2)
{
	function(param)
	{
		e1(param) * e2(param)
	}
})
divide <- function(e1, e2, ...){e1 / e2}
#removeMethods("divide")
setMethod("divide", signature(e1 = "numeric", e2 = "function"), function(e1, e2)
{
	function(param)
	{
		e1 / e2(param)
	}
})
setMethod("divide", signature(e1 = "function", e2 = "numeric"), function(e1, e2)
{
	multiply(e1 = 1 / e2, e2 = e1)
})
setMethod("divide", signature(e1 = "function", e2 = "function"), function(e1, e2)
{
	multiply(e1 = divide(e1 = 1, e2 = e2), e2 = e1)
})
add <- function(e1, e2, ...){e1 + e2}
#removeMethods("add")
setMethod("add", signature(e1 = "function", e2 = "function"), function(e1, e2)
{
	function(param)
	{
		e1(param) + e2(param)
	}
})
setMethod("Mean", signature(x = "list"), function(x)
{
	tot <- Sum(x, na.rm = FALSE)
	divide(e1 = tot, e2 = length(x))
})
setMethod("Mean", signature(x = "functions"), function(x, ...)
{
	av <- Mean(as(x, "list"), ...)
	assert.is(av, "function")
	av
})
setMethod("geoMean", signature(x = "list"), function(x, ...)
{
	ln.x <- ln(x)
	gm <- anti.ln(Mean(ln.x, ...))
	if(is(x[[1]], "function"))
		assert.is(gm, "function")
	gm
})
setMethod("geoMean", signature(x = "functions"), function(x)
{
	av <- geoMean(as(x, "list"))
	assert.is(av, "function")
	av
})

ln <- function(x, ...){printGeneric("ln"); browser()}
#removeMethods("ln")
setMethod("ln", signature(x = "numeric"), function(x, base = exp(1), translation = 0)
{
	stopifnot(base == exp(1))
	logb(x = x + translation, base = base)
})
setMethod("ln", signature(x = "function"), function(x, ...)
{
	function(param)
	{
		ln(x(param), ...)
	}
})
setMethod("ln", signature(x = "list"), function(x, ...)
{
	lapply(x, ln, ...)
})
setMethod("ln", signature(x = "matrix"), function(x, ...)
{
	mat <- t(sapply(1:nrow(x), function(i){ln(x[i, ], ...)}))
	stopifnot(all(dim(mat) == dim(x)))
	dimnames(mat) <- dimnames(x)
	mat
})
setMethod("ln", signature(x = "data.frame"), function(x, ...)
{
	ln(as.matrix(x), ...)
})


anti.ln <- function(x, ...){exp(x = x, ...)}
#removeMethods("anti.ln")
setMethod("anti.ln", signature(x = "function"), function(x, translation = 0, ...)
{
	function(param)
	{
		exp(x(param), ...) - translation
	}
})

delta.ln.from.datf <- function(x, samples1.j, samples2.j, ...)
{
	assert.is(x, "data.frame")
	stopifnot(length(samples1.j) == length(samples2.j))
	get.samples <- function(j){x[, j]}
	mat <- delta.ln(get.samples(samples1.j), get.samples(samples2.j), ...)
	oks <- try(is.matrix(mat) && all(sapply(1:ncol(mat), function(j){is.numeric(mat[, j])})))
	ok <- !is(oks, "try-error") && !any(is.na(oks)) && all(oks==T)
	if(!ok)
	{ message("bad mat!"); browser()}
	mat
}
delta.ln <- function(x, y, ...)
{
	cla <- if(class(x) == class(y))
		class(x)
	else
		stop("incompatible classes")
	if(is(x, "data.frame") || is(x, "matrix"))
		stopifnot(all(dim(x) == dim(y)))
	else if(is(x, "numeric"))
		stopifnot(length(x) == length(y))
	else
		stop("bad classes x, y")
	get.ln <- function(object){ln(object, ...)}
	get.ln(x) - get.ln(y)
}


#bind <- function(x, y, ...)
#{
#	printGeneric("bind")
#}
#removeMethods("merge")
setMethod("merge", signature(x = "numeric", y = "numeric"), function(x, y)
{
	vec <- c(x, y)
	if(is.character(names(vec)) && any(duplicated(names(vec))))
	{
		warning("making names unique")
		names(vec) <- make.names(names(vec), unique = TRUE)
	}
	vec
})
setMethod("merge", signature(x = "matrix", y = "matrix"), function(x, y, MARGIN, ...)
{
	fun <- if(MARGIN == 1)
		rbind
	else if(MARGIN == 2)
		cbind
	else
		stop("bad MARGIN")
	mat <- fun(x, y, ...)
	namfun <- if(MARGIN == 1)
		rownames
	else if(MARGIN == 2)
		colnames
	else
		stop("bad MARGIN 2")
	if(is.character(namfun(mat)) && any(duplicated(namfun(mat))))
	{
		warning("making dimnames unique")
		dimnames(mat)[[MARGIN]] <- make.names(namfun(mat), unique = TRUE)
	}
	mat
})

geomean <- function(x, ...){exp(mean(x = logb(x), ...))}
harmmean <- function(x, ...){rec <- function(y){1/y}; rec(mean(x = rec(x), ...))}

swap <- function(object, ...){printGeneric("swap"); browser()}
#removeMethods("swap")

in.interval <- function(x, lim, lower, upper, call.abs) # generalized 081002
{
	if(missing(call.abs))
	{
		call.abs <- FALSE
		message("in.interval call.abs = ", call.abs)
	}
	assert.is(call.abs, "logical")
	if(call.abs)
		x <- abs(x)
	if(missing(lim))
	{
		assert.is(lower, "numeric")
		assert.is(upper, "numeric")
		is.len.ok <- function(vec){length(vec) %in% c(1, length(x))}
		stopifnot(is.len.ok(lower) && is.len.ok(upper))
		x >= lower & x <= upper
	}
	else if(missing(lower) && missing(upper))
	{
		stopifnot(length(lim) == 2 && lim[1] <= lim[2])
		in.interval(x = x, lower = lim[1], upper = lim[2], call.abs = call.abs)
	}
	else
	{ message("bad argument combination"); browser()}
}

Is <- function(object, class2)
{
	any(sapply(class2, is, object = object))
}
assert.is <- function(object, class2, text = "")
{
	stopifnot(is.character(class2))
	if(missing(object))
		stop(paste(class2, "object missing in assert.is", text))
	stopifnot(length(class2) >= 1)
#	if(!is(object = object, class2 = class2))
#		stop(paste("got", class(object), "when", class2, "was required", text))
	if(!Is(object = object, class2 = class2))
	{
		warning(paste("got ", class(object), sep = ""))
		message("got ", class(object), " when one of these classes was required:")
		print(class2)
		stop(text)
	}
}
assert.are <- function(object, class2, ...)
{
	assert.is(object, "list")
	if(!is.nothing(object))
	{
		for(obj in object)
			assert.is(object = obj, class2 = class2, ...)
	}
}
are <- function(object, class2)
{
	all(sapply(object, Is, class2 = class2))
}
are.null <- function(object){are(object, "NULL")}

stopifnot.all <- function(...)
{
	ok <- all(...)
	arglis <- list(...)
	if(!ok)
	{ message("Not all ...; names(arglis):"); print(names(arglis)); browser()}
	stopifnot(ok)
}

new.argmax <- function(object, approx.max, from, to)
{
	new("argmax", scalar(object), approx.max = scalar(approx.max), from = scalar(from), to = scalar(to))
}
argmax <- function(object, from = stop("from missing"), to = stop("to missing"), length.out, length.multiplier, max.relevant, tolerance, max.iter, max.length)
{
	assert.is(object, "function")
	if(missing(length.out))
		length.out <- 100
	assert.is(length.out, "numeric")
	if(missing(length.multiplier))
		length.multiplier <- 5
	if(missing(max.length))
		max.length <- 10 ^ 12
	if(missing(max.iter))
		max.iter <- min(1000, floor(log(max.length / length.out) / log(length.multiplier)))
	if(missing(max.relevant) || length(max.relevant) == 0)
		max.relevant <- Inf
	if(missing(tolerance) || length(tolerance) == 0)
		stop("tolerance missing")
	arglis <- list(from = from, to = to, length.multiplier = length.multiplier, length.out = length.out, max.iter = max.iter, max.relevant = max.relevant, tolerance = tolerance)
	tr <- try(assert.are(arglis, "numeric"))
	if(is(tr, "try-error"))
	{ message("bad argmax arguments"); print(arglis); browser()}
	if(any(is.na(c(from, to))) || from > to)
	{ message("bad from, to"); browser()}
	get.argmax <- function(length.out)
	{
		if(length.out > max.length)
		{ message("length.out too high"); browser()}
		argvec <- seq(from = from, to = to, length.out = length.out)
		mapped <- object(argvec)
		if(any(is.na(mapped)))
		{
			message("bad mapped")
			browser()
		}
		boo <- mapped == max(mapped)
		if(!any(boo))
		{ message("no max."); browser()}
		if(sum(boo) > 1)
		{
			warning("no unique max.; taking median")
		}
		stopifnot(sum(boo) >= 1)
		median(argvec[boo])
	}
	oldArgmax <- newArgmax <- get.argmax(length.out = length.out)
	iter <- 1
	is.iter.low <- function(iter){iter <= max.iter}
	low.iter <- is.iter.low(iter = iter)
	while(low.iter && object(newArgmax) < max.relevant && (iter == 1 || abs(object(newArgmax) - object(oldArgmax)) > tolerance))
	{
		oldArgmax <- newArgmax
		newArgmax <- get.argmax(length.out = length.out * length.multiplier)
		iter <- iter + 1
		low.iter <- is.iter.low(iter = iter)
	}
	if(!low.iter)
	{	message("exited because exceeded ", max.iter, " iterations"); browser()}
	if(newArgmax[1] %in% c(from, to))
	{
		message("boundary value returned")
		browser()
	}
	new.argmax(newArgmax[1], approx.max = object(newArgmax[1]), from = from, to = to)
}
new.invert <- function(object, value, approx.value, fun)
{
	new("invert", object, value = scalar(value), approx.value = scalar(approx.value), fun = fun)
}
invert <- function(object, value, value.tolerance, ...) # from, to, length.out, length.multiplier, max.relevant, tolerance, max.iter, max.length
{
	assert.is(object, "function")
	assert.is(value, "numeric")
	assert.is(value.tolerance, "numeric")
	stopifnot(length(value) == 1)
	function.to.maximize <- function(x){-abs(value - object(x))} # patterned after cdf.quantile in "cdf.s"
	am <- argmax(object = function.to.maximize, ...)
	approx.value <- object(am)
	delta <- abs(value - approx.value)
	if(delta > value.tolerance)
	{
		message("value.tolerance exceeded")
		browser()
	}
	new.invert(am, value = value, approx.value = approx.value, fun = object)
}
Invert <- function(object, from = 0, to = 1, arg.tolerance = (to - from) / 1e4, ...)
{
	assert.is(object, "function")
	#before April2014: assert.is(value.tolerance, "numeric") #XXX|:no visible binding for global variable Ôvalue.toleranceÕ
	assert.is(arg.tolerance, "numeric")
	function(value)
	{
		assert.is(value, "numeric")
		y <- sapply(value, function(val)
		{
			invert(object = object, value = val, value.tolerance = arg.tolerance, from = from, to = to, ...)
		})
		names(y) <- names(value)
		y
	}
}

inversefun <- function(object, domain0, ...)
{
	assert.is(object, "function")
	tr <- try(assert.is(domain0, "numeric"))
	if(is(tr, "try-error"))
	{ message("error in domain0"); browser()}
	range0 <- sapply(domain0, object)
	assert.is(range0, "numeric")
	i <- 1
	af <- NULL
	while(!is.function(af) && i < 1000)
	{
		af <- try(approxfun(x = range0, y = domain0, ...))
		i <- i + 1
	}
	if(is(af, "try-error"))
	{ message("error in af; i:", i); browser()}
	af
}

overload <- function(f, ...)
{
	assert.is(f, "character")
	if(isGeneric(f, ...))
		stop(paste(f, "is generic and thus cannot be overloaded"))
	else
	{
		fun <- eval(parse(text = f))
		if(!is(fun, "function"))
			stop(paste(f, "does not exist and thus cannot be overloaded"))
		fun
	}
}

check.overload <- function(...)
{
	lis <- list(...)
	if(length(lis) >= 2)
	{
		if(any(sapply(lis[2:length(lis)], function(fun)
		{
			assert.is(fun, "function")
			identical(fun, lis[[1]])
		})))
			stop("overload error; try calling source(..., reload = TRUE)")
	}
	else
		stop("check.overload has less than 2 arguments")
}

impute <- function(object, ...) # for Bioinformatics revision of 15 Oct. 2008
{
	rule <- 2
	neighbors.equal <- TRUE
	if(neighbors.equal)
		warning("nonstandard linear imputation; consider changing impute function")
	if(is(object, "xprnSet"))
	{
		mat <- impute(exprs(object), ...)
		tr <- try(exprs(object@es) <- mat)
		if(is(tr, "try-error"))
		{ message("bad exprs assignment"); browser()}
		object
	}
	else if(is(object, "matrix"))
	{
		for(i in 1:nrow(object))
			object[i, ] <- impute(object[i, ], ...)
		object
	}
	else if(is(object, "numeric"))
	{
		x <- 1:length(object)
		y <- object
		boo <- is.finite(y)
		if(sum(boo) < 2)
		{
			replacement <- 0
			warning(paste("replacing all values with", replacement))
			return(rep(0, length(object)))
		}
		fun <- try(approxfun(x = x[boo], y = y[boo], rule = rule))
		if(is(fun, "try-error"))
		{ message("error"); browser()}
		imputed <- fun(x)
		is.ok <- function(vec){all(is.finite(vec) & (!is.finite(y) | vec == y))}
		ok <- is.ok(vec = imputed)
		if(!ok)
		{ message("imputation error"); browser()}
		if(neighbors.equal)
		{
#			z <- object
			adjusted <- imputed
			j <- 2
			while(j < length(y))
			{
				k <- j + 1
				while(!is.finite(y[j]) && !is.finite(y[k]) && k < length(y))
				{
					adjusted[j:k] <- mean(imputed[j:k])
					k <- k + 1
				}
				stopifnot(is.finite(y[k]) || is.finite(y[j]) || k == length(y))
				j <- k
			}
			ok <- is.ok(vec = adjusted)
			if(!ok)
			{ message("adjustment error"); browser()}
			adjusted
		}
		else
			imputed
	}
	else
		stop("bad object class")
}

default <- function(object, name, verbose, return.value = object)
{
	if(missing(verbose))
		verbose <- TRUE
	if(verbose)
	{
		if(missing(name))
			name <- "actual argument"
		if(is.function(name))
			name <- name(object)
		stopifnot(length(name) == 1 && is.character(name))
		prefix <- paste(name, "was set to default value of ")
		suffix <- paste(" on ", date(), ".", sep = "")
		object.name <- if(isS4(object) && any(c("annotation", "ann") %in% slotNames(object)))
		{
			object.nam <- try(paste("object of annotation '", annotation(object), "'", sep = ""))
			if(!is(object.nam, "character"))
			{ message("bad object.nam from object0"); browser()}
			object.nam
		}
		else
			object
		if(length(object.name) == 1 && (is(object.name, "character") || is(object.name, "numeric")))
			message(prefix, object.name, suffix)
		else
		{
			cat("\n", prefix)
			print(object.name)
			cat(paste(suffix, "\n\n"))
		}
	}
	return.value
}

defaultFUN <- function(object, ...)
{
	return.value <- if(is(object, "character"))
	{
		eval(parse(text = object))
	}
	else if(is(object, "function"))
		object
	if(is(object, "character"))
		object <- paste("function of name", object)
	if(is(return.value, "function") || is(return.value, "sampleFunction") || is(return.value, "sampleFunctions"))
		default(object = object, return.value = return.value, ...)
	else
		stop(paste("return.value is ", class(return.value), ", not a function", sep = ""))
}

export <- function(object, ...)
{
# printGeneric("export"); browser()
	datf <- as(object, "data.frame")
	if(is.data.frame(datf))
		export(object = datf, ...)
	else
	{ message("cannot coerce ", class(object), " to data frame"); browser()}
}
#removeMethods("export")
setMethod("export", signature(object = "data.frame"), function(object, file, ...)
{
	if(missing(file))
		file <- paste("exported ", Sys.Date(), ".csv", sep = "")
	write.csv(x = object, file = file, ...)
})

pmean <- function(x, y, ...)
{
	lx <- length(x)
	ly <- length(y)
	# stopifnot(min(lx, ly) %in% c(1, max(lx, ly)))
	stopifnot(lx == ly)
	vec <- numeric(lx)
	for(i in 1:lx)
		vec[i] <- mean(c(x[i], y[i]), ...)
	vec
}

Function.slot <- function(object, name) # see "Function.slot", signature(object = "dposterior", name = "character") in cdf.s
{
	assert.is(object, "Function")
	if(name %in% slotNames(object))
		slot(object, name = name)
	else
	{ message("bad Function.slot call; class(object): ", class(object)); browser()}
}
min_param <- function(object){Function.slot(object = object, name = "min_param")}
max_param <- function(object){Function.slot(object = object, name = "max_param")}
param.name <- function(object){Function.slot(object = object, name = "param.name")}
domain <- function(object, large.value = 3) # useful for xlim in plot: see fiducial.r
{
	xlim <- if(is(object, "Function"))
		c(min_param(object), max_param(object))
	else if(is(object, "Density"))
		range(xi(object))
	else
		stop("bad class for domain")
	ifelse(is.finite(xlim), xlim, sign(xlim) * abs(large.value))
}


Approx <- function(x, y, ...){printGeneric("Approx"); browser()}
#removeMethods("Approx")
setMethod("Approx", signature(x = "numeric", y = "numeric"), function(x, y, ...)
{
	stopifnot(length(x) == length(y))
	boo <- is.finite(x) && is.finite(y)
	approx(x = x[boo], y = y[boo], ...)
})
setMethod("Approx", signature(x = "numeric", y = "matrix"), function(x, y, ...)
{
	if(length(x) == nrow(y))
	{
		lis <- lapply(1:ncol(y), function(j){Approx(x = x, y = y[, j], ...)})
		app <- list(x = lis[[1]]$x, y = sapply(lis, function(elem){elem$y}))
		if(!all(sapply(lis, function(elem){all(app$x == elem$x)})))
		{ message("bad app"); browser()}
		stopifnot(ncol(y) == ncol(app$y))
		app
	}
	else
	{
		if(length(x) == ncol(y))
		{
			app <- Approx(x = x, y = t(y), ...)
			app$y <- t(app$y)
			app
		}
		else
			stop("bad dim")
	}
})

collapse <- function(...){paste(..., collapse = "-")}

rindex <- function(prob)
{
	assert.is(prob, "numeric")
	if(length(prob) == 1 && prob > 1 && floor(prob) == prob)
		prob <- rep(1 / prob, prob)
	stopifnot(sum(prob) == 1)
	ran <- which(as.numeric(rmultinom(n = 1, size = 1, prob = prob)) == 1)
	stopifnot(length(ran) == 1 && !is.na(ran) && ran >= 1 && ran <= length(prob) && ran == floor(ran))
	ran
}

print2file <- function(..., file = paste("print2file", Sys.Date(), rindex(999), "txt", sep = "."), FUN = function(...){cat(..., "\n", sep = " ")}, append = TRUE, add.date = TRUE)
{
#	text <- paste(, ..., sep = " ")
	sink(file = file, append = append)
	if(add.date)
		cat("\n", date(), "\n")
	FUN(...)
	sink(NULL)
}

figs <- function(text, file = text, call.pdf = TRUE, call.postscript = TRUE)
{

	if(is.null(file))
		file <- paste("figs", Sys.Date(), rindex(999), sep = " ")
	assert.is(text, "character")
	eval.arg <- parse(text = text)
	mes <- function(file){message("Finished plotting to ", file, " on ", date(), ".")}
	if(call.pdf)
	{
		pdf.file <- paste(file, "pdf", sep = ".")
		pdf(file = pdf.file)
		eval(eval.arg)
		dev.off()
		mes(pdf.file)
	}
	if(call.postscript)
	{
		ps.file <- paste(file, " PS.ps", sep = "")
		postscript(file = ps.file)
		eval(eval.arg)
		dev.off()
		mes(ps.file)
	}
}

zO.ECDF <- function(object, lowest = -Inf, highest = Inf, ...)
{
	if(is(object, "matrix"))
	{
		function(q, i)
		{
			if(missing(i))
			{
				if(length(q) == 1)
				{
					p <- sapply(1:nrow(object), function(ro)
					{
						pfun1 <- ECDF(object[ro, ], ...)
						p1 <- pfun1(q)
						if(length(p1) != length(q))
						{ message("p1 err: ", length(p1)); browser()}
						p1
					})
					if(!is(p, "numeric"))
					{ message("p err: ", class(p)); browser()}
					stopifnot(length(p) == nrow(object))
					names(p) <- rownames(object)
					p
				}
				else if(length(q) > 1 && length(q) == nrow(object))
				{
					p <- sapply(1:length(q), function(ro) #{pboot(point[i], i = i)}
					{
						pfun1 <- ECDF(object[ro, ], ...)
						p1 <- pfun1(q[ro], i = ro)
						if(length(p1) != 1)
						{ message("p1 error: ", c(p1 = length(p1), q = length(q))); browser()}
						p1
					})
					if(!is(p, "numeric"))
					{ message("p error: ", class(p)); browser()}
					stopifnot(length(p) == nrow(object))
					names(p) <- rownames(object)
					p
				}
				else
					stop("bad q")
			} # end if(missing(i))
			else if(length(i) == 1)
				ECDF(object[i, ], ...)
			else
				stop("bad i")
		} # end function
	}
	else if(is(object, "numeric"))
	{
		assert.is(object, "numeric")
		object <- object[is.finite(object)]
		if(length(object) == 0)
			function(q)
			{
				as.numeric(NA) # c(lowest, highest)
			}
		else
		{
			approxfun(x = sort(object), y = (1:length(object)) / (length(object) + 1), ...)
		}
	}
	else
		stop("fatal args")
}
ECDF <- function(...)
{
	EQF(..., inverse = TRUE)
}
EQF <- function(object, lowest = -Inf, highest = Inf, inverse = FALSE, ...) # equant qdata; cf. ecdf
{
	recall <- function(obj){EQF(object = obj, lowest = lowest, highest = highest, inverse = inverse, ...)}
	if(is(object, "matrix"))
	{
		function(p, i)
		{
			if(missing(i))
			{
				if(length(p) == 1)
				{
					quant <- sapply(1:nrow(object), function(ro)
					{
						qfun1 <- recall(object[ro, ])
						qfun1(p)
					})
					if(!is(quant, "numeric"))
					{ message("quant err: ", class(quant)); browser()}
					stopifnot(length(quant) == nrow(object))
					names(quant) <- rownames(object)
					quant
				}
				else if(length(p) > 1 && length(p) == nrow(object))
				{
					q <- sapply(1:length(p), function(ro) #{pboot(point[i], i = i)}
					{
						qfun1 <- recall(object[ro, ])
						q1 <- try(qfun1(p[ro]))
						if(!is.numeric(q1) || length(q1) != 1)
						{ message("q1 error: ", class(q1)); print(c(len.q1 = length(q1), len.p = length(p))); browser()}
						q1
					})
					if(!is(q, "numeric"))
					{ message("q error: ", class(q)); browser()}
					stopifnot(length(q) == nrow(object))
					names(q) <- rownames(object)
					q
				}
				else
					stop("bad p")
			} # end if(missing(i))
			else if(length(i) == 1)
				recall(object[i, ])
			else
				stop("bad i")
		} # end function
	}
	else if(is(object, "numeric"))
	{
		object <- object[is.finite(object)]
		if(length(object) == 0)
			function(p)
			{
				as.numeric(NA) # c(lowest, highest)
			}
		else
		{
			p <- (1:length(object)) / (length(object) + 1)
		#	fun <- ECDF(object)
			q <- sort(object)
			stopifnot(lowest <= q[1] && highest >= q[length(q)])
			if(is.finite(highest))
			{
				q <- c(q, highest)
				p <- c(p, 1)
			}
			if(is.finite(lowest))
			{
				q <- c(lowest, q)
				p <- c(0, p)
			}
			af <- try(if(inverse) # ECDF
				approxfun(x = q, y = p, ...)
			else
				approxfun(x = p, y = q, ...)) # see Likelihood.CDF
			if(is(af, "try-error"))
			{ message("af err"); print(list(...)); browser()}
			af
		}
	}
	else
		stop("fatal args")
} # end EQF

decile1 <- function(...){quantile(..., probs = 0.1)}
decile9 <- function(...){quantile(..., probs = 0.9)}
nx11<-function(height = 3.5, width = 3.5,pointsize=12,...){par(pin=c(width, height),ps=pointsize,...)} #par(new=T);fin,pin
Mfrow <- function(half = FALSE, quarter = FALSE, mfrow, ...)
{
	assert.is(half, "logical")
	if(missing(mfrow))
	{
		mfrow <- if(half)
		{
			nx11(height = 3.5)
			c(1, 2)
		}
		else if(quarter)
		{
			nx11(height = 3.5, width = 3.5)
#			c(1, 1)
		}
		else
			c(2, 2)
	}
	else
		stopifnot(!half && !quarter)
	par(mfrow = mfrow, ...)
}

halfAlphaRelativeEntropy <- function(mean1, mean2, sd1, sd2) # see http://eom.springer.de/B/b110490.htm and h090831
{
	assert.are(list(mean1, mean2, sd1, sd2), "numeric")
	sd <- sqrt((sd1 ^ 2 + sd2 ^ 2) / 2)
	BDistance <- (mean1 - mean2) ^ 2 / (8 * sd ^ 2) + logb(sd ^ 2 / (sd1 * sd2)) / 2 # Bhattacharyya coefficient
	stopifnot(BDistance >= 0)
	HellingerIntegral <- exp(-BDistance)
	-2 * logb(HellingerIntegral, base = 2) # verified by "Renyi entropy.nb" and "conditional-entropy.r"
}

Bonferroni <- function(object, ...)
{
	if(is(object, "cvalue"))
		Bonferroni(object = congruity(object), ...)
	else
	{
		x <- if(is(object, "numeric"))
			object
		else
			as(object, "numeric")
		x <- x[is.finite(x)]
		stopifnot(all(x >= 0 & x <= 1))
		min(x) * length(x)
	}
}

isInteger <- function(x){floor(x) == x}
Seq <- function(from, to, return.na.on.error = FALSE)
{
	stopifnot(isInteger(from) && isInteger(to))
	if(from > to)
	{
		if(return.na.on.error)
			NA
		else
			stop("from > to")
	}
	from:to
}

sup <- function(object, ...){printGeneric("sup"); browser()}
inf <- function(object, ...){printGeneric("inf"); browser()}

new.binomMetalevel <- function(valid, nonconservative, corrected, x, size, prob1, prob2)
{
	new("binomMetalevel", valid = Numeric(valid), nonconservative = Numeric(nonconservative), corrected = Numeric(corrected), x = Numeric(x), size = Numeric(size), prob1 = Numeric(prob1), prob2 = Numeric(prob2))
}
binomMetalevel <- function(x, size, prob1, prob2, true.prob, ...)
{
	assert.are(list(size, prob1, prob2), "numeric")
	get.x <- function(prob){ceiling(size * prob)}
	if(missing(x))
	{
		stopifnot(length(true.prob) == 1 && true.prob <= 1 && true.prob >= 0)
		x <- get.x(prob = true.prob)
	}
	else
		stopifnot(missing(true.prob) || x == get.x(prob = true.prob))
	if(length(x) == 1)
		x <- rep(x, length(size))
	else if(length(size) == 1)
		size <- rep(size, length(x))
	stopifnot(length(x) == length(size))
	if(all(1 == c(length(prob1), length(prob2))))
	{
		valid <- nonconservative <- corrected <- as.numeric(rep(NA, length(x)))
		for(i in 1:length(x))
		{
			lev <- binom.certainty(x = x[i], size = size[i], prob1 = prob1, prob2 = prob2, valid = 2, ...)
			stopifnot(length(lev) == 1)
			valid[i] <- lev@valid
			nonconservative [i] <- lev@nonconservative
			corrected[i] <- lev@corrected
		}
		if(any(is.na(valid)))
		{ message("NA in valid"); browser()}
		if(any(is.na(nonconservative)))
		{ message("NA in nonconservative"); browser()}
		if(any(is.na(corrected)))
		{ message("NA in corrected"); browser()}
		new.binomMetalevel(valid = valid, nonconservative = nonconservative, corrected = corrected, x = x, size = size, prob1 = prob1, prob2 = prob2)
	}
	else if(length(x) == 1)
	{
		stopifnot(length(prob1) == length(prob2))
		lis <- lapply(1:length(prob1), function(i)
		{
			binomMetalevel(x = x, size = size, prob1 = prob1[i], prob2 = prob2[i], true.prob = true.prob, ...)
		})
		extract.vec <- function(name){sapply(lis, function(elem){slot(elem, name = name)})}
		extract.sca <- function(name)
		{
			vec <- extract.vec(name = name)
			sca <- vec[1]
			stopifnot(all(sca == vec))
			sca
		}
		bml <- new.binomMetalevel(valid = extract.vec("valid"), nonconservative = extract.vec("nonconservative"), corrected = extract.vec("corrected"), x = extract.sca("x"), size = extract.sca("size"), prob1 = extract.vec("prob1"), prob2 = extract.vec("prob2"))
		stopifnot(all(prob1 == bml@prob1) && all(prob2 == bml@prob2))
		bml
	}
	else
		stop("bad args!!!!!!")
}
setMethod("sup", signature(object = "binomMetalevel"), function(object)
{
	v <- object@valid
	nc <- object@nonconservative
	ifelse(v > nc, v, nc)
})
setMethod("inf", signature(object = "binomMetalevel"), function(object)
{
	v <- object@valid
	nc <- object@nonconservative
	ifelse(v < nc, v, nc)
})

Pbinom <- function(x, size, prob, lower.tail = TRUE, correct = default(TRUE, "Pbinom correct"), correction = if(correct) 1/2 else 0, inclusive = TRUE, verbose = FALSE)
{
	stopifnot(length(x) == 1 && length(size) == 1 && length(prob) == 1)
	stopifnot(x >= 0 && size >= 1 && x <= size) # added 090916
	get.xvec <- function(from, to){Seq(from = from, to = to, return.na.on.error = TRUE)}
	xvec <- if(lower.tail)
		get.xvec(0, if(inclusive) x else x - 1)
	else if(x <= size)
		get.xvec(if(inclusive) x else x + 1, size)
	else
		stop("bad Pbinom args")
	xvec.ok <- !any(is.na(xvec)) && all(xvec %in% 0:size)
	mass <- function(x){dbinom(x = x, size = size, prob = prob)}
	get.delta.mass <- function()
	{
		stopifnot(correct)
		correction * mass(x = x)
	}
	if(!xvec.ok && !inclusive)
	{
		Mass <- 0
		p <- if(correct)
			Mass + get.delta.mass()
		else
			Mass
		if(verbose && correct && p == 0)
		{ message("cannot correct"); browser()}
		p
	}
	else if(xvec.ok) # begin nonzero mass
	{
		Mass <- sum(sapply(xvec, mass))
		if(correct)
		{
			delta.mass <- get.delta.mass()
			corrected.p <- if(inclusive)
				Mass - delta.mass
			else
			{
				if(correction == 1/2)
					warning("arguments canceled each other's effect")
				Mass + delta.mass
			}
			if(verbose && corrected.p == 0)
			{ message("corrected.p is ", corrected.p); browser()}
			corrected.p
		}
		else
			Mass
	} # end nonzero mass
	else
		stop("bad xvec.ok")
} # end Pbinom
binom.limit <- function(x, size, p, correct = default(TRUE, "binom.limit correct"), ...) # inverse Pbinom for fixed x
{
	stopifnot(is.prob(p))
	stopifnot(length(p) == 1)
	cdf <- function(prob)
	{
		Pbinom(x = x, size = size, prob = prob, correct = correct, ...)
	}
	probs <- binom.prob(x = x, size = size, p = p, alternative = 2)
	stopifnot(length(probs) >= 2)
	stopifnot(is.prob(probs))
	opt.prob <- optimize(f = function(prob){abs(cdf(prob) - p)}, lower = min(probs), upper = max(probs))$minimum
	stopifnot(is.prob(opt.prob))
	opt.prob
}


binom.certainty <- function(x, size, prob1 = default(0, "prob1"), prob2 = default(1, "prob2"), valid = 2, verbose = FALSE, correct = FALSE, use.correction)
{
	assert.are(list(prob1, prob2), "numeric")
	assert.is(use.correction, "logical")
	stopifnot(length(prob1) == 1 && length(prob2) == 1)
	stopifnot(prob1 >= 0)
	stopifnot(prob2 <= 1)
	stopifnot(prob1 <= prob2)
	if(is.numeric(valid) && length(valid) == 1 && valid == 2)
	{
		recall <- function(valid, correct){binom.certainty(x = x, size = size, prob1 = prob1, prob2 = prob2, valid = valid, verbose = verbose, correct = correct, use.correction = use.correction)}
		valid.certainty <- recall(valid = TRUE, correct = FALSE)
		nonconservative.certainty <- recall(valid = FALSE, correct = FALSE)
		corrected.certainty <- recall(valid = TRUE, correct = TRUE)
#		vec <- c(valid.certainty = valid.certainty, nonconservative.certainty = nonconservative.certainty)
#		return(if(vec[1] <= vec[2]) vec else rev(vec))
		if(all(c(valid.certainty, nonconservative.certainty, corrected.certainty) == 0))
		{ message("all 0"); browser()}
		if(verbose && corrected.certainty == 0)
		{ message("corrected.certainty is ", corrected.certainty); browser()}
		new.binomMetalevel(valid = valid.certainty, nonconservative = nonconservative.certainty, corrected = corrected.certainty, x = x, size = size, prob1 = prob1, prob2 = prob2)
	}
	else
	{
		certainty <- if(use.correction)
		{
			stopifnot(length(x) == 1)
			upper.tail.prob <- function(prob)
			{
				correction <- if(valid) 0 else 1
				adjustment <- correction * dbinom(x = x, size = size, prob = prob)
				Pbinom(x, size = size, prob = prob, lower.tail = FALSE, correct = correct, inclusive = FALSE, verbose = verbose) + adjustment
			}
			prob.higher.given.prob1 <- upper.tail.prob(prob = prob1)
			prob.higher.given.prob2 <- upper.tail.prob(prob = prob2)
			prob.higher.given.prob2 - prob.higher.given.prob1
		}
		else # original
		{
			prob.extreme <- function(prob, inclusive)
			{
				Pbinom(x, size = size, prob = prob, lower.tail = FALSE, correct = correct, inclusive = inclusive, verbose = verbose)
			}
			prob.higher.given.prob1 <- prob.extreme(prob = prob1, inclusive = valid)
			prob.higher.given.prob2 <- prob.extreme(prob = prob2, inclusive = !valid)
			prob.higher.given.prob2 - prob.higher.given.prob1
		}
		print.certainty <- function()
		{
			message("\n", if(valid) "valid" else "nonconservative")
			message("x: ", x, "; size: ", size, "; prob1: ", prob1, "; prob2: ", prob2)
			message(prob.higher.given.prob2, " - ", prob.higher.given.prob1, " = ", certainty)
		}
		if(certainty < 0 || certainty > 1)
		{
			print.certainty()
			certainty <- if(certainty < 0)
				0
			else if(certainty > 1)
				1
			message("certainty reset to ", certainty, "\n")
		}
		ok <- certainty >= 0 && certainty <= 1
		if(!ok)
		{ print.certainty(); message("bad certainty: ", certainty); browser()}
		certainty
	}
} # end binom.certainty
binom.interval <- function(x, size, conf.level = default(0.95, "conf.level"), alpha = default((1 - conf.level) / 2, "alpha"), p1, p2, valid = TRUE, by = 1e-4, verbose = FALSE) # slower than binom.ci
{
	stopifnot(length(x) == 1)
	stopifnot(length(size) == 1)
	stopifnot(x <= size)
	if(missing(p1))
		p1 <- alpha
	if(missing(p2))
		p2 <- p1 + conf.level
	a.ok <- function(a){length(a) == 1 && a >= 0 && a <= 1}
	stopifnot(a.ok(p1) && a.ok(p2))
	stopifnot(p2 >= p1)
	if(is.numeric(valid) && length(valid) == 1 && valid == 2)
	{
		recall <- function(valid){binom.interval(x = x, size = size, conf.level = conf.level, alpha = alpha, valid = valid, by = by, verbose = verbose)}
		valid.interval <- recall(valid = TRUE)
		nonconservative.interval <- recall(valid = FALSE)
		interval <- if(valid == 1)
		{
			stop("averaged interval not yet implemented")
		}
		else if(valid == 2)
			list(valid.interval = valid.interval, nonconservative.interval = nonconservative.interval)
		else
			stop("bad valid")
		return(interval)
	}
	get.theta <- if(is.logical(valid) && length(valid) == 1)
	{
		tail.mass <- function(theta, lower.tail)
		{
			Pbinom(x = x, size = size, prob = theta, lower.tail = lower.tail, correct = FALSE, inclusive = valid)
		}
		function(lower.tail)
		{
			target.mass <- if(lower.tail) p1 else 1 - p2
			old.theta <- new.theta <- if(lower.tail) 1 else 0
			add.mass <- TRUE
			mass <- 0
			if(verbose)
				message("\nlower.tail == ", lower.tail, "; target.mass == ", target.mass)
			while(mass <= target.mass && a.ok(new.theta))
			{
				if(!a.ok(new.theta))
				{ message("bad new.theta: ", new.theta); browser()}
				old.theta <- new.theta
				new.theta <- new.theta + by * (if(lower.tail) -1 else +1)
				mass <- tail.mass(theta = new.theta, lower.tail = lower.tail)
				if(verbose && mass > 1 - 1e-3)
				{ message("mass is ", mass); browser()}
				if(verbose)
					print(c(old.theta = old.theta, new.theta = new.theta, mass = mass))
#				if(call.browser)
#				{ message("call.browser is ", call.browser); browser()}
			}
	#		stopifnot(mass == target.mass)
			if(!a.ok(old.theta))
			{ message("bad old.theta: ", old.theta); browser()}
			old.theta
		}
	}
	else
	{
		function(lower.tail)
		{
			stop("not ready to use qbinom")
		}
	}
	theta.minus <- get.theta(lower.tail = FALSE)
	theta.plus <- get.theta(lower.tail = TRUE)
	if(theta.minus > theta.plus)
	{
		distance <- function(theta){abs(theta - 1/2)}
		theta.minus <- theta.plus <- if(distance(theta.minus) > distance(theta.plus))
			theta.plus
		else
			theta.minus
	}
	stopifnot(theta.minus <= theta.plus)
	c(theta.minus, theta.plus)
} # end binom.interval
binom.prob <- function(x, size, p, alternative = character(0), correct = TRUE, ...)
{
	stopifnot(is.prob(p))
	if(length(p) > 1)
		return(sapply(p, function(component){binom.prob(x = x, size = size, p = component, alternative = alternative, correct = correct, ...)})) # correction = correction, call.qbinom = call.qbinom,
	stopifnot(length(p) == 1)
	if(all(p %in% c(0, 1)))
		return(p)
	if(TRUE) # missing(correction))
	{
		get.prob <- function(alternative)
		{
			conf <- if(alternative == "less")
				p
			else if(alternative == "greater")
				1 - p
			else
				stop("very bad alternative")
			i <- if(alternative == "less")
				2
			else if(alternative == "greater")
				1
			else
				stop("very bad alternative")
			limit <- sapply(conf, function(conf.level){binom.test(x = x, n = size, alternative = alternative, conf.level = conf.level, ...)$conf.int[i]})
			if(alternative == "less")
				limit
			else if(alternative == "greater")
				limit
			else
				stop("very bad alternative")
		}
		if(is.nothing(alternative))
		{
			binom.limit(x = x, size = size, p = p, correct = correct, ...)
		}
		else if(is.numeric(alternative) && alternative == 2)
			c(less = get.prob("less"), greater = get.prob("greater")) # mean
		else
			get.prob(alternative)
	}
#	else
#	if(is.numeric(correction) && length(correction) == 1 && is.nothing(alternative))
#	{
#		if(missing(correction))
#			correction <- 1 / 2
#		correction <- Scalar(correction)
#		stopifnot(is.prob(correction))
#		prob <- 1 - binom.limit(x = x, size = size, p = p, correction = correction, lower.tail = TRUE, ...) # confidence.lyx
#		if(!is.prob(prob))
#		{ message("bad prob:"); print(prob); browser()}
#		prob
#	}
	else
		stop("bad binom.prob arguments")
} # end binom.prob
binom.ci <- function(x, size, p1 = 1 / 40, p2 = 1 - p1, valid = TRUE, ...) # faster than binom.interval
{
	stopifnot(p1 <= p2)
	get.limit <- function(p, alternative){binom.prob(x = x, size = size, p = p, alternative = alternative, ...)}
	ci <- c(get.limit(p = p1, alternative = if(valid) "greater" else "less"), get.limit(p = p2, alternative = if(!valid) "greater" else "less"))
	stopifnot(ci[1] <= ci[2])
	ci
}
binom.rprob <- function(x, size, n, FUN = binom.prob, ...)
{
	sapply(runif(n = n), function(p){FUN(x = x, size = size, p = p, ...)})
}

maxEntropy <- function(object)
{
	assert.is(object, "binomMetalevel")
	maxent <- object@prob2 - object@prob1
	dis <- function(vec){abs(vec - maxent)}
	v <- object@valid
	nc <- object@nonconservative
	if(length(maxent) == 1)
		maxent <- rep(maxent, length(v))
	stopifnot(all(c(length(v), length(nc)) == length(maxent)))
	ifelse((v < maxent & nc > maxent) | (v > maxent & nc < maxent), maxent, ifelse(dis(v) < dis(nc), v, nc))
}
meanOverConvexSet <- function(object)
{
	assert.is(object, "binomMetalevel")
	v <- object@valid
	nc <- object@nonconservative
	stopifnot(length(v) == length(nc))
	return((v + nc) / 2)
}

new.ttest <- function(pvalue, stat, df, alternative, level1, level2)
{
	new("ttest", pvalue = Numeric(pvalue), stat = stat, df = Numeric(df), alternative = alternative, level1 = level1, level2 = level2)
}
ttest <- function(x, y, factor.name, level1, level2, alternative = "greater", ...)
{
	if(is(x, "XprnSet"))
		x <- logb(x)
	if(!missing(y) && is(y, "XprnSet"))
		y <- logb(y)
	assert.is(x, "xprnSet")
	if(missing(y) && !missing(factor.name))
	{
		get.sub <- function(level)
		{
			xprnSubset(object = x, level = level, factor.name = factor.name)
		}
		stopifnot(length(level2) == 1)
		get.ttest <- function(lev1)
		{
			stopifnot(length(lev1) == 1)
			ttest(x = get.sub(level = lev1), y = get.sub(level = level2), level1 = lev1, level2 = level2, alternative = alternative, ...)
		}
		if(length(level1) == 1)
			get.ttest(lev1 = level1)
		else if(length(level1) == 2)
		{
			t1 <- get.ttest(lev1 = level1[1])
			t2 <- get.ttest(lev1 = level1[2])
			com <- Combine(t1, t2, xlab = level1[1], ylab = level1[2])
			if(length(com) != length(t1) + length(t2))
			{ message("length(com) != length(t1) + length(t2)"); browser()}
			com
		}
		else
			stop("bad length(level1)")
	}
	else
	{
		nam <- featureNames(x)
		x <- exprs(x)
		y <- if(missing(y))
		{
			NULL
		}
		else if(missing(factor.name))
		{
			assert.is(y, "xprnSet")
			stopifnot(all(nam == featureNames(y)))
			exprs(y)
		}
		else
			stop("bad args!!!!!!!!!!!!!")
		if(is.null(y))
		{
			level1 <- "x"
			level2 <- "missing"
		}
		else
			stopifnot(nrow(x) == nrow(y))
		get.s3 <- function(i)
		{
			test <- function(...)
			{
				t.test(x = x[i, ], alternative = alternative, ...)
			}
			if(is.null(y))
				test(...)
			else
				test(..., y = y[i, ])
		}
		pvalue <- stat <- df <- as.numeric(rep(NA, nrow(x)))
		for(i in 1:length(stat))
		{
			s3 <- get.s3(i)
			pvalue[i] <- s3$p.value
			stat[i] <- s3$statistic
			df[i] <- s3$parameter
		}
		stopifnot(length(pvalue) == length(nam))
		names(pvalue) <- names(stat) <- names(df) <- nam
		new.ttest(pvalue = pvalue, stat = stat, df = df, alternative = alternative, level1 = level1, level2 = level2)
	}
}

# begin moved from cdf.s on 27 February 2010
new.rParameter <- function(funs, ann)
{
	rp <- new("rParameter", functions(funs), ann = ann)
	names(rp) <- names(funs)
	rp
}
rParameter <- function(rfunctions, ann, ...)
{
	assert.is(rfunctions, "list")
	assert.is(ann, "character")
	stop("this internal function rParameter seems to fail due to strange scoping (using the last rfunction for all?)")
	funs <- lapply(rfunctions, function(rfunction)
	{
		function(ndraw)
		{
			rfunction(ndraw = ndraw, ...)
		}
	})
	new.rParameter(funs = funs, ann = ann)
}
rMean <- function(size = default(2, "size"), sample.mean = default(0, "sample.mean"), sample.var = default(1, "sample.var"), ndraw = default(1e3, "ndraw"), ...) # ... not used but included for rParameter
{
	args <- list(...)
	if(length(args) > 0)
		warning(paste("unused args:", names(args)))
	z <- rt(n = ndraw, df = size - 1)
	se.sample.mean <- sqrt(sample.var / size)
	z * se.sample.mean + sample.mean
}
rVar <- function(size = default(2, "size"), sample.var = default(1, "sample.var"), ndraw = default(1e3, "ndraw"), type = default("CD", "type"), known.mean = default(FALSE, "known.mean"), verbose = FALSE, ...) # tested by compare.rVar in logical.r # ... not used but included for rParameter
{
	args <- list(...)
	if(length(args) > 0)
		warning(paste("unused args:", names(args)))
	if(verbose)
		message("\ngenerating ", ndraw, " draws of type ", type, " and known.mean ", known.mean, " on ", date())
	if(type == "CD") # h090415b, h090416b
	{
		if(known.mean)
			size * sample.var / rchisq(n = ndraw, df = size)
		else
			(size - 1) * sample.var / rchisq(n = ndraw, df = size - 1)
	}
	else if(type == "default") # Bernardo & Smith, p. 440
	{
		prec <- if(known.mean)
			rgamma(n = ndraw, shape = size / 2, rate = size * sample.var / 2)
		else
			rgamma(n = ndraw, shape = (size - 1) / 2, rate = (size - 1) * sample.var / 2)
		1 / prec
	}
	else
		stop("type not recognized")
}
rSd <- function(...){sqrt(rVar(...))}
rNormalParameter <- function(x, type = default("CD", "rNormalParameter type"), ...)
{
	assert.is(x, "numeric")
	sample.mean <- mean(x, na.rm = TRUE)
	sample.var <- var(x, na.rm = TRUE)
	size <- length(x)
#	rParameter(rfunctions = list(mean = rMean, var = rVar), ann = "default normal", sample.mean = sample.mean, sample.var = sample.var, size = size, type = type, known.mean = FALSE)
	mean.fun <- function(ndraw)
	{
		rMean(ndraw = ndraw, sample.mean = sample.mean, sample.var = sample.var, size = size)
	}
	var.fun <- function(ndraw)
	{
		rVar(ndraw = ndraw, sample.var = sample.var, size = size, type = type, known.mean = FALSE)
	}
	funs <- list(mean = mean.fun, var = var.fun)
	new.rParameter(funs = funs, ann = "default normal")
}
Ran <- function(object, name, ndraw)
{
	assert.is(object, "rParameter")
	assert.is(ndraw, "numeric")
	assert.is(name, "character")
	stopifnot(name %in% names(object))
	fun <- object[name][[1]]
	assert.is(fun, "function")
	vec <- fun(ndraw = ndraw)
	stopifnot(length(vec) == ndraw)
	vec
}
Ran.diff <- function(x, y, name, ndraw)
{
	R <- function(object){Ran(object = object, name = name, ndraw = ndraw)}
	R(x) - R(y)
}
# end moved from cdf.s on 27 February 2010

uncertain <- function(object, FUN, ...)
{
	recall <- function(x){uncertain(object = x, FUN = FUN, ...)}
	if(is(object, "XprnSet"))
		recall(logb(object))
	else if(is(object, "xprnSet"))
		recall(exprs(object))
	else if(is(object, "matrix"))
	{
		vec <- sapply(1:nrow(object), function(i)
		{
			recall(object[i, ])
		})
		stopifnot(length(vec) == nrow(object))
		names(vec) <- rownames(object)
		vec
	}
	else if(is(object, "numeric"))
	{
		sca <- FUN(object, ...)
		scalar(sca)
	}
}
new.nonparametricBootstrap <- function(object, original)
{
	new("nonparametricBootstrap", object, original = original)
}
nonparametricBootstrap <- function(object, FUN, nsample, factor.name = NULL, parametric = FALSE, bootstrap = TRUE)
{
	get.sample <- if(is.function(bootstrap) && is.null(factor.name)) # good way to use "empiricalNull" (a hierarchical resampling model)
		bootstrap
	else if(bootstrap) # non-hierarchical resampling model
		function(x){Sample(x, replace = TRUE, factor.name = factor.name, parametric = parametric)}
	else if(is.null(factor.name)) # not bootstrap
		NULL
	else
		stop("argument currently incompatible with non-bootstrap")
	get.col <- if(is.function(get.sample)) # bootstrap
		function(dummy)
		{
			sampled <- get.sample(object)
			if(!is(sampled, class(object)))
			{ message(class(object), " changed to ", class(sampled)); browser()}
			FUN(sampled)
		}
	else if(is.null(get.sample)) # not bootstrap
		function(dummy)
		{
			uncertain(object = object, FUN = FUN)
		}
	else
		stop("weird")
	mat <- sapply(1:nsample, get.col)
	if(is(mat, "numeric"))
		mat <- matrix(mat, ncol = nsample)
	assert.is(mat, "matrix")
	if(any(is.na(mat)))
	{
		stopifnot(nsample == ncol(mat))
		nsample.na <- try(sum(sapply(1:nsample, function(j)
		{
			boo <- try(all(is.na(mat[, j])))
			if(is.logical(boo))
				boo
			else
			{ message("boo error"); browser()}
		})))
		if(!is.numeric(nsample.na))
		{ message("nsample.na error"); browser()}
		nel.na <- sum(is.na(mat))
		message(nsample.na, " of ", nsample, " samples (", 100 * signif(nsample.na / nsample, 2), "%) & ", nel.na - nsample.na * nrow(mat), " of ", length(mat), " more elements NA.")
		if(nsample.na == nsample)
		{ message("nsample.na == nsample"); browser()}
	}
	original <- if(is.function(get.sample)) # bootstrap
		FUN(object)
	else if(is.null(get.sample)) # not bootstrap
		numeric(0) # as.numeric(rep(NA, nrow(mat)))
	else
		stop("strange")
	if(!is.matrix(mat))
	{ message("mat err"); browser()}
	if(length(original) > 0)
	{
		ok <- nrow(mat) == length(original)
		if(length(ok) != 1 || !ok)
		{ message("length(ok) != 1 || !ok"); print(list(nrow(mat), length(original))); browser()}
		tr <- try(rownames(mat) <- names(original))
		if(is(tr, "try-error"))
		{ message("try-error"); browser()}
	}
	new.nonparametricBootstrap(mat, original = original)
} # end nonparametricBootstrap

new.CI <- function(object, expectation, zeroCI, point)
{
	tr <- try(colnames(object) <- c("lo", "hi"))
	if(is(tr, "try-error"))
	{ message("cannot assign column names for confidence interval"); browser()}
	rownames(object) <- names(expectation)
	ci <- new("CI", object, expectation = expectation, zeroCI = zeroCI, point = point) # Numeric(expectation)
	stopifnot(validObject(ci))
	ci
}
BC.percentile <- function(xstar, point, lowest = 0, highest = 1, ptrans = pnorm, qtrans = qnorm, qboot, nsample = length(xstar), rule = 2, ...) # shao:tu:1995, p. 131
{
	assert.are(list(xstar, point), "numeric")
	if(!(length(xstar) > 1 && length(point) == 1))
	{ message("bad length in BC.percentile"); browser()}
	assert.are(list(ptrans, qtrans), "function")
	stopifnot(!xor(identical(ptrans, pnorm), identical(qtrans, qnorm)))
	EF <- function(fun){fun(xstar, rule = rule, lowest = lowest, highest = highest, ...)}
	if(missing(qboot))
		qboot <- try(EF(fun = EQF)) # EQF(xstar))
	if(is(qboot, "try-error"))
	{ message("bad qboot in BC.percentile"); browser()}
	pboot <- try(EF(fun = ECDF)) # ECDF(xstar, lowest = lowest, highest = highest, ...))
	if(is(pboot, "try-error"))
	{ message("bad pboot"); browser()}
	prob <- pboot(point) # sapply(1:length(point), function(i){pboot(point[i], i = i)}) # sapply(1:length(point), function(i){pboot[[i]](point[i])})
	assert.is(prob, "numeric")
	if(length(prob) != length(point))
	{ message("length(prob) != length(point)"); browser()}
	bias <- qtrans(prob)
	stopifnot(length(prob) == 1)
	get.param <- function(p)
	{
		if(length(p) != 1)
		{ message("length(p) != 1"); browser()}
		if(is.finite(bias))
			qboot(ptrans(qtrans(p) + 2 * bias))
		else
		{ message("!is.finite(bias); bias:", bias); browser()}
	}
	params <- sapply(1:nsample, function(i)
	{
		get.param(p = (i + 1) / (nsample + 1))
	})
	if(any(is.na(params)))
	{ message("bad params:"); print(params); browser()}
	params
}
CI <- function(object, ...){printGeneric("CI"); browser()}
setMethod("CI", signature(object = "ANY"), function(object, FUN, nsample, factor.name = NULL, parametric = FALSE, p1 = default(0.025, "p1"), p2 = default(1 - p1, "p2"), ptrans, qtrans, bootstrap = TRUE, ...) # Invert.arglis,
{
	xstar <- nonparametricBootstrap(object = object, FUN = FUN, nsample = nsample, factor.name = factor.name, parametric = parametric, bootstrap = bootstrap)
	point <- xstar@original
	qboot <- try(EQF(xstar, ...))
	if(is(qboot, "try-error"))
	{ message("bad qboot in CI"); browser()}
	if(missing(ptrans) && missing(qtrans)) # percentile
	{
		get.param <- qboot
		params <- xstar
	}
	else if(is.function(ptrans) && bootstrap) # BC percentile, shao:tu:1995
	{
		warning("changed 110420 to call BC.percentile; without testing")
		BC.percentile(xstar = xstar, point = point, ptrans = ptrans, qtrans = qtrans, qboot = qboot, nsample = nsample, ...)
	}
	else
		stop("bad arg")
	ok <- any(is.finite(as(params, "numeric")))
	if(!ok)
	{ message("!ok params"); browser()}
	expectation <- Mean(params)
	if(!(length(point) %in% c(0, length(expectation))))
	{ message("cannot make CI"); browser()}
	Param <- function(p, name)
	{
		param <- get.param(p)
		if(any(is.nan(param)))
		{ message("bad param of name ", name, "; consider increasing the number of samples or looking at intervals of lower confidence"); browser()}
		if(length(param) != length(expectation))
		{ message("length(param) != length(expectation)"); browser()}
		param
	}
	zeroCI <- Param(p = 0.5, name = "zeroCI") # [median,median] is 0% CI
	lo <- Param(p = p1, name = "lo")
	hi <- Param(p = p2, name = "hi")
	stopifnot(length(lo) == length(hi))
	stopifnot(length(lo) == length(zeroCI))
	nam <- names(lo)
	if(!is.null(nam))
		stopifnot(all(nam == names(hi)))
	mat <- cbind(lo, hi)
	names(zeroCI) <- rownames(mat) <- nam
	stopifnot(all(names(expectation) == rownames(mat)))
	new.CI(mat, expectation = expectation, zeroCI = zeroCI, point = point) # new.CI(mat, point = FUN(object))
})
new.Coverage <- function(cover, interest.value, covered.value, nominal.rate, base = numeric(0))
{
	if(all(covered.value == interest.value))
		covered.value <- matrix(numeric(0))
	if(!is.matrix(covered.value))
	{
		assert.is(covered.value, "numeric")
		matrix(covered.value, ncol = 1)
	}
	new("Coverage", cover, interest.value = interest.value, covered.value = covered.value, nominal.rate = Scalar(nominal.rate), base = Scalar(base))
}

new.pointEstimate <- function(expectation, zeroCI, point)
{
	pe <- new("pointEstimate", expectation = expectation, zeroCI = zeroCI, point = point)
	stopifnot(validObject(pe))
	pe
}

subtract <- function(x, y, ...){printGeneric("subtract"); browser()}
#removeMethods("subtract")
setMethod("subtract", signature(x = "CI", y = "CI"), function(x, y, point = TRUE, ...)
{
	get.pair <- function(name)
	{
		if(missing(name))
		{
			get.mat <- function(object)
			{
				mat <- as(object, "matrix")
				if(length(mat) == 0)
				{ message("empty matrix"); browser()}
				mat
			}
			list(x = get.mat(x), y = get.mat(y))
		}
		else
			list(x = slot(x, name = name), y = slot(y, name = name))
	}
	get.diff <- function(...)
	{
#		if(missing(pair))
#			pair <- get.pair(name = name)
#		else
#			stopifnot(missing(name))
		pair <- get.pair(...)
		ok <- if(is.matrix(pair$x) && is.matrix(pair$y))
			(nrow(pair$y) == 1 && ncol(pair$x) == ncol(pair$y)) || all(dim(pair$x) == dim(pair$y))
		else if(is.numeric(pair$x) && is.numeric(pair$y))
			length(pair$y) %in% c(1, length(pair$x))
		else
			FALSE
		if(!ok)
		{ message("cannot subtract"); browser()}
		dif <- if((is.numeric(pair$x) && is.numeric(pair$y)) || nrow(pair$x) == nrow(pair$y))
			pair$x - pair$y
		else
		{ message("subtract is not yet implemented for matrices of differing numbers of rows"); browser()}
		long <- length(dif) > 0 || (is.numeric(dif) && "point" == list(...)[[1]])
		if(is.na(long) || !long)
		{ message("dif is empty"); browser()}
		dif
	}
	if(point)
	{
		expectation <- get.diff("expectation")
		zeroCI <- get.diff("zeroCI")
		point <- get.diff("point")
		if(sameNames(x@expectation, y@expectation))
			names(expectation) <- names(x@expectation)
		else
			stop("x, y expectation names conflict with each other")
		new.pointEstimate(expectation = expectation, zeroCI = zeroCI, point = point)
	}
	else
		get.diff() #	matrix
})

is.bootstrap <- function(object)
{
	estimate <- if(is(object, "nonparametricBootstrap"))
		object@original
	else if(is(object, "posteriorEstimate"))
		object@point
	else
		stop("wrong question")
	length(estimate) > 0
}

Lim <- function(x){range(x[is.finite(x)])}
percent <- function(object, digits = 3, return.numeric = TRUE)
{
	x <- 100 * signif(object, digits = digits)
	if(return.numeric)
		x
	else
		paste(x, "%", sep = "")
}
big.width  <- 8.17
big.height <- 8.17
big.pointsize <- 14
y11 <- function(width = big.width, height = big.height, pointsize = big.pointsize, ...){nx11(width = width, height = width, pointsize = pointsize, ...)}
Pdf <- function(width = big.width, height = big.height, pointsize = big.pointsize, ...)
{
	print(list(width = width, height = width, pointsize = pointsize))
	pdf(width = width, height = width, pointsize = pointsize, ...)
}

Stopifnot <- function(x, ...)
{
	if(!x)
	{
		message(...)
		stopifnot(x)
	}
}

is.nothing <- function(object){length(object) == 0}
are.nothing <- function(object)
{
	assert.is(object, "list")
	all(sapply(object, is.nothing))
}

argnames <- function(object)
{
	if(is(object, "function"))
		names(as.list(args(object)))
	else if(is(object, "list"))
		names(object)
	else
		stop("cannot get the argument names")
}

max.finite <- function(x, ...){max(x[is.finite(x)], ...)}
min.finite <- function(x, ...){min(x[is.finite(x)], ...)}

finite.optimize <- function(f, domain, maximum = FALSE, ...)
{
	assert.is(f, "function")
	assert.is(domain, "numeric")
	extreme <- if(maximum) max else min
	domain <- unique(domain)
	codomain <- sapply(domain, f, ...)
	objective <- extreme(codomain)
	boo <- codomain == objective
	extremum <- if(any(boo))
		domain[boo]
	else
		numeric(0)
	lis <- list(extremum = extremum, objective = objective)
	names(lis)[1] <- if(maximum) "maximum" else "minimum"
	lis
}

is.err <- function(object){is(object, "try-error")}

odds <- function(object, ...)
{
	if(is.numeric(object) && all(object < 1.001) && all(object > -0.001))
		object / (1 - object)
	else if(is(object, "empiricalNull")) # congruity.s
		odds(1 - lfdr(object), ...)
	else
	{ message("unsupported object"); browser(); stop("unsupported object")}
}

ve.t.test <- function(...){t.test(..., var.equal = TRUE)}

vectorize <- function(FUN)
{
	if(is(FUN, "function"))
	{
		function(x)
		{
			sapply(x, FUN)
		}
	}
	else if(is(FUN, "list") && are(FUN, "function"))
	{
		stop("what?")
		function(x)
		{
			lapply(FUN, FUN=function(y){y(x)})#lapply(FUN,functi)#XXX|:no visible binding for global variable ÔfunctiÕ
		}

	}
	else
		stop("cannot vectorize")
}

attribute <- function(x)
{
	assert.is(x, "function")
	so <- attributes(x)$source
	if(length(so) == 0)
		character(0)
	else if(is(so, "character"))
		so
	else
		stop("cannot coerce this function to character")
}
sameFunction <- function(x, y)
{
	assert.are(list(x, y), "function")
	ax <- attribute(x)
	ay <- attribute(y)
	boo <- identical(x, y) || (length(ax) > 0 && length(ax) == length(ay) && all(ax == ay))
	!is.na(boo) && boo
}

slots <- function(object, name = ".Data")
{
	assert.is(object, "list")
	vec <- sapply(object, slot, name = name)
	names(vec) <- names(object)
	if(is.list(vec) && all(sapply(vec, length) == 0 & sapply(vec, is, "numeric")))
		numeric(0)
	else
		vec
}

Optimize <- function(f, lower, upper, boundary.tol = numeric(0), maximum = FALSE, call.browser = FALSE, f.bounds = 0, excuse.infinities = FALSE, length.out = 100, lower.reliable = default(FALSE, "lower.reliable"), ...)
{
	stopifnot(lower <= upper)
	if(!is.nothing(length.out))
	{
		X <- seq(lower, upper, length.out = length.out)
		Y <- sapply(X, f)
		if(lower.reliable)
		{
			stopifnot(maximum)
			xi.max <- if(length(Y) > 1)
			{
				dif <- diff(Y)
				i <- 0
				inc <- TRUE
				while(inc && i < length(dif))
				{
					i <- i + 1
					inc <- dif[i] > 0
				}
				dec <- TRUE
				while(dec && i < length(dif))
				{
					i <- i + 1
					dec <- dif[i] < 0
				}
				stopifnot(length(dif) == length(X) - 1)
				stopifnot(i >= 1 && i <= length(X))
				if(dec) length(X) else i
			}
			else
				stop("bad X")
#			xi <- 1:xi.max
			lower2 <- X[1]
			upper2 <- X[xi.max]
		} # end if(lower.reliable)
		else
		{
			extreme <- function(vec){if(maximum) max(vec) else min(vec)}
			get.op <- function(i)
			{
				x <- X[i]
				y <- Y[i]
				stopifnot(length(x) == length(y) && length(x) > 0)
				extremum <- extreme(y)
				boo <- y == extremum
				k <- 1:length(boo)
				stopifnot(length(x) == length(k))
				j <- k[boo]
				if(min(j) > 1)
					j <- c(min(j) - 1, j)
				if(max(j) < length(k))
					j <- c(j, max(j) + 1)
				stopifnot(all(j >= 1 & j <= length(x)))
				op <- x[j]
				stopifnot(extreme(sapply(op, f)) == extremum)
				op
			}
			opX <- get.op(i = 1:length(X))
			lower2 <- min(opX)
			upper2 <- max(opX)
		} # end else
#		eq <- X == opX
#		if(!all(eq))
#			opX <- c(opX, get.op(i = !eq))
#		if(maximum)
#		{
#			rang <- range(opX)
#
#		}
#		else
#			stop("not implemented for maximum = FALSE")
		stopifnot(lower2 <= upper2)
	}
	optimum.name <- if(maximum)
		"maximum"
	else
		"minimum"
	get.lis <- function(optimum)
	{
		li <- list(objective = f(optimum))
		li[optimum.name] <- optimum
		li
	}
	if(lower2 == upper2)
	{
		optimum <- lower
		lis <- get.lis(optimum = optimum)
	}
	else if(lower2 < upper2)
	{
		if(is.nothing(boundary.tol))
			boundary.tol <- abs(lower - upper) / 1e2
		lis <- optimize(f = f, lower = lower2, upper = upper2, maximum = maximum, ...)
		optimum <- lis[optimum.name][[1]]
		too.close <- function(boundary)
		{
			dif <- boundary - optimum
			assert.is(dif, "numeric")
			assert.is(boundary.tol, "numeric")
			at.f.bound <- !is.nothing(f.bounds) && boundary %in% f.bounds
			tc <- try(!at.f.bound && (abs(dif) <= boundary.tol))
			if(is.err(tc))
			{ message("bad tc"); browser()}
			tc
		}
		if(call.browser)
		{ message("call.browser == TRUE"); browser()}
		if(too.close(lower) || too.close(upper))
		{
			objective.lower <- f(lower)
			objective.upper <- f(upper)
			excused <- !all(is.finite(c(objective.lower, objective.upper))) && (is.nothing(excuse.infinities) || excuse.infinities)
			optimum0 <- optimum
			optimum <- if(excused && objective.lower < objective.upper)
			{
				if(maximum) upper else lower
			}
			else if(excused && objective.lower > objective.upper)
			{
				if(maximum) lower else upper
			}
			else if(excused && objective.lower == objective.upper)
			{
				if(lower %in% f.bounds && !(upper %in% f.bounds))
					lower
				else if(upper %in% f.bounds && !(lower %in% f.bounds))
					upper
				else
					NULL
			}
			else
				NULL
			plot.f <- function(length.out = 100, ...)
			{
				X <- seq(lower, upper, length.out = length.out)
				Y <- sapply(X, f)
				plot(x = X, y = Y, xlab = "argument", ylab = "f", ...)
			}
			if(is.null(optimum))
			{ message("\noptimum is too close to lower or upper:"); print(list(optimum0 = optimum0, optimum = optimum, lower = lower, upper = upper)); plot.f(); stop("optimum is too close to lower or upper")}
			lis <- get.lis(optimum = optimum)
			if(is.nothing(excuse.infinities))
			{ message("\nOptimum is too close to lower or upper!"); print(list(lis = lis, lower = lower, upper = upper)); plot.f(); browser()}
		}
	} # end if(lower < upper)
	else
		stop("bad arguments to Optimize")
	oks <- c(length(lis) >= 2, optimum.name %in% names(lis), optimum == lis[optimum.name][[1]], "objective" %in% names(lis))
	ok <- all(oks==T)
	if(!ok)
	{ message("bad return value"); browser()}
	lis
}

is.prob <- function(P, tolerance = 1e-3, ...)
{
#	any(is.finite(P)) &&
	boo <- is(P, "numeric") && all(P >= 0 - tolerance & P <= 1 + tolerance, ...)
	!is.na(boo) && boo
}
are.prob <- function(object, ...)
{
	is.list(object) && all(sapply(object, is.prob, ...))
}

setGeneric("Median", function(x,...) standardGeneric("Median"))#Median <- function(x, ...){} # generic function for methods: marta Nov2013
setMethod("Median", signature(x = "numeric"), function(x, ...){ median(x = x, ...)})#marta Nov2013
Add <- function(x, ...)
{
	sapply(1:length(x), function(i){sum(x[1:i], ...)})
}

plot.polygon <- function(x, y, col = "gray", border = col, in.ylim = numeric(0), add = FALSE, ylim, ...) # plot.polygon(1:3, list(function(x){x ^ 2}, function(x){x ^ 2 - 10}))
{
	assert.is(x, "numeric")
	if(is(y, "list") && are(y, "function") && length(y) == 2)
	{
		Y <- sapply(x, y[[1]])
		x2 <- rev(x)
		Y2 <- sapply(x2, y[[2]])
		X <- c(x, x2)
		Y <- c(Y, Y2)
	}
	else
	{
		Y <- if(is(y, "function"))
			sapply(x, y)
		else if(is(y, "list") && are(y, "function"))
		{
			yy <- numeric(0)
			for(fun in y)
				yy <- c(yy, sapply(x, fun))
			yy
		}
		else
			y
		fac <- length(Y) / length(x)
		stopifnot(floor(fac) == fac && fac >= 1)
		X <- rep(x, fac)
	}
	assert.is(Y, "numeric")
	assert.is(X, "numeric")
	if(length(Y) < 50)
		print(data.frame(X = X, Y = Y))
	if(missing(ylim))
		ylim <- range(c(Y, in.ylim))
	if(!add)
		plot(x = X, y = Y, col = col, type = "l", ylim = ylim, ...)
	polygon(x = X, y = Y, col = col, border = border)
}

zeta <- function(object, ...)
{
	assert.is(object, "numeric")
	stopifnot(is.prob(object))
	vec <- qnorm(object, ...)
	names(vec) <- names(object)
	vec
}

divergence.Bernoulli <- function(P1, P2, base = exp(1), nsample = numeric(0), independent = !is.nothing(nsample), na.rm = logical(0), verbose = FALSE) # see expRelativeEntropy of odd.s, divergence of MDL.s, and relative.entropy of verisimilitude.s
{
	stopifnot(length(base) == 1)
	stopifnot(length(P1) == length(P2))
	assert.are(list(P1, P2), "numeric")
	stopifnot(is.prob(P1) && is.prob(P2))
	lg <- function(vec, ...)
	{
		vec <- as(vec, "numeric")
		logb(vec, ..., base = base)
	}
	if(!independent && is.nothing(nsample) && is.nothing(na.rm))
	{
		if(verbose)
		{
			message("no independence assumption")
			# browser()
		}
		term <- function(p1, p2)
		{
			p1 <- as(p1, "numeric")
			ifelse(p1 == 0, 0, p1 * (lg(p1) - lg(p2)))
		}
		di <- try(term(P1, P2) + term(1 - P1, 1 - P2))
		if(is.err(di))
		{ message("bad divergence!")}#; browser()
		if(any(di < 0))
		{ message("negative divergence!")}#; browser()
		if(verbose) message("returning divergence")
		di
	}
	else if(independent)
	{
		if(verbose) message("independent case for Bernoulli")
		if(is.nothing(na.rm))
			na.rm <- FALSE
		summary.stat <- function(vec, FUN)
		{
			sca <- FUN(vec, na.rm = na.rm)
			if(sca < 0)
			{
				mess <- paste("divergence of ", sca, " truncated at 0", sep = "")
				if(verbose) message(mess)
				warning(mess)
				sca <- 0
			}
			Scalar(sca)
		}
		if(is.nothing(nsample))
		{
			if(verbose) message("making recursive call")
			divergences <- divergence.Bernoulli(P1 = P1, P2 = P2, base = base, nsample = nsample, independent = FALSE, na.rm = logical(0), verbose = verbose)
			if(any(divergences < 0))
			{ message("negative divergence?")}#; browser()
			if(verbose) message("adding up")
			summary.stat(divergences, FUN = sum) # chain rule of relative entropy
		}
		else if(!is.nothing(nsample) && nsample >= 1 && base == exp(1)) # begin Monte Carlo
		{
			warning("Monte Carlo is probably not needed; try independent = TRUE instead of specifying nsample")
			get.lik <- function(p, success)
			{
				stopifnot(is.logical(success) && length(success) == length(p))
				lik <- sum(lg(ifelse(success, p, 1 - p)))
	#			if(!is.finite(lik))
	#			{ message("bad lik:"); print(stats(lik))}
				lik
			}
			get.success <- function(p)
			{
				ran <- runif(n = length(p))
				stopifnot(length(ran) == length(p) && is.prob(ran) && is.prob(p))
				ran <= p
			}
			lik.ratio <- sapply(1:nsample, function(i)
			{
				success <- get.success(P1)
				get.lik(p = P1, success = success) - get.lik(p = P2, success = success)
			})
			if(any(is.na(lik.ratio)))
			{ message("bad lik.ratio"); print(stats(lik.ratio))}#; browser()
			summary.stat(lik.ratio, FUN = mean)
		} # end Monte Carlo
		else
			stop("bad nsample, independent case")
	} # end if(independent)
	else
		stop("bad nsample")
}

classes <- function(...)
{
	lis <- list(...)
	cla <- sapply(lis, class)
	names(cla) <- names(lis)
	cla
}

lowerBayesFactor0 <- function(p){ifelse(p < exp(-1), - exp(1) * p * log(p), 1)}

#zplot <- function(x, y, ...){printGeneric("zplot"); browser()} # new 25 April 2008
#removeMethods("zplot")
#setMethod("zplot", signature(x = "Fdr", y = "missing"), function(x, y, ...)#XXX| class"Fdr" not defined
#{
#	get.lab <- function(prefix){paste(prefix, "(z space)")}
#	plot(zvalue(x, type = "pvalue"), zvalue(x, type = "fdr"), xlab = get.lab(prefix = "p-value"), ylab = get.lab(prefix = "local false discovery rate"), ...)
#	abcol <- "gray"
#	abline(v = 0, col = abcol)
#	abline(h = 0, col = abcol)
#})
#setMethod("zplot", signature(x = "Prob0", y = "missing"), function(x, y, ...)#XXX| class"Prob0" not defined
#{
#	zplot(x = as(x, "Fdr"), ...)
#})
Character <- function(object, ...){message("generic Character"); browser()}
#removeMethods("Character")
#setMethod("Character", signature(object = "xprnScore"), function(object, size, ...)
#{
#  if(missing(size)) stop("size missing")
#  arglis <- list(...)
#  fun <- function(xS, ...)
#  {
#		stop <- length(xS)
#		start <- stop - size + 1
#		stopifnot(start <= stop)
#  	nam <- names(do.call("sort", c(arglis, list(xS, ...)))[start:stop])
#  	if(!(is.character(nam) && length(nam) == size))
#  	{ message("bad nam in Union"); browser()}
#  	nam
#  }
#	fun(object)
#})

Union <- function(x, y, ...){message("generic Union"); browser()}
#removeMethods("Union")
#setMethod("Union", signature(x = "xprnScore", y = "xprnScore"), function(x, y, size, ...)
#{
#  if(missing(size)) stop("size missing")
#  stopifnot(length(x) == length(y))
#  fun <- function(xS){Character(object = xS, size = size, ...)}
#  union(fun(x), fun(y))
#})

Intersect <- function(x, y, ...){message("generic Intersect"); browser()}
#removeMethods("Intersect")
#setMethod("Intersect", signature(x = "xprnScore", y = "xprnScore"), function(x, y, size, ...)
#{
#  if(missing(size)) stop("size missing")
#  stopifnot(length(x) == length(y))
#  fun <- function(xS){Character(object = xS, size = size, ...)}
#  intersect(fun(x), fun(y))
#})

# Sellke et al. (RefWorks:1218) according to "blended inference.nb"
#From our Rprofile.site April2014
refresh <- function(file = getwd()){setwd(dir = file)} # needed if disk is disconnected from Wintel PC
loadi <- function(file.base = "", envir = .GlobalEnv, ...)
{
	refresh()
	file <- paste(file.base, "RData", sep = ".")
	message("loading ", file, " on ", date())
	load(file = file, envir = envir, ...)
}
cd <- function(dir = stop("file name not specified"), ...)
{
	setwd(dir = dir)
	try(loadi(...))
}
# EOF

