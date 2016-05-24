
nx11<-function(height = 2.5, width = 1.5,pointsize=8,mai=c(0.7,0.7,0.6,0.3),...){par(pin=c(width, height),ps=pointsize,mai=mai,...)}
#par(new=T);fin,pin
Mfrow<-function(half = FALSE, quarter = FALSE, mfrow, height = 2.5,width= 1.5,mai=c(0.7,0.7,0.6,0.3),...)
{
	assert.is(half, "logical")
	
	if(missing(mfrow))
	{
		mfrow <- if(half)
		{
			nx11(height = height, width = width ,mai=mai,...)#3.5
			c(1, 2)
		}
		else if(quarter)
		{
			mai[1]<-mai[1]+0.1
			mai[2]<-mai[2]+0.1
			width<-height
			nx11(height = height, width = width ,mai=mai,...)
			c(1, 1)
		}
		else {
			nx11(height = height, width = width ,mai=mai,...)
			c(2, 2)}
			
	}
	else
		stopifnot(!half && !quarter)
	par(mfrow = mfrow, ...)
}

blank.plot <- function(legend, ...)
{
	plot(x = 0, type = "n", col.axis = "white", xlab = "", ylab = "", tcl = 0, xaxt = "n", yaxt = "n", bty = "n", ...)
	if(!missing(legend))
		legend(x = "center", legend = legend, col = "white", bty = "n")
}
#--------------------

setGeneric("Rank", function(object,...) standardGeneric("Rank"))#removeMethods("Rank")
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
		new_Numeric(ra)
	}
})




#
s3_recursive <- function(object, name, first.call)
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
		s3_recursive(x, first.call = FALSE)
	else
		x
}

setGeneric("s3", function(object, name) standardGeneric("s3")) #removeMethods("s3")
setMethod("s3", signature(object = "ANY", name = "missing"), function(object, name)
{
	s3_recursive(object = object, name = name)
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
#
silence <- function(z, nsilence, silent = NULL)
{
	assert.is(z, "numeric")
	az <- abs(z)
	get.rank <- function(x){rank(x, ties.method = "random")}
	raz <- get.rank(-az)
	maz <- max(az, na.rm = TRUE)
	if(az[raz == 1] != maz)
	{ message("bad raz"); browser()}
	stopifnot(nsilence >= 1 && floor(nsilence) == nsilence && nsilence < length(z))
	boo <- raz <= nsilence
	stopifnot(sum(boo) == nsilence)
	stopifnot(maz %in% az[boo])
	z[boo] <- if(is.null(silent))
	{
		rz <- get.rank(z)
#		oz <- order(z)
		expected.order.stats <- rnorm(ppoints(length(z)))[rz]
		stopifnot(all(length(expected.order.stats) %in% c(length(rz), length(boo))))
		new_z <- expected.order.stats[rz[boo]]
		stopifnot(length(new_z) == sum(boo))
		new_z
	}
	else if(length(silent) == 0)
		rnorm(n = sum(boo), mean = 0, sd = 1)
	else if((is.numeric(silent) || is.na(silent)) && length(silent) == 1)
		silent
	else
		stop("bad silent")
	z
}
#---------
geomean <- function(x, ...){exp(mean(x = logb(x), ...))}
harmmean <- function(x, ...){rec <- function(y){1/y}; rec(mean(x = rec(x), ...))}

setGeneric("Mean", function(x,...) standardGeneric("Mean"))#removeMethods("Mean") 
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


setGeneric("Median", function(x,...) standardGeneric("Median"))
setMethod("Median", signature(x = "numeric"), function(x, ...){ median(x = x, ...)})
setMethod("Median", signature(x = "matrix"), function(x, ...)
{
	vec <- sapply(1:nrow(x), function(i){Median(x[i, ], ...)})
	stopifnot(nrow(x) == length(vec) && is.numeric(vec))
	names(vec) <- rownames(x)
	vec
})

setGeneric("Sd", function(object,...) standardGeneric("Sd"))
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
	nScalar(Sca)
})




setGeneric("geoMean",function(x, ...) standardGeneric("geoMean"))
#try(removeMethods("geoMean"), silent = TRUE)
setMethod("geoMean", signature(x = "numeric"), function(x, ...) # na.rm corrected 18 March 2010
{
	Mean(x = x, geometric = TRUE, ...)
})

se.mean <- function(...)
{
	Sd(..., se.of.mean = TRUE)
}
Var <- function(...){Sd(...) ^ 2}
#
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

