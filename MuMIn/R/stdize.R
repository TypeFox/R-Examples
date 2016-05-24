isFALSE <- function(x) identical(FALSE, x)

stdize <-
function(x, ...) UseMethod("stdize")

rootmeansq <- function(v) {
	v <- as.numeric(v[!is.na(v)])
	sqrt(sum(v^2) / max(1, length(v) - 1L))
}

stdize.default <-
stdize.numeric <-
function(x, center = TRUE, scale = TRUE, ...) {
	if(is.function(scale)) {
		scaleFunc <- scale
		scale <- TRUE
	} else scaleFunc <- function(x) sd(x, na.rm = TRUE)
	
	#if(!missing(...)) warning("additional arguments ignored")
	
	for(i in c("scale", "center")) {
		if(length(v <- get(i, inherits = FALSE)) > 1L)
			cry(, "only first element of '%s' is used", i)
		assign(i, v[1L])
	}

	if(is.logical(scale)) scale <- if(scale) scaleFunc(x) else 1
	if(is.logical(center)) center <- if(center) mean(x, na.rm = TRUE) else 0
	
	if(scale == 0) scale <- 1
	if(!is.na(scale)) {
		x <- (x - center) / scale
		attr(x, "scaled:center") <-  center
		attr(x, "scaled:scale") <- scale
	}
	x
}

stdize.matrix <-
function(x, center = TRUE, scale = TRUE, ...) {
	if(!is.numeric(x)) return(x)
	#if(!missing(...)) warning("additional arguments ignored")
	if(is.function(scale)) {
		scaleFunc <- scale
		scale <- TRUE
	} else scaleFunc <- function(x) sd(x, na.rm = TRUE)
	
	for(i in c("scale", "center"))
		if((nv <- length(v <- get(i, inherits = FALSE))) > 1L && nv != ncol(x))
			cry(, "length of '%s' (%d) not equal to number of columns in 'x' (%d)", i, nv, ncol(x))
	
	if(is.logical(scale)) scale <- if(scale) apply(x, 2L, scaleFunc) else 1
	if(is.logical(center)) center <- if(center) colMeans(x, na.rm = TRUE) else 0
	nc <- ncol(x)
	center <- rep(center, length.out = nc)
	scale <- rep(scale, length.out = nc)
	ok <- which(scale != 0)
	scale[-ok] <- center[-ok]  <- NA
	for(i in ok) x[, i] <- (x[, i] - center[i]) / scale[i]
	attr(x, "scaled:center") <- center
	attr(x, "scaled:scale") <- scale
	x
}


stdize.factor <-
function(x, binary = c("center", "scale", "binary", "half", "omit"),
center = TRUE, scale = FALSE, 
...) {
	#if(!missing(...)) warning("additional arguments ignored")
	if(nlevels(x) == 2L) {
		stdize.logical(as.numeric(x) - 1, binary, center, scale)
	} else x
}

stdize.logical <-
function(x, binary = c("center", "scale", "binary", "half", "omit"),
center = TRUE, scale = FALSE, 
...) {
	#if(!missing(...)) warning("additional arguments ignored")
	binary <- if(is.null(binary) || is.na(binary)) "" else match.arg(binary)
	switch(binary,
			   center = stdize.numeric(x, center = TRUE, scale = 1),
			   scale = stdize.numeric(x, center = TRUE, scale = TRUE),
			   half = stdize.numeric(x, center = 0.5, scale = 1),
			   binary = stdize.numeric(x, center = 0, scale = 1),
			   omit = x,
			   stdize.numeric(as.numeric(x), center = center, scale = scale)
			   )
}

stdize.data.frame <-
function(x, binary = c("center", "scale", "binary", "half", "omit"),
	center = TRUE, scale = TRUE,
	omit.cols = NULL,
	source = NULL, prefix = TRUE,
	append = FALSE, ...) {

	if(is.function(scale)) {
		scaleFunc <- scale
		scale <- TRUE
	} else scaleFunc <- function(x) sd(x, na.rm = TRUE)
	
	for(i in c("scale", "center"))
		if((nv <- length(v <- get(i, inherits = FALSE))) > 1L && nv != ncol(x))
			cry(, "length of '%s' (%d) not equal to number of columns in 'x' (%d)", i, nv, ncol(x))
	
	if(!is.null(source)) {
		if(!missing(center) || !missing(scale) || !missing(binary))
			warning("arguments 'center', 'scale' and 'binary' ignored if 'source' is given")
		j <- match(colnames(x), attr(source, "orig.names"), nomatch = 0L)

		if(all(j == 0L)) stop("no columns in 'source' match columns in 'x'")
		
		x <- x[, j != 0L, drop = FALSE]
		center <- attr(source, "scaled:center")[j]
		scale <- attr(source, "scaled:scale")[j]
		binary <- ""
		if(is.null(center) || is.null(scale)) stop("invalid 'source' object")
	} else
		binary <- if(is.null(binary) || is.na(binary)) "" else match.arg(binary)
	
	dataClasses <- vapply(x, function(x) {
		if (is.logical(x)) return("logical")
		if (is.factor(x)) if(nlevels(x) == 2L) return("factor2") else 
			return("other")
		if (is.matrix(x) && is.numeric(x)) return("nmatrix")
		if (is.numeric(x)) return("numeric")
		 return("other")
	}, "")
	
	origx <- x
	
	if(is.character(omit.cols)) dataClasses[colnames(x) %in% omit.cols] <- "omit"
		else if(is.numeric(omit.cols)) dataClasses[omit.cols] <- "omit"
		
	numData <- dataClasses == "numeric"

	if(binary == "omit")  {
		binaryData <- FALSE
	} else {
		binaryData <- dataClasses == "factor2" | dataClasses == "logical"
		for (i in which(binaryData))
			x[, i] <- as.numeric(x[, i]) - if(dataClasses[i] == "factor2") 1 else 0
	}
		
	nc <- ncol(x)

	f <- function(x, bin) {
		if(is.numeric(x)) {
			calc <- isTRUE(bin) & binaryData
			x <- rep(x, length.out = nc)
			do <- !is.na(x) & (numData | calc | (binaryData & (is.na(bin) | !isFALSE(bin))))
			x[!calc & do & binaryData & !is.na(bin)] <-  bin
			return(list(num = x, calc = calc, do = do))
		} else {
			calc <- (numData & x) | (binaryData & ((is.na(bin) & x) | isTRUE(bin)))
			do <- calc | (binaryData & (!is.na(bin) & !isFALSE(bin)))
			num <- numeric(nc)
			num[!calc & do & binaryData & !is.na(bin)] <-  bin
			return(list(num = num, calc = calc, do = do))
		}
	}
	
	
	ctr <- f(center, switch(binary, center = TRUE, scale = TRUE, half = .5, 
		binary = 0, omit = FALSE, NA))
	scl <- f(scale, switch(binary, center = FALSE, scale = TRUE, half = 1, 
		binary = FALSE, omit = FALSE, NA))
	
	center <- ctr$num
	scale <- scl$num
	center[ctr$calc] <- colMeans(x[, ctr$calc, drop = FALSE], na.rm = TRUE)
	scale[scl$calc] <- apply(x[, scl$calc, drop = FALSE], 2L, scaleFunc)
	
	#dp(center)
	#dp(scale)
	
	scl$do[scl$do & scale == 0] <- FALSE
	jTransformed <- ctr$do | scl$do
	center[jTransformed & !ctr$do] <- 0
	scale[jTransformed & !scl$do] <- 1
	for (i in which(jTransformed)) x[, i] <- (x[, i] - center[i]) / scale[i]
	
	doprefix <-  FALSE
	if(is.character(prefix) ||
	   (doprefix <- (is.logical(prefix) && isTRUE(prefix)))) {
		prefix <- if(doprefix) c("z.", "c.") else rep(prefix, length.out = 2L)
		
		Dcenter<- as.data.frame(ctr)
		Dscale <- as.data.frame(scl)
		#dp(Dcenter)
		#dp(Dscale)
		#dp(jTransformed)
		#dp(paste0(prefix[jTransformed + (ctr$do & !scl$do)], 
		#		colnames(x)[jTransformed]))
		#dp(jTransformed + (ctr$do & !scl$do))
		
		colnames(x)[jTransformed] <-
			paste0(prefix[jTransformed + (ctr$do & !scl$do)], 
				colnames(x)[jTransformed])
	}
	if(append) x <- cbind(origx, x[,!(colnames(x) %in% names(origx)), drop = FALSE], deparse.level = 0L)
	attr(x, "scaled:center") <- ifelse(jTransformed, center, NA)
	attr(x, "scaled:scale") <- ifelse(jTransformed, scale, NA)		
	attr(x, "orig.names") <- colnames(origx)
	x
}

stdize.formula <-
function(x, data = NULL, response = FALSE,
binary = c("center", "scale", "binary", "half", "omit"),
center = TRUE, scale = TRUE,
omit.cols = NULL,
prefix = TRUE,
append = FALSE, ...) {
	mf <- model.frame(x, data = data, drop.unused.levels = TRUE, ...)
	if(!is.null(omit.cols)) 
		omit.cols <- if(is.character(omit.cols)) 
			which(colnames(mf) %in% omit.cols) else
			stop("'omit.cols' must be a character vector")
	if(!response) omit.cols <- unique(c(omit.cols, 1L))
	attr(mf, "terms") <- NULL
	stdize.data.frame(mf, center = center, scale = scale, 
		omit.cols = omit.cols, binary = binary, prefix = prefix,
		append = append)
}

stdizeFit <-
function(object, data, which = c("formula", "subset", "offset", "weights"),
		 evaluate = TRUE, quote = NA) {
	thiscall <- match.call()
	if(is.na(quote)) quote <- is.call(thiscall$object) && !is.primitive(match.fun(thiscall$object[[1L]]))
	cl <- if(quote) thiscall$object
		else if(is.expression(object)) object[[1L]]
		else if(is.call(object)) object
		else get_call(object)
	
	cl <- match.call(Fun <- match.fun(cl[[1L]]), cl)
	
	if(!("data" %in% (formalnames <- names(formals(Fun)))))
		warning(gettextf("%s does not have a formal argument 'data'", as.character(cl[[1L]])))
	
	which <- which[which %in% formalnames]
	
	i <- names(data) != attr(data, 'orig.names')
	env <- structure(lapply(names(data)[i], as.name), names = attr(data, 'orig.names')[i])
	if(isTRUE(which)) cl <- subst(cl, env)
		else
		for(i in which) if(!is.null(cl[[i]])) cl[[i]] <- subst(cl[[i]], env)
	cl[['data']] <- thiscall[['data']] #substitute(data)
	if(evaluate) eval.parent(cl) else cl
}

