"adjamat.mathgraph" <-
function(x, general = FALSE, ...)
{
	if(general != FALSE) .NotYetUsed("general != FALSE")
	x <- unclass(x)
	xdir <- attr(x, "directed")
	ischar <- is.character(x)
	if(ischar) {
		nnam <- unique(x)
		nnode <- length(nnam)
		has.names <- TRUE
	}
	else {
		if(!is.numeric(x))
			stop("nodes must be character or numeric")
		nnode <- max(x)
		nnam <- paste("node", 1:nnode)
		has.names <- FALSE
	}
	ans <- array(0, c(nnode, nnode), list(nnam, nnam))
	if(ischar) {
		dx <- dim(x)
		x <- match(x, nnam)
		dim(x) <- dx
	}
	ans[x] <- 1
	if(any(!xdir))
		ans[x[!xdir, 2:1]] <- 1
	attr(ans, "call") <- match.call()
	attr(ans, "has.names") <- has.names
	class(ans) <- "adjamat"
	ans
}

"adjamat" <- function(x, ...)
UseMethod("adjamat")

"alldirected.default" <-
function(x, ...)
{
	stop("do not know how to handle")
}

"alldirected.mathgraph" <-
function(x, ...)
{
	dir <- attr(x, "directed")
	if(all(dir))
		return(x)
	x <- unclass(x)
	ans <- rbind(x, x[!dir, 2:1])
	attr(ans, "directed") <- rep(TRUE, nrow(ans))
	class(ans) <- "mathgraph"
	ans
}

"alldirected" <- function(x, ...)
UseMethod("alldirected")

"build.mathgraph" <-
function(formula, data)
{
	cterm <- paste(collapse = "", deparse(formula[[2]]))
	eyes <- .find.I.of(cterm)
	if(length(eyes)) {
		# calls to I() need to be evaluated
		nI <- nrow(eyes)
		Inames <- paste("Build.mathgraphI", sys.nframe(), 1:nI, sep = 
			".")
		Iexpr <- substring(cterm, eyes[, 1] + 2, eyes[, 2] - 1)
		for(i in 1:nI) {
			this.val <- eval(parse(text = Iexpr[i]), data)
			assign(Inames[i], this.val, pos = 1)
		}
		eye.ext <- matrix(c(0, t(eyes), nchar(cterm) + 2), nrow = 2)
		allI.strings <- character(2 * nI + 1)
		allI.strings[2 * (1:nI)] <- Inames
		allI.strings[seq(1, 2 * nI + 1, by = 2)] <- substring(cterm, 
			eye.ext[1,  ] + 1, eye.ext[2,  ] - 1)
		cterm <- paste(allI.strings, collapse = "")
	}
	allterms <- strsplit(cterm, split = "\\+")#changed by Claus Dethlefsen
#	allterms <- strsplit(cterm, split = "+") 
	ans <- NULL
	adir <- logical(0)
	for(i in seq(along = allterms)) {
		raw <- allterms[[i]]
		this.parse <- parse(file = "", 
			text = paste(raw, sep = "", collapse = "") )
		if(length(this.parse[[1]]) == 1) {
			this.op <- " "
		}
		else {
			this.op <- as.character(this.parse[[1]][[1]])
		}
		switch(this.op,
			"/" = {
				e1 <- eval(this.parse[[1]][[2]], data)
				e2 <- eval(this.parse[[1]][[3]], data)
				this.ev <- cbind(e1, e2)
				ans <- rbind(ans, this.ev)
				adir <- c(adir, rep(NA, dim(this.ev)[1]))
			}
			,
			"*" = {
				e1 <- eval(this.parse[[1]][[2]], data)
				e2 <- eval(this.parse[[1]][[3]], data)
				e1 <- unique(e1)
				e2 <- unique(e2)
				le1 <- length(e1)
				le2 <- length(e2)
				e1 <- e1[rep(1:le1, le2)]
				e2 <- e2[rep(1:le2, rep(le1, le2))]
				this.ev <- cbind(e1, e2)
				ans <- rbind(ans, this.ev)
				adir <- c(adir, rep(NA, dim(this.ev)[1]))
			}
			,
			{
				this.ev <- eval(this.parse, data)
				if(inherits(this.ev, "mathgraph")) {
				  ans <- rbind(ans, this.ev)
				  adir <- c(adir, attr(this.ev, "directed"))
				}
				else stop(paste(
				    "do not know how to handle term:", raw))
			}
		)
	}
	attr(ans, "directed") <- adir
	ans
}

"c.mathgraph" <-
function(...)
{
	dots <- list(...)
	the.class <- commontail(lapply(dots, class))
	if(!match("mathgraph", the.class, nomatch = 0))
		stop("not all mathgraph objects")
	amat <- do.call("rbind", dots)
	adir <- unlist(lapply(dots, function(x)
	attr(x, "directed")))
	if(dim(amat)[1] != length(adir))
		stop("garbled object")
	attr(amat, "directed") <- adir
	class(amat) <- the.class
	amat
}

"commontail" <-
function(x)
{
	the.lens <- unlist(lapply(x, length))
	min.len <- min(the.lens)
	if(min.len == 0)
		return(NULL)
	ans <- NULL
	x <- lapply(x, rev)
	for(i in 1:min.len) {
		this.ans <- unique(unlist(lapply(x, function(y, i)
		y[i], i = i)))
		if(length(this.ans) == 1) {
			ans <- c(this.ans, ans)
		}
		else {
			break
		}
	}
	ans
}

".find.I.of" <-
function(string, nesting.ok = FALSE)
{
	if(nesting.ok != FALSE) .NotYetUsed("nesting.ok != FALSE")
#	xx <- regexpr(pattern = "I\([^)]*\)", text = string)
	xx <- regexpr(pattern = "I([^)]*)", text = string)
	if(xx[1] == -1)
		return(NULL)
	# Test the ends of the pattern found
	yy <- substring(string, xx[1], last = xx[1] + attr(xx, "match.length"))
	if(substring(yy, first = 2, last = 2) != "(") {
		return(NULL)
	}
	ans <- array(data = 0, dim = c(1, 2))
	# For `nesting.ok', loop over length(Istart.num) rows
	if(substring(yy, first = nchar(yy), last = nchar(yy)) == ")") {
		ans[1, ] <- c(xx[1], xx[1] + attr(xx, "match.length"))
	} else {
		warning("no closing parenthesis for I")
		ans <- NULL
	}
	ans
}

"getpath.adjamat" <-
function(x, start, end, ...)
{
	if(start == end) {
# distinguish this from no path possible
		return(mathgraph())
	}
	has.names <- attr(x, "has.names")
	if(has.names)
		node.names <- dimnames(x)[[1]]
	else node.names <- NULL
	if(is.character(c(end, start))) {
		start <- match(start, node.names, nomatch = NA)
		end <- match(end, node.names, nomatch = NA)
		bad.in <- c("start", "end")[is.na(c(start, end))]
		if(length(bad.in))
			stop(paste(paste(bad.in, collapse = " and "), 
				"not right"))
	}
	tset <- start
	prev <- 0
	unchecked <- TRUE
	nseq <- 1:dim(x)[1]
	repeat {
		this.index <- (1:length(unchecked))[unchecked][1]
		newn <- nseq[x[tset[this.index],  ] > 0]
		newn <- setdiff(newn, tset)
		unchecked[this.index] <- FALSE
		if(n <- length(newn)) {
			if(match(end, newn, nomatch = 0)) {
				tset <- c(tset[!unchecked], end)
				prev <- c(prev[!unchecked], tset[this.index])
				pseq <- 1:length(tset)
				path <- this.index <- length(tset)
				this.node <- end
				while(prev[this.index] != start) {
				  this.index <- pseq[tset == prev[this.index]]
				  path <- c(this.index, path)
				}
				if(has.names) {
				  return(mathgraph( ~ node.names[prev[path]]/
				    node.names[tset[path]], directed = TRUE))
				}
				else {
				  return(mathgraph( ~ prev[path]/tset[path], 
				    directed = TRUE))
				}
			}
			tset <- c(tset, newn)
			prev <- c(prev, rep(tset[this.index], n))
			unchecked <- c(unchecked, rep(TRUE, n))
		}
		if(!any(unchecked))
			return(NULL)
	}
}

"getpath.default" <-
function(x, start, end, ...)
{
	getpath(incidmat(x), start, end, ...)
}

"getpath.incidmat" <-
function(x, start, end, ...)
{
	if(start == end)
		return(mathgraph())
	x <- unclass(x)
	dircheck <- rep(1, dim(x)[1]) %*% x
	if(any(dircheck))
		stop("need matrix for directed graph")
	has.names <- attr(x, "has.names")
	if(has.names["edges"])
		enames <- dimnames(x)[[2]]
	else enames <- NULL
	if(has.names["nodes"])
		node.names <- dimnames(x)[[1]]
	else node.names <- NULL
	if(is.character(c(end, start))) {
		start <- match(start, node.names, nomatch = NA)
		end <- match(end, node.names, nomatch = NA)
		bad.in <- c("start", "end")[is.na(c(start, end))]
		if(length(bad.in))
			stop(paste(paste(bad.in, collapse = " and "), 
				"not right"))
	}
	tset <- start
	prev <- 0
	edges <- 0
	unchecked <- TRUE
	nseq <- 1:dim(x)[1]
	eseq <- 1:dim(x)[2]
	repeat {
		this.index <- (1:length(unchecked))[unchecked][1]
		newe <- eseq[x[tset[this.index],  ] > 0.5]
		if(length(newe)) {
			newn <- row(x[, newe, drop = FALSE])[as.vector(x[, newe] < 
				-0.5
				)]
			newt <- !duplicated(newn)
			newn <- newn[newt]
			newe <- newe[newt]
			newt <- match(newn, tset, nomatch = 0) == 0
			newn <- newn[newt]
			newe <- newe[newt]
		}
		else newn <- NULL
		unchecked[this.index] <- FALSE
		if(n <- length(newn)) {
			if(endind <- match(end, newn, nomatch = 0)) {
# have a path
				tset <- c(tset[!unchecked], end)
				prev <- c(prev[!unchecked], tset[this.index])
				edges <- c(edges[!unchecked], newe[endind])
				pseq <- 1:length(tset)
				path <- this.index <- length(tset)
				this.node <- end
				while(prev[this.index] != start) {
				  this.index <- pseq[tset == prev[this.index]]
				  path <- c(this.index, path)
				}
				if(has.names["nodes"]) {
				  ans <- mathgraph( ~ node.names[prev[path]]/
				    node.names[tset[path]], directed = TRUE)
				}
				else {
				  ans <- mathgraph( ~ prev[path]/tset[path], 
				    directed = TRUE)
				}
				if(has.names["edges"]) {
				  names(ans) <- enames[edges[path]]
				}
				return(ans)
			}
			tset <- c(tset, newn)
			edges <- c(edges, newe)
			prev <- c(prev, rep(tset[this.index], n))
			unchecked <- c(unchecked, rep(TRUE, n))
		}
		if(!any(unchecked))
			return(NULL)
	}
}

"getpath.mathgraph" <-
function(x, start, end, ...)
{
	if(start == end)
		return(mathgraph())
	getpath(incidmat(x), start, end, ...)
}

"getpath" <- function(x, start, end, ...)
UseMethod("getpath")

"incidmat.mathgraph" <-
function(x, expand = TRUE, general = FALSE, ...)
{
	x <- unclass(x)
	xdir <- attr(x, "directed")
	nedge <- dim(x)[1]
	eseq <- 1:nedge
	cnam <- dimnames(x)[[1]]
	has.names <- c(nodes = is.character(x), edges = TRUE)
	if(!length(cnam)) {
		has.names["edges"] <- FALSE
		cnam <- paste(ifelse(xdir, "arc", "edge"), eseq)
	}
	ischar <- is.character(x)
	if(ischar) {
		nnam <- unique(x)
		nnode <- length(nnam)
		dx <- dim(x)
		x <- match(x, nnam)
		dim(x) <- dx
	}
	else {
		if(!is.numeric(x))
			stop("nodes must be character or numeric")
		nnode <- max(x)
		nnam <- paste("node", 1:nnode)
	}
	loops <- x[, 1] == x[, 2]
	if(expand && !all(xdir)) {
		cnam <- rep(cnam, 2 - xdir)
		ans <- array(0, c(nnode, length(cnam)), list(nnam, cnam))
		reseq <- match(eseq, rep(eseq, 2 - xdir))
		ans[cbind(x[, 2], reseq)] <- -1
		ans[cbind(x[, 1], reseq)] <- 1
		ans[cbind(x[!xdir, 1], reseq[!xdir] + 1)] <- -1
		ans[cbind(x[!xdir, 2], reseq[!xdir] + 1)] <- 1
		if(any(loops)) {
			rloops <- rep(loops, 2 - xdir)
			loop.node <- rep(x[loops, 1], 2 - xdir[loops])
			if(general) {
				ans[cbind(loop.node, seq(along = rloops)[rloops
				  ])] <- 1
			}
			else {
				ans[, rloops] <- 0
			}
		}
	}
	else {
		ans <- array(0, c(nnode, nedge), list(nnam, cnam))
		ans[cbind(x[, 1], eseq)] <- 1
		ans[cbind(x[, 2], eseq)] <- ifelse(xdir, -1, 1)
		if(any(loops)) {
			if(general) {
				ans[cbind(x[loops, 1], eseq[loops])] <- 2 - 
				  xdir[loops]
			}
			else {
				ans[, loops] <- 0
			}
		}
	}
	attr(ans, "has.names") <- has.names
	attr(ans, "call") <- match.call()
	class(ans) <- "incidmat"
	ans
}

"incidmat" <- function(x, ...)
UseMethod("incidmat")

"is.adjamat" <- function(x) inherits(x, "adjamat")

"is.incidmat" <- function(x) inherits(x, "incidmat")

"is.mathgraph" <- function(x) inherits(x, "mathgraph")

"justify" <-
function(x, type = "r")
{
	type <- match.arg(type, c("right", "center", "left"))
	x <- as.character(x)
	ncx <- nchar(x)
	blanks <- paste(rep(" ", max(ncx)), collapse = "")
	blanks <- substring(blanks, 1, max(ncx) - ncx)
	switch(type,
		right = paste(blanks, x, sep = ""),
		left = paste(x, blanks, sep = ""),
		center = {
			blank.half <- nchar(blanks) %/% 2
			paste(substring(blanks, 1, blank.half), x, substring(
				blanks, blank.half + 1), sep = "")
		}
	)
}

"length.mathgraph" <-
function(x)
{
	x <- unclass(x)
	if(length(x))
		dim(x)[1]
	else 0
}

"mathgraph" <-
function(formula, directed = FALSE, data = sys.parent())
{
	if(missing(formula)) {
		ans <- matrix(NA, nrow=0, ncol=2)
	}
	else {
		ans <- build.mathgraph(formula, data = data)
		adir <- attr(ans, "directed")
		adir[is.na(adir)] <- directed
		attr(ans, "directed") <- adir
	}
	class(ans) <- "mathgraph"
	ans
}

"names.mathgraph" <-
function(x)
{
	dimnames(unclass(x))[[1]]
}

"names<-.mathgraph" <-
function(x, value)
{
	cl <- class(x)
	x <- unclass(x)
	dimnames(x)[[1]] <- value
	class(x) <- cl
	x
}

"plot.mathgraph" <-
function(x, ...)
{
	x <- unclass(x)
	xdir <- attr(x, "directed")
	if(ischar <- is.character(x)) {
		node.names <- unique(x)
		maxx <- length(node.names)
		dx <- dim(x)
		x <- match(x, node.names)
		dim(x) <- dx
	}
	else {
		maxx <- max(x)
		node.names <- as.character(1:maxx)
	}
	px <- cos((2 * 0:(maxx - 1) * pi)/maxx)
	py <- sin((2 * 0:(maxx - 1) * pi)/maxx)
	plot(px, py, axes = FALSE, xlab = "", ylab = "", xlim = c(-1.04, 1.04), 
		ylim = c(-1.04, 1.04), ...)
	box()
	px <- 0.97999999999999998 * px
	py <- 0.97999999999999998 * py
	if(!all(xdir))
		segments(px[x[!xdir, 1]], py[x[!xdir, 1]], px[x[!xdir, 2]], py[
			x[!xdir, 2]])
	if(any(xdir))
		arrows(px[x[xdir, 1]], py[x[xdir, 1]], px[x[xdir, 2]], py[x[
			xdir, 2]])
	text(px * 1.0700000000000001, py * 1.0700000000000001, node.names)
	invisible()
}

"print.mathgraph" <-
function(x, prefix.node = if(is.character(xu)) "" else "node", ...)
{
	if(length(unclass(x))) {
		xu <- unclass(x)
		if(length(the.nams <- names(x))) {
			the.nams <- paste(justify(the.nams, "r"), " ", sep = ""
				)
		}
		else {
			the.nams <- paste("[", format(1:length(x)), "] ", sep
				 = "")
		}
		out <- paste(the.nams, prefix.node, xu[, 1], ifelse(attr(x, 
			"directed"), "->", "--"), prefix.node, xu[, 2])
		cat(out, sep = "\n")
		cat("\nclass:", class(x), "\n")
	}
	else {
		cat("mathgraph()\n")
	}
	invisible(x)
}

"sortmathgraph" <-
function(x, nodes = TRUE, edges = TRUE)
{
	dir <- attr(x, "directed")
	cl <- class(x)
	x <- unclass(x)
	if(nodes && !all(dir)) {
		x[!dir,  ] <- stable.apply(x[!dir,  , drop = FALSE], 1, sort)
	}
	if(edges) {
		ord <- order(x[, 1], x[, 2])
		x <- x[ord,  ]
		attr(x, "directed") <- dir[ord]
	}
	class(x) <- cl
	x
}
"stable.apply"<-
function(X, MARGIN, FUN, ...)
{
	ldx <- length(dim(X))
	if(length(MARGIN) != ldx - 1) {
		warning("stability not performed")
		return(apply(X, MARGIN, FUN, ...))
	}
	ans <- apply(X, MARGIN, FUN, ...)
	if(length(dim(ans)) != ldx)
		ans
	else aperm(ans, order(c((1:ldx)[ - MARGIN], MARGIN)))
}

"[.mathgraph" <-
function(x, i)
{
	cl <- class(x)
	x <- unclass(x)
	xdir <- attr(x, "directed")
	if(is.character(i))
		names(xdir) <- dimnames(x)[[1]]
	if(is.logical(i)){
		if(ncol(i)==2)
			i <- i[ ,1] & i[ ,2]
	}
	x <- x[i,  , drop = FALSE]
	attr(x, "directed") <- xdir[i]
	class(x) <- cl
	x
}

"[<-.mathgraph" <-
function(x, i, value)
{
	if(!inherits(value, "mathgraph"))
		stop("need mathgraph on right-hand side")
	cl <- class(x)
	x <- unclass(x)
	value <- unclass(value)
	ilen <- length(i)
	if(is.logical(i))
		ilen <- sum(rep(i, length = nrow(x)))
	if(nrow(value) != ilen)
		stop("replacement value not correct length")
	xdir <- attr(x, "directed")
	if(is.character(i))
		names(xdir) <- dimnames(x)[[1]]
	x[i,  ] <- value
	xdir[i] <- attr(value, "directed")
	attr(x, "directed") <- xdir
	class(x) <- cl
	x
}

"unique.mathgraph" <-
function(x, incomparables = FALSE, ...)
{
	dir <- attr(x, "directed")
	cl <- class(x)
	x <- unclass(x)
	xwork <- x
	if(!all(dir)) {
# worry about order of nodes in undirected edges
		xwork[!dir,  ] <- stable.apply(x[!dir,  , drop = FALSE], 1, sort)
	}
	xdup <- duplicated(paste(xwork[, 1], xwork[, 2], dir))
	x <- x[!xdup,  , drop = FALSE]
	dir <- dir[!xdup]
	attr(x, "directed") <- dir
	class(x) <- cl
	x
}
