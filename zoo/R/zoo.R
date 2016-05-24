zoo <- function (x = NULL, order.by = index(x), frequency = NULL) 
{
    ## process index "order.by"    
    if(length(unique(MATCH(order.by, order.by))) < length(order.by))
      warning(paste("some methods for", dQuote("zoo"),
      "objects do not work if the index entries in", sQuote("order.by"), "are not unique"))
    index <- ORDER(order.by)
    order.by <- order.by[index]

    if(is.matrix(x) || is.data.frame(x)) x <- as.matrix(x)
    if(is.matrix(x) && sum(dim(x)) < 1L) x <- NULL
    if(missing(x) || is.null(x))
      x <- numeric()
    else if(is.factor(x))         
      x <- factor(rep(as.character(x), length.out = length(index))[index],
        levels = levels(x), ordered = is.ordered(x))
    else if(is.matrix(x) || is.data.frame(x)) 
      x <- (x[rep(1:NROW(x), length.out = length(index)), , 
        drop = FALSE])[index, , drop = FALSE]
    else if(is.atomic(x)) 
      x <- rep(x, length.out = length(index))[index]
    else stop(paste(dQuote("x"), ": attempt to define invalid zoo object"))

    if(!is.null(frequency)) {
      delta <- suppressWarnings(try(diff(as.numeric(order.by)), silent = TRUE))
      freqOK <- if(class(delta) == "try-error" || any(is.na(delta))) FALSE
        else if(length(delta) < 1) TRUE
        else identical(all.equal(delta*frequency, round(delta*frequency)), TRUE)
      if(!freqOK) {
        warning(paste(dQuote("order.by"), "and", dQuote("frequency"),
        	"do not match:", dQuote("frequency"), "ignored"))
        frequency <- NULL
      } else {
        if(frequency > 1 && identical(all.equal(frequency, round(frequency)), TRUE))
	  frequency <- round(frequency)
      }
      if(!is.null(frequency) && identical(class(order.by), "numeric") | identical(class(order.by), "integer")) {
        orig.order.by <- order.by
        order.by <- floor(frequency * order.by + .0001)/frequency
        if(!isTRUE(all.equal(order.by, orig.order.by))) order.by <- orig.order.by
      }
    }

    attr(x, "oclass") <- attr(x, "class")
    attr(x, "index") <- order.by
    attr(x, "frequency") <- frequency
    class(x) <- if(is.null(frequency)) "zoo" else c("zooreg", "zoo")
    return(x)
}

print.zoo <- function (x, style = ifelse(length(dim(x)) == 0,
    "horizontal", "vertical"), quote = FALSE, ...) 
{
    style <- match.arg(style, c("horizontal", "vertical", "plain"))
    if (is.null(dim(x)) && length(x) == 0) style <- "plain"
    if (length(dim(x)) > 0 && style == "horizontal") style <- "plain"
    if (style == "vertical") {
	y <- as.matrix(coredata(x))
        if (length(colnames(x)) < 1) {
            colnames(y) <- rep("", NCOL(x))
        }
		if (NROW(y) > 0) {
			rownames(y) <- index2char(index(x), frequency = attr(x, "frequency"))
		}
        print(y, quote = quote, ...)
    }
    else if (style == "horizontal") {
        y <- as.vector(x)
        names(y) <- index2char(index(x), frequency = attr(x, "frequency"))
        print(y, quote = quote, ...)
    }
    else {
        cat("Data:\n")
        print(coredata(x))
        cat("\nIndex:\n")
        print(index(x))
    }
    invisible(x)
}

summary.zoo <- function(object, ...) 
{
	y <- as.data.frame(object)
	if (length(colnames(object)) < 1) {
		lab <- deparse(substitute(object))
		colnames(y) <- if (NCOL(object) == 1) lab
		  else paste(lab, 1:NCOL(object), sep=".")
	}
	if (NROW(y) > 0) {
		summary(cbind(data.frame(Index = index(object)), y), ...)
	} else summary(data.frame(Index = index(object)), ...)
}


is.zoo <- function(object)
  inherits(object, "zoo")

str.zoo <- function(object, ...)
{
  cls <- if(inherits(object, "zooreg")) "zooreg" else "zoo"
  if(NROW(object) < 1) cat(paste(sQuote(cls), "series (without observations)\n")) else {
    cat(paste(sQuote(cls), " series from ", start(object), " to ", end(object), "\n", sep = ""))
    cat("  Data:")
    str(coredata(object), ...)
    cat("  Index: ")
    str(index(object), ...)
    if(cls == "zooreg") cat(paste("  Frequency:", attr(object, "frequency"), "\n"))
  }
}

"[.zoo" <- function(x, i, j, drop = TRUE, ...)
{
  if(!is.zoo(x)) stop("method is only for zoo objects")
  rval <- coredata(x)
  n <- NROW(rval)
  n2 <- if(nargs() == 1) length(as.vector(rval)) else n
  if(missing(i)) i <- 1:n

  if (all(class(i) == "matrix")) i <- as.vector(i)
  ## also support that i can be index:
  ## if i is not numeric/integer/logical, it is interpreted to be the index
  if (all(class(i) == "logical"))
    i <- which(rep(i, length.out = n2))
  else if (inherits(i, "zoo") && all(class(coredata(i)) == "logical")) {
    i <- which(coredata(merge(zoo(,time(x)), i)))
  } else if(!((all(class(i) == "numeric") || all(class(i) == "integer")))) 
    i <- which(MATCH(index(x), i, nomatch = 0L) > 0L)
  
  if(length(dim(rval)) == 2) {
    drop. <- if (length(i) == 1) FALSE else drop
    rval <- if (missing(j)) rval[i, , drop = drop., ...]
      else rval[i, j, drop = drop., ...]
    if (drop && length(rval) == 1) rval <- c(rval)
    rval <- zoo(rval, index(x)[i])
  } else
    rval <- zoo(rval[i], index(x)[i])

  attr(rval, "oclass") <- attr(x, "oclass")
  attr(rval, "levels") <- attr(x, "levels")
  attr(rval, "frequency") <- attr(x, "frequency")
  if(!is.null(attr(rval, "frequency"))) class(rval) <- c("zooreg", class(rval))

  return(rval)
}

"[<-.zoo" <- function (x, i, j, value) 
{
  ## x[,j] <- value and x[] <- value can be handled by default method
  if(missing(i)) return(NextMethod("[<-"))

  ## otherwise do the necessary processing on i
  n <- NROW(coredata(x))
  n2 <- if(nargs() == 1) length(as.vector(coredata(x))) else n
  n.ok <- TRUE
  value2 <- NULL
  
  if (all(class(i) == "matrix")) i <- as.vector(i)
  if (all(class(i) == "logical")) {
    if (length(i) == n) {
		i <- which(i)
		n.ok <- TRUE
	} else {
		i <- which(rep(i, length.out = n))
		n.ok <- all(i <= n2)
	}
  } else if (inherits(i, "zoo") && all(class(coredata(i)) == "logical")) {
    i <- which(coredata(merge(zoo(,time(x)), i)))
    n.ok <- all(i <= n2)
  } else if(!((all(class(i) == "numeric") || all(class(i) == "integer")))) {
    ## all time indexes in index(x)?
    i.ok <- MATCH(i, index(x), nomatch = 0L) > 0L
    if(any(!i.ok)) {
      if(is.null(dim(value))) {
        value2 <- value[!i.ok]
        value <- value[i.ok]
      } else {
        value2 <- value[!i.ok,, drop = FALSE]
        value <- value[i.ok,, drop = FALSE]      
      }
      i2 <- i[!i.ok]
      i <- i[i.ok]
    }
    i <- which(MATCH(index(x), i, nomatch = 0L) > 0L)
    n.ok <- all(i <= n2)
  }
  if(!n.ok | any(i < 1)) stop("Out-of-range assignment not possible.")
  rval <- NextMethod("[<-")

  if(!is.null(value2)) {
    rval2 <- if(missing(j)) zoo(value2, i2) else {
      value2a <- matrix(NA, nrow = length(i2), ncol = NCOL(rval))
      colnames(value2a) <- colnames(rval)
      value2a[, j] <- value2
      zoo(value2a, i2)
    }
    rval <- c(rval, rval2)
  }
  return(rval)
}

"$.zoo" <- function(object, x) {
  if(length(dim(object)) != 2) stop("not possible for univariate zoo series")
  if(is.null(colnames(object))) stop("only possible for zoo series with column names")
  wi <- pmatch(x, colnames(object))
  if(is.na(wi)) NULL else object[, wi]
}

"$<-.zoo" <- function(object, x, value) {
  if(length(object) == 0L) {
    is.plain <- function(x)
      all(class(x) %in% c("array", "integer", "numeric", "factor", "matrix", "logical"))
    if(is.plain(value)) value <- zoo(value,
      if(length(index(object))) index(object) else seq_along(value), attr(object, "frequency"))
    return(setNames(merge(object, value, drop = FALSE), x))
  }
  if(length(dim(object)) != 2) stop("not possible for univariate zoo series")
  if(NCOL(object) > 0L & is.null(colnames(object))) stop("only possible for zoo series with column names")
  wi <- match(x, colnames(object))
  if(is.na(wi)) {
    if(!is.null(value)) {
      object <- cbind(object, value)
      if(is.null(dim(object))) dim(object) <- c(length(object), 1)
      colnames(object)[NCOL(object)] <- x  
    }
  } else {
    if(is.null(value)) {
      object <- object[, -wi, drop = FALSE]
    } else {   
      object[, wi] <- value
    }
  }
  object
}

head.zoo <- function(x, n = 6, ...) {
	stopifnot(length(n) == 1L)
	xlen <- NROW(x)
    n <- if (n < 0L) 
        max(NROW(x) + n, 0L)
    else min(n, xlen)
	if (length(dim(x)) == 0) x[seq_len(n)]
	else x[seq_len(n),, drop = FALSE]
}
 
tail.zoo <- function(x, n = 6, ...) {
    stopifnot(length(n) == 1L)
    xlen <- NROW(x)
    n <- if (n < 0L) 
        max(xlen + n, 0L)
    else min(n, xlen)
	if (length(dim(x)) == 0) x[seq.int(to = xlen, length.out = n)]
	else x[seq.int(to = xlen, length.out = n),, drop = FALSE]
}

range.zoo <- function(..., na.rm = FALSE)
    range(sapply(list(...), coredata), na.rm = na.rm)


scale.zoo <- function (x, center = TRUE, scale = TRUE) {
	x[] <- xs <- scale(coredata(x), center = center, scale = scale)
	attributes(x) <- c(attributes(x), attributes(xs))
	x
}

with.zoo <- function(data, expr, ...) {
    stopifnot(length(dim(data)) == 2)
    eval(substitute(expr), as.list(data), enclos = parent.frame())
}

xtfrm.zoo <- function(x) coredata(x)

subset.zoo <- function (x, subset, select, drop = FALSE, ...) 
{
    if (missing(select)) 
        vars <- TRUE
    else {
        nl <- as.list(1:ncol(x))
        names(nl) <- colnames(x)
        vars <- eval(substitute(select), nl, parent.frame())
    }
    if (missing(subset)) {
        subset <- rep(TRUE, NROW(x))
    } else {
        e <- substitute(subset)
	if("time" %in% colnames(x)) {
	  xdf <- as.data.frame(x)
          subset <- eval(e, xdf, parent.frame())
          xdf$time <- time(x)
          subset2 <- eval(e, xdf, parent.frame())
	  if(!identical(subset, subset2))
  	      warning("'time' is a column in 'x' (not the time index)")
	} else {
          subset <- eval(e, cbind(as.data.frame(x), time = time(x)), parent.frame())
	}
        if (!is.logical(subset)) stop("'subset' must be logical")
    }
    x[subset & !is.na(subset), vars, drop = drop]
}

names.zoo <- function(x) {
  cx <- coredata(x)
  if(is.matrix(cx)) colnames(cx) else names(cx)
}

"names<-.zoo" <- function(x, value) {
  if(is.matrix(coredata(x))) {
    colnames(x) <- value
  } else {
    names(coredata(x)) <- value
  }
  x
}

rev.zoo <- function(x) {
  ix <- rev(ORDER(time(x)))
  zoo(coredata(x), time(x)[ix])
}

ifelse.zoo <- function(test, yes, no) {
	if(!is.zoo(test)) test <- zoo(test, index(yes))
	merge(test, yes, no, retclass = NULL)
	ifelse(test, yes, no)
}

mean.zoo <- function(x, ...)  mean(coredata(x), ...)

median.zoo <- function(x, na.rm = FALSE)  median(coredata(x), na.rm = na.rm)

quantile.zoo <- function(x, ...) quantile(coredata(x), ...)

transform.zoo <- function(`_data`, ...)
{
  ## zoo series elements
  if(is.null(dim(coredata(`_data`)))) warning("transform() is only useful for matrix-based zoo series")
  ix <- index(`_data`)
  freq <- attr(`_data`, "frequency")
  `_data` <- as.data.frame(`_data`)

  ## get transformations
  e <- eval(substitute(list(...)), `_data`, parent.frame())
  trafo <- names(e)
  inx <- match(trafo, names(`_data`))
  matched <- !is.na(inx)

  ## evaluate transformations
  if(any(matched)) {
    `_data`[inx[matched]] <- e[matched]
    `_data` <- data.frame(`_data`)
  }
  if (!all(matched)) 
    `_data` <- do.call("data.frame", c(list(`_data`), e[!matched]))

  ## return zoo
  zoo(`_data`, ix, freq)
}

`dim<-.zoo` <- function(x, value) {
  d <- dim(x)
  l <- length(x)
  ok <- isTRUE(all.equal(d, value)) ||                                  ## no change
    (is.null(d) && l == 0L && all(value == c(length(index(x)), 0L))) || ## zero-length vector -> 0-column matrix
    (is.null(d) && l > 0L && all(value == c(l, 1L))) ||                 ## positive-length vector -> 1-column matrix
    (!is.null(d) && d[2L] <= 1L && is.null(value))                      ## 0- or 1-column matrix -> vector
  if(!ok) warning("setting this dimension may lead to an invalid zoo object")
  NextMethod()
}
