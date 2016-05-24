read.zoo <- function(file, format = "", tz = "", FUN = NULL,
  regular = FALSE, index.column = 1, drop = TRUE, FUN2 = NULL, 
  split = NULL, aggregate = FALSE, ..., text)
{

  if (missing(file) && !missing(text)) {
        file <- textConnection(text)
        on.exit(close(file))
  }

  ## if file is a vector of file names
  if (is.character(file) && length(file) > 1) {
	mc <- match.call()
	pf <- parent.frame()
	L <- sapply(file, function(file) eval(replace(mc, 2, file), pf), 
			simplify = FALSE)
	return(do.call("merge", L))
  }

  ## read data
  rval <- if (is.data.frame(file)) file else read.table(file, ...)

  ## if time index appears to be already processed, use FUN = identity
  if (is.data.frame(file) && 
      length(index.column) == 1 && 
      !is.character(rval[[index.column]]) &&
      !is.factor(rval[[index.column]]) &&
      missing(tz) &&
      missing(format) &&
      missing(FUN)) FUN <- identity

  ## if time index is POSIXlt it is coerced to POSIXct
  if (is.data.frame(file) && 
      length(index.column) == 1 && 
      inherits(rval[[index.column]], "POSIXlt")) rval[[index.column]] <- as.POSIXct(rval[[index.column]])

  # returns TRUE if a formal argument x has no default
  no.default <- function(x) typeof(x) %in% c("symbol", "language")

  if (is.null(FUN) && is.null(FUN2)) {
	index.column <- as.list(index.column)
  } else if (identical(FUN, paste)) {
	index.column <- as.list(index.column)
  } else if (is.null(FUN) && identical(FUN2, paste)) {
	index.column <- as.list(index.column)
  } else if (!is.null(FUN) && !is.list(index.column) && length(index.column) <=
		length(sapply(formals(match.fun(FUN)), no.default))) {
	index.column <- as.list(index.column)
  } else if (is.null(FUN) && !is.null(FUN2) && length(index.column) <= 
		length(sapply(formals(match.fun(FUN2)), no.default))) {
	index.column <- as.list(index.column)
  }

  if (is.list(index.column) && length(index.column) == 1 && 
	index.column[[1]] == 1) index.column <- unlist(index.column)

  is.index.column <- seq_along(rval) %in% unname(unlist(index.column)) |
	names(rval) %in% unname(unlist(index.column))

  name.to.num <- function(x) if (is.character(x))
		match(x, names(rval), nomatch = 0) else x
  index.column <- if (is.list(index.column)) index.column <- 
	lapply(index.column, name.to.num)
  else name.to.num(index.column)

  ## convert factor columns in index to character
  is.fac <- sapply(rval, is.factor)
  is.fac.index <- is.fac & is.index.column
  if (any(is.fac.index)) rval[is.fac.index] <- 
	lapply(rval[is.fac.index], as.character)

  ## if file does not contain index or data
  if(NROW(rval) < 1) {
    if(is.data.frame(rval)) rval <- as.matrix(rval)
    if(NCOL(rval) > 1) rval <- rval[, ! is.index.column, drop = drop]
    rval2 <- zoo(rval)
    return(rval2)
  }

  ## extract index
  if(NCOL(rval) < 1) stop("data file must specify at least one column")
  
  ## extract index, retain rest of the data
  ix <- if (identical(index.column, 0) || identical(index.column, list(0)) ||
	identical(index.column, 0L) || identical(index.column, list(0L))) {
	attributes(rval)$row.names
  } else if (is.list(index.column)) {
	sapply(index.column, function(j) rval[,j], simplify = FALSE)
  } else rval[,index.column]

  # split. is col number of split column (or Inf, -Inf or NULL)
  split. <- if (is.character(split)) match(split, colnames(rval), nomatch = 0)
  else split

  rval2 <- if (is.null(split.)) {
    rval[ , ! is.index.column, drop = FALSE]
  } else {

     split.values <- if (is.character(split) || is.finite(split)) rval[, split]
	 else {
		# Inf: first value in each run is first series, etc.
	    # -Inf: last value in each run is first series, etc.
		if (identical(split, Inf)) ave(ix, ix, FUN = seq_along)
	    else if (identical(split, -Inf)) ave(ix, ix, FUN = function(x) rev(seq_along(x)))
	    else ix
	 }

	 if (0 %in% split.) {
		stop(paste("split:", split, "not found in colnames:", colnames(rval)))
	 }
	 rval[,-c(if (is.finite(split.)) split. else 0, which(is.index.column)), drop = FALSE]
  }

  if(is.factor(ix)) ix <- as.character(ix)
  rval3 <- if(is.data.frame(rval2)) as.matrix(rval2) else  if(is.list(rval2)) t(rval2) else rval2
  
  if(NCOL(rval3) == 1 && drop) rval3 <- drop(rval3)

    
  ## index transformation functions

  toDate <- if(missing(format) || is.null(format)) {
     function(x, ...) as.Date(format(x, scientific = FALSE))
  } else {
     function(x, format) {
       if(any(sapply(c("%H", "%M", "%S"), function(y) grepl(y, format, fixed = TRUE)))) {
         warning("the 'format' appears to be for a date/time, please specify 'tz' if you want to create a POSIXct time index")
       }
       as.Date(format(x, scientific = FALSE), format = format)
     }
  }

  toPOSIXct <- if (missing(format) || is.null(format)) {
        function(x, tz) as.POSIXct(format(x, scientific = FALSE), tz = tz)
  } else function(x, format, tz) {
        as.POSIXct(strptime(format(x, scientific = FALSE), tz = tz, format = format))
  }

  toDefault <- function(x, ...) {
    rval. <- try(toPOSIXct(x, tz = ""), silent = TRUE)
    if(inherits(rval., "try-error"))
      rval. <- try(toDate(x), silent = TRUE)
    else {
      hms <- as.POSIXlt(rval.)
      hms <- hms$sec + 60 * hms$min + 3600 * hms$hour
      if(isTRUE(all.equal(hms, rep.int(hms[1], length(hms))))) {
        rval2. <- try(toDate(x), silent = TRUE)
        if(!inherits(rval2., "try-error")) rval. <- rval2.
      }
    }
    if(inherits(rval., "try-error")) rval. <- rep(NA, length(x))
    return(rval.)
  }

  toNumeric <- function(x, ...) x
  
  ## setup default FUN

  if ((missing(FUN) || is.null(FUN)) && !missing(FUN2) && !is.null(FUN2)) {
	FUN <- FUN2
	FUN2 <- NULL
  }

  FUN0 <- NULL
  if(is.null(FUN)) {
	if (is.list(index.column)) FUN0 <- paste
    FUN <- if (!missing(tz) && !is.null(tz)) toPOSIXct
        else if (!missing(format) && !is.null(format)) toDate
        else if (is.numeric(ix)) toNumeric
        else toDefault
  }

  FUN <- match.fun(FUN)

 processFUN <- function(...) {
	if (is.data.frame(..1)) FUN(...)
	else if (is.list(..1)) {
		if (is.null(FUN0)) do.call(FUN, c(...))
		else {
			L <- list(...)
			L[[1]] <- do.call(FUN0, L[[1]])
			do.call(FUN, L)
		}
	} else FUN(...)
  }

  ## compute index from (former) first column
  ix <- if (missing(format) || is.null(format)) {
    if (missing(tz) || is.null(tz)) processFUN(ix) else processFUN(ix, tz = tz)
  } else {
    if (missing(tz) || is.null(tz)) processFUN(ix, format = format) 
    else processFUN(ix, format = format, tz = tz)
  }

  if (!is.null(FUN2)) {
	FUN2 <- match.fun(FUN2)
	ix <- FUN2(ix)
  }
  
  ## sanity checking
  if(any(is.na(ix))) {
    idx <- which(is.na(ix))
	msg <- if (length(idx) == 1)
		paste("index has bad entry at data row", idx)
	else if (length(idx) <= 100)
		paste("index has bad entries at data rows:", paste(idx, collapse = " "))
	else paste("index has", length(idx), "bad entries at data rows:", 
		paste(head(idx, 100), collapse = " "), "...")
	stop(msg)
  }
  if(length(ix) != NROW(rval3)) stop("index does not match data")
  
  ## setup zoo object and return 
  ## Suppress duplicates warning if aggregate specified
  if(identical(aggregate, TRUE)) {
    agg.fun <- mean
  } else if(identical(aggregate, FALSE)) {
    agg.fun <- NULL
  } else {
    agg.fun <- match.fun(aggregate)
    if(!is.function(agg.fun)) stop(paste("invalid specification of", sQuote("aggregate")))
  }
  remove(list = "aggregate")

  if (is.null(split)) {

	rval4 <- if (!is.null(agg.fun)) aggregate(zoo(rval3), ix, agg.fun)
	else zoo(rval3, ix)
    rval8 <- if(regular && is.regular(rval4)) as.zooreg(rval4) else rval4

  } else {

	split.matrix <- split.data.frame
	rval4 <- split(rval3, split.values)
	ix <- split(ix, split.values)
	# rval5 <- mapply(zoo, rval4, ix)
    rval5 <- if (!is.null(agg.fun)) {
		lapply(seq_along(rval4), function(i) {
			aggregate(zoo(rval4[[i]]), ix[[i]], agg.fun)
		})
	} else lapply(seq_along(rval4), function(i) zoo(rval4[[i]], ix[[i]]))
	names(rval5) <- names(rval4)
    rval6 <- if(regular) {
		lapply(rval5, function(x) if (is.regular(x)) as.zooreg(x) else x)
	} else rval5

	rval8 <- do.call(merge, rval6)

  }
	
  return(rval8)
}

write.zoo <- function(x, file = "", index.name = "Index",
  row.names = FALSE, col.names = NULL, ...)
{
  if(is.null(col.names)) col.names <- !is.null(colnames(x))
  dx <- as.data.frame(x)  
  stopifnot(all(names(dx) != index.name))
  dx[[index.name]] <- index(x)
  dx <- dx[, c(ncol(dx), 1:(ncol(dx)-1))]
  write.table(dx, file = file, row.names = row.names, col.names = col.names, ...)
}
