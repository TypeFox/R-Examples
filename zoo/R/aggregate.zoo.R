aggregate.zoo <- function(x, by, FUN = sum, ..., regular = NULL, frequency = NULL)
{
  ## index processing
  my.unique <- function(x) {
    ix <- MATCH(x, x) == seq_len(length(x))
    x[ix]
  }
  if(is.function(by)) by <- by(index(x))
  if(!is.list(by)) by <- list(by)

  ## sanity checks and option processing
  stopifnot(length(time(x)) == length(by[[1]]))
  if(is.null(frequency)) {
    if(is.null(regular)) regular <- inherits(x, "zooreg")
  } else {
    if(identical(regular, FALSE)) warning(paste(sQuote("regular"), "is ignored"))
    regular <- TRUE
  }

  ## aggregate data
  by_integer <- list(MATCH(by[[1]], by[[1]]))
  df <- aggregate(coredata(x), by_integer, match.fun(FUN), ...)
  if(length(unique(as.character(df[,1]))) == length(df[,1]))
      row.names(df) <- df[, 1]
  df <- df[, -1]
  if(is.matrix(x)) df <- as.matrix(df)
  
  ## regularity processing, set up return value
  ix <- my.unique(by[[1]])
  rval <- zoo(df, ix[!is.na(ix)])
  
  if(regular) {
    freq <- ifelse(is.null(frequency), frequency(rval), frequency)
    rval <- zoo(df, ix, freq)
  }
  
  return(rval)
}

# works even if zoo series has duplicates among its times
split.zoo <- function(x, f, drop = FALSE, ...) {
    ix <- time(x)
	xc <- coredata(x)
	if (length(dim(xc)) < 2) {
		lapply(split(seq_along(xc), f, drop = drop, ...), 
			function(ind) zoo(xc[ind], ix[ind]))
	} else {
		lapply(split(seq_len(nrow(xc)), f, drop = drop, ...), 
			function(ind) zoo(xc[ind, , drop = drop], ix[ind]))
	}
}


