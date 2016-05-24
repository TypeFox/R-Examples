window.zoo <- function(x, index. = index(x), start = NULL, end = NULL, ...)
{
  all.indexes <- index(x)
  in.index <- MATCH(all.indexes, index., nomatch = 0) > 0

  if(length(start) == 2 && !is.null(attr(x, "frequency")) && is.numeric(all.indexes)) {
    freq <- attr(x, "frequency")
    start <- floor(start[1]*freq + (start[2] - 1) + .0001)/freq
  }
  if(length(end) == 2 && !is.null(attr(x, "frequency")) && is.numeric(all.indexes)) {
    freq <- attr(x, "frequency")
    end <- floor(end[1]*freq + (end[2] - 1) + .0001)/freq
  }

  if(is.null(start)) {
    if(is.null(end)) {
      wi <- which(MATCH(all.indexes, index., nomatch = 0) > 0)
      return(x[wi, , drop = FALSE])
    } else {
      wi <- which(in.index & all.indexes <= end)
      return(x[wi, , drop = FALSE])
    }
  } else {
    if(is.null(end)) {
      wi <- which(in.index & all.indexes >= start)
    } else {
      wi <- which(in.index & all.indexes >= start & all.indexes <= end)
    }
    return(x[wi, , drop = FALSE])
  }
}

"window<-.zoo" <- function(x, index. = index(x), start = NULL, end = NULL, ..., value)
{
  ix <- index(x)
  stopifnot(all(MATCH(index., ix, nomatch = 0) > 0))
  
  if(length(start) == 2 && !is.null(attr(x, "frequency")) && is.numeric(ix)) {
    freq <- attr(x, "frequency")
    start <- floor(start[1]*freq + (start[2] - 1) + .0001)/freq
  }
  if(length(end) == 2 && !is.null(attr(x, "frequency")) && is.numeric(ix)) {
    freq <- attr(x, "frequency")
    end <- floor(end[1]*freq + (end[2] - 1) + .0001)/freq
  }
  
  if (!is.null(start)) index. <- index.[index. >= start]
  if (!is.null(end)) index. <- index.[index. <= end]

  wi <- which(MATCH(ix, index., nomatch = 0) > 0)
  if (length(dim(x)) == 0)
  	  x[wi] <- value
  else
  	  x[wi,] <- value
  return(x)
}
 
lag.zoo <- function(x, k = 1, na.pad = FALSE, ...)
{
   if (length(k) > 1) {
	if (is.null(names(k))) names(k) <- paste("lag", k, sep = "")
	return(do.call("merge.zoo", lapply(k, lag.zoo, x = x, na.pad = na.pad, ...)))
   }
   nr <- NROW(x)
   if (k != round(k)) {
	k <- round(k)
	warning("k is not an integer")
   }
   if (k == 0) return(x)
   if (abs(k) > nr) k <- nr
   if (k > 0)  {
	   xx <- x[-seq(1, length.out = k),, drop = FALSE]
	   attr(xx, "index") <- index(x)[-seq(to = nr, length.out = k)]
   } else {
	   xx <- x[-seq(to = nr, length.out = -k),, drop = FALSE]
	   attr(xx, "index") <- index(x)[-seq(1, length.out = -k)]
   }
   if (na.pad) merge(zoo(,time(x)), xx, all = c(TRUE, FALSE)) else xx
}



lag.zooreg <- function(x, k = 1, na.pad = FALSE, ...)
{
   if (length(k) > 1) {
	if (is.null(names(k))) names(k) <- paste("lag", k, sep = "")
	return(do.call("merge.zoo", lapply(k, lag.zooreg, x = x, na.pad = na.pad, ...)))
   }
   x0 <- x
   nr <- NROW(x)
   freq <- attr(x, "frequency")
   
   if (k != round(k)) warning("k is not an integer")
   k <- round(k)

   ix <- index(x)
   ix <- if(identical(class(ix), "numeric") | identical(class(ix), "integer"))
     floor(freq*ix - k + .0001)/freq else ix - k/freq
   index(x) <- ix

   if (na.pad) merge(x, zoo(, time(x0))) else x
}

diff.zoo <- function(x, lag = 1, differences = 1, arithmetic = TRUE, na.pad = FALSE, ...)
{
    ix <- index(x)
    stopifnot(differences >= 1)
    if (!arithmetic) x <- log(x)
	if (lag > 0) for(i in 1:differences) x <- x - lag(x, k = -lag, ...)
	else for(i in 1:differences) x <- lag(x, k = -lag, ...) - x
    if (!arithmetic) x <- exp(x)
    if (na.pad) merge(zoo(,ix), x, all = c(TRUE, FALSE)) else x
}
