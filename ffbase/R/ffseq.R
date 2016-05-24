#' Sequence Generation of \code{ff} vectors.
#'
#' Similar as \code{seq_len} in the base package, generating an \code{ff} vector.
#'
#' @export
#' @example ../examples/ffseq.R
#' @param length.out desired length of the sequence. Only non-negative numbers larger than 0 are allowed.
#' @return An ff vector of integers with range from 1 to length.out
#' @seealso \code{\link[base]{seq_len}}
ffseq_len <- function(length.out){
	BATCHBYTES <- getOption("ffbatchbytes")
	length.out <- as.integer(length.out)
	if(length.out < 1){
		stop("length.out needs to be positive")
	}
	bysize <- floor(BATCHBYTES / .rambytes["integer"])
	x <- ff(NA, length=length.out, vmode = "integer")
	for (i in chunk(1, length.out, by=bysize)){
    Log$chunk(i)
		idx <- as.integer(hi(from=min(i), to=max(i), by = 1L, maxindex = max(i)))
    x[idx] <- idx
  }
	x
}


#' Sequence Generation of \code{ff} vectors.
#'
#' Similar as \code{seq} in the base package, generating an \code{ff} vector.
#'
#' @export
#' @example ../examples/ffseq.R
#' @param from the starting value of the sequence
#' @param to the end (maximal) value of the sequence
#' @param by number, increment of the sequence
#' @param length.out desired length of the sequence. Only non-negative numbers larger than 0 are allowed.
#' @param along.with take the length from the length of this argument
#' @param ... arguments passed to or from methods
#' @return An ff vector with the generated sequence, similar as what \code{seq} generates but as an ff vector. \cr
#' Mark: in case this would generate a sequence of length 0, will return integer().
#' @seealso \code{\link[base]{seq}}
ffseq <- function(from = 1, to = 1, by = ((to - from)/(length.out - 1)), 
                  length.out = NULL, along.with = NULL, ...){
  if ((One <- nargs() == 1L) && !missing(from)) {
    lf <- length(from)
    return(if (mode(from) == "numeric" && lf == 1L) 1L%ff:%from else if (lf) 1L%ff:%lf else integer())
  }
  if (!missing(along.with)) {
    length.out <- length(along.with)
    if (One) 
      return(if (length.out) ffseq_len(length.out) else integer())
  }else if (!missing(length.out)) {
    len <- length(length.out)
    if (!len) 
      stop("argument 'length.out' must be of length 1")
    if (len > 1L) {
      warning("first element used of 'length.out' argument")
      length.out <- length.out[1L]
    }
    length.out <- ceiling(length.out)
  }
  if (length(list(...))) 
    warning(gettextf("extra argument(s) %s will be disregarded", 
                     paste(sQuote(names(list(...))), collapse = ", ")), 
            domain = NA)
  if (!missing(from) && length(from) != 1L) 
    stop("'from' must be of length 1")
  if (!missing(to) && length(to) != 1L) 
    stop("'to' must be of length 1")
  if (is.null(length.out)) {
    if (missing(by)){
      from%ff:%to
    }else {
      del <- to - from
      if (del == 0 && to == 0) 
        return(ff(to))
      n <- del/by
      if (!(length(n) && is.finite(n))) {
        if (length(by) && by == 0 && length(del) && del == 0) 
          return(ff(from))
        stop("invalid (to - from)/by in seq(.)")
      }
      if (n < 0L) 
        stop("wrong sign in 'by' argument")
      if (n > .Machine$integer.max) 
        stop("'by' argument is much too small")
      dd <- abs(del)/max(abs(to), abs(from))
      if (dd < 100 * .Machine$double.eps) 
        return(ff(from))
      if (is.integer(del) && is.integer(by)) {
        n <- as.integer(n)
        from + (0L%ff:%n) * by
      }else {
        n <- as.integer(n + 1e-10)
        x <- from + (0L%ff:%n) * by
        if (by > 0) {
          idx <- to < x
          idx <- ffwhich(idx, idx==TRUE)
          if(length(idx) != 0){
            x[idx] <- to  
          }
        }else{
          idx <- to > x
          idx <- ffwhich(idx, idx==TRUE)
          if(length(idx) != 0){
            x[idx] <- to  
          }
        }
        return(x)
        #if (by > 0) 
        #  pmin(x, to)
        #else pmax(x, to)
      }
    }
  }else if (!is.finite(length.out) || length.out < 0L){
    stop("length must be non-negative number")
  }else if (length.out == 0L) {
    integer()
  }else if (One) {
    ffseq_len(length.out)
  }else if (missing(by)) {
    if (missing(to)) 
      to <- from + length.out - 1L
    if (missing(from)) 
      from <- to - length.out + 1L
    if (length.out > 2L){
      if (from == to){
        rep.int(from, length.out)
      }else {
        c(ff(from), from + ffseq_len(length.out - 2L) * by, ff(to)) 
      }
    }else{
      c(ff(from), ff(to))[ffseq_len(length.out)]
    } 
  }
  else if (missing(to)) 
    from + (0L%ff:%(length.out - 1L)) * by
  else if (missing(from)) 
    to - ((length.out - 1L)%ff:%0L) * by
  else stop("too many arguments")
}
"%ff:%" <- function(from, to){
  ## mimics the colon operator but creates an ff integer vector
  n <- (to - from)+1
  cumsum(ff(1, length = n, vmode = "integer")) + (from-1)
}



