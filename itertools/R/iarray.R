iarray <- function(X, MARGIN, ..., chunks, chunkSize, drop,
                   idx=lapply(dim(X), function(i) TRUE)) {
  dimx <- dim(X)

  # Verify that X has the dim attribute set and length > 0
  if (length(dimx) == 0)
    stop('dim(X) must have a positive length')

  # Check for unknown arguments
  if (length(list(...)) > 0) {
    nms <- names(list(...))
    if (is.null(nms) || '' %in% nms)
      stop('arguments other than X and MARGIN must be named')
    else
      stop('unused argument(s) ', paste(nms, collapse=', '))
  }

  # Don't allow both chunks and chunkSize
  if (! missing(chunks) && ! missing(chunkSize)) 
    stop('chunks and chunkSize cannot both be specified')

  # Get the number of value this iterator will return
  i <- 0
  mlen <- length(MARGIN)
  n <- dimx[MARGIN[mlen]]

  # Create an iterator based on chunking
  if (! missing(chunks)) {
    if (length(chunks) != 1 && length(chunks) != mlen)
      stop('length of chunks must be 1 or the same as MARGIN')
    if (missing(drop))
      drop <- FALSE
    if (length(chunks) == 1)
      chunks <- rep(chunks, mlen)
    it <- idiv(n, chunks=chunks[mlen])
  } else if (! missing(chunkSize)) {
    if (length(chunkSize) != 1 && length(chunkSize) != mlen)
      stop('length of chunkSize must be 1 or the same as MARGIN')
    if (missing(drop))
      drop <- FALSE
    if (length(chunkSize) == 1)
      chunkSize <- rep(chunkSize, mlen)
    it <- idiv(n, chunkSize=chunkSize[mlen])
  } else {
    if (missing(drop))
      drop <- TRUE
    it <- irep(1, times=n)
  }

  # Create a call object if this is the final dimension
  if (mlen == 1) {
    q <- as.call(c(list(as.name('['), as.name('X')), idx, list(drop=drop)))
    iq <- MARGIN + 2L
  }

  # Define the "nextElem" function
  nextEl <- if (mlen == 1) {
    function() {
      m <- nextElem(it)
      j <- i + m
      q[[iq]] <- if (m > 1) call(':', i + 1, j) else j
      i <<- j
      eval(q)
    }
  } else if (! missing(chunks)) {
    function() {
      m <- nextElem(it)
      j <- i + m
      idx[[MARGIN[mlen]]] <- if (m > 1) call(':', i + 1, j) else j
      i <<- j
      iarray(X, MARGIN[-mlen], chunks=chunks[-mlen], drop=drop, idx=idx)
    }
  } else if (! missing(chunkSize)) {
    function() {
      m <- nextElem(it)
      j <- i + m
      idx[[MARGIN[mlen]]] <- if (m > 1) call(':', i + 1, j) else j
      i <<- j
      iarray(X, MARGIN[-mlen], chunkSize=chunkSize[-mlen], drop=drop, idx=idx)
    }
  } else {
    function() {
      nextElem(it)  # returns 1 or throws 'StopIteration'
      i <<- i + 1
      idx[[MARGIN[mlen]]] <- i
      iarray(X, MARGIN[-mlen], drop=drop, idx=idx)
    }
  }

  obj <- list(nextElem=nextEl)
  class(obj) <- c('abstractiter', 'iter')
  obj
}
