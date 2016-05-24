# ========================================================================
# rvmapply  -  apply a function to multiple rv objects
# ========================================================================
# Note. Won't work with functions allowing "blank" arguments
# such as "[" (e.g. x[y,,]). The functions "[" and "[<-" use
# modified versions of simmapply.
#

rvmapply <- function (FUN, ..., MoreArgs=NULL, SIMPLIFY=FALSE, USE.NAMES=TRUE, SAMPLESIZE=NULL) {
  a <- list(...)
  dim.a.names <- dimnames(a)
  a.names <- names(a)
  dimnames(a) <- dim.a.names
  names(a) <- a.names
  a <- lapply(a, .sims.as.list)
  if (!is.null(SAMPLESIZE)) {
    m <- max(sapply(a, length))
    s <- (sample(1:m, size=SAMPLESIZE, replace=TRUE)-1)
    a <- lapply(a, function (x) x[(s %% length(x))+1])
  }
  a <- .Primitive("c")(FUN = FUN, a, SIMPLIFY = FALSE, USE.NAMES=USE.NAMES)
  a$MoreArgs <- MoreArgs
  S <- do.call(mapply, args = a)
  S <- lapply(S, function (x) if (is.null(x)) NA else x) ## DEBUG:: OK??
  r <- rvsims(S)
  ## DEBUG: match the largest-dimensional param in list(...) and 
  ## set the dimnames if they match
  ## if (isTRUE(all.equal(dim(r), dim(x)))) {
  ##  dimnames(r) <- dimnames(x)
  ##}
  return(r)
}



