# ========================================================================
# rvattr  -  return the rvsim attribute for each component
# ========================================================================
# returns a list.

rvattr <- function(x, attrib=NULL)
{
  a <- lapply(unclass(x), "attr", "rvsim")
  if (!is.null(attrib)) {
    a <- lapply(a, "[[", attrib)
    nulls <- sapply(a, is.null)
    a[nulls] <- NA
  }
  return(a)
}

# rvattr(x, 'name') <- vector of values for each component or x; or,
# rvattr(x) <- list of attributes: list(attributename=vector of values, attributename2=...)
# e.g. list(Rhat=c(1.01,1.03,1.23), n.eff=c(200,100,16)).
#



"rvattr<-" <- function(x, attrib=NULL, by.name=FALSE, value) {
  "Set attributes of each component of a random vector"
  ## Note: 'value' can be a list or a vector - components are extracted by [[...]]
  nv <- names(value)
  .stop <- function (...) {
    stop("[rvattr<-] ", ..., call.=FALSE)
  }
  if (is.null(attrib)) {
    if (length(nv) > 0) {
      if (any(is.na(nv))) {
        .stop("If 'attrib' is NULL, then value must be a *named* list with no missing values")
      }
      for (a in nv) {
        rvattr(x, attrib=a) <- value[[a]]
      }
    } else {
      .stop("If 'attrib' is NULL, then value must be a named list mapping attributes to values")
    }
  } else if (length(value) > length(x)) {
    .stop("Too many values")
  } else if (by.name) {
    if (length(nv) == 0) {
      .stop("No names - cannot match by name")
    }
    if (any(duplicated(nv))) {
      .stop("Duplicate names in the list of values")
    }
    bad.nv <- (is.na(nv) | nchar(nv) == 0) ## Cannot match by NA or ""
    ## Attempt to match by name of component - but only in strict conditions
    if (any(bad.nv)) {
      .stop("Missing values or empty strings in the names of 'value' => cannot assign values by name")
    }
    nx <- names(x)
    bad.nx <- (is.na(nx) | nchar(nx) == 0) ## Cannot match by NA or ""
    good.names <- nx[! bad.nx]
    if (all(bad.nx)) {
      .stop("Only NAs or empty strings in the names of the rv => cannot assign values by name")
    }
    if (! all(nv %in% good.names)) {
      .stop("Some names in value do not exist in names of the rv")
    }
    remaining <- rep(TRUE, length(x))
    for (i in seq_along(value)) {
      name <- nv[i]
      w <- which(nx %in% name)
      if (length(w) == 0L) {
        .stop("Internal error - all names were supposed to be there in names(x)!!!")
      }
      for (i in w) {
        A <- attr(x[[i]], "rvsim")
        if (is.null(A)) {
          A <- list()
        }
        A[[attrib]] <- value[[i]]
        attr(x[[i]], "rvsim") <- A
      }
    }
  } else if (length(x) != length(value)) {
    .stop("The value must be of the same length as that of the rv, or the matching must be done by name.")
  } else {
    ## Match by position
    for (i in seq_along(x)) {
      A <- attr(x[[i]], "rvsim")
      if (is.null(A)) { A <- list() }
      A[[attrib]] <- value[[i]]
      attr(x[[i]], "rvsim") <- A
    }
  }
  return(x)
}


