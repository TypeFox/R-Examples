##' Create a mutaframe, a mutable data.frame
##' @export
##' @param ... Objects to coerce to a mutaframe and combine column-wise
##' @param row.names optional, the character vector of row names
##' @return a mutaframe
##' @exportClass mutaframe
##' @aliases mutaframe-class
mutaframe <- function(..., row.names = NULL) {
  listData <- list(...)

  # If no data, return NULL mutaframe
  if (length(listData) == 0) return(.mutaframe())

  # Work out names
  orig_names <- names(listData) %||% rep(NA_character_, length(listData))
  # Don't touch the child elements that already have colnames
  null_names <- sapply(sapply(listData, colnames, simplify = FALSE), is.null)
  names(listData)[null_names] <- variable_names(orig_names[null_names])
  names(listData)[!null_names] <- orig_names[!null_names] <- NA_character_
  
  varnames <- as.list(names(listData))
  varlist <- vector("list", length(listData))
  nrows <- ncols <- integer(length(varnames))
  for (i in seq_along(listData)) {
    element <- listData[[i]]

    if (is.function(element)) {
      element <- element()
    }
    if (!is.mutaframe(element)) {
      element <- as.data.frame(element)
    }
    
    nrows[i] <- nrow(element)
    ncols[i] <- ncol(element)
    
    varlist[[i]] <-
      if (is.environment(listData[[i]]))
        proxy_bindings(element, names(element))
      else if (is.function(listData[[i]]))
        listData[[i]]
      else as.list(element)
    if ((length(dim(listData[[i]])) > 1) ||
        (ncol(element) > 1)) {
          
      if (is.na(orig_names[i])) {
        varnames[[i]] <- colnames(element)
      } else {
        varnames[[i]] <- paste(varnames[[i]], colnames(element), sep = ".")
      }
    }
  }
  nr <- max(nrows)
  for (i in which((nrows > 0L) & (nrows < nr) & (nr %% nrows == 0L))) {
    recycle <- rep(seq_len(nrows[i]), length.out = nr)
    varlist[[i]] <- lapply(varlist[[i]], function(x) {
      if (is.function(x))
        function(v) {
          if (!missing(v))
            .irreversible("replication")
          x()[recycle, drop=FALSE]
        }
      else x[recycle, drop=FALSE]
    })
    nrows[i] <- nr
  }
  if (!all(nrows == nr))
    stop("different row counts implied by arguments")
  varlist <- unlist(varlist, recursive = FALSE, use.names = FALSE)
  names(varlist) <- make.names(unlist(varnames[ncols > 0L]), unique = TRUE)
  
  if (!is.null(row.names)) {
    if (any(is.na(row.names)))
      stop("missing values in 'row.names'")
    if (length(varlist) && length(row.names) != nr)
      stop("invalid length of row names")
    if (any(duplicated(row.names)))
      stop("duplicate row names")
    row.names <- as.character(row.names)
  } else row.names <- as.character(seq(max(nr)))


  env <- .mutaframe(varlist, row.names)
  # provenance(env) <- sys.call()
  env
}

##' Raw constructor.
##' Constructs a mutaframe from a list of variables/bindings.
##' @noRd
.mutaframe <- function(varlist = list(), row.names = NULL) {
  mf <- new.env(parent = emptyenv())
  
  # Ensure all atomic vectors converted to binding functions
  binders <- as.list(varlist)
  fun <- sapply(varlist, is.function)
  binders[!fun] <- raw_bindings(mf, binders[!fun])
  
  # Activate bindings
  for(name in names(binders)) {
    makeActiveBinding(name, binders[[name]], mf)
  }

  structure(mf,
    col.names = names(varlist), # keep variable order
    row.names = row.names,
    changed = Signal(i, j)$accumulator(combine_data_events),
    class = c("mutaframe", class(mf))
  )
}

