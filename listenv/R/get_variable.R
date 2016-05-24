#' Get name of variable for a specific element of list environment
#'
#' @param x A list environment.
#' @param name The name or index of element of interest.
#' @param mustExist If TRUE, an error is generated if \code{name}
#'       does not exist.
#' @param create If TRUE, element \code{name} is created if missing.
#'
#' @return The name of the underlying variable
#'
#' @aliases get_variable.listenv
#' @export
#' @keywords internal
get_variable <- function(...) UseMethod("get_variable")

#' @export
get_variable.listenv <- function(x, name, mustExist=FALSE, create=!mustExist, ...) {
  if (is.character(name)) {
  } else if (is.numeric(name)) {
  } else {
    stop("Subscript must be a name or an index: ", mode(name), call.=FALSE)
  }

  dim <- dim(x)
  if (is.null(dim)) {
    if (length(name) != 1L) {
      stop("Subscript must be a scalar: ", length(name), call.=FALSE)
    }
  } else {
    ndim <- length(dim)
    if (length(name) != 1L && length(name) != ndim) {
      stop(sprintf("Subscript must be a scalar or of equal length as the number of dimension (%d): %d", ndim, length(name)), call.=FALSE)
    }

    ## Map multi-dimensional index to scalar index
    if (length(name) > 1L) {
      stopifnot(is.numeric(name))
      idxs <- name
      if (anyNA(idxs)) stop("Unknown index detected")

      for (kk in seq_len(ndim)) {
        if (idxs[kk] < 1 || idxs[kk] > dim[kk]) {
          stop(sprintf("Index #%d out of range [1,%d]: %s", kk, dim[kk], idxs[kk]))
        }
      }
      bases <- rev(c(cumprod(dim[-ndim]), 1))
      idx <- sum(bases * (idxs-1)) + 1
      name <- idx
    }
  }

  map <- map(x)

  ## Existing variable?
  var <- map[name]
  if (length(var) == 1L && !is.na(var)) return(var)

  if (mustExist) {
    stop(sprintf("No such %s element: %s", sQuote(class(x)[1]), name))
  }

  ## Create new variable
  if (is.character(name)) {
    var <- name

    ## Non-existing name?
    if (!is.element(name, names(map))) {
      map <- c(map, var)
      names(map)[length(map)] <- var
    }
  } else if (is.numeric(name)) {
    i <- name
    ## Expand map?
    if (i > length(map)) {
      extra <- rep(NA_character_, times=i-length(map))
      map <- c(map, extra)
    }
    ## Create internal variable
    var <- new_variable(x, value=NULL, create=create)
    map[i] <- var
  }

  ## Update map?
  if (create) map(x) <- map

  var
}
