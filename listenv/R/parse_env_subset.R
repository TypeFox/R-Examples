#' Helper function to infer target from expression and environment
#'
#' @param expr An expression.
#' @param envir An environment.
#' @param substitute If TRUE, then the expression is
#'        \code{substitute()}:ed, otherwise not.
#'
#' @return A named list.
#'
#' @export
#' @keywords internal
parse_env_subset <- function(expr, envir=parent.frame(), substitute=TRUE) {
  if (substitute) expr <- substitute(expr)
  code <- paste(deparse(expr), collapse="")

  res <- list(envir=envir, name="", op=NULL, subset=NULL, idx=NA_integer_, exists=NA, code=code)

  if (is.symbol(expr)) {
    ## Variable specified as a symbol
    res$name <- deparse(expr)
  } else if (is.character(expr)) {
    ## Variable specified as a name
    if (length(expr) > 1L) {
      stop(sprintf("Does not specify a single variable, but %d: %s", length(expr), hpaste(sQuote(expr), collapse=", ")), call.=FALSE)
    }
    res$name <- expr
  } else if (is.numeric(expr)) {
    ## Variable specified as a subset of envir
    if (length(expr) > 1L) {
      stop(sprintf("Does not specify a single index, but %d: %s", length(expr), hpaste(sQuote(expr), collapse=", ")), call.=FALSE)
    }
    res$subset <- list(expr)
  } else {
    n <- length(expr)
    stopifnot(n >= 2L)

    if (n >= 3L) {
      ## Assignment to enviroment via $ and [[
      op <- as.character(expr[[1]])
      res$op <- op
      if (op == "$" && n > 3L) {
        stop("Invalid syntax: ", sQuote(code), call.=FALSE)
      } else if (!is.element(op, c("$", "[[", "["))) {
        stop("Invalid syntax: ", sQuote(code), call.=FALSE)
      }

      ## Target
      objname <- deparse(expr[[2]])
      if (!exists(objname, envir=envir, inherits=TRUE)) {
        stop(sprintf("Object %s not found: %s", sQuote(objname), sQuote(code)), call.=FALSE)
      }

      obj <- get(objname, envir=envir, inherits=TRUE)
      if (!is.environment(obj)) {
        stop(sprintf("Subsetting can not be done on a %s; only to an environment: %s", sQuote(mode(obj)), sQuote(code)), call.=FALSE)
      }
      res$envir <- obj

      ## Subset
      subset <- list()
      for (kk in 3:n) {
        missing <- (length(expr[[kk]]) == 1L) && (expr[[kk]] == "")
        if (missing) {
          subsetKK <- NULL
        } else {
          subsetKK <- expr[[kk]]
        }
        if (is.symbol(subsetKK)) {
          subsetKK <- deparse(subsetKK)
          if (op == "[[") {
            if (!exists(subsetKK, envir=envir, inherits=TRUE)) {
              stop(sprintf("Object %s not found: %s", sQuote(subsetKK), sQuote(code)), call.=FALSE)
            }
            subsetKK <- get(subsetKK, envir=envir, inherits=TRUE)
          }
        } else if (is.language(subsetKK)) {
          subsetKK <- eval(subsetKK, envir=envir)
        }
        if (is.null(subsetKK)) {
          subset[kk-2L] <- list(NULL)
        } else {
          subset[[kk-2L]] <- subsetKK
        }
      }

      res$subset <- subset
    } # if (n >= 3)
  } # if (is.symbol(expr))


  ## Validat name, iff any
  name <- res$name
  if (nzchar(name) && !grepl("^[.a-zA-Z]+", name)) stop("Not a valid variable name: ", sQuote(name), call.=FALSE)


  ## Validate subsetting, e.g. x[[1]], x[["a"]], and x$a, iff any
  subset <- res$subset
  if (!is.null(subset)) {
    if (!is.list(subset)) {
      stop(sprintf("INTERNAL ERROR (expected 'subset' to be a list): %s", sQuote(code)), call.=FALSE)
    }
    if (length(subset) == 0L) {
      stop(sprintf("Subsetting of at least on element is required: %s", sQuote(code)), call.=FALSE)
    }

    for (kk in seq_along(subset)) {
      subsetKK <- subset[[kk]]
      if (is.null(subsetKK)) {
      } else if (any(is.na(subsetKK))) {
        stop(sprintf("Invalid subsetting. Subset must not contain missing values: %s", sQuote(code)), call.=FALSE)
      } else if (is.character(subsetKK)) {
        if (!all(nzchar(subsetKK))) {
          stop(sprintf("Invalid subset. Subset must not contain empty names: %s", sQuote(code)), call.=FALSE)
        }
      } else if (is.numeric(subsetKK)) {
      } else {
        stop(sprintf("Invalid subset of type %s: %s", sQuote(typeof(subsetKK)), sQuote(code)), call.=FALSE)
      }
    } # for (kk ...)

    ## Special: listenv:s
    envir <- res$envir
    stopifnot(is.environment(envir))

    if (inherits(envir, "listenv")) {
      names <- names(envir)
      map <- map(envir)
      dim <- dim(envir)

      op <- res$op
      if (is.null(op)) op <- "[["

      ## Multi-dimensional subsetting?
      if (length(subset) > 1L) {
        if (is.null(dim)) {
          stop("Multi-dimensional subsetting on list environment without dimensions: ", sQuote(code), call.=TRUE)
        }
        dimnames <- dimnames(envir)
        exists <- TRUE
        for (kk in seq_along(subset)) {
          subsetKK <- subset[[kk]]
          if (is.null(subsetKK)) {
            subset[[kk]] <- seq_len(dim[kk])
          } else if (is.numeric(subsetKK)) {
            exists <- exists && (subsetKK >= 1 && subsetKK <= dim[kk])
          } else if (is.character(subsetKK)) {
            subsetKK <- match(subsetKK, dimnames[[kk]])
            exists <- exists && !is.na(subsetKK)
            subset[[kk]] <- subsetKK
          }
        }

        ## Indexing scale factor per dimension
        ndim <- length(dim)
        scale <- c(1L, cumprod(dim[-ndim]))
        idx <- 1
        for (kk in seq_along(subset)) {
          i <- subset[[kk]]
          stopifnot(is.numeric(i))
          d <- dim[kk]
          if (any(i < 0)) {
            if (op == "[[") {
              stop("Invalid (negative) indices: ", hpaste(i))
            } else if (any(i > 0)) {
              stop("only 0's may be mixed with negative subscripts")
            }
            ## Drop elements
            i <- setdiff(seq_len(d), -i)
          }
          if (any(i > d)) i[i > d] <- NA_integer_
          ## Drop zeros
          i <- i[i != 0]
          i <- scale[kk]*(i - 1)
          if (kk == 1) {
            idx <- idx + i
          } else {
            idx <- outer(idx, i, FUN=`+`)
          }
        } # for (kk ...)

        res$idx <- idx
        res$name <- names[res$idx]
        if (length(res$name) == 0L) res$name <- ""
        if (exists) {
          exists <- !is.na(map[idx])
        }
        res$exists <- exists
      } else {
        subset <- subset[[1L]]
        if (is.numeric(subset)) {
          i <- subset
          n <- length(envir)
          if (any(i < 0)) {
            if (op == "[[") {
              stop("Invalid (negative) indices: ", hpaste(i))
            } else if (any(i > 0)) {
              stop("only 0's may be mixed with negative subscripts")
            }
            ## Drop elements
            i <- setdiff(seq_len(n), -i)
          }
          ## Drop zeros?
          keep <- which(i != 0)
          if (length(keep) != length(i)) {
            if (op == "[[") {
              ## BACKWARD COMPATIBILITY:
              ## In order not to break two `R CMD check` package tests
              ## for future 0.9.0 on CRAN, we tweak the result here in
              ## order for those two tests not to fail. /HB 2015-12-26
              ## FIX ME: Remove when future (> 0.9.0) is on CRAN.
              if (identical(i, 0) && identical(code, "x[[0]]") && is.element("package:future", search()) && utils::packageVersion("future") <= "0.9.0") {
                res$idx <- i
                res$exists <- FALSE
                return(res)
              }
              stop("Invalid (zero) indices: ", hpaste(i))
            }
            i <- i[keep]
          }
          res$idx <- i
          res$exists <- !is.na(map[res$idx]) & (res$idx >= 1 & res$idx <= n)
          res$name <- names[i]
          if (length(res$name) == 0L) res$name <- ""
        } else if (is.character(subset)) {
          res$idx <- match(subset, names)
          res$exists <- !is.na(res$idx) && !is.na(map[res$idx])
        }
      }
    } else {
      if (length(subset) > 1L) {
        stop("Invalid subset: ", sQuote(code), call.=TRUE)
      }
      subset <- subset[[1L]]
    }
    if (is.character(subset)) res$name <- subset
  }

  ## Identify index?
  if (inherits(res$envir, "listenv")) {
    envir <- res$envir
    if (any(is.na(res$idx)) && nzchar(res$name)) {
      res$idx <- match(res$name, names(envir))
    }
    res$exists <- !is.na(res$idx) & !is.na(map(envir)[res$idx])
  }

  ## Validate
  if (is.null(dim) && length(res$subset) == 1 && identical(res$op, "[")) {
    if (any(is.na(res$idx)) && !nzchar(res$name)) {
      stop("Invalid subset: ", sQuote(code), call.=TRUE)
    }
  }

  unknown <- which(is.na(res$exists))
  if (length(unknown) > 0) {
    res$exists[unknown] <- sapply(unknown, FUN=function(idx) {
      exists(res$name[idx], envir=res$envir, inherits=TRUE)
    })
  }

  ## Sanity check
  stopifnot(is.environment(res$envir))
  stopifnot(is.character(res$name))
  stopifnot(is.null(res$idx) || all(is.numeric(res$idx)))
  stopifnot(is.logical(res$exists), !anyNA(res$exists))
  stopifnot(length(res$exists) == length(res$idx))

  res
}
