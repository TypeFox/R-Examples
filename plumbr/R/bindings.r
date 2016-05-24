#' Generate binding for proxies.
#'
#' @param mf mutaframe to inherit from
#' @param j columns to generate bindings for
proxy_bindings <- function(mf, j = names(mf)) {
  force(mf)
  binder <- function(sym) {
    function(v) {
      if (missing(v)) {
        get(sym, mf)
      } else {
        assign(sym, v, mf)
      } 
    }
  }
  names(j) <- j
  lapply(j, function(name) {
    force(name)
    binder(name)
  })
}

#' @param i rows to filter
#' @return named list of binding functions
#' @noRd
filter_bindings <- function(mf, j = names(mf), i) {
  binder <- function(sym) {
    function(v) {
      if (is.character(i))
        i <- pmatch(i, rownames(mf), duplicates.ok = TRUE)
      xval <- get(sym, mf)
      if (missing(v))
        xval[i]
      else {
        xval[i] <- v
        assign(sym, xval, mf)
      }
    }
  }    
  names(j) <- j
  lapply(j, function(name) {
    force(name)
    binder(name)
  })
}

##' Generate binding for raw values
##' 
##' @param mf mutaframe
##' @param name name
##' @param data vector to store
##' @return named list of binding functions
raw_binding <- function(mf, name, data) {
  force(name)
  force(data)
  function(new) {
    if (missing(new)) {
      data
    } else {
      rows_changed <- which(data != new)
      data <<- new
      notify_listeners(mf, rows_changed, name)
    }
  }
}

##' Generate binding for raw values
##' 
##' @param mf mutaframe
##' @param data list of values
##' @return named list of binding functions
raw_bindings <- function(mf, data) {
  if (length(data) && is.null(names(data)))
    stop("'data' must have names")
  mapply(raw_binding, names(data), data, MoreArgs = list(mf = mf))
}
