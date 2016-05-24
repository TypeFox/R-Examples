equal <- function(x,y) {
  all(x==y)
}

dots <- function(...) {
  eval(substitute(alist(...)))
}

getnames <- function(x, default = character(sign(length(x)))) {
  if(length(x) == 0L) character(0L)
  else if(is.null(names(x))) character(1L)
  else names(x)
}

checktype <- function(x, type, obj) {
  if(!type(x))
    stop(obj," fails type checking with ", deparse(type),
      call. = FALSE)
  invisible(NULL)
}

checkcond <- function(x, cond, obj, envir = parent.frame(), ...) {
  Map(function(name,value) {
    if(nchar(name) == 0L) {
      name <- value
      value <- TRUE
    }
    pcond <- is.logical(value)
    actual <- do.call(name,list(x),envir = envir)
    valid <- equal(actual,value)
    if(!valid) {
      if(pcond) {
        stop(obj, " violates condition [",deparse(name,...),"]",
          call. = FALSE)
      } else {
        stop(obj," [",
          name," = ",deparse(actual,
            width.cutoff = 20L,nlines = 1L, ...),
          "] violates condition [",
          name," = ",deparse(value,
            width.cutoff = 20L,nlines = 1L, ...),"]",call. = FALSE)
      }
    }
  },getnames(cond),cond)
  invisible(NULL)
}

check <- function(x, value, type, ...) {
  cond <- list(...)
  if(!missing(type)) {
    type <- match.fun(type)
    if(!is.null(x))
      checktype(x, type, "symbol")
    checktype(value, type, "value")
  }
  checkcond(value, cond, "value")
  value
}
