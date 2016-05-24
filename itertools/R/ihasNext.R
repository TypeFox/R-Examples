# This is based on code that was contributed by Hadley Wickham

hasNext.ihasNext <- function(obj, ...) {
  obj$hasNext()
}

ihasNext <- function(iterable) {
  it <- iter(iterable)

  if (inherits(it, 'ihasNext')) {
    it
  } else {
    cache <- NULL
    hasnext <- NA

    nextEl <- function() {
      if (! hasNx()) {
        stop('StopIteration', call.=FALSE)
      }

      hasnext <<- NA
      cache
    }

    hasNx <- function() {
      if (is.na(hasnext)) {
        tryCatch({
          cache <<- nextElem(it)
          hasnext <<- TRUE
        },
        error=function(e) {
          if (identical(conditionMessage(e), 'StopIteration')) {
            hasnext <<- FALSE
          } else {
            stop(e)
          }
        })
      }

      hasnext
    }

    obj <- list(nextElem=nextEl, hasNext=hasNx)
    class(obj) <- c('ihasNext', 'abstractiter', 'iter')
    obj
  }
}
