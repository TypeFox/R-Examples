setGeneric('mergesolaR', signature='...', function(...){standardGeneric('mergesolaR')})

fooMeteo <- function(object, var){yY <- getData(object)[,var]}

fooG0 <- function(object, var){yY <- as.zooD(object)[,var]}

mergeFunction <- function(..., foo, var){
  dots <- list(...)
  dots <- lapply(dots, as, class(dots[[1]])) ##el primer elemento es el que dicta la clase a todos
  nms0 <- substitute(list(...))
  if (!is.null(names(nms0))){ ##estamos dentro de do.call
    nms <- names(nms0[-1])
  } else { ##llamada convencional
    nms <- as.character(nms0[-1])
  }
  cdata <- sapply(dots, FUN=foo, var, simplify=FALSE)
  names(cdata) <- nms
  z <- do.call(merge, cdata)
  z
}

setMethod('mergesolaR',
          signature='Meteo',
          definition=function(...){
            res <- mergeFunction(..., foo=fooMeteo, var='G0')
            res
          }
          )

setMethod('mergesolaR',
          signature='G0',
          definition=function(...){
            res <- mergeFunction(..., foo=fooG0, var='G0d')
            res
          }
          )

setMethod('mergesolaR',
          signature='Gef',
          definition=function(...){
            res <- mergeFunction(..., foo=fooG0, var='Gefd')
            res
          }
          )

setMethod('mergesolaR',
          signature='ProdGCPV',
          definition=function(...){
            res <- mergeFunction(..., foo=fooG0, var='Yf')
            res
          }
          )

setMethod('mergesolaR',
          signature='ProdPVPS',
          definition=function(...){
            res <- mergeFunction(..., foo=fooG0, var='Yf')
            res
          }
          )
