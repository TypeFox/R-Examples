## compareFunction: no visible binding for global variable ‘name’
## compareFunction: no visible binding for global variable ‘x’
## compareFunction: no visible binding for global variable ‘y’
## compareFunction: no visible binding for global variable ‘group.value’

if(getRversion() >= "2.15.1") globalVariables(c('name', 'x', 'y', 'group.value'))

setGeneric('compare', signature='...', function(...){standardGeneric('compare')})

compareFunction <- function(..., vars){
  dots <- list(...)
  nms0 <- substitute(list(...))
  if (!is.null(names(nms0))){ ##estamos dentro de do.call
    nms <- names(nms0[-1])
  } else {
    nms <- as.character(nms0[-1])
  }
  foo <- function(object, label){
    yY <- colMeans(as.data.frameY(object, complete=TRUE)[vars])
    yY <- cbind(stack(yY), name=label)
    yY
  }
  cdata <- mapply(FUN=foo, dots, nms, SIMPLIFY=FALSE)
  z <- do.call(rbind, cdata)
  z$ind <- ordered(z$ind, levels=vars)
  p <- dotplot(ind~values, groups=name, data=z, type='b',
               par.settings=solaR.theme)
  print(p+glayer(panel.text(x[length(x)], y[length(x)],
                            label=group.value, cex=0.7, pos=3, srt=45)))
  return(z)
}


setMethod('compare',
          signature='G0',
          definition=function(...){
            vars <- c('D0d', 'B0d', 'G0d')
            res <- compareFunction(..., vars=vars)
            return(res)
          }
          )

setMethod('compare',
          signature='Gef',
          definition=function(...){
            vars <- c('Defd', 'Befd', 'Gefd')
            res <- compareFunction(..., vars=vars)
            return(res)
          }
          )

setMethod('compare',
          signature='ProdGCPV',
          definition=function(...){
            vars <- c('G0d', 'Gefd', 'Yf')
            res <- compareFunction(..., vars=vars)
            return(res)
          }
          )
