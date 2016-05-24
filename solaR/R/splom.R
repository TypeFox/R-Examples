## splom
## splomsolaR <- function(x, ...){
##   splom(x,
##         panel=panel.hexbinplot,
##         diag.panel = function(x, ...){
##           yrng <- current.panel.limits()$ylim
##           d <- density(x, na.rm=TRUE)
##           d$y <- with(d, yrng[1] + 0.95 * diff(yrng) * y / max(y) )
##           panel.lines(d)
##           diag.panel.splom(x,...)
##         },
##         lower.panel = function(x, y, ...){
##           panel.hexbinplot(x, y, ...)
##           panel.loess(x, y, ..., col = 'red')
##         },
##         pscale=0, varname.cex=0.7
##         )
##   }

## setMethod('splom',
##           signature='Meteo',
##           definition=function(x, ...){
##             df <- as.data.frame(getData(x))
##             splomsolaR(df)
##           }
##           )



## setGeneric('compareSplom', signature='...', function(...){standardGeneric('compareSplom')})

## setMethod('compareSplom',
##           signature='Meteo',
##           definition=function(...){
##             z <- mergesolaR(..., var='G0')
##             df <- as.data.frame(z)
##             splomsolaR(df)
##           }
##           )

## setMethod('compareSplom',
##           signature='G0',
##           definition=function(...){
##             z <- mergesolaR(..., var='G0d')
##             df <- as.data.frame(z)
##             splomsolaR(df)
##           }
##           )

## setMethod('compareSplom',
##           signature='Gef',
##           definition=function(...){
##             z <- mergesolaR(..., var='Gefd')
##             df <- as.data.frame(z)
##             splomsolaR(df)
##           }
##           )

## setMethod('compareSplom',
##           signature='ProdGCPV',
##           definition=function(...){
##             z <- mergesolaR(..., var='Yf')
##             df <- as.data.frame(z)
##             splomsolaR(df)
##           }
##           )

