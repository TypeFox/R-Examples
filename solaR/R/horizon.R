## setGeneric('horizonsolaR', signature='...', function(...){standardGeneric('horizonsolaR')})

## setMethod('horizonsolaR',
##           signature='Meteo',
##           definition=function(...){
##             z <- mergesolaR(..., var='G0')
##             z <- z-rowMeans(z, na.rm=1)
##             horizonplot(z, colorkey=TRUE)
##           }
##           )

## setMethod('horizonsolaR',
##           signature='G0',
##           definition=function(...){
##             z <- mergesolaR(..., var='G0d')
##             z <- z-rowMeans(z, na.rm=1)
##             horizonplot(z, colorkey=TRUE)
##           })


## setMethod('horizonsolaR',
##           signature='Gef',
##           definition=function(...){
##             z <- mergesolaR(..., var='Gefd')
##             z <- z-rowMeans(z, na.rm=1)
##             horizonplot(z, colorkey=TRUE)
##           })

## setMethod('horizonsolaR',
##           signature='ProdGCPV',
##           definition=function(...){
##             z <- mergesolaR(..., var='Yf')
##             z <- z-rowMeans(z, na.rm=1)
##             horizonplot(z, colorkey=TRUE)
##           })
