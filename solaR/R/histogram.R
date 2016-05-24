setGeneric('histogram')

setMethod('histogram',
          signature=c(x='formula', data='zoo'),
          definition=function(x, data, par.settings=solaR.theme,
            strip=strip.custom(strip.levels=c(TRUE, TRUE)),breaks=50,...){            
            data0=as.data.frame(data)
            ind=index(data)
            data0$day=doy(ind) ##Incorporo dia, mes y año para facilitar la formula.
            data0$month=month(ind)
            data0$year=year(ind)
            if (!('w' %in% names(data0))){
              data0$w=h2r(hms(ind)-12) ##hora solar en radianes
            }
            histogram(x, data0,
                      scales=list(x=list(relation='free'),
                        y=list(relation='free',
                          draw=FALSE)),
                      breaks=breaks, col='gray',
                      xlab='',
                      strip.names=c(TRUE, TRUE), bg='gray', fg='darkblue',
                      ...)
          }
          )

setMethod('histogram',
          signature=c(x='formula', data='Meteo'),
          definition=function(x, data, ...){
            data0=getData(data)
            histogram(x, data0, ...) ##es un zoo, luego ahora aplica el método data='zoo'
          }
          )

setMethod('histogram',
          signature=c(x='formula', data='Sol'),
          definition=function(x, data, ...){
            data0=as.zooI(data, complete=TRUE, day=TRUE)
            histogram(x, data0, ...)
          }
          )

setMethod('histogram',
          signature=c(x='formula', data='G0'),
          definition=function(x, data, ...){
            data0=as.zooI(data, complete=TRUE, day=TRUE)
            histogram(x, data0, ...)
          }
          )


setMethod('histogram',
          signature=c(x='Meteo', data='missing'),
          definition=function(x, data, par.settings=solaR.theme,
            strip=TRUE, strip.left=FALSE, nbins=50,...){
            x0=getData(x)
            nms <- names(x0)
            N=ncol(x0)
            form <- paste(nms, collapse='+')
            form <- as.formula(paste('~', form, sep=''))
            histogram(form, x0, par.settings=par.settings,
##                   layout=c(N, 1), ##xlab='',
##                   scales=list(cex=0.6, rot= 0),
                   strip=strip, strip.left=strip.left,
                   par.strip.text=list(cex=0.6), nbins=nbins,
                   ...)
          }
          )
