##################################################################
## THEMES
##################################################################
xscale.solar <- function(...){ans <- xscale.components.default(...); ans$top=FALSE; ans}
yscale.solar <- function(...){ans <- yscale.components.default(...); ans$right=FALSE; ans}

solaR.theme <- function(pch=19, cex=0.7, region=rev(brewer.pal(9, 'YlOrRd')), ...) {
  theme <- custom.theme.2(pch=pch, cex=cex, region=region, ...)
  theme$strip.background$col='transparent'
  theme$strip.shingle$col='transparent'
  theme$strip.border$col='transparent'
  theme
}

solaR.theme.2 <- function(pch=19, cex=0.7, region=rev(brewer.pal(9, 'YlOrRd')), ...) {
  theme <- custom.theme.2(pch=pch, cex=cex, region=region, ...)
  theme$strip.background$col='lightgray'
  theme$strip.shingle$col='lightgray'
  theme
}

##################################################################
## XYPLOT
##################################################################
setGeneric('xyplot')

setMethod('xyplot',
          signature=c(x='formula', data='zoo'),
          definition=function(x, data,
            par.settings=solaR.theme,
            xscale.components=xscale.solar,
            yscale.components=yscale.solar,
            ...){
            data0=as.data.frame(data)
            ind=index(data)
            data0$day=doy(ind) ##Incorporo dia, mes y año para facilitar la formula.
            data0$month=month(ind)
            data0$year=year(ind)
            if (!('w' %in% names(data0))){
              data0$w=h2r(hms(ind)-12) ##hora solar en radianes
            }
            xyplot(x, data0, par.settings=par.settings,
                   xscale.components=xscale.components,
                   yscale.components=yscale.components,
                   strip=strip.custom(strip.levels=c(TRUE, TRUE)),...)
          }
          )

setMethod('xyplot',
          signature=c(x='formula', data='Meteo'),
          definition=function(x, data, ...){
            data0=getData(data)
            xyplot(x, data0, ...) ##es un zoo, luego ahora aplica el método data='zoo'
          }
          )

setMethod('xyplot',
          signature=c(x='formula', data='Sol'),
          definition=function(x, data, ...){
            data0=as.zooI(data, complete=TRUE, day=TRUE)
            xyplot(x, data0, ...)
          }
          )

setMethod('xyplot',
          signature=c(x='formula', data='G0'),
          definition=function(x, data, ...){
            data0=as.zooI(data, complete=TRUE, day=TRUE)
            xyplot(x, data0, ...)
          }
          )


setMethod('xyplot',
          signature=c(x='Meteo', data='missing'),
          definition=function(x, data,
            par.settings=solaR.theme.2,
            ## xscale.components=xscale.solar,
            ## yscale.components=yscale.solar,
            strip=FALSE, strip.left=TRUE,...){
            x0=getData(x)
            N=ncol(x0)
            xyplot(x0, par.settings=par.settings,
                   ## xscale.components=xscale.components,
                   ## yscale.components=yscale.components,
                   layout=c(1, N),
                   scales=list(cex=0.6, rot= 0),
                   strip=strip, strip.left=TRUE,
                   par.strip.text=list(cex=0.6),
                   ...)
          }
          )

setMethod('xyplot',
          signature=c(x='G0', data='missing'),
          definition=function(x, data,
            par.settings=solaR.theme.2,
            ## xscale.components=xscale.solar,
            ## yscale.components=yscale.solar,
            ...){
            x0=as.zooD(x, complete=FALSE)
            xyplot(x0, par.settings=par.settings,
                   ## xscale.components=xscale.components,
                   ## yscale.components=yscale.components,
                   superpose=TRUE,
                   auto.key=list(space='right'),
                   ylab='Wh/m\u00b2',
                   ...)
          }
          )

setMethod('xyplot',
          signature=c(x='ProdGCPV', data='missing'),
          definition=function(x, data,
            par.settings=solaR.theme.2,
            ## xscale.components=xscale.solar,
            ## yscale.components=yscale.solar,
            ...){
            x0=as.zooD(x, complete=FALSE)
            xyplot(x0, layout=c(1, 3),
                   par.settings=par.settings,
                   ## xscale.components=xscale.components,
                   ## yscale.components=yscale.components,
                   strip=FALSE,
                   strip.left=TRUE,
                   ...)
          }
          )

setMethod('xyplot',
          signature=c(x='ProdPVPS', data='missing'),
          definition=function(x, data,
            par.settings=solaR.theme.2,
            ## xscale.components=xscale.solar,
            ## yscale.components=yscale.solar,
            ...){
            x0=as.zooD(x, complete=FALSE)
            xyplot(x0, layout=c(1, 3),
                   par.settings=par.settings,
                   ## xscale.components=xscale.components,
                   ## yscale.components=yscale.components,
                   strip=FALSE,
                   strip.left=TRUE,
                   ...)
          }
          )
