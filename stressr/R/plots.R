#' @title Financial stress index component data xyplot
#' @description Provides a convenience function for passing an \code{cfsi} object to \code{xyplot}.
#' @param x an object of class \code{stress} as returned by \code{\link[stressr]{getStressIndex}}.
#' @param ... other parameters passed to \code{\link[lattice]{xyplot}}.
#' @importFrom lattice xyplot
#' @importFrom lattice panel.xyplot
#' @importFrom lattice panel.grid
#' @export
#' @seealso stressLineChart stressAreaChart getStressIndex xyplot.stress
#' @examples
#' \dontrun{
#' ci = getStressIndex()
#' xyplot(ci)
#' }
xyplot.cfsi <- function(x,...) {
  xyplot(x$df,
         main=x$main,
         col=x$colors,
         sub="Source: Federal Reserve Bank of Cleveland",
         ylab=x$ylab,
         xlab=NULL,
         panel=function(x,...) {
           panel.xyplot(x,...)
           panel.grid(-1,0,lty=2,...)
         },
         ...)
}


#' @title Financial stress index component data xyplot
#' @description Provides a convenience function for passing an \code{stress} object to \code{xyplot}.
#' @param x an object of class \code{stress} as returned by \code{\link[stressr]{getStressComponents}} and its many offspring.
#' @param ... other parameters passed to \code{\link[lattice]{xyplot}}.
#' @importFrom lattice xyplot
#' @importFrom lattice panel.xyplot
#' @importFrom lattice panel.abline
#' @export
#' @seealso stressLineChart stressAreaChart getStressComponents xyplot.cfsi
#' @examples
#' \dontrun{
#' require(lattice)
#' fs <- getFundingStress()
#' xyplot(fs)
#' }
xyplot.stress <- function(x,...) {
  xyplot(x$df[,x$columns],
         main=x$main,
         col=x$colors,
         sub="Source: Federal Reserve Bank of Cleveland",
         ylab=x$ylab,
         xlab=NULL,
         panel=function(x,...){
           panel.xyplot(x,...)
           panel.abline(h=0,col="gray",lty=2)
         },
         ...)
}

#' @title Financial stress index data line chart with regions.
#' @description Provides a convenience function for passing a \code{cfsi} object to \code{xyplot} with attributes as presented by the source.
#' @details Provides several assumptions about the display of the \code{cfsi} data to correspond to similar presentations at the Cleveland Fed's data site.
#' @param e an object of class \code{cfsi} as returned by \code{\link[stressr]{getStressIndex}}.
#' @param range a range string as used by \code{xts} to subset time series dates, e.g. "1996/1997".  Defaults to NA for full range. 
#' @param showGradeRegions whether to show the stress grade regions and labels
#' @importFrom lattice xyplot
#' @importFrom lattice panel.xyplot
#' @importFrom lattice panel.rect
#' @importFrom lattice panel.text
#' @export
#' @seealso xyplot.cfsi getStressIndex
#' @examples
#' \dontrun{
#' idx <- getStressIndex()
#' stressIndexChart(idx)
#' }
stressIndexChart <- function(e,range=NA,showGradeRegions=TRUE) {
  stopifnot(class(e) == "cfsi")
  
  # grade thresholds according to CFSI model
  g1 = -0.733
  g2 = 0.544
  g3 = 1.82
  g4 = 2.30 # convenience
  
  # region and annotation colors
  darker="#aec4d1"
  lighter="#bed2df"
  grade="black"
  
  df <- e$df
  if ( ! is.na(range) )
    df <- df[range,]
  
  xyplot(df,
         par.settings = list(superpose.symbol = list(pch=15, 
                                                     col=e$colors, 
                                                     cex=1.2),
                             superpose.line = list(lwd=2, 
                                                   lty=1, 
                                                   col=e$colors)
         ),
         main=e$main,
         sub="Source: Federal Reserve Bank of Cleveland",
         ylab=e$ylab,
         xlab=NULL,
         auto.key=FALSE,
         panel=function(x,y,...){
           
           if ( showGradeRegions ) {
             # plot rectangle thresholds
             xmin=1
             xmax=1e6
             ymin=-20
             ymax=20
             
             panel.rect(xleft=xmin,ybottom=ymin,
                        xright=xmax,ytop=g1,
                        col=darker,border=NA)
             panel.rect(xleft=xmin,ybottom=g1,
                        xright=xmax,ytop=g2,
                        col=lighter,border=NA)
             panel.rect(xleft=xmin,ybottom=g2,
                        xright=xmax,ytop=g3,
                        col=darker,border=NA)
             panel.rect(xleft=xmin,ybottom=g3,
                        xright=xmax,ytop=ymax,
                        col=lighter,border=NA)
           }
           
           panel.xyplot(x,y,...)
           
           if ( showGradeRegions ) {
             position = c(4,1)
             yadj = -0.2
             panel.text(x[1],y=g1+yadj,labels="Low",
                        pos=position,cex=1.2,col=grade)
             panel.text(x[1],y=g2+yadj,labels="Normal",
                        pos=position,cex=1.2,col=grade)
             panel.text(x[1],y=g3+yadj,labels="Moderate",
                        pos=position,cex=1.2,col=grade)
             panel.text(x[1],y=g4+yadj,labels="Significant",
                        pos=position,cex=1.2,col=grade)
           }
         }
  )
}

#' @title Financial stress component data as an unstacked line chart.
#' @description Provides a convenience function for passing a \code{stress} object to \code{xyplot}.
#' @details Provides several assumptions about the display of the \code{stress} data to correspond to similar presentations at the Cleveland Fed's data site.
#' @param e an object of class \code{stress} as returned by \code{\link[stressr]{getStressComponents}} and its many offspring.
#' @param range a range string as used by \code{xts} to subset time series dates, e.g. "1996/1997".  Defaults to NA for full range. 
#' @importFrom lattice xyplot
#' @importFrom lattice panel.xyplot
#' @importFrom lattice panel.grid
#' @export
#' @seealso xyplot.stress stressAreaChart getStressComponents getComponentSummary getEquityMarkets getFundingMarkets getCreditMarkets getForeignExchangeMarkets getRealEstateMarkets getSecuritizationMarkets
#' @examples
#' \dontrun{
#' es <- getEqityStress()
#' stressLineChart(es,"2007/2009)
#' }
stressLineChart <- function(e,range=NA) {
  stopifnot(class(e) == "stress")
  
  df <- e$df[,e$columns]
  
  if ( ! is.na(range) )
    df <- df[range,]
  
  xyplot(df,
         par.settings = list(superpose.symbol = list(pch=15, 
                                                     col=e$colors, 
                                                     cex=1.2),
                             superpose.line = list(lwd=2, 
                                                   lty=1, 
                                                   col=e$colors)),
         main=e$main,
         sub="Source: Federal Reserve Bank of Cleveland",
         ylab=e$ylab,
         xlab=NULL,
         superpose=TRUE,
         scales=list(y=list(limits=c(0,max(df)))),
         auto.key=
           list(columns=min(ncol(df),3),
                text=colnames(df),
                cex=0.8),
         panel=function(x,...){
           panel.xyplot(x,...)
           panel.grid(-1,0,...)
         }
  )
}

#' @title Financial stress component data as a stacked area chart.
#' @description Provides a convenience function for passing a \code{stress} object to \code{xyplot} to render a sand (stacked area) chart.
#' @details Provides several assumptions about the display of the \code{stress} data to correspond to similar presentations at the Cleveland Fed's data site.  To implement the stacked area chart the function first computes the column-wise value accumulations, then passes these values to the \code{latticeExtra} \code{xyarea} polygon rendering tools.  Plots the columns in reverse stacking order to show the desired overlaps.
#' @param e an object of class \code{stress} as returned by \code{\link[stressr]{getStressComponents}} and its many offspring.
#' @param range a range string as used by \code{xts} to subset time series dates, e.g. "1996/1997".  Defaults to NA for full range. 
#' @importFrom lattice xyplot
#' @importFrom latticeExtra panel.xyarea
#' @importFrom lattice panel.grid
#' @export
#' @seealso xyplot.stress stressLineChart getStressComponents getComponentSummary getEquityMarkets getFundingMarkets getCreditMarkets getForeignExchangeMarkets getRealEstateMarkets getSecuritizationMarkets
#' @examples
#' \dontrun{
#' es <- getEquityStress()
#' stressAreaChart(es)
#' }
stressAreaChart <- function(e,range=NA) {
  stopifnot(class(e) == "stress")
  
  # we don't modify the class's list value
  df <- e$df[,e$columns]
  
  if ( ! is.na(range) )
    df <- df[range,]
  
  # convert to cumulative across columns for sand chart
  if ( ncol(df) > 1 ) {
    df[is.na(df)] <- 0
    for (i in 2:ncol(df)) {
      df[,i] <- df[,i]+df[,i-1]
    }
  }
  
  # use area polygons to implement the stacked area chart
  # plot in reverse column order to show polygon layers
  xyplot(df[,ncol(df):1],
         par.settings = list(superpose.symbol = list(pch=15, 
                                                     col=rev(e$colors), 
                                                     cex=1.2),
                             superpose.polygon = list(lwd=1, 
                                                      lty=1, 
                                                      col=rev(e$colors))),
         main=paste(e$main,ifelse(ncol(df)>1,"(Stacked)","")),
         sub="Source: Federal Reserve Bank of Cleveland",
         ylab=e$ylab,
         xlab=NULL,
         superpose=ifelse(ncol(df)>1,TRUE,FALSE),
         scales=list(y=list(limits=c(0,max(df)))),
         auto.key=
           list(points=TRUE,
                lines=FALSE,
                columns=min(ncol(df),3),
                text=rev(colnames(df)),
                cex=0.8),
         panel=function(x,...) {
           panel.xyarea(x,origin=0,...)
           panel.grid(-1,0,...)
         }
  )
}



