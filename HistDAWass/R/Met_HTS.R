## ---- Show overridding for Time series of TdistributionH and TMatH ----
setMethod("show",
          signature(object="TdistributionH"),
          definition=function(object){
            newdist=distributionH(object@x, object@p)
            show(newdist)
          }
)
setMethod("show",
          signature(object="TMatH"),
          definition=function(object){
            cat("to be implemented")
          }
)
setMethod("show",
          signature(object="HTS"),
          definition=function(object){
            x=summary.HTS(object)
                show(x[,1:7])
          }
)

## --- Plot overloading for Timed TdistributionH  TMatH and HTS----
#' plot for a TdistributionH object
#' @name plot-TdistributionH
#' @docType methods
#' @aliases plot,TdistributionH-method
#' @description A plot function for a \code{TdistributionH} object. The function returns a representation 
#' of the histogram.
#' @param x  a \code{TdistributionH} object
#' @param type (optional) a string describing the type of plot, default="HISTO".\cr Other allowed types are 
#' \cr"CDF"=Cumulative distribution function, \cr"QF"= quantile function, \cr"DENS"=a density approximation, 
#' \cr"HBOXPLOT"=horizontal boxplot, \cr"VBOXPLOT"= vertical boxplot,
#' @param col (optional) a string the color of the plot, default="green".
#' @param border (optional) a string the color of the border of the plot, default="black".
#' @export
setMethod("plot",
          signature(x = "TdistributionH"),
          function (x,  type="HISTO",col="green",border="black") 
          {
            newdist=distributionH(x@x, x@p)
            plot(newdist,  type=type,col=col,border=border)
          }
)
# setMethod("plot",
#           signature(x = "TMatH"),
#           function (x, y="missing", type="HISTO",col="green",border="black") 
#           {
#             stop("to be implemented")            
#           }
# )
## --- Plot overloading
#' Method plot for a histogram time series
#' @name plot-HTS
#' @docType methods
#' @rdname plot-HTS
#' @aliases plot,HTS-method
#' @description An overloading plot function for a \code{HTS} object. The method returns a graphical representation 
#' of a histogram time series. 
#' @param x a \code{distributionH} object
#' @param y not used in this implementation
#' @param type (optional) a string describing the type of plot, default="BOXPLOT".\cr
#'  Other allowed types are \cr
#'  "VIOLIN"=a violin-plot representation, 
#' @param border (optional) a string the color of the border of the plot, default="black".
#' @param maxno.perplot An integer (DEFAULT=30). Maximum number of timestamps per row. 
#'  It allows a plot organized by rows, each row of the plot contains a max number of time stamps
#'   indicated by maxno.perplot.  
#' @examples
#'  plot(subsetHTS(RetHTS,from=1,to=40)) #plots RetHTS dataset
#' \dontrun{
#'  plot(RetHTS, type="BOXPLOT",  border="blue", maxno.perplot=20) 
#'  plot(RetHTS, type="VIOLIN",  border="blue", maxno.perplot=20)
#'  plot(RetHTS, type="VIOLIN",  border="blue", maxno.perplot=10)
#'  }
#' @export
 setMethod("plot",
           signature(x = "HTS"),
           function (x, y="missing", type="VIOLIN",border="black",maxno.perplot=30) 
           {
            plot.HTS.1v(x, type=type, border=border, maxno.perplot=maxno.perplot)
           }
 )

summary.HTS=function (x)
{
  df=data.frame(Tstamp=numeric(), 
                T_period_start=numeric(),
                T_period_end=numeric(),
                mean=numeric(),
                std=numeric(),
                skew=numeric(),
                kurt=numeric(),
                min=numeric(),
                Q1=numeric(),
                MED=numeric(),
                Q3=numeric(),
                Max=numeric())
  nr=length(x@data)
  for (i in 1:nr){
    tmpo=x@data[[i]]
    newrow=data.frame(Tstamp=tmpo@tstamp, 
                      T_period_start=tmpo@period$start,
                      T_period_end=tmpo@period$end,
                      mean=tmpo@m,
                      std=tmpo@s,
                      skew=skewH(tmpo),
                      kurt=kurtH(tmpo),
                      min=tmpo@x[1],
                      Q1=compQ(tmpo,0.25),
                      MED=compQ(tmpo,0.5),
                      Q3=compQ(tmpo,0.75),
                      Max=tmpo@x[length(tmpo@x)])
    df=rbind(df,newrow)
    
  }
  return(df)
}

#' Method \code{subsetHTS}: extract a subset of a histogram time series
#' @name subsetHTS
#' @rdname subsetHTS-methods
#' @exportMethod subsetHTS
setGeneric("subsetHTS", function(object,from,to) standardGeneric("subsetHTS"))
#' @rdname subsetHTS-methods
#' @aliases subsetHTS,HTS,numeric,numeric-method
#' @description This functon return the mean of a \code{distributionH} object.
#' @param object a \code{HTS} object. A histogram 1d time series
#' @param from an integer, the initioal timepont
#' @param to an integer, a final timepoint
#' @return a \code{HTS} object. A histogram 1d time series
#' @examples
#' SUB_RetHTS=subsetHTS(RetHTS,from=1, to=20)# the first 20 elements

setMethod(f="subsetHTS",
          signature=c(object = "HTS",from="numeric",to="numeric"),
          function(object,from,to){
            nel=to-from+1
            SubHTS=new("HTS")
            cc=0
            for (i in from:to){
              cc=cc+1
              SubHTS@data[[cc]]=object@data[[i]]
            }
            return(SubHTS)
          }
)
