########################################################
## FUNCTIONS THAT CREATE AND INTERACT WITH THE STdata ##
########################################################
##Functions in this file:
## createSTdata          EX:ok
## print.STdata          EX:NO
## summary.STdata        EX:NO
## print.summarySTdata   EX:NO
## plot.STdata           EX:ok
## qqnorm.STdata         EX:ok
## scatterPlot.STdata    EX:ok
## scatterPlot           EX:NO, S3 compliance methods

##' Creates a \code{STdata} object that can be used as input for
##' \code{\link{createSTmodel}}. Names and dates are derived from the input data,
##' either using predefined fields or \code{rownames} / \code{colnames}; for
##' details see the sub-functions linked under the relevant Arguments.
##' 
##' @title Construct STdata Object
##' 
##' @param obs Either a data.frame with fields \code{obs}, \code{date},
##'   \code{ID} giving obsevations, time-points and location names; or a matrix,
##'   e.g. output from \code{\link{createDataMatrix}}.
##' @param covars matrix/data.frame of covariates; should include both
##'   geographic covariates and coordinates of all locations, see
##'   \code{\link{stCheckCovars}}.
##' @param SpatioTemporal possible spatio-temporal covariate, see
##'   \code{\link{stCheckSTcovars}}.
##' @param transform.obs function to apply to the observations,
##'   defaults to an identity transform. Possible options are
##'   \code{\link[base:log]{log}}, \code{\link[base:sqrt]{sqrt}}, and
##'   \code{\link[base:exp]{exp}}.
##' @param mean.0.ST Call \code{\link{removeSTcovarMean}} to produce a mean-zero
##'   spatio-temporal covariate?
##' @param n.basis Number of temporal components in the smooth trends computed
##'   by \code{\link{updateTrend.STdata}}, if \code{NULL} no trend is computed 
##'   (implies only a constant).
##' @param extra.dates Additional dates for which smooth trends should be
##'   computed, used by \code{\link{updateTrend.STdata}}. If \code{n.basis=NULL}
##'   this will force n.basis=0; since the dates are stored in the trend..
##' @param ... Additional parameters passed to \code{\link{updateTrend.STdata}}.
##' @param detrend Use \code{\link{detrendSTdata}} to remove a termporal trend
##'   from the observations; requires \code{n.basis!=NULL}.
##' @param region,method Additional parameters passed to
##'   \code{\link{detrendSTdata}}.
##' @return A \code{STdata} object with, some or all of, the following elements:
##'   \item{covars}{Geographic covariates, locations and names
##'     of the observation locations (the later in \code{covars$ID}),
##'     \code{\link{createSTmodel}} will extract covariates (land use regressors),
##'     observations locations, etc from this data.frame when constructing the
##'     model specification.}
##'   \item{trend}{The temporal trends with \emph{one of the} columns
##'     being named \code{date}, preferably of class \code{\link[base:Date]{Date}}
##'     providing the time alignment for the temporal trends.}
##'   \item{obs}{A data.frame with columns:
##'     \describe{
##'       \item{obs}{The value of each observation.}
##'       \item{date}{The observations time, preferably of class
##'                   \code{\link[base:Date]{Date}}.}
##'       \item{ID}{A \code{character}-class giving observation locations;
##'                 should match elements in \code{locations$ID}.}
##'     }
##'   }
##'   \item{SpatioTemporal}{A 3D-array of spatio-temporal covariates, or \code{NULL}
##'                     if no covariates exist. The array should be
##'                     (number of timepoints) - by - (number of locations) -
##'                     by - (number of covariates) and provide spatio-temporal
##'                     covariates for \emph{all} space-time locations, even
##'                     unobserved ones (needed for prediction).
##'                     The \code{rownames} of the array should represent dates/times
##'                     and \code{colnames} should match the observation location
##'                     names in \code{covars$ID}.}
##'   \item{old.trend,fit.trend}{Additional components added if the observations
##'                              have been detrended, see
##'                              \code{\link{detrendSTdata}}.}
##' 
##' @example Rd_examples/Ex_createSTdata.R
##' 
##' @author Johan Lindström and Assaf P. Oron
##' @family STdata methods
##' @family STdata functions
##' @export
createSTdata <- function(obs, covars, SpatioTemporal=NULL,
                         transform.obs=function(x){return(x)}, mean.0.ST=FALSE,
                         n.basis=0, extra.dates=NULL, ...,
                         detrend=FALSE, region=NULL, method=NULL){

  ##n.basis should default to 0 if not given
  if( is.null(n.basis) ){
    n.basis <- 0
  }
  ##Check if observations are given as a matrix or as a vector
  if( length(obs)==0 ){
    ##empty observation vector
    obs <- data.frame(obs=double(0), date=integer(0), ID=character(0),
                      stringsAsFactors=FALSE)
    if( missing(extra.dates) || is.null(extra.dates) ){
      stop("If 'obs' is empty, 'extra.dates' is needed to define time-points")
    }
  }else if( is.data.frame(obs) ){
    if( dim(obs)[2]<3 ){
      stop("'obs' needs at least 3 columns.")
    }
    ##try to locate "obs", "date" and "ID" in data.frame obs.
    i.obs <- which("obs"==names(obs))
    i.date <- which("date"==names(obs))
    i.ID <- which("ID"==names(obs))
    if( length(i.obs)==0 ){
      warning("Unable to find column 'obs$obs', using 'obs[,1]")
      i.obs <- 1
    }
    if( length(i.date)==0 ){
      warning("Unable to find column 'obs$date', using 'obs[,2]")
      i.date <- 2
    }
    if( length(i.ID)==0 ){
      warning("Unable to find column 'obs$ID', using 'obs[,3]")
      i.ID <- 3
    }
    obs <- obs[,c(i.obs[1], i.date[1], i.ID[1])]
    ##Force ID to character
    if( is.factor(obs$ID) ){
      obs$ID <- as.character(obs$ID)
    }
    ##add names to the observation matrix.
    names(obs) <- c("obs","date","ID")
  }else if( is.matrix(obs) ){
    ID <- colnames(obs)
    if( is.null(ID) ){
      warning("No colnames(obs), using ID <- as.character(1:dim(obs)[2])")
      ID <- as.character(1:dim(obs)[2])
    }
    date <- convertCharToDate( rownames(obs) )
    if( is.null(date) ){
      warning("Unable to coerce rownames(obs) to Date/double, using date <- 1:dim(obs)[1]")
      date <- 1:dim(obs)[1]
    }
    obs <- data.frame(obs=c(obs), date=rep(date, dim(obs)[2]),
                      ID=rep(ID, each=dim(obs)[1]), stringsAsFactors=FALSE)
  }else{
    stop("'obs' has to be either a data.frame or a matrix.")
  }
  ##drop missing observations
  obs <- obs[!is.na(obs$obs),]
  ##and rownames
  rownames(obs) <- NULL
  ##transform observations
  if( is.function(transform.obs) ){
    obs$obs <- transform.obs(obs$obs)
  }else{
    stop("'transform.obs' should be a function")
  }
  ##check that the resulting obs data.frame is valid
  stCheckObs(obs)

  ##extract unique dates and locations from the observation vector
  ID.unique <- unique(obs$ID)
  date.unique <- unique(obs$date)

  ##check LUR-covariates
  covars <- stCheckCovars(covars, ID.unique)

  ##check SpatioTemporal covariates
  SpatioTemporal <- stCheckSTcovars(SpatioTemporal, unique(c(ID.unique, covars$ID)),
                                    date.unique)

  ##create output object
  out <- list(obs=obs, covars=covars, SpatioTemporal=SpatioTemporal)
  ##set class of return object
  class(out) <- "STdata"

  ##now attempt to modify output object
  ##compute smooth trends?
  out <- updateTrend(out, n.basis=n.basis, extra.dates=extra.dates, ...)

  ##detrend data? (will throw an internal warning if n.basis=NULL since this
  ##implies no trend)
  if( detrend ){
    out <- detrendSTdata(out, region, method)
  }
  ##mean seperable STcovariate?
  if( mean.0.ST ){
    out <- removeSTcovarMean(out)
  }
  
  return(out)
}##function createSTdata

###########################
## S3-METHODS FOR STdata ##
###########################

##' \code{\link[base:print]{print}} method for class \code{STdata}.
##'
##' @title Print details for \code{STdata} object
##' @param x \code{STdata} object to print information for.
##' @param type Factorial of \code{length(x$covars$ID)}, if not \code{NULL}
##'   the output also presents summaries of number of sites and observations
##'   as well as time periods per type of site.
##' @param ... Ignored additional arguments.
##' @return Nothing
##' 
##' @author Johan Lindström
##' 
##' @family STdata methods
##' @method print STdata
##' @export
print.STdata <- function(x, type=x$covars$type, ...){
  ##check class belonging
  stCheckClass(x, "STdata", name="x")
  
  ##general information regarding number of observations and trends
  commonPrintST(x, "STdata", 1)

  ##covariates
  cat( sprintf("%d covariate(s):\n", dim(x$covars)[2]))
  print( names(x$covars) )
  cat("\n")
  ##spatio-temporal covariates
  if( is.null(x$SpatioTemporal) ){
    cat("No spatio-temporal covariates.\n")
  }else{
    cat( sprintf("%d spatio-temporal covariate(s):\n",
                 dim(x$SpatioTemporal)[3]))
    print( dimnames(x$SpatioTemporal)[[3]] )
  }
  cat("\n")

  ##group observations by type.
  commonPrintST(x, "STdata", 2, type)

  return(invisible())
}##function print.STdata

##' \code{\link[base:summary]{summary}} method for class \code{STdata}.
##'
##' @title Computes summary details for \code{STdata} object
##' @param object \code{STdata} object to compute summary information for.
##' @param type Factorial of \code{length(x$covars$ID)}, if not \code{NULL}
##'   summaries for the observations are computed per type of site.
##' @param ... Ignored additional arguments. 
##' @return A \code{summary.STdata} object.
##' 
##' @author Johan Lindström
##' 
##' @family STdata methods
##' @method summary STdata
##' @export
summary.STdata <- function(object, type=object$covars$type, ...){
  ##check class belonging
  stCheckClass(object, "STdata", name="object")
  
  ##allocate output object
  out <- list()
  class(out) <- "summary.STdata"
  ##covariates
  out$covars <- summary(object$covars)
  ##common summary computations observations and trend
  out <- c(out, commonSummaryST(object, type))
  ##ST-covariates
  if( !is.null(object$SpatioTemporal) ){
    out$SpatioTemporal <- summary( matrix(object$SpatioTemporal,
                                          prod(dim(object$SpatioTemporal)[1:2]),
                                          dim(object$SpatioTemporal)[3]) )
    colnames(out$SpatioTemporal) <- dimnames(object$SpatioTemporal)[[3]]
  }

  ##return the object
  class(out) <- "summary.STdata"
  return(out)
}##function summary.STdata

##' \code{\link[base:print]{print}} method for class \code{summary.STdata}.
##'
##' @title Print details for \code{summary.STdata} object
##' @param x \code{summary.STdata} object to print information for.
##' @param ... Ignored additional arguments.
##' @return Nothing
##'
##' @author Johan Lindström
##' 
##' @family STdata methods
##' @method print summary.STdata
##' @export
print.summary.STdata <- function(x, ...){
  ##check class belonging
  stCheckClass(x, "summary.STdata", name="x")

  ##print data
  if( is.null(x[["obs"]]) ){
    cat("No observations.\n");
  }else{
    cat("Summary of observations:\n");
    print(x$obs)
  }
  cat("\nSummary of covariates:\n");
  print(x$covars)
  if( is.null(x$trend) ){
    cat("\nNo trend specified\n")
  }else{
    cat("\nSummary of smooth trends:\n");
    print(x$trend)
  }
  if( is.null(x$SpatioTemporal) ){
    cat("\nNo spatio-temporal covariates.\n")
  }else{
    cat("\nSummary of spatio-temporal covariates.\n")
    print(x$SpatioTemporal)
  }
  if( !is.null(x$obs.by.type) ){
    for(i in 1:length(x$obs.by.type)){
      if( is.null(x$obs.by.type[[i]]) ){
        cat( paste("\nNo observations of type", names(x$obs.by.type)[i],"\n") )
      }else{
        cat( paste("\nSummary for observations of type",
                   names(x$obs.by.type)[i],"\n") )
        print(x$obs.by.type[[i]])
      }
    }
  }##if( !is.null(x$obs.by.type) )
  
  return(invisible())
}##function print.summary.STdata


##' \code{\link[graphics:plot]{plot}} method for class \code{STdata} or
##' \code{STmodel}. Provides several different plots of the data.
##' When calles for \code{STmodel}, \code{STmodel$locations} acts as
##' \code{STdata$covars}.
##'
##' Performs a variety of different plots determined by \code{y}:
##' \describe{
##'   \item{"obs"}{Plot observations for location \code{ID}, along with the
##'     fitted temporal trend.}
##'   \item{"res"}{Plot residuals for the fitted temporal trend at location
##'     \code{ID}; adds the \code{y=0} line for reference.}
##'   \item{"acf"}{Plot autocorrelation function for the residuals from the
##'     fitted temporal trend at location \code{ID}.}
##'   \item{"pacf"}{Plot partial autocorrelation function for the residuals
##'     from the fitted temporal trend at location \code{ID}.}
##'   \item{"loc"}{Plot the observation location index number as a function
##'     of the observation date, for all observations. Possibly coded by the
##'     \code{type} of observations locations.}
##'   \item{"loc.obs"}{Plot the observation value as a function
##'     of the observation date, for all observations. Possibly coded by the
##'     \code{type} of observations locations.}
##' }
##'
##' For \code{y=c("obs","res")} the first element of \code{col,pch,cex,lty} is
##' used to specify plotting of the observations, and the second element is used
##' to  specify plotting of the fitted temporal trend, or 0-line for
##' \code{"res"}. Defaults: \code{col=1}, \code{pch=c(1,NA)}, \code{cex=1},
##' \code{lty=c(NA,1)}. Elements of length one are repeated.
##'
##' For \code{y=c("acf","pacf")} \code{col,pch,cex,lty} are ignored.
##'
##' For \code{y=c("loc","loc.obs")} \code{col,pch,cex} are used to specify the
##' points for each of the different levels in \code{type} and should be of
##' length 1 or \code{length(levels(type))}. \code{lty} is ignored.
##' Default: \code{col=1:length(levels(type))}, \code{pch=19}, \code{cex=.1}
##' 
##' For \code{y=c("loc","loc.obs")} a legend is added if \code{legend.loc!=NULL}.
##' The vector \code{legend.names} should have length equal to the number of
##' unique location types.  The default legend names are \code{levels(type)}.
##' 
##' @title Different Plots for \code{STdata}/\code{STmodel} object
##' @param x \code{STdata}/\code{STmodel} object to plot.
##' @param y Type of plot, options are \code{"obs"}, \code{"res"}, \code{"acf"},
##'   \code{"pacf"}, \code{"loc"}, or \code{"loc.obs"}, see details below.
##' @param ID The location for which we want to plot observations. Either a
##'   string matching the names in \code{x$covars$ID} or an integer; if
##'   an integer the functions will plot data from
##'   \code{ID=x$covars$ID[ID]}.
##' @param type Factorial of \code{length(x$covars$type)}, used by
##'   \code{"loc"} and \code{"loc.obs"} to determine how many groups should be
##'   plotted and colour/type coded.
##' @param col,pch,cex,lty Colour, type of points, size of points, and type of
##'   lines. Exact meaning depends on value of \code{y}, see Details.
##' @param legend.loc The location of the legend, for \code{"loc"} and
##'   \code{"loc.obs"}. See \code{\link[graphics:legend]{legend}}.
##' @param legend.names A vector of character strings to be used in the legend,
##'   for \code{"loc"} and for \code{"loc.obs"}
##' @param add Add to existing plot, only relevant if \code{y} is
##'   \code{"obs"}, \code{"res"}, \code{"loc"}, or \code{"loc.obs"}.
##' @param ... Additional parameters passed to \code{\link[graphics:plot]{plot}}
##'   or \code{\link[stats:plot.acf]{plot.acf}}.
##' @return Nothing
##'
##' @example Rd_examples/Ex_plot_STdata.R
##' 
##' @author Johan Lindström and Assaf P. Oron
##' 
##' @family STdata methods
##' @method plot STdata
##' @export
plot.STdata <- function(x, y=c("obs", "res", "acf", "pacf", "loc", "loc.obs"),
                        ID=x$covars$ID[1], type=x$covars$type,
                        col=NULL, pch=NULL, cex=NULL, lty=NULL,
                        legend.loc="topleft", legend.names=NULL,
                        add=FALSE, ...){
  ##check class belonging
  stCheckClass(x, "STdata", name="x")

  ##we have to use y, cast to resonable name
  plot.type <- match.arg(y)
  
  ##defualt values of col,pch,cex,lty
  if( (plot.type %in% c("obs", "res")) ){
    if( missing(col) || is.null(col) ){ col <- 1 }
    if( missing(pch) || is.null(pch) ){ pch <- c(1,NA) }
    if( missing(cex) || is.null(cex) ){ cex <- 1 }
    if( missing(lty) || is.null(lty) ){ lty <- c(NA,1) }
    ##expand to length=2
    if( length(col)==1 ){ col <- rep(col,2) }
    if( length(pch)==1 ){ pch <- rep(pch,2) }
    if( length(cex)==1 ){ cex <- rep(cex,2) }
    if( length(lty)==1 ){ lty <- rep(lty,2) }
    
  }else if( (plot.type %in% c("loc", "loc.obs")) ){
    if( missing(pch) || is.null(pch) ){ pch <- 19 }
    if( missing(cex) || is.null(cex) ){ cex <- 0.1 }
    ##lty ignored, and defaults for col=NULL handled below
  }
  
  if( (plot.type %in% c("obs", "res", "acf", "pacf")) ){
    ##previous plotMesaData function
    if( length(ID)!=1 ){
      stop("length(ID)!=1")
    }
    ##ID not character, assume it's a number and pick this component
    if( !is.character(ID) ){
      ID <- as.character(x$covars$ID)[ID]
    }
    ##pick out the site to use
    IND <- x$obs$ID==ID
    if( sum(IND)==0 )
      stop( paste("No observations for ID =", ID) )
    ##extract times and observations
    date <- x$obs$date[IND]
    y <- x$obs$obs[IND]
    ##extract trend
    if( is.null(x$trend) ){
      warning("No trend provided, assuming constant temporal trend.")
      trend <- NULL
      date.trend <- sort(unique(x$obs$date))
    }else{
      trend <- x$trend[ match(date, x$trend$date),,drop=FALSE]
      ##do we have basis vectors or only the dates -> i.e. constant trend
      if(dim(trend)[2]==1){
        trend <- NULL
      }else{
        trend <- trend[, -which(names(trend)=="date"),drop=FALSE]
      }
      date.trend <- x$trend$date
    }##if( is.null(x$trend) ){...}else{...}
    N.trend <- length(date.trend)
    
    ##linear regression
    if( is.null(trend) ){
      y.fit <- lm(y~1)
    }else{
      y.fit <- lm(paste("y~", paste(colnames(trend), collapse="+")),
                  data=cbind(y,trend))
    }
    if( plot.type!="obs" ){
      ##if not observations we only need the residuals
      y <- y.fit$residuals
      y.p <- rep(0, length(date.trend))
    }else{
      ##prediction of y for all times (works also for empty trend)
      suppressWarnings( y.p <- predict(y.fit, x$trend) )
    }
    if( y.fit$df.residual==0 ){
      if(plot.type %in% c("acf","pacf")){
        warning("Site ",ID," has too few observations, acf/pacf may FAIL.")
      }else{
        warning("Site ",ID," has too few observations, fitted smooths may be inacurate.")
      }
    }
    
    if(plot.type %in% c("acf","pacf")){
      ##plot correlation functions
      ##expand the dates to equidistant locations
      date.exp <- seq(min(date), max(date), min(diff(date.trend)))
      y.exp <- rep(NA, length(date.exp))
      ##and match the dates suitably.
      y.exp[match(date,date.exp)] <- y
      ##construct list of arguments
      args <- internalPlotFixArgs(list(...),
                                  default=list(main=switch(plot.type,
                                                 "acf"=paste("ACF for",ID),
                                                 "pacf"=paste("PACF for",ID))),
                                  add=list(x=y.exp, na.action=na.pass))
      if(plot.type == "acf"){
        do.call(acf, args)
      }else{
        do.call(pacf, args)
      }
    }else{
      ##plot time series or residuals
      if(!add){
        args <- internalPlotFixArgs(list(...),
                                    default=list(xlab="date", ylab="y",
                                      main=switch(plot.type, "obs"=ID,
                                        "res"=paste("Residuals",ID))),
                                    add=list(x=date, y=y, xlim =
                                      range(c(date,x$trend$date),na.rm=TRUE),
                                      ylim = range(c(y,y.p),na.rm=TRUE),
                                      type="n"))
        do.call(plot, args)
      }
      ##add points and lines
      if( !is.na(pch[1]) ){
        points(date, y, pch=pch[1], cex=cex[1], col=col[1])
      }
      if( !is.na(lty[1]) ){ lines(date, y, lty=lty[1], col=col[1]) }
      ##plot observations or residuals (different fitted lines)
      if( !is.na(pch[2]) ){
        points(date.trend, y.p, pch=pch[2], col=col[2], cex=cex[2])
      }
      if( !is.na(lty[2]) ){ lines(date.trend, y.p, lty=lty[2], col=col[2]) }
    }##if(plot.type %in% c("acf","pacf")){...}else{...}
  }else{
    ##previous plotMonitoringLoc function
    if( is.null(type) ){
      ##add a constant type argument to enable plotting
      type <- factor( rep("Observations", dim(x$covars)[1]) )
    }
    if( !is.factor(type) ){
      warning("Attempting to coerc 'type' to factor")
      type <- as.factor(type)
    }
    ##find matching locations
    idx <- match( x$obs$ID, x$covars$ID)
    
    ##set up colours
    if( is.null(col) ){ col <- 1:length(levels(type)) }
    if( length(col)==1 ){ col <- rep(col, length(levels(type))) }
    ##point-types
    if( length(pch)==1 ){ pch <- rep(pch, length(levels(type))) }
    ##...and cex-types
    if( length(cex)==1 ){ cex <- rep(cex, length(levels(type))) }
    ##set up legend
    if( !is.null(legend.loc) ){
      if( is.null(legend.names) ){
        legend.names <- as.character(levels(type))
      }
      if( length(legend.names) != length(levels(type)) ){
        stop( sprintf("Needs one 'legend.names' for every one of %d types",
                      length(levels(type))) )
      }
    }##if( !is.null(legend.loc) ){
    ##start plotting - select y-values
    if( plot.type=="loc.obs" ){
      ##Plot observations rather than location ID:s
      y.vals <- x$obs$obs
      y.text <- "Observations"
    }else{
      ##Plot location ID:s
      y.vals <- idx
      y.text <- "Location ID"
    }
    if(!add){
      args <- internalPlotFixArgs(list(...),
                                  default=list(xlab="Date", ylab=y.text),
                                  add=list(x=x$obs$date, y=y.vals, type="n"))
      do.call(plot, args)
    }
    ##loop over all possible types
    j=1
    for(i in levels(type) ){
      Ind <- (idx %in% which(type==i))
      if( sum(Ind)!=0 ){
        points(x$obs$date[Ind], y.vals[Ind], col=col[j], pch=pch[j], cex=cex[j])
      }
      j=j+1
    }
    ##possibly add a legend
    if( !is.null(legend.loc) ){
      legend(legend.loc, legend.names, pch=pch, col=col, pt.cex=cex)
    }
  }##if( (plot.type %in% c("obs", "res", "acf", "pacf")) ){...}else{...}
  
  return(invisible())
}##function plot.STdata

##' \code{\link[stats:qqnorm]{qqnorm}} method for classes
##' \code{STdata}/\code{STmodel}/\code{predCVSTmodel}. 
##' Used for data and residual analysis of the cross validation.
##' 
##' @title QQ-norm for \code{STdata}/\code{STmodel}/\code{predCVSTmodel} objects
##' @param y \code{STdata}/\code{STmodel}/\code{predCVSTmodel} object for the
##'   qqnorm.
##' @param ID The location for which we want to norm-plot observations/residuals
##'   or \code{"all"} to plot for all locations.
##' @param main Title of the plot
##' @param group Do the norm-plot both for all data and then for each subset
##'   defined by the factor/levels in group variable.
##' @param col Colour of points in the plot, either a scalar or a vector
##'   with length matching the number of observations/residuals.
##' @param line If non-zero add a \code{\link[stats:qqline]{qqline}} with
##'   \code{lty=line}, to the plot; if 0 \emph{do not} add a line.
##' @param ...  Arguments passed on to the plotting function,
##'   \code{\link[stats:qqnorm]{qqnorm}}. 
##' @return Nothing
##'
##' @author Johan Lindström
##' 
##' @example Rd_examples/Ex_qqnorm_STdata.R
##' 
##' @family STdata methods
##' @importFrom stats qqnorm
##' @method qqnorm STdata
##' @export
qqnorm.STdata <- function(y, ID="all", main="Q-Q plot for observations",
                          group=NULL, col=1, line=0, ...){ 
  ##check class belonging
  stCheckClass(y, "STdata", name="y")

  Y <- y$obs[, c("ID","obs")]

  internalQQnormPlot(Y, ID, main, group, col, FALSE, line, ...)
  
  return(invisible())
}##function qqnorm.STdata

##' Does a scatterPlot of observations/residuals against covariates (either
##' geographic or temporal trends), adding a spline fit (similar to
##' \code{\link[stats:scatter.smooth]{scatter.smooth}}.
##'
##' @title scatterPlot for \code{STdata}/\code{STmodel}/\code{predCVSTmodel} objects
##' 
##' @param x \code{STdata}/\code{STmodel}/\code{predCVSTmodel} object to plot.
##' @param covar,trend Plot observations as a function of? Only \emph{one} of
##'   these should be not \code{NULL}. \code{covar} uses location covariates,
##'   and \code{trend} uses temporal trend (or dates); \code{trend=0} uses a
##'   temporal intercept (i.e. a constant).
##' @param pch,cex Point and point size for the plot, a single value or
##'   \code{nlevels(group)} 
##' @param col,lty Color of points and smooth lines. A single value or
##'  \code{nlevels(group)+1} values; the last value is used for fitting a line
##'  to \emph{all} data. Use \code{lty=NA} to supress smooth lines.
##' @param subset A subset of locations for which to plot observations as a
##'   function of covariates.
##' @param group A vector of factors of the same length as the number of
##'   observations (typically \code{length(x$obs$obs)}, or
##'   \code{length(x$pred.obs$obs)}) used to group data and fit different
##'   smooths to each group.
##' @param add Add to existing plot
##' @param smooth.args List of arguments for
##'   \code{\link[stats:loess.smooth]{loess.smooth}}.
##' @param ... Additional parameters passed to \code{\link[graphics:plot]{plot}}.
##' @return Nothing
##'
##' @example Rd_examples/Ex_scatterPlot_STdata.R
##' 
##' @author Johan Lindström
##' @family STdata methods
##' @method scatterPlot STdata
##' @export
scatterPlot.STdata <- function(x, covar=NULL, trend=NULL, pch=1, col=1, cex=1,
                               lty=1, subset=NULL, group=NULL, add=FALSE,
                               smooth.args=NULL, ...){
  ##check class belonging
  stCheckClass(x, "STdata", name="x")

  ##add trend to avoid future problems
  if( is.null(x$trend) ){
    x <- updateTrend(x, n.basis=0)
  }
  
  ##pass data to internalScatterPlot function
  internalScatterPlot(obs=x$obs[, c("obs","ID","date")],
                      covar=covar, trend=trend, subset=subset,
                      data=list(covars=x$covars, trend=x$trend),
                      group=group, pch=pch, col=col, cex=cex, lty=lty,
                      add=add, smooth.args=smooth.args, ...) 
}##function scatterPlot.STdata

################################
## S3 methods for scatterPlot ##
################################

##' Scatter plot of data in x
##'
##' @title Scatter plot
##' @param x object to scatter plot
##' @param ... additional parameters
##' @return Nothing
##' 
##' @author Johan Lindström
##' @export
scatterPlot <- function(x, ...){
  UseMethod("scatterPlot")
}
