#########################################################
## FUNCTIONS THAT CREATE AND INTERACT WITH THE STmodel ##
#########################################################
##Functions in this file:
## createSTmodel         EX:ok
## print.STmodel         EX:ok
## summary.STmodel       EX:ok
## print.summary.STmodel EX:implicit in summary.STmodel
## plot.STmodel          EX:aliased as plot.STdata
## qqnorm.STmodel        EX:with qqnorm.STdata
## scatterPlot.STmodel   Ex:with scatterPlot.STdata

##' Creates a \code{STmodel} object that can be for estimation and prediction.
##' For details see the sub-functions linked under the relevant Arguments.
##'
##' The object holds observations, trends, geographic, and spatio-temporal
##' covariates, as well as a number of precomputed fields that speed up
##' log-likelihood evaluations. To improve performance the locations are also
##' \strong{reorder} so that observed locations come before unobserved.
##' 
##' @title Construct STmodel Object
##' 
##' @param STdata \code{STdata} object with observations, covariates, trends, etc;
##'   see \code{\link{createSTdata}} or \code{\link{mesa.data.raw}} for an example.
##' @param LUR Specification of covariates for the beta-fields,
##'   see \code{\link{processLUR}}.
##' @param ST Specification of spatio-temporal covariates,
##'   see \code{\link{processST}}.
##' @param cov.beta,cov.nu Specification of the covariance functions,
##'   see \code{\link{updateCovf}}.
##' @param locations Specification of the sites (both monitored and un-monitored),
##'   see \code{\link{processLocation}}.
##' @param strip Should unobserved locations be dropped?
##' @param scale Scale the covariates? If \code{TRUE} all non-factor covariates
##'   are scaled \emph{after} the locations have been extracted but before
##'   constructing the covariate matrix for the beta-fields. (NOTE: If set to
##'   \code{TRUE} this scales the \code{LUR.all} elements to mean=0, sd=1).
##' @param scale.covars list with elements \code{mean} and \code{sd} giving the
##'   mean and standard deviation to use when scaling the covariates. Computed
##'   from \code{STdata$covars} if not given.
##' 
##' @return A \code{STmodel} object, see \code{\link{mesa.model}} for an example.
##' 
##' @example Rd_examples/Ex_createSTmodel.R
##' 
##' @author Johan Lindström
##' @family STmodel methods
##' @family STmodel functions
##' @family STdata functions
##' @export
createSTmodel <- function(STdata, LUR=NULL, ST=NULL,
                          cov.beta=list(covf="exp", nugget=FALSE),
                          cov.nu=list(covf="exp", nugget=TRUE,
                            random.effect=FALSE),
                          locations=list(coords=c("x","y"), long.lat=NULL,
                            coords.beta=NULL, coords.nu=NULL, others=NULL),
                          strip=FALSE, scale=FALSE, scale.covars=NULL){
  ##check class belonging
  stCheckClass(STdata, "STdata", name="STdata")

  ##Default values for several of the inputs
  cov.beta <- defaultList(cov.beta, eval(formals(createSTmodel)$cov.beta) )
  cov.nu <- defaultList(cov.nu, eval(formals(createSTmodel)$cov.nu) )
  locations <- defaultList(locations, eval(formals(createSTmodel)$locations ))
  
  ##add any nu-field nugget covariates to the locations
  if( class(cov.nu$nugget)=="formula" ){
    nugget.names <- all.vars(cov.nu$nugget)
  }else if( is.character(cov.nu$nugget) ){
    nugget.names <- cov.nu$nugget
  }else{
    nugget.names <- NULL
  }
  ##add the names needed by nugget to locations$others
  locations$others <- unique(c(locations$others, nugget.names))

  ##create output object, start with a copy
  STmodel <- STdata
  ##and assign both classes
  class(STmodel) <- c("STmodel","STdata")

  ##first add trend if missing - mainly backwards compability
  if( is.null(STmodel$trend) ){
    warning("No smooth trend in 'STdata', assuming only an intercept.")
    STmodel <- updateTrend.STdata(STmodel, n.basis=0)
  }else if( is.null(STmodel$trend.fnc) && dim(STmodel$trend)[2]!=1 ){
    ##trend function missing (and more than constant), add a trend function (backwards comp.)
    message("No trend $trend.fnc object detected, STdata probably from old ",
            "version of the package.\n$trend.fnc has been added based on ",
            "spline fit to elements in STmodel$trend.")
    STmodel$trend.fnc <- internalCreateTrendFnc(STmodel$trend)
  }
  ##sort temporal trend
  STmodel$trend <- STmodel$trend[order(STmodel$trend$date),,drop=FALSE]
  ##sort SpatioTemporal covariate to match time-points
  if( !is.null(STmodel$SpatioTemporal) ){
    STmodel$SpatioTemporal <-
      STmodel$SpatioTemporal[as.character(STmodel$trend$date),,,drop=FALSE]
  }
  
  ##should we drop some covariates?
  ID.unique <- sort(unique(STmodel$obs$ID))
  IND <- STmodel$covars$ID %in% ID.unique
  if( strip ){
    STmodel$covars <- STmodel$covars[IND,,drop=FALSE]
  }else{
    ##if not dropping at least order sites with observed locations first
    STmodel$covars <- rbind(STmodel$covars[IND,,drop=FALSE],
                           STmodel$covars[!IND,,drop=FALSE])
  }
  ##sort SpatioTemporal covariate to match locations
  if( !is.null(STmodel$SpatioTemporal) ){
    STmodel$SpatioTemporal <- STmodel$SpatioTemporal[,STmodel$covars$ID,,drop=FALSE]
  }
  
  ##create a seperate locations data.frame
  STmodel$locations.list <- locations
  STmodel$locations <- processLocation(STmodel, locations)
  
  ##add idx variables for the observations
  STmodel$obs$idx <- match(STmodel$obs$ID, STmodel$locations$ID)
  ##sort the observations
  STmodel$obs <- STmodel$obs[order(STmodel$obs$date,
                                   STmodel$obs$idx),,drop=FALSE]
  ##and remove rownames
  rownames(STmodel$obs) <- NULL

  ##Covariance specification, default values
  ##process covariance specification
  STmodel <- updateCovf(STmodel, cov.beta, cov.nu)
  
  ##should covariates be scaled
  if( scale ){
    if( is.null(scale.covars) ){
      scale.covars <- unlist(lapply(STmodel$covars,is.factor))
      ##mean and sd
      suppressWarnings( mean.covars <- unlist(lapply(STmodel$covars,mean)) )
      suppressWarnings( sd.covars <- unlist(lapply(STmodel$covars,sd)) )
      ##which fields do not have mean/sd
      scale.covars <- (is.na(mean.covars) | is.na(sd.covars) | scale.covars)
      mean.covars[scale.covars] <- NA
      sd.covars[scale.covars] <- NA
      ##save for future use
      STmodel$scale.covars <- list(mean=mean.covars, sd=sd.covars)
    }else{
      STmodel$scale.covars <- scale.covars
    }
  }
  
  ##process LUR
  STmodel$LUR.list <- processLUR(STmodel, LUR)
  ##...and ST specification
  STmodel$ST.list <- processST(STmodel, ST)
  ##create covariate matrices
  STmodel <- createLUR(STmodel, STmodel$LUR.list)
  ##create spatio-temporal covariate(s) for observations, ST=M and ST.all
  STmodel <- createST(STmodel, STmodel$ST.list)
  
  ##create temporal trends for observations, F
  STmodel$F <- internalSTmodelCreateF(STmodel)
  
  ##test for colocated monitoring sites, beta-fields
  I.idx <- unique(STmodel$obs$idx)
  if( anyDuplicated(STmodel$locations[I.idx,c("x.beta","y.beta")]) ){
    ##if so, require nugget
    if( any(!STmodel$cov.beta$nugget) ){
      stop("Colocated monitoring sites, nugget is required for beta-fields.")
    }
  }
  if( anyDuplicated(STmodel$locations[I.idx,c("x.nu","y.nu")]) ){
    ##if so, require nugget
    if( dim(STmodel$cov.nu$nugget.matrix)[2]==0 ){
      stop("Colocated monitoring sites, nugget is required for nu-fields.")
    }
  }
  ##compute distance matrix and nt
  STmodel <- createSTmodelInternalDistance(STmodel)
  
  ##drop SpatioTemporal, has been replaced by other things
  STmodel$SpatioTemporal <- NULL
  STmodel$covars <- NULL
  ##and asign one class
  class(STmodel) <- "STmodel"
  ##return
  return( STmodel )
}##function createSTmodel

############################
## S3-METHODS FOR STmodel ##
############################
##Additional (long) S3 methods in:
## STmodel_simulate.R
## STmodel_predict.R
## STmodel_combine.R

##' \code{\link[base:print]{print}} method for class \code{STmodel}.
##'
##' @title Print details for \code{STmodel} object
##' @param x \code{STmodel} object to print information for.
##' @param type Factorial of \code{length(x$locations$ID)}, if not \code{NULL}
##'   the output also presents summaries of number of sites and observations
##'   as well as time periods per type of site.
##' @param ... Ignored additional arguments.
##' @return Nothing
##' 
##' @examples
##' ##load some data
##' data(mesa.model)
##' ##print basic information regarding obs, locations, dates, etc
##' print(mesa.model)
##'
##' @author Johan Lindström
##' 
##' @family STmodel methods
##' @method print STmodel
##' @export
print.STmodel <- function(x, type=x$locations$type, ...){
  ##check class belonging
  stCheckClass(x, "STmodel", name="x")

  ##add covars field so we can use similarity with STdata objects
  x$covars <- x$locations
  
  ##general information regarding number of observations and trends
  commonPrintST(x, "STmodel", 1)

  ##covariate models
  cat( "Models for the beta-fields are:\n")
  print( x$LUR.list )
  if( !is.null(x$scale.covars) ){
    cat("Covariates have been scaled.\n\n")
  }
  
  ##spatio-temporal models
  if( length(x$ST.list)==0 ){
    cat("No spatio-temporal covariates.\n")
  }else{
    cat( sprintf("%d spatio-temporal covariate(s):\n", length(x$ST.list)) )
    print( x$ST.list )
  }
  cat("\n")

  ##covariance models
  cat("Covariance model for the beta-field(s):\n")
  cat( paste("\tCovariance type(s):", paste(x$cov.beta$covf, collapse=", "), "\n"))
  cat( paste("\tNugget:", paste(c("No","Yes")[x$cov.beta$nugget+1],
                                collapse=", "), "\n"))
  
  cat("Covariance model for the nu-field(s):\n")
  cat( paste("\tCovariance type:", x$cov.nu$covf, "\n") )
  if( dim(x$cov.nu$nugget.matrix)[2]==0 ){
    cat( "\tNugget: No\n" )
  }else{
    cat( paste("\tNugget:", paste(x$cov.nu$nugget,collapse=""), "\n") )
  }
  cat( paste("\tRandom effect:", c("No","Yes")[x$cov.nu$random.effect+1], "\n"))

  ##group observations by type.
  commonPrintST(x, "STmodel", 2, type)

  return(invisible())
}##function print.STmodel

##' \code{\link[base:summary]{summary}} method for class \code{STmodel}.
##'
##' @title Computes summary details for \code{STmodel} object
##' @param object \code{STmodel} object to compute summary information for.
##' @param type Factorial of \code{length(x$locations$ID)}, if not \code{NULL}
##'   the output also presents summaries of number of sites and observations
##'   as well as time periods per type of site.
##' @param ... Ignored additional arguments. 
##' @return A \code{summary.STmodel} object.
##' 
##' @examples
##' ##load some data
##' data(mesa.model)
##' ##Summary of data fields.
##' summary(mesa.model)
##'
##' @author Johan Lindström
##' 
##' @family STmodel methods
##' @method summary STmodel
##' @export
summary.STmodel <- function(object, type=object$covars$type, ...){
  ##check class belonging
  stCheckClass(object, "STmodel", name="object")
  
  ##allocate output object
  out <- list()
  ##compute various summaries, locations
  out$locations <- summary(object$locations)
  ##LUR:s

  ##Are all locations observed? if not compute also only obs summary
  if( dim(object$locations)[1]!=length(unique(object$obs$ID)) ){
    out$LUR <- lapply(object$LUR,summary)
  }
  out$LUR.all <- lapply(object$LUR.all,summary)
  
  ##ST-covariates
  if( !is.null(object$ST.all) ){
    out$SpatioTemporal <- summary( matrix(object$ST.all, prod(dim(object$ST.all)[1:2]),
                                          dim(object$ST.all)[3]) )
    colnames(out$SpatioTemporal) <- dimnames(object$ST.all)[[3]]
  }

  ##add covars field so we can use similarity with STdata objects
  object$covars <- object$locations
  ##common summary computations observations and trend
  out <- c(out, commonSummaryST(object, type))

  ##return the object
  class(out) <- "summary.STmodel"
  return(out)
}##function summary.STmodel

##' \code{\link[base:print]{print}} method for class \code{summary.STmodel}.
##'
##' @title Print details for \code{summary.STmodel} object
##' @param x \code{summary.STmodel} object to print information for.
##' @param ... Ignored additional arguments.
##' @return Nothing
##'
##' @author Johan Lindström
##' 
##' @family STmodel methods
##' @method print summary.STmodel
##' @export
print.summary.STmodel <- function(x, ...){
  ##check class belonging
  stCheckClass(x, "summary.STmodel", name="x")

  ##print data
  if( is.null(x[["obs"]]) ){
    cat("No observations.\n");
  }else{
    cat("Summary of observations:\n");
    print(x$obs)
  }
  cat("\nSummary of locations:\n");
  print(x$locations)
  ##Are all locations observed?
  if( !is.null(x[["LUR"]]) ){
    if( !is.null(x$LUR) ){
      cat("\nSummary of geographic covariates (for observation locations):\n");
      print(x$LUR)
    }
    cat("\nSummary of geographic covariates (for all locations):\n");
    print(x$LUR.all)
  }else{
    cat("\nSummary of geographic covariates:\n");
    print(x$LUR.all)
  }
  cat("\nSummary of smooth trends:\n");
  print(x$trend)

  if( is.null(x$SpatioTemporal) ){
    cat("\nNo spatio-temporal covariates.\n")
  }else{
    cat("\nSummary of spatio-temporal covariates.\n")
    print(x$SpatioTemporal)
  }
  if( !is.null(x$obs.by.type) ){
    for(i in 1:length(x$obs.by.type)){
      cat( paste("\nSummary for observations of type",
                 names(x$obs.by.type)[i],"\n") )
      print(x$obs.by.type[[i]])
    }
  }##if( !is.null(x$obs.by.type) )
  
  return(invisible())
}##function print.summary.STmodel

##' @rdname plot.STdata
##' @family STmodel methods
##' @method plot STmodel
##' @export
plot.STmodel <- function(x, y="obs", ID=x$locations$ID[1],
                         type=x$locations$type, ...){
  ##check class belonging
  stCheckClass(x, "STmodel", name="x")

  ##recast and introduce covars
  x$covars <- x$locations
  class(x) <- "STdata"

  plot(x, y=y, ID=ID, type=type, ...)
  
  return(invisible())
}##function plot.STmodel


##' @rdname qqnorm.STdata
##' @family STmodel methods
##' @importFrom stats qqnorm
##' @method qqnorm STmodel
##' @export
qqnorm.STmodel <- function(y, ID="all", main="Q-Q plot for observations",
                           group=NULL, col=1, line=0, ...){ 
  ##check class belonging
  stCheckClass(y, "STmodel", name="y")

  Y <- y$obs[, c("ID","obs")]

  internalQQnormPlot(Y, ID, main, group, col, FALSE, line, ...)
  
  return(invisible())
}##function qqnorm.STmodel

##' @rdname scatterPlot.STdata
##' @family STmodel methods
##' @method scatterPlot STmodel
##' @export
scatterPlot.STmodel <- function(x, covar=NULL, trend=NULL, pch=1, col=1, cex=1,
                                lty=1, subset=NULL, group=NULL, add=FALSE,
                               smooth.args=NULL, ...){
  ##check class belonging
  stCheckClass(x, "STmodel", name="x")
  ##collect location & LUR information
  covars <- cbind(x$locations, do.call(cbind, x$LUR))
  covars <- covars[,unique(colnames(covars))]
  
  ##pass data to internalScatterPlot function
  internalScatterPlot(obs=x$obs[, c("obs","ID","date")],
                      covar=covar, trend=trend, subset=subset,
                      data=list(covars=covars, trend=x$trend),
                      group=group, pch=pch, col=col, cex=cex,
                      lty=lty, add=add, smooth.args=smooth.args, ...)
}##function scatterPlot.STmodel
