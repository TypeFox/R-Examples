############################################
## S3-METHOD THAT PREDICTS FROM A STmodel ##
############################################
##Functions in this file:
## c.STmodel              Ex:ok
## combineSTmodel         Ex:INTERNAL
## areConsistent          Ex:INTERNAL
## areSTmodelsConsistent  Ex:INTERNAL


##' Combines several locations and covariates for several STmodel/STdata objects.
##' Temporal trend, observations and covariance model (both spatial and
##' spatio-temporal) are taken from the first object in the call. Any additional
##' covariates/trends/observations not present in the first argument are dropped
##' from the additional arguments \emph{without warning}. 
##' Locations and covariates (both spatial and spatio-temporal) from
##' additional objects are added to those in the first object.
##' 
##' For additional \code{STdata} objects the covariates are transformed according to
##' \code{STmodel$scale.covars} of the first object, see
##' \code{\link{createSTmodel}}.
##' 
##' For \code{STmodel} objects \strong{can not} be combined if either has scaled
##' covariates.
##' 
##' @title Combine Several STmodel/STdata Objects
##' 
##' @param ... \code{STmodel} and \code{STdata} objects to combine, the first
##'   object has to be a \code{STmodel}.
##' @param recursive For S3 compatibility; the function will ALWAYS run
##'   recursively
##'
##' @return An updated \code{STmodel} object.
##' 
##' @example Rd_examples/Ex_c_STmodel.R
##'
##' @author Johan Lindström
##' 
##' @family STmodel methods
##' @family STdata functions
##' @method c STmodel
##' @export
c.STmodel <- function(...,recursive=FALSE){
  input <- list(...)
  ##check class belongings
  stCheckClass(input[[1]], "STmodel", name="First argument")
  if( length(input)==1 ){
    ##only one input, return it
    return( input[[1]] )
  }
  ##combine the objects, one by one
  out <- input[[1]]
  for(i in 2:length(input)){
    out <- combineSTmodel(out, input[[i]], i)
  }
  ##update covariance functions
  out <- updateCovf(out)
  
  ##return combined object.
  return(out)
}

##internal function that combined STmodel with a second STdata or STmodel object
combineSTmodel <- function(STmodel, STdata, i.arg){
  ##check second input, first already checked by c.STmodel
  stCheckClass(STdata, c("STmodel","STdata"), name=paste("argument no.",i.arg))
    
  ##is second object of type STdata, convert to STmodel
  if( inherits(STdata,"STdata") ){
    ##set a trend of the right size, use fnc in model-object
    suppressMessages( STdata <- updateTrend(STdata, fnc=STmodel$trend.fnc,
                                            extra.dates=STdata$trend$date) )
    ##since we're updating the covariance fucntions once everything has been
    ##added together, pick a simple covariance structure
    cov.nu <- STmodel$cov.nu
    cov.nu$nugget <- TRUE
    STdata <- createSTmodel(STdata, LUR=STmodel$LUR.list, ST=STmodel$ST.list,
                            cov.beta=STmodel$cov.beta, cov.nu=cov.nu,
                            locations=STmodel$locations.list,
                            scale=!is.null(STmodel$scale.covars),
                            scale.covars=STmodel$scale.covars)
  }else{
    ##trend is ignored, no need to check, but check everything else
    areSTmodelsConsistent(STmodel, STdata, i.arg)
  }##if( inherits(STdata,"STdata") ){...}else{...}

  ##combine the two datasets, trend of STdata is ignored.
  ##locations
  I.new <- !(STdata$locations$ID %in% STmodel$locations$ID)
  if( sum(I.new)==0 ){
    ##no new locations
    return( STmodel )
  }
  STmodel$locations <- rbind(STmodel$locations,
                             STdata$locations[I.new,,drop=FALSE])
  ##LURs
  for(i in 1:length(STmodel$LUR.all) ){
    STmodel$LUR.all[[i]] <- rbind(STmodel$LUR.all[[i]],
                                  STdata$LUR.all[[i]][I.new,,drop=FALSE])
  }
  ##ST
  if( !is.null(STmodel$ST.all) ){
    tmp <- array(NA, c(dim(STmodel$ST.all)[1], dim(STmodel$locations)[1],
                       dim(STmodel$ST.all)[3]) )
    dimnames(tmp) <- list(rownames(STmodel$ST.all), STmodel$locations$ID,
                          dimnames(STmodel$ST.all)[[3]])
    tmp[, colnames(STmodel$ST.all), ] <- STmodel$ST.all
    tmp[, STdata$locations$ID[I.new], ] <- STdata$ST.all[rownames(tmp), I.new,]
    if( any(is.na(tmp)) )
      stop( paste("Some ...$ST.all values missing, probably from argument no.",
                    i.arg) )
    STmodel$ST.all <- tmp
  }
  
  ##return object
  return( STmodel )
}

##internal function that compares two matrices to determine if they have
##the same columns (number and names)
areConsistent <- function(X, Y, names="X, Y"){
  if( dim(X)[2]!=dim(Y)[2] || any(colnames(X)!=colnames(Y)) ){
    stop( paste("Matrices", names, "do not have the same columns (no./names).") )
  }
  return(invisible())
}

##internal function that compares two STmodel objects to determine if there
##locations and LUR columns are consistent (and if ST-covars in the second
##model, matches those of the first.)
areSTmodelsConsistent <- function(model1, model2, i.arg){
  if( !isTRUE( all.equal(model1$scale.covars, model2$scale.covars) )){
    stop( "Un-able to combine STmodel objects with different scaling." )
  }

  ##check locations
  areConsistent(model1$locations, model2$locations, "...$locations")
  ##check LUR
  if( length(model1$LUR.all) != length(model2$LUR.all) ){
    stop( paste("Unequal number of LURs for argument no.", i.arg) )
  }
  for(i in 1:length(model1$LUR.all) ){
    areConsistent(model1$LUR.all[[i]], model2$LUR.all[[i]],
                  paste("...$LUR[[", i, "]]", sep="") )
  }
  ##check ST
  if( !is.null(model1$ST.all) ){
    if( is.null(model2$ST.all) ){
      stop( paste("spatio-temporal covariate missing from argument no.",
                  i.arg) )
    }
    if( dim(model1$ST.all)[3]!=dim(model2$ST.all)[3] ||
       dimnames(model1$ST.all)[[3]]!=dimnames(model2$ST.all)[[3]] ){
      stop( paste("spatio-temporal covariates inconsistent, argument no.",
                  i.arg) )
    }
    if( any(!(rownames(model1$ST.all) %in% rownames(model2$ST.all))) ){
      stop( paste("spatio-temporal dates missing, argument no.",
                  i.arg) )
    }
  }##if( !is.null(model1$ST.all) ){
}##function areSTmodelsConsistent
