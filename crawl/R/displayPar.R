#'Display the order of parameters along with fixed values and starting values
#'
#'This function takes the model spesification arguments to the \code{\link{crwMLE}} function and displays a table
#'with the parameter names in the order that \code{crwMLE} will use during model fitting. This is useful for specifying 
#'values for the \code{fixPar} or \code{theta} (starting values for free parameters) arguments. 
#'
#'@param mov.model formula object specifying the time indexed covariates for
#'movement parameters.
#'@param err.model A 2-element list of formula objects specifying the time
#'indexed covariates for location error parameters.
#'@param activity formula object giving the covariate for the stopping
#'portion of the model.
#'@param drift logical indicating whether or not to include a random
#'drift component.
#'@param data data.frame object containg telemetry and covariate data. A
#'\code{SpatialPointsDataFrame} object from the package 'sp' will also be accepted.
#'@param theta starting values for parameter optimization.
#'@param fixPar Values of parameters which are held fixed to the given value.
#'@param ... Additional arguments (probably for testing new features.)
#' 
#'@return A data frame with the following columns
#' 
#'\item{ParNames}{The names of the parameters specified by the arguments.}
#'
#'\item{fixPar}{The values specified by the \code{fixPar} argument for fixed values of the parameters. In model fitting, 
#'these values will remain fixed and will not be estimated.}
#'
#'\item{thetaIndex}{This column provides the index of each element of the theta argument and to which parameter it corresponds.}
#'
#'\item{thetaStart}{If a value is given for the \code{theta} argument it will be placed in this column and its elements will 
#'correspond to the \code{thetaIdx} column.}
#' 
#'@author Devin S. Johnson
#'@seealso \code{demo(northernFurSealDemo)} for example.
#'  
#'@export 

displayPar <- function(mov.model=~1, err.model=NULL, activity=NULL, drift=FALSE, data, theta, fixPar, ...){
  if(inherits(data, "trip")){
    Time.name <- data@TOR.columns[1]
  }
  if(inherits(data, "SpatialPoints")) {  
    if("+proj=longlat" %in% strsplit(sp::proj4string(data), " ")[[1]]) stop("Location data must be projected.")	
    coordVals <- as.data.frame(sp::coordinates(data))	
    coord <- names(coordVals)	
    data <- cbind(slot(data,"data"), coordVals)    
  }
#   if(inherits(data[,Time.name],"POSIXct")){
#     data$TimeNum <- as.numeric(data[,Time.name])#/3600
#     Time.name <- "TimeNum"
#   }
  
  
  ### Check for duplicate time records ###
  #if(any(diff(data[,Time.name])==0)) stop("There are duplicate time records for some data entries! Please remove before proceeding.")
  
  
  ## SET UP MODEL MATRICES AND PARAMETERS ##
  errMod <- !is.null(err.model)
  activeMod <- !is.null(activity)
  driftMod <- drift
  mov.mf <- model.matrix(mov.model, model.frame(mov.model, data, na.action=na.pass))
  if (any(is.na(mov.mf))) stop("Missing values are not allowed in movement covariates!")
  n.mov <- ncol(mov.mf)
  if (errMod) {
    err.mfX <- model.matrix(err.model$x,model.frame(err.model$x, data, na.action=na.pass))
    err.mfX <- ifelse(is.na(err.mfX), 0, err.mfX)
    n.errX <- ncol(err.mfX)
    if (!is.null(err.model$y)) {
      err.mfY <- model.matrix(err.model$y,model.frame(err.model$y, data, na.action=na.pass))
      err.mfY <- ifelse(is.na(err.mfY), 0, err.mfY)
      n.errY <- ncol(err.mfY)
    } else {
      err.mfY <- NULL
      n.errY <- 0
    }
    if(!is.null(err.model$rho)){
      rho = model.matrix(err.model$rho,model.frame(err.model$rho, data, na.action=na.pass))[,-1]
      if(any(rho > 1 | rho < -1, na.rm=TRUE)) stop("Error model correlation outside of the range (-1, 1).")
    } else rho = NULL
  } else {
    n.errY <- n.errX <- 0
    err.mfX <- err.mfY <- rho <- NULL
  }
  if (activeMod) {
    #stop.model
    activity <- model.matrix(activity, model.frame(activity, data, na.action=na.pass))
    if (ncol(activity) > 2) stop("There can only be one activity variable.")
    activity <- as.double(activity[,2])
    if (any(activity < 0) | any(activity > 1)) stop("'activity' variable must be >=0 and <=1.")
    if (any(is.na(activity))) stop("Missing values are not allowed in the activity variable.")
  } else activity <- NULL
  n.drift <- as.integer(driftMod)
  n.activ <- as.integer(activeMod)
  b.nms <- paste("ln beta ", colnames(mov.mf), sep="")
  sig.nms <- paste("ln sigma ", colnames(mov.mf), sep="")
  if (errMod) {
    if (!is.null(err.model$y)) {
      tau.nms <- c(paste("ln tau.x ", colnames(err.mfX), sep=""),
                   paste("ln tau.y ", colnames(err.mfY), sep=""))
    } else tau.nms <- paste("ln tau ", colnames(err.mfX), sep="")
  } else tau.nms <- NULL
  if (activeMod){
    active.nms <- "ln phi"
  } else active.nms <- NULL
  if (driftMod) {
    drift.nms <- c("ln sigma.drift/sigma", "ln psi-1")
  } else drift.nms <- NULL
  nms <- c(tau.nms, sig.nms, b.nms, active.nms, drift.nms)
  n.par <- length(nms)
  if (missing(fixPar)) fixPar <- rep(NA, n.par)
  if (length(fixPar)!=n.par) stop("'fixPar' argument is not the right length! The number of parameters in the model is ", n.par, "\n")
  if (!missing(theta)) if(length(theta) != sum(is.na(fixPar))) stop("\nWrong number of parameters specified in start value.\n")
  thetaIdx <- fixPar
  thetaIdx[is.na(fixPar)] <- 1:sum(is.na(fixPar))
  thetaIdx[!is.na(fixPar)] <- NA
  out <- data.frame(ParNames=nms, fixPar=fixPar, thetaIdx=thetaIdx)
  if(!missing(theta)){
    thetaStart <- thetaIdx
    thetaStart[!is.na(thetaStart)] <- theta
    out <- cbind(out, thetaStart=thetaStart)
  }
  return(out)
}