# Purpose        : Fit a 2D or 3D regression-kriging model;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : linear models with normally distributed residuals;


## Fit a 'simple' 2D RK model:
setMethod("fit.gstatModel", signature(observations = "SpatialPointsDataFrame", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), function(observations, formulaString, covariates, method = list("GLM", "rpart", "randomForest", "quantregForest", "xgboost", "ranger"), dimensions = list("2D", "3D", "2D+T", "3D+T"), fit.family = gaussian(), stepwise = TRUE, vgmFun = "Exp", subsample = 5000, subsample.reg = 10000, ...){
 
  method <- method[[1]]
  dimensions <- dimensions[[1]]
  ## TH: the function only works with 2D maps at the moment:
  if(length(attr(coordinates(observations), "dimnames")[[2]])>2){
    warning("This method uses only 2D coordinates of the points. For 3D data consider using the 'geosamples-class'.")
  }
  ## TH: the model has to include at least one covariate:
  if(length(all.vars(formulaString)[-1])==0){
    stop("No covariates have been specified in the 'formulaString'")
  }
  ## check argument 'fit.family'
  if(!missing(fit.family) & !method=="GLM"){ warning("'fit.family' argument will be ignored. Using 'method=\"GLM\"' instead.")  }

  ## overlay:
  ov <- over(observations, covariates)
  ## all variables of interest:
  tv <- all.vars(formulaString)[1]
  seln <- names(covariates) %in% all.vars(formulaString)[-1]
  ## check if all covariates are available: 
  if(sum(!is.na(seln))==0){
      stop("None of the covariates in the 'formulaString' matches the names in the 'covariates' object")
  }
  ov <- cbind(data.frame(observations[,tv]), ov)  
  
  ## copy coordinate column names for consistency:
  xyn <- which(names(ov) %in% attr(observations@bbox, "dimnames")[[1]])
  if(length(xyn)==2) { 
    names(ov)[xyn] <- attr(covariates@bbox, "dimnames")[[1]][1:2] 
    } else {
    names(ov)[xyn] <- attr(covariates@bbox, "dimnames")[[1]]     
  }

  ## check the size of output:
  if(nrow(ov)==0|is.null(ov[,tv])) {
    stop("The 'over' operations resulted in an empty set.")
  }
  
  ## fit/filter the regression model:
  m <- fit.regModel(formulaString = formulaString, rmatrix = ov, predictionDomain = covariates[seln], method = method, dimensions = dimensions, fit.family = fit.family, stepwise = stepwise, vgmFun = vgmFun, subsample = subsample, subsample.reg = subsample.reg, ...)
  return(m)
  
})


## Fit a RK model using geosamples class:
setMethod("fit.gstatModel", signature(observations = "geosamples", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), function(observations, formulaString, covariates, method = list("GLM", "rpart", "randomForest", "quantregForest", "xgboost", "ranger"), dimensions = list("2D", "3D", "2D+T", "3D+T"), fit.family = gaussian(), stepwise = TRUE, vgmFun = "Exp", subsample = 5000, subsample.reg = 10000, ...){
  
  method <- method[[1]]
  dimensions <- dimensions[[2]]
  ## TH: the model has to include at least one covariate:
  if(length(all.vars(formulaString)[-1])==0){
    stop("No covariates have been specified in the 'formulaString'")
  } 
  
  ## all columns of interest:
  methodid <- all.vars(formulaString)[1]
  seln <- names(covariates) %in% all.vars(formulaString)[-1]
  xyn <- attr(covariates@bbox, "dimnames")[[1]]
  ## check if all covariates are available: 
  if(sum(!is.na(seln))==0){
      stop("None of the covariates in the 'formulaString' matches the names in the 'covariates' object")
  }
  
  ## prepare regression matrix:
  ov <- over(x=covariates, y=observations, method=methodid, var.type = "numeric")
  if(nrow(ov)==0|is.null(ov$observedValue)) {
    warning("The 'over' operations resulted in an empty set. Check 'methodid' column.")
  }
  ## geostats only possible with numeric variables:
  ov[,methodid] = as.numeric(ov$observedValue)
  ov$observedValue = NULL
   
  ## subset to columns of interest:  
  if(dimensions == "3D"){ ov <- ov[, c(methodid, names(covariates)[seln], xyn, "altitude")] }
  if(dimensions == "2D"){ ov <- ov[, c(methodid, names(covariates)[seln], xyn)] }
  
  ## fit/filter the regression model:
  m <- fit.regModel(formulaString = formulaString, rmatrix = ov, predictionDomain = covariates[seln], method = method, dimensions = dimensions, fit.family = fit.family, stepwise = stepwise, vgmFun = vgmFun, subsample = subsample, subsample.reg = subsample.reg, ...)
  return(m)  

})


## Fit a RK model to a list of covariates / formulas:
setMethod("fit.gstatModel", signature(observations = "geosamples", formulaString = "list", covariates = "list"), function(observations, formulaString, covariates, method = list("GLM", "rpart", "randomForest", "quantregForest", "xgboost", "ranger"), dimensions = list("2D", "3D", "2D+T", "3D+T"), fit.family = gaussian(), stepwise = TRUE, vgmFun = "Exp", subsample = 5000, subsample.reg = 10000, ...){

  method <- method[[1]]
  dimensions <- dimensions[[2]]
  if(!length(formulaString)==length(covariates)){
    stop("'formulaString' and 'covariates' lists of same size expected")
  }
  
  rkm.l <- list(NULL)
  for(l in 1:length(covariates)){   
    rkm.l[[l]] <- fit.gstatModel(observations, formulaString[[l]], covariates[[l]], method = method, fit.family = fit.family, stepwise = stepwise, subsample = subsample, subsample.reg = subsample.reg, dimensions = dimensions, ...)
  }
  return(rkm.l)  

})


## Fit a RK model and return an object of class "gstatModel" for a list of multiscale grids:
setMethod("fit.gstatModel", signature(observations = "geosamples", formulaString = "formula", covariates = "list"), function(observations, formulaString, covariates, method = list("GLM", "rpart", "randomForest", "quantregForest", "xgboost", "ranger"), dimensions = list("2D", "3D", "2D+T", "3D+T"), fit.family = gaussian(), stepwise = TRUE, vgmFun = "Exp", subsample = 5000, subsample.reg = 10000, ...){

  method <- method[[1]]                                                              
  dimensions <- dimensions[[2]]
  if(!any(sapply(covariates, class)=="SpatialPixelsDataFrame")){
    stop("List of covariates of class 'SpatialPixelsDataFrame' expected")
  }

  if(!missing(fit.family) & !method=="GLM"){ 
    warning("'fit.family' argument will be ignored. Use 'method=\"GLM\"' instead.")
  }

  ## covariate names:
  covs <- unlist(sapply(covariates, FUN=function(x){names(x)}))
  if(!length(unique(covs))==length(covs)){ stop("'Covariates' column names must be unique") }
  
  ## target variable:
  methodid <- all.vars(formulaString)[1]
     
  ## prepare regression matrix:
  ov <- list(NULL)
  for(j in 1:length(covariates)){
    ov[[j]] <- over(x=covariates[[j]], y=observations, method=methodid, var.type="numeric")
      if(nrow(ov[[j]])==0|is.null(ov[[j]]$observedValue)) {
      warning("The 'over' operations resulted in an empty set. Check 'methodid' column.")
      }
  }
  xyn = attr(covariates[[1]]@bbox, "dimnames")[[1]]
  ## coordinates in the local space:
  xy <- ov[[1]][,xyn]
  
  gn <- names(observations@data)[!(names(observations@data) %in% c("observationid", xyn))]
  ## merge all regression matrices:
  ov <- Reduce(function(x,y) {merge(x[!(names(x) %in% gn)],y[!(names(y) %in% gn)], by="observationid")}, ov)
  ov <- merge(observations@data, ov[,c("observationid", covs)], by="observationid")
  ov <- cbind(ov, xy) 

  ## geostats only possible with numeric variables:  
  ov[,methodid] <- as.numeric(ov$observedValue)
  ov$observedValue <- NULL
  
  ## check the size of the output object:
  if(nrow(ov)==0|is.null(ov[,methodid])) {
    warning("The 'over' operations resulted in an empty set. Check 'methodid' column.")
  }
  
  ## fit/filter the regression model:
  m <- fit.regModel(formulaString = formulaString, rmatrix = ov, predictionDomain = covariates[[1]], method = method, dimensions = dimensions, fit.family = fit.family, stepwise = stepwise, vgmFun = vgmFun, subsample = subsample, subsample.reg = subsample.reg, ...)
  return(m) 

})


"print.gstatModel" <- function(x, ...){
  print(x@regModel)
  print(x@vgmModel)
  summary(x@sp)
}

## BK: plot functionality; see also ?plotKML::plot.SpatialPredictions
"plot.gstatModel" <- function(x, ...){
  dev.new(width=9, height=5)
  par(mfrow=c(1,2))
  if(any(class(x@regModel) == "lm")|any(class(x@regModel)=="quantregForest")|any(class(x@regModel)=="randomForest")|any(class(x@regModel)=="rpart")|any(class(x@regModel)=="ranger")|any(class(x@regModel)=="xgboost")){
    if(any(class(x@regModel) == "lm")){
      plot(y=fitted(x@regModel), x=x@regModel$y, col = alpha("grey18", 0.6), pch=16, xlab='observed', ylab='predicted', main='Goodness of fit', asp=1, xlim=range(x@regModel$y), ylim=range(x@regModel$y), ...)
    }
    if(class(x@regModel)[1]=="quantregForest"){
      plot(y=x@regModel$predicted, x=x@regModel$y, col = alpha("grey18", 0.6), pch=16, xlab='observed', ylab='predicted', main='Goodness of fit', asp=1, xlim=range(x@regModel$y), ylim=range(x@regModel$y), ...)
    } else {
      if(any(class(x@regModel)=="randomForest")|any(class(x@regModel)=="rpart")|any(class(x@regModel)=="xgboost")){
        plot(y=predict(x@regModel), x=x@regModel$y, col = alpha("grey18", 0.6), pch=16, xlab='observed', ylab='predicted', main='Goodness of fit', asp=1, xlim=range(x@regModel$y), ylim=range(x@regModel$y), ...)
      }
      if(any(class(x@regModel)=="ranger")){
        plot(y=predict(x@regModel, x@sp@data)$predictions, x=x@sp@data[,1], col = alpha("grey18", 0.6), pch=16, xlab='observed', ylab='predicted', main='Goodness of fit', asp=1, xlim=range(x@sp@data[,1]), ylim=range(x@sp@data[,1]), ...)
      }
    }
  } else {   
    pred <- fitted(x@regModel)
    obs <- pred+resid(x@regModel)
    plot(y=pred, x=obs, col = alpha("grey18", 0.6), pch=16, xlab='observed', ylab='predicted', main='Goodness of fit', asp=1, xlim=range(obs), ylim=range(obs), ...)
  }
  abline(a=0, b=1, lwd=2, lty=2, col="green")
  vgmmodel <- x@vgmModel
  class(vgmmodel) <- c("variogramModel","data.frame")
  plot(x=x@svgmModel$dist, y=x@svgmModel$gamma, pch="+", col = "grey18", xlab='distance', cex=1.1, ylab='gamma', ylim = c(0, max(x@svgmModel$gamma)), main='Residual variogram')
  vline <- variogramLine(vgmmodel, maxdist=max(x@svgmModel$dist), n=length(x@svgmModel$dist))
  lines(x=vline$dist, y=vline$gamma, col="green", lwd=2)
}

setMethod("plot", signature("gstatModel"), plot.gstatModel)

# end of script;