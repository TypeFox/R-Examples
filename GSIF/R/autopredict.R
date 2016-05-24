# Purpose        : Automate prediction of soil properties / soil types;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ;
# Dev Status     : Pre-Alpha
# Note           : ;

setMethod("autopredict", signature(target = "SpatialPointsDataFrame", covariates = "SpatialPixelsDataFrame"), function(target, covariates, auto.plot=TRUE, ...){
  
  ## parent call:
  parent_call <- as.list(substitute(list(...)))[-1]
  ## TH: TO-DO estimate processing time
  
  ## predictive components:
  spc.fm <- as.formula(paste("~", paste(names(covariates), collapse = "+")))
  covariates <- spc(covariates, spc.fm)
  ## Model:
  fm <- as.formula(paste(names(target)[1], "~", paste(names(covariates@predicted), collapse = "+")))
  if(is.factor(target@data[,1])){
    m <- spmultinom(fm, target, covariates@predicted, ...)
    if(auto.plot==TRUE){ plotKML(m@predicted, folder.name=names(target)[1], file.name=paste0(names(target)[1], "_predicted.kml")) }
    return(m)
  }
  if(is.numeric(target@data[,1])){
    if(!any(names(parent_call) %in% "method")){
      method <- "ranger"
    }
    ## TH: TO-DO add ensemble predictions
    m <- fit.gstatModel(target, fm, covariates@predicted, method=method, ...)
    ## predict:
    p <- predict(m, covariates@predicted)
    if(auto.plot==TRUE){ plotKML(p, folder.name=names(target)[1], file.name=paste0(names(target)[1], "_predicted.kml")) }
    return(p)
  }

})