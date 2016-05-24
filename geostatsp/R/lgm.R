setGeneric('lgm', function(
        formula,data,grid=NULL,covariates=NULL, ...) 
			standardGeneric("lgm")
)

# sort out formula
# null formula
setMethod("lgm", 
		signature("NULL"), 
    function(formula=NULL, data, grid, 
        covariates=NULL, ...) {
    formula =  1 
  callGeneric(formula, data, grid, covariates, ...)
  }
  )


setMethod("lgm", 
		signature("numeric"),
    function(formula, data, grid, 
        covariates=NULL, ...) {
      
      formula = names(data)[formula]
      callGeneric(formula, data, grid, covariates, ...)
    }
)

# change character to formula
setMethod("lgm", 
		signature("character"),  
		function(formula, data, grid, 
        covariates=NULL, ...) {
      
      if(length(names(covariates)))
        names(covariates) = gsub("[[:punct:]]|[[:space:]]","_", names(covariates))
      if(length(covariates) & !length(names(covariates))) 
        names(covariates) = paste("c", 1:length(covariates),sep="")			
      
      if(length(formula)==1)
        formula = unique(c(formula, names(covariates)))
      if(length(formula)==1)
        formula = c(formula, '1')
      
      formula = paste(formula[1] , "~",
          paste(formula[-1], collapse=" + ")
      )
      formula = as.formula(formula)
      
      callGeneric(formula, data, grid, covariates, ...)
    }
)





# numeric cells, create raster from data bounding box

setMethod("lgm", 
		signature("formula", "Spatial", "numeric", "ANY"),
    function(formula, data, grid, covariates=NULL, ...) {
      grid = squareRaster(data, grid)
      callGeneric(formula, data, grid, covariates, ...)
    }
)

# missing covariates, create empty list
setMethod("lgm", 
    signature("formula", "Spatial", "Raster", "missing"),
    function(formula, data, grid=NULL, covariates=NULL, ...) {
      callGeneric(formula, data, grid, 
          covariates=list(), 
          ...)
    }
)



setMethod("lgm",
			signature("formula", "Spatial", "Raster", "list"),
      function(formula, 
          data, grid, 
          covariates=NULL, 
          buffer=0,...) {

        dataCov = gm.dataSpatial(
            formula, data, 
            grid, covariates, buffer)
        
        callGeneric(formula, 
            data=dataCov$data, 
            grid=dataCov$grid, 
            covariates=dataCov$covariates, ...)
      }
)

setMethod("lgm",
    signature("formula", "Spatial", "Raster", "Raster"),
    function(formula, 
        data, grid, 
        covariates=NULL, 
        buffer=0,...) {
      
      dataCov = gm.dataSpatial(
          formula, data, 
          grid, covariates, buffer)
      
      callGeneric(formula, 
          data=dataCov$data, 
          grid=dataCov$grid, 
          covariates=dataCov$covariates, ...)
    }
)

# the real work
setMethod("lgm", 
		signature("formula", "Spatial", "Raster","data.frame"), 
		function(formula, 
				data,
				grid,
		covariates=NULL, 
		shape=1, boxcox=1, nugget = 0, 
		expPred=FALSE, nuggetInPrediction=TRUE,
		reml=TRUE,mc.cores=1,
		aniso=FALSE,
		fixShape=TRUE,
		fixBoxcox=TRUE,
		fixNugget = FALSE,
		...){
  
	locations = grid
	
	
	dots <- list(...)  
  param=dots$param	
  if(!length(param)) {
		param=c()
	}
	

	paramToEstimate	= c(
      "variance", "range", "shape","nugget","boxcox"
          )[c(
                  TRUE, TRUE, !fixShape, !fixNugget, !fixBoxcox
          )]		
	range=NA
	Spar = c(shape=as.numeric(shape),
      nugget=as.numeric(nugget),
      range=NA,
      boxcox=as.numeric(boxcox))
	
	if(aniso) {
		Spar = c(Spar, anisoAngleRadians=NA,anisoRatio=NA)
		paramToEstimate = c(paramToEstimate,
				"anisoAngleDegrees","anisoRatio")		
	}

	Spar = Spar[!names(Spar) %in% names(param)]
	param = c(param, Spar)
	
# to do: make sure factors in rasters are set up correctly
	# have baseline as first entry in cov@data@attributes,
	# NA's for levels without data
	# have most common level the baseline
	
# call likfit
	
	dots$param = param
	dots$formula=formula
	dots$data=data
	dots$paramToEstimate=paramToEstimate
  dots$reml = reml

 	likRes = do.call(likfitLgm, dots)
 

  

# call krige	
	krigeRes =  krigeLgm(
			formula=formula,data=data,
			grid=grid,
			covariates=covariates, param=likRes$param, 
			expPred=expPred,
			nuggetInPrediction=nuggetInPrediction
			)
			

	res = c(predict=krigeRes, likRes)
#  stuff <<- res
  
  # add confidence intervals for covariance parameters
	theInf=informationLgm(res)
	res$varBetaHat = list(beta=res$varBetaHat)
	names(res) = gsub("varBetaHat", "varParam", names(res))
	res$varParam$information = theInf$information

  
  res$summary = 	theInf$summary
  
  if(is.na(res$summary ['anisoAngleDegrees','ci0.05']) &
      !is.na(res$summary ['anisoAngleRadians','ci0.05']) ){
    ciCols = grep("^ci0\\.[[:digit:]]+$", colnames(res$summary ))
    res$summary ['anisoAngleDegrees',ciCols] =
        (360/(2*pi))*res$summary ['anisoAngleRadians',ciCols]
    res$summary ['anisoAngleDegrees','Estimated'] = 
        res$summary ['anisoAngleRadians','Estimated']
  }
  
  
  
	if(FALSE){
	for(Dvar in names(covariates)) {
		theLevels =levels(covariates[[Dvar]])[[1]]
		if(!is.null(nrow(theLevels))){
			for(D in 1:nrow(theLevels)) {
				rownames(res$summary) = gsub(
						paste("(factor)?(\\()?", Dvar, "(\\))?:?", 
								theLevels[D,1],"$",sep=""),
						paste(Dvar, ":",theLevels[D,2],sep=""), 
						rownames(res$summary))
			}
		}
	}
	}
  
  theOrder = c('sdNugget','sdSpatial', 'range', 'shape','anisoRatio',
      'anisoAngleRadians','anisoAngleDegrees', 'boxcox')  
  theOrder = na.omit(match(theOrder, rownames(res$summary)))
  notInOrder = (1:nrow(res$summary))[-theOrder]
  res$summary = res$summary[c(notInOrder,theOrder),]
  
  
  
	# if range is very big, it's probably in metres, convert to km
  	if(res$summary['range','estimate']>1000) {
  		logicalCol = names(res$summary) == "Estimated"
  		res$summary["range",!logicalCol] = 
  				res$summary["range",!logicalCol] /1000
  		rownames(res$summary) = gsub("^range$", "range/1000", 
  				rownames(res$summary))
  	}
	
  class(res) = c('lgm',class(res))    
	return(res)
}

)
AIC.lgm = function(object, ..., k = 2) {
  
  AIC(logLik(object, ..., k))
  
}


logLik.lgm = function(object, ...){
  res = object$opt$logL
  res = res[grep('^logL', names(res))]
  
  srows = rownames(object$summary)
  df = sum(
      object$summary[
          grep("^anisoAngle", srows, invert=TRUE),
          'Estimated']
      ) + any(
          object$summary[
              grep("^anisoAngle", srows, invert=FALSE),
              'Estimated']
  )
  attributes(res)$df = df
  attributes(res)$nobs = nrow(data.frame(object$data))
  class(res)= 'logLik'
  res
}

summary.lgm = function(object, ...) {
  object$summary
}