
# data is a raster.  grid is ignored
setMethod("lgm", 
    signature("formula", "Raster", "ANY", "ANY"),
    function(
        formula, 
        data,  
        grid=NULL,
        covariates=NULL, ...) {
      
      dataCov = gm.dataRaster(
          formula, data,
          grid=raster(data),
          covariates=covariates,
          buffer=0)
      
      callGeneric(formula, 
          data=dataCov$data, 
          grid=dataCov$grid, 
          covariates=dataCov$covariates, ...)
    }
)


setMethod("lgm", 
    signature("formula", "data.frame", "Raster", "data.frame"), 
    function(formula,
        data, grid,
				covariates=NULL,
        shape=1, boxcox=1, nugget=0,
				  expPred=FALSE, nuggetInPrediction=TRUE,
          reml = TRUE, mc.cores=1,
          fixBoxcox=TRUE,
          fixNugget = FALSE,
          ...)
      {
  NN=NNmat(grid)
  csiz = xres(grid)
        
	Yvar = all.vars(formula)[1]
	if(!length(grep('[[:digit:]]$', Yvar))) {
	  allYvar = grep(paste("^", Yvar, "[[:digit:]]*$",sep=""), names(data), value=TRUE)
  } else {
		allYvar = Yvar
	}
  Yvec = as.matrix(data[,allYvar, drop=FALSE])
	
	if(!Yvar %in% colnames(data)) {
		data[,Yvar] = 0
	}
  Xmat = model.matrix(formula, data)		  

	if(nrow(Xmat) != ncell(grid))
		warning("dimensions of data and grid are not compatible")
  
  if(!fixNugget & (length(nugget)<2))
    nugget = NULL
  
  thel = loglikGmrf(Yvec=Yvec,Xmat=Xmat,
                    NN=NN, 
                    propNugget=nugget,
                    boxcox=boxcox, fixBoxcox=fixBoxcox,
                    shape=shape,mc.cores=mc.cores,
                    reml=reml, ...)
  mle = thel$mle   
  lArray = thel$mlArray
  lMat = thel$mlMat
  

  if (reml){
    chooseLike = 'logLreml'
  }else{
    chooseLike = 'logLml'
  }
  

  res = list(param = drop(mle))
  if(!is.null(lArray)){
    res$array = lArray
    res$profL = list()
		if(ncol(Yvec)>1){
			# multiple datasets
     		forProfRange = apply(lArray[,,chooseLike,,drop=FALSE],
         		c(1,4), which.max) 
				forProfNugget =  apply(lArray[,,chooseLike,,drop=FALSE],
         		c(1,2), which.max) 
				
				res$profL$range = lArray[,1,,,drop=FALSE]
				res$profL$nugget = lArray[,,,1,drop=FALSE]
				
				for(Dy in 1:dim(lArray)[1]){
					for(Drange in 1:dim(lArray)[4]){
						res$profL$range[Dy,1,,Drange] = 
							lArray[Dy,forProfRange[Dy,Drange],,Drange]
					}
					for(Dnugget in 1:dim(lArray)[2]){
						res$profL$nugget[Dy,Dnugget,,1] = 
							lArray[Dy,Dnugget,,forProfNugget[Dy,Dnugget] ]
					}
				}
				
				# subtract best L, range
				maxL = apply(res$profL$range[,,c('logLreml', 'logLml'),,drop=FALSE],1:3,max)
				deviance = array(maxL, c(dim(maxL), dim(lArray)[4]),
					dimnames = c(
							dimnames(maxL),
							dimnames(lArray)[4]
							)
					)
				dimnames(deviance)[[3]] = gsub("^logL","deviance",dimnames(deviance)[[3]])
				deviance = - deviance + res$profL$range[,,c('logLreml','logLml'),,drop=FALSE]

				res$profL$range = abind(res$profL$range, deviance, along=3)				
				res$profL$range = drop(res$profL$range)
				
				# subtract best L, nugget	
				# aperm to put the nugget variable last
		res$profL$nugget = aperm(
				res$profL$nugget[,,c('logLreml', 'logLml'),,drop=FALSE],
				c(1,4,3,2))
				maxL = apply(res$profL$nugget, 1:3, max)
		
		deviance = array(maxL, c(dim(maxL), dim(lArray)[2]),
				dimnames = c(
						dimnames(maxL),
						dimnames(lArray)[2]
				)
		)
		dimnames(deviance)[[3]] = gsub("^logL","deviance",dimnames(deviance)[[3]])
		
		deviance = - deviance + res$profL$nugget[,,c('logLreml','logLml'),,drop=FALSE]
		res$profL$nugget = drop(res$profL$nugget)				
				
				
		} else {
			# dimension 1 of the array is box-cox
    # nugget
    if(dim(lArray)[2]>1){ # have nugget
			# find best range and boxcox for this nugget
      best = apply(lArray[,,chooseLike,,drop=FALSE],
          2, which.max) 
      best=arrayInd(best, dim(lArray)[-(2:3)])
      res$profL$propNugget = NULL
      for(D in 1:nrow(best)){
        res$profL$propNugget = rbind(
            res$profL$propNugget,
            lArray[best[1], D, 
                c('propNugget',chooseLike),
                best[2]])
      }
    } # end have nugget
     
    if(dim(lArray)[4]>1){ # have oneminusar
      best = apply(lArray[,,chooseLike,,drop=FALSE],
          4, which.max) 
      best=arrayInd(best, dim(lArray)[-(3:4)])
      
      res$profL$oneminusar = NULL
      for(D in 1:nrow(best)){
        res$profL$oneminusar = rbind(
            res$profL$oneminusar,
            lArray[best[1], best[2], 
                c('oneminusar',chooseLike,'range'),
                D])
      }
    } # end have oneminusar
    res$profL$range = res$profL$oneminusar[,c(3,2)]
    res$profL$oneminusar = res$profL$oneminusar[,c(1,2)]
    
    if(all(dim(lArray)[c(2,4)]>1)){ # have both
      x=lArray[1,,'propNugget',1]
      orderx = order(x)
      y = lArray[1,1,'range',]
      ordery=order(y)
      
      res$profL$twoDim = list(
          x=lArray[1,orderx,'propNugget',1],
          y=lArray[1,1,'range',ordery],
          z=apply(
          lArray[,,chooseLike,,drop=FALSE],
          c(2,4), max, na.rm=TRUE)[orderx,ordery],
      oneminusar=lArray[1,1,'oneminusar',ordery]      
    )
 
    } # end have both
		} # end box-cox
  } # end lArray not null
	
	if(!is.null(thel$profileBoxCox)){
		res$profL$boxcox = thel$profileBoxCox
	}		
	
  
  res$data = data
  res$model$reml = reml
  res$model$trend = formula
 
  # summary table
	covInMle = grep("Se$", rownames(mle), value=TRUE)

  scovariates = gsub(
      'Se$','', covInMle
  )
  
  srownames = c('sdNugget','sdSpatial','range','shape')

  scolnames = c("estimate", "stdErr", "ci0.005", "ci0.995", "ci0.025", "ci0.975", 
      "ci0.05", "ci0.95", "ci0.1", "ci0.9", "pval", "Estimated")
	ress = list()
	for(D in 1:ncol(mle)){
    ress[[D]] = as.data.frame(
      matrix(
          NA,
          length(scovariates) + length(srownames),
          length(scolnames),
          dimnames = list(
              c(scovariates, srownames),
              scolnames
              )
          )
      )
   ress[[D]][c('sdNugget','sdSpatial','range','shape'),'estimate'] = 
       c(sqrt(mle[c('tausq','sigmasq'),D]),
        mle[c('range','optimalShape'),D]   
       )
   ress[[D]][c('sdNugget','sdSpatial','range','shape'),'Estimated'] =
       c(fixNugget, TRUE, TRUE, FALSE)
   ress[[D]][scovariates,'Estimated']  = TRUE   
   ress[[D]][scovariates,'estimate']  = mle[covInMle,D]   
   ress[[D]][scovariates,'stdErr']  = mle[covInMle,D]   
 } # for D
       
  if(length(ress)==1) ress = ress[[1]]
	  res$summary = ress

	
  return(res)
  

}
)

