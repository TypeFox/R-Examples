loglikLgm = function(param,
		data, formula, coordinates=data,
		reml=TRUE, 
		minustwotimes=TRUE,
		moreParams=NULL) {

	# create 'covariates', and 'observations'
 	
	trend = formula
	
	if(class(trend)=="formula") {
		observations = all.vars(trend)[1]
		if(!any(names(data)==observations))
			warning("can't find observations ", observations, "in data")
		# frame first, so we see which row names have been omitted due to NA's
		# sometimes model.matrix removes row names
		covariates = model.frame(trend, data.frame(data))
		observations = covariates[,	observations]
		theRowNames = rownames(covariates)
		covariates = model.matrix(trend,covariates)
		rownames(covariates) = theRowNames
		theNA = which(!rownames(data.frame(data)) %in% theRowNames)
	} else if(!is.null(trend)){
		# observations must be a vector and covariates a matrix
		covariates= as.matrix(trend)
		theNA = apply(covariates, 1, function(qq) any(is.na(qq)))
		observations=data[!theNA]
		theNA = which(theNA)
	} else {
		theNA = NULL
	}
	
		if(length(grep("SpatialPoints", class(coordinates)))) {
			if(length(theNA))
				coordinates = coordinates[-theNA,]
		} else if(	class(coordinates) == "dist")	{
			if(length(theNA)) {
				coordinates = 
					as.matrix(coordinates)[-theNA,-theNA]
      } else {
        coordinates= as.matrix(coordinates)
      } 
    } else if(class(coordinates) == "matrix" |
        class(coordinates) == 'dsyMatrix')	{
      if(length(theNA)) {
        coordinates = coordinates[-theNA,-theNA]
      }
    } else {
				warning("coordinates must be a SpatialPoints object\n or matrix a dist object.  It's being assumed \n it's a matrix of coordinates")
				coordinates = SpatialPoints(coordinates)
				if(length(theNA))				
					coordinates = coordinates[-theNA,]
		}

		# sort out the parameters
		param = c(param, moreParams)
		
		names(param) = gsub("^var$", "variance", names(param))
		
		param=param[!is.na(param)]
	
		haveVariance = any(names(param)=="variance")
		haveNugget = any(names(param)=="nugget")
		if(haveNugget){
			if(param["nugget"] == 0) {
				haveNugget=FALSE
			}
		} else {
			param["nugget"] = 0
		}

		if(!haveVariance) {
			# variance will be estimated.  
		  # nugget is assumed to be nugget/variance
			param["variance"] = 1
		}
			
		# box cox transform
		jacobian = 0
		if(any(names(param)=="boxcox")) {
			
			
			if(abs(param["boxcox"]-  1 ) < 0.001) {
			 # boxcox close to 1, don't transform
				jacobian=0
				
			} else { # box cox is not one.
				jacobian = -2*(param["boxcox"]-1)* 
					sum(log(observations))	
			
				if(is.nan(jacobian))
					warning("boxcox shouldnt be used with negative data")

				if(abs(param["boxcox"])<0.001) {
					observations = log(observations) 
				} else if(abs(param["boxcox"]-1)>0.001) {
  				observations <- 
						((observations^param["boxcox"]) - 1)/
							param["boxcox"]
				}
			} # end boxcox param far from 1
		} # end have box cox
		
		
    obsCov = cbind(observations, covariates)
  
  Nobs = nrow(obsCov)
  Ncov = ncol(obsCov)-1
  Nrep = 1
  
  paramFull = fillParam(param)
  Ltype = c(ml=0, reml=1, mlFixed=2, remlFixed=3)
  Ltype = reml + 2*haveVariance

  if(class(coordinates)=='matrix'|
      class(coordinates) == 'dsyMatrix'){
    xcoord = as.vector(coordinates)
    ycoord = -99
    aniso=FALSE
  } else if(length(grep("^SpatialPoints", class(coordinates)))){
    xcoord=coordinates@coords[,1] 
    ycoord=coordinates@coords[,2]
    aniso=TRUE
  } else {
    warning('coordinates should be SpatialPoints or matrix')
    xcoord=ycoord=aniso=NULL
  }
  
  resultC = .C("maternLogL",
      xcoord=as.double(xcoord), 
      ycoord=as.double(ycoord),
      param=as.double(paramFull[
              c('nugget','variance','range',
                  'shape','anisoRatio', 'anisoAngleRadians')]),
      aniso=as.integer(aniso),
      obsCov = as.double(obsCov),
      N= as.integer(c(Nobs,Nrep,Ncov)),
      boxcox=as.double(-9.9),
      boxcoxType=as.integer(0),
      logL=as.double(rep(-9.9, Nrep)),
      totalVarHat=as.double(rep(-9.9, Nrep)),
      betaHat = as.double(rep(-9.9, Ncov)), 
      varBetaHat = as.double(rep(-9.9, Ncov* Ncov)),
      Ltype=as.integer(Ltype)
  )

  totalVarHat = resultC$totalVarHat
  betaHat = resultC$betaHat
  names(betaHat) = colnames(obsCov)[-1]
  
  varBetaHat = new("dsyMatrix", 
      Dim = as.integer(c(Ncov, Ncov)), 
      uplo="L",
      x=resultC$varBetaHat)
  dimnames(varBetaHat) = list(names(betaHat),names(betaHat))
  

    if(!haveVariance) {
    varBetaHat = varBetaHat*totalVarHat
    param[c("variance","nugget")] = 
        totalVarHat * param[c("variance","nugget")]
  }

  result = resultC$logL[1] + jacobian
  
  

if(minustwotimes) {
  names(result) = "minusTwoLogLik"
} else {
  result = -0.5*result
  
  names(result)="logLik"
} 

if(reml)
  names(result) = gsub("Lik$", "RestrictedLik", names(result))

attributes(result)$param = param
attributes(result)$totalVarHat = resultC$totalVarHat
attributes(result)$betaHat = betaHat
attributes(result)$varBetaHat = varBetaHat
attributes(result)$reml=reml
attributes(result)$jacobian = jacobian
attributes(result)$Ltype = as.integer(Ltype)
attributes(result)$Lorig = resultC$logL
attributes(result)$aniso = aniso
attributes(result)$determinants = resultC$obsCov[1:2]
result
}

 


likfitLgm = function(
		formula, data,
    paramToEstimate = c("range","nugget"),
    reml=TRUE, 
    coordinates=data,
		param=NULL,
		upper=NULL,
    lower=NULL, 
    parscale=NULL,
    verbose=FALSE) {

  param = param[!is.na(param)]
  coordinatesOrig = coordinates
  
  # check if model is isotropic
  # if it is coordinates will be a distance matrix
  # if not coordinates will be SpatialPoints
  aniso = as.logical(length(grep("^aniso", c(names(param), paramToEstimate)))) |
      any(abs(param['anisoRatio']-1) > 0.00001,na.rm=TRUE)
  if(aniso){
    if(is.matrix(coordinates)){
      if(ncol(coordinates)!= 2 | nrow(coordinates) != nrow(data))
        stop("anisotropic model requested but coordinates appears to be a distance matrix")
    }
    coordinates = SpatialPoints(coordinates)
    maxDist = dist(t(bbox(coordinates)))
  } else { # isotropic
    if(is.matrix(coordinates)){
      if(ncol(coordinates)== 2) {# assume the columns are x and y coordinates
        coordinates = dist(coordinates)
      } else {
        coordinates = as(coordinates, 'dsyMatrix')
      }
    }
    if(class(coordinates)=='dist')
      coordinates = as(as.matrix(coordinates), 'dsyMatrix')
    if(length(grep("^Spatial", class(coordinates))))
      coordinates = as(spDists(coordinates), 'dsyMatrix')
    maxDist = max(coordinates,na.rm=TRUE)
  }

  trend = formula
  theFactors = NULL
	if(class(trend)=="formula") {
    # convert input data to a model matrix
		data = data.frame(data)
		theNA = apply(
				data[,all.vars(trend),drop=FALSE],
				1, function(qq) any(is.na(qq)))
		noNA = !theNA
	
		theFactors = model.frame(trend, data[noNA,])
		whichFactors = unlist(lapply(theFactors, is.factor))
		theFactors = theFactors[,whichFactors,drop=FALSE]
		theFactors = lapply(theFactors, levels)
		theFactors = unlist(lapply(theFactors, function(qq) qq[1]))
		
		covariates = model.matrix(trend, data[noNA,])
		observations = all.vars(trend)[1]
		
		if(!any(names(data)==observations))
			warning("can't find observations ", observations, "in data")
		observations = data[noNA,observations]
	} else {
		# observations must be a vector and covariates a matrix
		trend = as.matrix(trend)
		theNA = is.na(data) | apply(trend, 1, 
					function(qq) any(is.na(qq))
				)
		noNA = !theNA
				
		observations=data[noNA]
		covariates=trend[noNA,,drop=FALSE]
	}

	if(any(theNA)) {
		if(length(grep("^SpatialPoints", class(coordinates)))) {
			coordinates = coordinates[noNA]	
		} else {
			coordinates = coordinates[noNA,noNA]
		}
	}
  

  # check for the variance parameter
  estimateVariance = TRUE
  if(any(paramToEstimate=="variance")) {
    # remove varinace, it's estimated by profile likelihood
    paramToEstimate = paramToEstimate[paramToEstimate != "variance"]
    param = param[names(param)!="variance"]
  } else {
    if(any(names(param)=="variance")){
      estimateVariance = FALSE
    }		
  }
  
  
  # make sure angle is radians
  paramToEstimate = gsub("anisoAngleDegrees","anisoAngleRadians", paramToEstimate)
  
  degToRad = function(par) {
  if(length(grep("anisoAngleDegrees", names(par)))){
    if(!length(grep("anisoAngleRadians", names(par))))
      par['anisoAngleRadians'] = 2*pi*par['anisoAngleDegrees']/360
    par = par[grep("anisoAngleDegrees",names(par),invert=TRUE)]
  }
   par
  }
  param = degToRad(param)
  lower = degToRad(lower)
  upper = degToRad(upper)
  parscale = degToRad(parscale)
  
  # parameter defaults
  lowerDefaults = c(
      nugget=0,
      range=maxDist/1000,
      anisoRatio=0.01,
      anisoAngleRadians=-pi/2,
      shape=0.1,boxcox=-1.5,variance=0)
  
  upperDefaults= c(
      nugget=Inf,
      range=10*maxDist,
      anisoRatio=100,
      anisoAngleRadians=pi/2,
      shape=4,boxcox=2.5,variance=Inf)
  
  paramDefaults = c(
      nugget=0,
      anisoRatio=1, 
      anisoAngleRadians=0,
      shape=1.5, boxcox=1,
      range=maxDist/10
  )
  
  if(any(names(paramToEstimate)=='nugget')) {
    paramDefaults['nugget'] = 1
  }
  
  parscaleDefaults = c(
      range=maxDist/5,
      nugget=0.05,
      boxcox=0.5,
      anisoAngleRadians=0.2,
      anisoRatio=1,
      variance=1,
      shape=0.2)
  
  ndepsDefault = c(
      range=0.01,
      nugget=0.05,
      boxcox=0.005,
      anisoAngleRadians=0.01,
      anisoRatio=0.01,
      variance=0.01,
      shape=0.01
      )
  
  # replace defaults with user supplied values
  paramDefaults[names(param)] = param
  parscaleDefaults[names(parscale)] = parscale
  lowerDefaults[names(lower)]=lower
  upperDefaults[names(upper)] = upper
  
  
  # don't let nugget start on the boundary
  if(any(paramToEstimate=='nugget') & paramDefaults['nugget'] == lowerDefaults['nugget']){
    paramDefaults['nugget'] = min(c(0.5, upperDefaults['nugget']))
  }
  
  startingParam = paramDefaults[paramToEstimate]
  names(startingParam) = paramToEstimate # fixes names lost when no starting value provided
  
  naStarting = is.na(startingParam)
  startingParam[naStarting]= paramDefaults[names(startingParam)[naStarting]]
  
  moreParams = paramDefaults[!names(paramDefaults) %in% paramToEstimate]

  allParams = c(startingParam, moreParams)
  allParams = fillParam(allParams)
  
  paramsForC = allParams[c('nugget','variance','range','shape',
          'anisoRatio','anisoAngleRadians','boxcox')]
  
  Sparam = names(paramsForC) %in% paramToEstimate
  names(Sparam) = names(paramsForC)
  paramToEstimate = names(Sparam)[Sparam]

  
  
  parOptions = cbind(
      lower=lowerDefaults[paramToEstimate], 
      upper=upperDefaults[paramToEstimate],
      parscale = parscaleDefaults[paramToEstimate],
      ndeps=ndepsDefault[paramToEstimate] # for derivatives
  )

  # parameters for l-bfgs-b
  forO = list(
      scalarF = c(
          fnscale=-1, 
          abstol=-1,
          reltol = -1,
          alpha=-1, beta=-1,gamma=-1,
          factr= 1e6,
          pgtol=0),
      scalarInt=c(
          trace=0,
          maxit=200,
          REPORT=1,
          type=-1,lmm=25,
          tmax=-1,temp=-1
      ),
      pars = parOptions[,c('lower','upper','parscale','ndeps'),drop=FALSE]
  )
  
  if(verbose){
    forO$scalarInt['trace']=6
    forO$scalarInt['REPORT']=200
  }
  
  forO$parsInt = rep(0, nrow(forO$pars))
  names(forO$parsInt) = rownames(forO$pars)
  forO$parsInt[
      forO$pars[,'lower'] != -Inf &
          forO$pars[,'upper'] != Inf
      ] = 2

  forO$parsInt[
          forO$pars[,'lower'] != -Inf & 
              forO$pars[,'upper'] == Inf
  ] = 1

  forO$parsInt[
      forO$pars[,'lower'] == -Inf & 
          forO$pars[,'upper'] != Inf
  ] = 3
  
  forOparsOrig = forO$pars
  
  forO$pars[!is.finite(forO$pars)]=-1
  forO$pars = c(forO$pars, rep(0.0, ncol(covariates)+ncol(covariates)^2))
  

  if(aniso){
    xcoord=coordinates@coords[,1] 
    ycoord=coordinates@coords[,2]
  } else {
    xcoord = as.vector(coordinates)
    ycoord = -99
  }
  
  obsCov = cbind(y1=observations, y2=0, y3=0, covariates)
  
  if(all(paramToEstimate=='variance') &
      param['nugget']==0){
    
    theL = loglikLgm(param,
        data=observations, 
        formula=covariates, 
        coordinates=coordinates,
        reml=reml)
    fromOptim = attributes(theL)
    
    result = list(
        optim = list(
            mle=fillParam(fromOptim$param),
            logL = c(m2logL = as.numeric(theL),
                logL = - as.numeric(theL)/2),
            totalVarHat = fromOptim$totalVarHat,
            message = 'numerical optimization not needed',
            options=NULL,
            detail = NULL
        ),
        betaHat = fromOptim$betaHat,
        varBetaHat =  
            fromOptim$varBetaHat
    )

} else {
  
  fromOptim = .C(
          "maternLogLOpt",
      start=paramsForC,    
      Sparam=as.integer(Sparam),
      obsCov=as.double(obsCov), 
      as.double(xcoord),
      as.double(ycoord),
      as.integer(aniso),
      N=as.integer(c(nrow(obsCov), 3, ncol(covariates))),
    Ltype=as.integer(reml+2*!estimateVariance),
    optInt = as.integer(forO$scalarInt),
    optF = as.double(forO$scalarF),
    betas=as.double(forO$pars),
    limType = as.integer(forO$parsInt),
    message=format(" ",width=80)
  )


  result = list(
      optim = list(
        mle=fromOptim$start,
        logL = c(m2logL = fromOptim$optF[1],
          logL = - fromOptim$optF[1]/2),
        totalVarHat = fromOptim$optF[2],
        boxcox = fromOptim$optF[3:5],
        determinants=fromOptim$optF[6:7],
        message = fromOptim$message,
        options=cbind(
            start=paramsForC[Sparam],
            opt = fromOptim$start[Sparam],
            parOptions[,c('parscale','lower','upper','ndeps')]),
        detail = fromOptim$optInt[1:3]      ),
    betaHat = fromOptim$betas[1:ncol(covariates)],
    varBetaHat =  
        new("dsyMatrix", 
            Dim = as.integer(rep(ncol(covariates),2)), 
            uplo="L",
            x=fromOptim$betas[
                seq(1+ncol(covariates), len=ncol(covariates)^2)]
        )
  )
  
# names(result$optim$boxcox) = c('param','sumLogY','twoLogJacobian')
# names(result$optim$determinants) = c('variance','reml')
 
 result$optim$options = cbind(result$optim$options,
      gradient=fromOptim$betas[
      seq(ncol(covariates)^2+ncol(covariates)+1, 
          len=sum(Sparam))
  ])
  
  names(result$optim$detail)  = c(
      "fail","fncount","grcount"
        )
 names(result$betaHat) = colnames(covariates)
 dimnames(result$varBetaHat) = list(
     names(result$betaHat),names(result$betaHat)
     )
     

   } # end not only variance to estimate

   result$parameters = fillParam(
       c(
       result$optim$mle, result$betaHat
   )
)
   
if(estimateVariance) {
   result$parameters[c('nugget', 'variance')] = 
       result$parameters[c('nugget', 'variance')] * 
       result$optim$totalVarHat  
 }
   
   names(result$optim$logL) = paste(
       names(result$optim$logL),
       c('.ml', '.reml')[reml+1],
       sep=''
   )
   
   result$data = cbind(
       data.frame(
       observations = observations,
       fitted=
           covariates %*% result$parameters[colnames(covariates)]
        ),
       data.frame(data)[noNA,]
   )
   


   
   if(abs(result$parameters["boxcox"]-1)>0.0001){ # boxcox not 1
       
     if(abs(result$parameters["boxcox"])<0.001) {
       result$data$obsBC = log(observations) 
     } else  { #boxcox far from 0 and 1
       result$data$obsBC <- 
           ((observations^result$parameters["boxcox"]) - 1)/
           result$parameters["boxcox"]
     }
   } else {
     result$data$obsBC <- 
         result$data$observations 
   }
   result$data$resid = result$data$obsBC - result$data$fitted
   

 if(length(grep("^Spatial", class(coordinatesOrig)))){
	 
	forDf = rep(NA, length(noNA))
	forDf[noNA] = seq(1, sum(noNA))

	theDf = result$data[forDf,] 

	result$data = SpatialPointsDataFrame(
       coords=SpatialPoints(coordinatesOrig),
       data=theDf)
   
   projection(result$data) = projection(coordinatesOrig)
 }
   
   
   result$model = list(reml=reml, baseline=theFactors)
	if(class(trend)=="formula") {
		result$model$formula = trend
    result$data[[all.vars(formula)[1]]] =
        result$data$observations
	} else {
		result$model$formula= names(trend)
	}

	
	if(length(result$model$baseline)){
		
		rparams = result$parameters
		
		baseParams =data.frame(
				var = names(result$model$baseline),
				base = result$model$baseline,
				pasted = paste(names(result$model$baseline), result$model$baseline, sep='')
				)
				
		for(D in which(! baseParams$pasted %in% names(rparams)) ) {
			
			sameFac = grep(paste("^",baseParams[D,'var'],sep=""),
					names(rparams))
			pseq = 1:length(rparams)
			if(length(sameFac)){
				
				minFac = min(sameFac)
				toAdd = 0
				names(toAdd) = baseParams[D,'pasted']
				rparams = c(
						rparams[pseq < minFac],
						toAdd,
						rparams[pseq >= minFac]
						)
			}
		}
		result$parameters = rparams
	}
	
	parameterTable = data.frame(estimate=result$parameters)
	rownames(parameterTable) =  names(result$parameters)

	parameterTable$stdErr = NA

  stdErr = sqrt(Matrix::diag(result$varBetaHat))
	# sometimes varBetaHat doesn't have names
	parameterTable[rownames(result$varBetaHat), "stdErr"] =
			stdErr

	thelims = c(0.005, 0.025, 0.05, 0.1)
	thelims = c(rbind(thelims, 1-thelims))

  theQ = qnorm(thelims)
	toadd = outer(parameterTable$stdErr, theQ)
	toadd = toadd + matrix(parameterTable$estimate, 
			ncol=length(thelims), nrow=dim(parameterTable)[1])
	colnames(toadd)= paste("ci", thelims, sep="")
	parameterTable = cbind(parameterTable, toadd)

  parameterTable[,"pval"] = pchisq(
			parameterTable$estimate^2  / parameterTable$stdErr^2,
			df=1,lower.tail=FALSE)
	
	parameterTable[,"Estimated"] = FALSE
	parameterTable[paramToEstimate,"Estimated"] = TRUE
	parameterTable[rownames(result$varBetaHat),"Estimated"] = TRUE
	if(estimateVariance)
		parameterTable["variance","Estimated"] = TRUE
	
	rownames(parameterTable)=gsub("^variance$", "sdSpatial", 
			rownames(parameterTable))	
	rownames(parameterTable)=gsub("^nugget$", "sdNugget", 
			rownames(parameterTable))	
	parameterTable[c("sdSpatial", "sdNugget"),"estimate"] = 
			sqrt(parameterTable[c("sdSpatial", "sdNugget"),"estimate"])

  #	dimnames(parameterTable) = unlist(lapply(dimnames(parameterTable),
#			function(qq) {
#				qq=gsub("_", "\\\\textunderscore ", qq)
#				qq=gsub("\\$", "\\\\textdollar ", qq)
#				qq=gsub("<", "\\\\textless ", qq)
#				qq=gsub(">", "\\\\textgreater ", qq)
#				qq
#			}
#	))

#  if(is.na(result$summary['anisoAngleDegrees','ci0.05']) &
#      !is.na(result$summary['anisoAngleRadians','ci0.05']) ){
#    ciCols = grep(colnames("^ci0\\.[[:digit:]]+$", result$summary))
#    result$summary['anisoAngleDegrees',ciCols] =
#      (360/(2*pi))*result$summary['anisoAngleRadians',ciCols]
#  }
  

	result$summary = as.data.frame(parameterTable)
  result$summary$Estimated = as.logical(result$summary$Estimated)

  result
}
