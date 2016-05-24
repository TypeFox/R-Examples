

informationLgm = function(fit, ...) {
	nonLinearParams = c('boxcox','shape','nugget','variance',
			'anisoAngleRadians','anisoRatio','range')
	
	reEstimate = rownames(fit$summary)[
			fit$summary[,"Estimated"]
	]
	reEstimate = gsub("sdNugget", "nugget", reEstimate)
	reEstimate = gsub("sdSpatial", "variance", reEstimate)
	reEstimate = gsub("range \\(km\\)", "range", reEstimate)
	reEstimate = intersect(reEstimate, nonLinearParams)
	
	baseParam = fit$parameters[reEstimate]
	moreParams = fit$parameters[
			!names(fit$parameters) %in% reEstimate &
					names(fit$parameters) %in% nonLinearParams]
	
	parToLog = c("nugget","variance","anisoRatio","range")
	parToLog = intersect(reEstimate, parToLog)
	
	if(!all(baseParam[parToLog]>0))
		return(list(summary=fit$summary,information=NULL))

  # get rid of NA's
  fit$data = na.omit(fit$data)
  
	aniso = length(grep("^aniso", reEstimate)) |
      any(abs(moreParams['anisoRatio']-1) > 0.00001,na.rm=TRUE)
  if(!aniso) {
    coordinates = as(spDists(fit$data), 'dsyMatrix')  
  } else{
    coordinates=fit$data
  }
  
	oneL = function(param) {
#    parToExp = grep("^log\\(", names(param))
#		param[parToExp] = exp(param[parToExp])
#    names(param) = gsub("^log\\(|\\)$", "", names(param))
    param[parToLog] = exp(param[parToLog])
		loglikLgm(param, 
        data=fit$data,
        formula=fit$model$formula,
        coordinates=coordinates,
        reml=fit$model$reml,
        moreParams=moreParams,
        minustwotimes=FALSE)
	}
	
	baseParam[parToLog] = log(baseParam[parToLog])
	

	hess = numDeriv::hessian(oneL, baseParam, ...)
	
	whichLogged = which(names(baseParam)%in% parToLog)
	names(baseParam)[whichLogged] = paste("log(", 
			names(baseParam)[whichLogged], ")",sep="")
	
	dimnames(hess) = list(names(baseParam),names(baseParam))
  
  infmat = -hess
  infmat = try(solve(infmat), silent=TRUE)
  if(class(infmat)=='try-error') {
#    stuff <<- fit
    return(list(summary=fit$summary,information=NULL, error=infmat))
  } 

  
  if(length(grep("anisoAngleRadians", colnames(infmat))) ) {
  anisoAngleDegrees = (360/(2*pi))*infmat[,'anisoAngleRadians']
  infmat = rbind(infmat, anisoAngleDegrees=anisoAngleDegrees)
  anisoAngleDegrees = c(anisoAngleDegrees,
      anisoAngleDegrees = (360/(2*pi))*
          as.numeric(anisoAngleDegrees['anisoAngleRadians'])
  )
  infmat = cbind(infmat, anisoAngleDegrees=anisoAngleDegrees)  
  }
  
	pvec = grep("^ci([[:digit:]]|\\.)+$", colnames(fit$summary),
			value=TRUE)
	pvec = as.numeric(gsub("^ci","", pvec))
	qvec = qnorm(pvec)
	names(qvec) = paste("ci", pvec, sep="")
	
	stdErr = diag(infmat)
  if(any(is.na(stdErr)))
    return(list(summary=fit$summary,information=infmat))
    
	if(!all(stdErr>0))
		return(list(summary=fit$summary,information=infmat))
	stdErr = sqrt(stdErr)
	
	toAdd = outer(stdErr, qvec, FUN="*")

	forSummary = baseParam[rownames(toAdd)] + toAdd
	expAdd = exp(forSummary[
					grep("^log\\(",rownames(forSummary))
					,])
	rownames(expAdd) = gsub("^log\\(|\\)$","",rownames(expAdd))
	forSummary = rbind(forSummary, expAdd)	

	summary = fit$summary
	
	if(any(rownames(forSummary)=='nugget'))
		forSummary = rbind(forSummary,
				sdNugget = sqrt(
						pmax(0,forSummary["nugget",]))
		)
	if(any(rownames(forSummary)=='variance'))
		forSummary = rbind(forSummary,
				sdSpatial = sqrt(
						pmax(0,forSummary["variance",]))
		)

  if(any(rownames(forSummary)=='anisoAngleRadians'))
    forSummary = rbind(forSummary,
        anisoAngleDegrees = forSummary['anisoAngleRadians',]*360/(2*pi)
    )
	
	inBoth = intersect(rownames(summary), rownames(forSummary))
	
	summary[inBoth,colnames(forSummary)] = 
			forSummary[inBoth,]
	
	list(summary=summary,information=infmat)
}
