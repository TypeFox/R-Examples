profLlgm = function(fit,mc.cores=1, ...) {
	
	dots = list(...)
  fit$parameters = fillParam(fit$parameters)
	varying = intersect(names(dots), names(fit$parameters))

  for(D in varying)
    dots[[D]] = sort(dots[[D]])
  
  
	nonLinearParams = c('boxcox','shape','nugget','variance',
			'anisoAngleRadians','anisoAngleDegrees','anisoRatio','range')
	
	reEstimate = rownames(fit$summary)[
			fit$summary[,"Estimated"]
	]
	reEstimate = gsub("sdNugget", "nugget", reEstimate)
	reEstimate = gsub("sdSpatial", "variance", reEstimate)
	reEstimate = gsub("range \\(km\\)", "range", reEstimate)
  reEstimate = gsub("/1000", "", reEstimate)
  reEstimate = intersect(reEstimate, nonLinearParams)
	reEstimate = reEstimate[!reEstimate %in% varying]
  
  
	baseParams = fit$parameters
	baseParams = baseParams[names(baseParams)%in%
					nonLinearParams]
	
	baseParams=baseParams[!names(baseParams)%in% varying]
  if(any(reEstimate=='variance')){
    baseParams['nugget'] = baseParams['nugget'] / baseParams['variance']
  }
	baseParams=baseParams[names(baseParams) != 'variance']
	
	
  if(length(grep("^anisoAngle", varying))){
    reEstimate = grep("^anisoAngle", reEstimate,value=TRUE,invert=TRUE)
    baseParams = baseParams[
        grep("^anisoAngle", names(baseParams),value=TRUE,invert=TRUE)
    ]
  } 
  
  parValues = do.call(expand.grid, dots[varying])
  
  parList = apply(parValues,1,list)
  parList = lapply(parList, function(qq) c(unlist(qq), baseParams))
  
 
	if(mc.cores>1) {
		resL = parallel::mcmapply(
        likfitLgm, 
        param=parList,
        MoreArgs=list(
        data=fit$data, 
        formula=fit$model$formula,
        paramToEstimate=reEstimate,
        reml=fit$model$reml,
        coordinates=fit$data),
				mc.cores=mc.cores,
        SIMPLIFY=FALSE
		)
	} else {
    resL = mapply(
        likfitLgm, 
        param=parList,
        MoreArgs=list(
            data=fit$data, 
            formula=fit$model$formula,
            paramToEstimate=reEstimate,
            reml=fit$model$reml,
            coordinates=fit$data),
        SIMPLIFY=FALSE
    )
  }
  resL2 = lapply(resL, function(qq) qq$optim$logL)
  resL = simplify2array(resL2)

	if(length(varying)==1) {
    forNames =  c(
				dots[[varying]],
				fit$param[varying]
		)
		theorder = order(forNames)
		L = c(
        resL[grep("^m2logL",rownames(resL)),],
        fit$optim$logL[
            grep("^m2logL",
            names(fit$optim$logL))]
    )[theorder]
    dots[[varying]] = forNames = forNames[theorder]
    names(L) = paste(varying, forNames,sep="")
	} else {
		thedimnames=dots[varying]
		for(D in names(thedimnames)) {
			thedimnames[[D]] = paste(
					D, thedimnames[[D]],sep='_')
    }
		L = array(resL[grep("^m2logL",rownames(resL)),], 
				unlist(lapply(thedimnames, length)),
				dimnames=thedimnames)	
	} 
	

  Sprob = c(1, 0.999, 0.99, 0.95, 0.9, 0.8, 0.5, 0)
  Squant = qchisq(Sprob, df=length(varying))

  # get Scol from RColorBrewer
#	dput( rev(
#			RColorBrewer::brewer.pal(
#					length(Scontour)-1,
#					'Spectral')	))
  Scol=c("#3288BD", "#99D594", 
		  "#E6F598", "#FFFFBF", "#FEE08B", "#FC8D59", 
		  "#D53E4F")
  names(Scol) = as.character(Sprob[-length(Sprob)])
  
  res = list(
	  logL=-L/2,
	  legend = list(breaks=Sprob, col=Scol),
	  prob = Sprob,
	  col=Scol,
	  MLE=fit$param[varying],
	  basepars=baseParams
  )
  
  res$maxLogL = max(res$logL)
  res$breaks= res$maxLogL - Squant/2 	
  res$breaks[1] = min(c(res$breaks[2],min(res$logL)))-1

  res = c(dots[varying],res)
	
	if(length(varying)==1) {
    Skeep = seq(2, length(Sprob)-1)
    res$ci = cbind(
        prob=Sprob[Skeep],
        lower=NA, upper=NA)
    
    
    smaller = seq(1,which.max(res$logL))
		bigger = seq(which.max(res$logL), length(res$logL))
    
    # make the likelihood unimodal
    monotoneLik = res$logL
    
    if(length(smaller)>1) {
      for(D in 2:length(smaller)){
        monotoneLik[D] = max(monotoneLik[c(D-1,D)])
      }
      
      resCi = try(stats::spline(
          x=monotoneLik[smaller], 
          y=res[[1]][smaller],
          xout=res$breaks[Skeep],
          method= "hyman"), silent=TRUE)
      
      if(class(resCi)!='try-error')
  			res$ci[,'lower'] = resCi$y
  }

  if(sum(bigger)>1) {
    for(D in seq(length(smaller)+1, length(monotoneLik)-1)){
      monotoneLik[D] = min(monotoneLik[c(D-1,D)])
    }

    resCi = try(stats::spline(
        x=monotoneLik[bigger], 
        y=res[[1]][bigger], 
        xout=res$breaks[Skeep], method='hyman'),
    silent=TRUE)
    
    if(class(resCi)!='try-error')
      res$ci[,'upper']= resCi$y     
  }
		
  res$ci = cbind(
      res$ci, 
      logL = res$breaks[match(res$ci[,'prob'], res$prob)]
  )
		res$ciLong = na.omit(
				reshape(as.data.frame(res$ci[,c('prob','upper','lower')]), 
						direction="long",
						varying=list(par=c('upper','lower')),
						v.names='par',
						times = c('upper','lower'),
						timevar=c('direction'),
						idvar='prob')	
		)
		res$ciLong = rbind(
				res$ciLong,
				data.frame(prob=0,direction='upper',
						par=res$MLE)[
						colnames(res$ciLong)
				])
		res$ciLong$quantile = (1-res$ciLong$prob)/2
		res$ciLong[res$ciLong$direction=='upper','quantile'] =
				1 - res$ciLong[res$ciLong$direction=='upper','quantile'] 
		
		
		res$ciLong = res$ciLong[order(res$ciLong$par),]
	}

	res
}
