setClass("CommunityClimateStatistics",slots=c(obsStats="list",nullStats="list",deviations="list"))
setClass("CommunityClimateInput",slots=c(species_list_tinf="character", regional_pool_tinf="character", regional_pool_weights_tinf="numeric", climate_niches_tinf="data.frame", observed_climate_tobs="numeric"))

getnichesforspecieslist <- function(speciesnames, obs, climateaxes,verbose)
{
	if (verbose)
	{
		cat('\tBeginning niche sample...')
	}
	result <- vector("list",length=length(speciesnames))
	for (i in 1:length(speciesnames))
	{
		if (verbose) 
		{
			cat(sprintf(" %.2f ", (i/length(speciesnames))))
		}
		
		os <- subset(obs, obs$taxon==speciesnames[i])
		
		tr <- os[,climateaxes]
		
		result[[i]] <- tr
	}
	if (verbose)
	{
		cat('\t...done.\n')
	}
	
	return(result)
}

samplepointsfromnicheofspecieslist <- function(specieslist, obs, climateaxes, numsamplesperspecies,verbose)
{
	niches <- getnichesforspecieslist(speciesnames=specieslist, obs=obs, climateaxes=climateaxes,verbose)
	sampledpoints <- lapply(niches, function(x) { 
		indices <- sample(1:nrow(x),numsamplesperspecies,replace=TRUE) 
		return(x[indices,])
	})
	
	return(sampledpoints)
}

climatevolume <- function(data)
{
	volume <- NA
	data <- as.data.frame(na.omit(data))
	
	if (nrow(data) > 0)
	{
		inferredclimate <- apply(data, 2, median, na.rm=TRUE)
		
		distances <- rep(NA, nrow(data))
		for (i in 1:nrow(data))
		{
			vec <- data[i,] - inferredclimate
			distances[i] <- sqrt(sum(vec^2))
		}
		
		mediandistance <- median(distances, na.rm=TRUE)
	}
	else
	{
		mediandistance <- NA
	}
	
	return(mediandistance)
}

climatemismatch <- function(data, climate)
{	
	n <- names(climate)
	climate <- as.numeric(climate); names(climate) <- n
	
	inferredclimate <- apply(data, 2, median, na.rm=TRUE)
	
	mismatch <- inferredclimate - climate
	
	mismatchmagnitude <- sqrt(sum(mismatch^2))
	
	result <- list(inferredclimate=inferredclimate,observedclimate=climate, mismatch=mismatch, mismatchmagnitude=mismatchmagnitude)

	return(result)
}



computestatisticsfromsamples <- function(sampledpoints, observedclimate)
{
	nreplicates <- nrow(sampledpoints[[1]])
	nspecies <- length(sampledpoints)
	
	result <- vector("list", length=nreplicates)
	for (i in 1:nreplicates)
	{
		nichesample <- do.call("rbind",lapply(sampledpoints, function(x) { return(x[i,]) }))

		volume <- climatevolume(nichesample)
		
		climatemismatch <- climatemismatch(nichesample, observedclimate)
		
		result[[i]] <- list(nichesample=nichesample, volume=volume, climatemismatch=climatemismatch)
	}
	
	return(result)
}

averagestatisticsacrossreplicates <- function(stats)
{
	sampled <- do.call("rbind",lapply(stats, function(x) {z <- as.data.frame(x$nichesample); z$taxon <- paste(1:nrow(z)); return(z) }))
	meanNiches <- do.call("rbind",by(sampled[,1:(ncol(sampled)-1)],sampled[,ncol(sampled)], colMeans, na.rm=TRUE))

	
	volumeMagnitude <- mean(sapply(stats, function(x) { x$volume }), na.rm=TRUE)

	mismatchMagnitude <- mean(sapply(stats, function(x) { x$climatemismatch$mismatchmagnitude }), na.rm=TRUE)
	
	mismatchDirections <- rowMeans(sapply(stats, function(x) { x$climatemismatch$mismatch }), na.rm=TRUE)

	inferredClimate <- rowMeans(sapply(stats, function(x) { x$climatemismatch$inferredclimate }), na.rm=TRUE)
	
	observedClimate <- rowMeans(sapply(stats, function(x) { x$climatemismatch$observedclimate }), na.rm=TRUE)
	
	return(list(meanNiches= meanNiches, inferredClimate=inferredClimate, observedClimate=observedClimate, volumeMagnitude=volumeMagnitude, mismatchMagnitude=mismatchMagnitude, mismatchDirections=mismatchDirections))
}

pvalue_twotailed <- function(observed, nulldistribution)
{
	observed_rescaled <- observed - mean(nulldistribution, na.rm=TRUE)
	null_rescaled <- nulldistribution - mean(nulldistribution, na.rm=TRUE)

	result <- NULL
  if (!is.na(observed_rescaled))
  {
  	if (observed_rescaled < 0)
  	{
  		return(2*length(which(observed_rescaled > null_rescaled))/length(null_rescaled))
  	}
  	else
  	{
  		return(2*length(which(observed_rescaled < null_rescaled))/length(null_rescaled))
  	}
  }
  else
  {
    return(NA)
  }
}


reportdeviations <- function(observed, nulldistribution)
{
	ses <- as.numeric((observed - quantile(nulldistribution, 0.5, na.rm=TRUE) ) / ( quantile(nulldistribution,0.75, na.rm=TRUE) - quantile(nulldistribution,0.25, na.rm=TRUE) ))
	pvalue <- pvalue_twotailed(observed, nulldistribution)
	
	return(c(ses=ses, pvalue=pvalue))
}

assigndeviations <- function(observed, null)
{
	obs_volumeMagnitude = observed$volumeMagnitude
	null_volumeMagnitude <- sapply(null, function(x) { x$volumeMagnitude })
	deviation_volumeMagnitude <- reportdeviations(observed=obs_volumeMagnitude, nulldistribution=null_volumeMagnitude)
	
	obs_mismatchMagnitude = observed$mismatchMagnitude
	null_mismatchMagnitude = sapply(null, function(x) { x$mismatchMagnitude })
	deviation_mismatchMagnitude <- reportdeviations(observed= obs_mismatchMagnitude, nulldistribution=null_mismatchMagnitude)

	null_mismatchDirections <- sapply(null, function(x) {x$mismatchDirections})
	obs_mismatchDirections <- observed$mismatchDirections
	deviation_mismatchDirections <- matrix(NA, nrow=length(obs_mismatchDirections),ncol=2)
	dimnames(deviation_mismatchDirections) <- list(names(obs_mismatchDirections),c("ses","pvalue"))
	for (i in 1:nrow(null_mismatchDirections))
	{
		deviation_mismatchDirections[i,] <- reportdeviations(observed=obs_mismatchDirections[i], null_mismatchDirections[i,])
	}
	
	return(list(deviation_volumeMagnitude=deviation_volumeMagnitude, deviation_mismatchMagnitude=deviation_mismatchMagnitude, deviation_mismatchDirections=deviation_mismatchDirections))
}

communityclimate <- function(object, climateaxes=NULL, numreplicates=100, numsamplesperspecies=100,verbose=TRUE)
{		
	species_list_tinf = object@species_list_tinf
	regional_pool_tinf = object@ regional_pool_tinf
	regional_pool_weights_tinf=object@regional_pool_weights_tinf
	climate_niches_tinf = object@climate_niches_tinf
	observed_climate_tobs = object@observed_climate_tobs
	
	# use default set of axes from the observed climate
	if (is.null(climateaxes))
	{
		climateaxes = names(observed_climate_tobs)
	}
	# equiprobable draws from the regional pool as default
	if (length(regional_pool_weights_tinf) == 0)
	{
		regional_pool_weights_tinf = rep(1, length(regional_pool_tinf))
	}

	# check input
	if (!all(species_list_tinf %in% climate_niches_tinf$taxon))
	{
		warning(sprintf('Not all species in community found in climate niche data. Missing species will not be sampled:\n%s', paste(species_list_tinf[(!species_list_tinf %in% climate_niches_tinf$taxon)],collapse="\n")))
	}
	if (!all(regional_pool_tinf %in% climate_niches_tinf$taxon))
	{
		warning(sprintf('Not all species in regional pool found in climate niche data. Missing species will not be sampled:\n\t%s', paste(species_list_tinf[(! regional_pool_tinf %in% climate_niches_tinf$taxon)],collapse="\n\t")))
	}
	if (!all(climateaxes %in% names(climate_niches_tinf)))
	{
		stop('Requested axes not found in climate niche data')
	}
	if (!all(climateaxes %in% names(observed_climate_tobs)))
	{
		stop('Requested axes not found in observed climate data')
	}
	
	# subset the climate niches for speed
	climate_niches_tinf <- subset(climate_niches_tinf, climate_niches_tinf$taxon %in% union(species_list_tinf, regional_pool_tinf))

	# get observed statistics
	if (verbose)
	{
		cat('*** Observed statistics ***\n')
	}
	samples <- samplepointsfromnicheofspecieslist(specieslist=species_list_tinf, obs=climate_niches_tinf,
		climateaxes=climateaxes, numsamplesperspecies=numsamplesperspecies,verbose=verbose)

	stats <- computestatisticsfromsamples(samples, observed_climate_tobs[climateaxes])
	
	avgstats <- averagestatisticsacrossreplicates(stats)

	# get null distribution
	null_distribution <- vector("list",length=numreplicates)
	for (i in 1: numreplicates)
	{
		if (verbose)
		{
			cat(sprintf('*** Null statistics (%d / %d) ***\n', i, numreplicates))
		}
		# sample from the regional pool preserving richness
		nulllist <- regional_pool_tinf[sample(length(regional_pool_tinf),size=length(species_list_tinf),prob=regional_pool_weights_tinf,replace=TRUE)]
		
		samples_null <- samplepointsfromnicheofspecieslist(specieslist=nulllist, obs=climate_niches_tinf, 
		climateaxes=climateaxes, numsamplesperspecies=numsamplesperspecies,verbose=verbose)
		stats_null <- computestatisticsfromsamples(samples_null, observed_climate_tobs[climateaxes])
		avgstats_null <- averagestatisticsacrossreplicates(stats_null)
		
		null_distribution[[i]] <- avgstats_null
	}

	deviations <- assigndeviations(avgstats, null_distribution)

	result <- new("CommunityClimateStatistics")
	result@obsStats <- avgstats
	result@nullStats <- null_distribution
	result@deviations <- deviations

	return(result)
}

plotdensity <- function(obs, null,label)
{
	
	dev <- reportdeviations(obs, null)
	

	plot(density(null),xlab='',ylab='Probability',main=label,xlim=c(min(null, obs),max(null, obs)))
	rect(quantile(null,0.25),-1,quantile(null,0.75),1,col=gray(0.75,alpha=0.8),border=NA)
	abline(v=obs,col='red',lwd=2)	
	abline(v=quantile(null,0.5),lty=1)
	abline(v=quantile(null,0.25),lty=2)	
	abline(v=quantile(null,0.75),lty=2)
	legend('topleft',bty='n',legend=sprintf('SES=%.2f, p=%f',dev["ses"], dev["pvalue"]))
	

}

climatestatistics <- function(object)
{
	return(object@obsStats)
}


climatedeviations <- function(object)
{
	return(object@deviations)
}

summary.CommunityClimateStatistics <- function(object, ...)
{
	cat("*******************************************************************************\n")
	cat(sprintf('*** Community climate statistics - observed\n'))
	cat("*******************************************************************************\n")
	print(object@obsStats)
	
	cat("*******************************************************************************\n")
	cat(sprintf('*** Community climate - deviations (n=%d nulls)\n', length(object@nullStats)))
	cat("*******************************************************************************\n")
	print(climatedeviations(object))
}

summary.CommunityClimateInput <- function(object, ...)
{
	print(str(object))
}

plotcircle <- function(xcenter, ycenter, radius, lwd, ...)
{
	np <- 72
	
	theta <- seq(0,2*pi,length.out=np)
	
	xv <- xcenter + radius * cos (theta)
	yv <- ycenter + radius * sin (theta)
	
	lines(xv,yv, lwd=lwd, ...)
}

plotDeviation <- function(x, whichdeviation)
{
	if (whichdeviation=="Volume")
	{
		obs = x@obsStats$volumeMagnitude
		null <- sapply(x@nullStats, function(x) { x$volumeMagnitude })
		name <- expression(paste(Delta))
	}
	else if (whichdeviation=="Mismatch")
	{
		obs = x@obsStats$mismatchMagnitude
		null <- sapply(x@nullStats, function(x) { x$mismatchMagnitude })
		name <- expression(paste(Lambda))
	}
	else if (whichdeviation %in% row.names(x@deviations$deviation_mismatchDirections))
	{
		obs = x@obsStats$mismatchDirections[whichdeviation]
		null <- sapply(x@nullStats, function(x) { x$mismatchDirections[whichdeviation] })
		name <- paste("Lambda",whichdeviation)
	}
	else
	{
		stop('Not a valid axis name')
	}

	plotdensity(obs, null,name)
}

plot.CommunityClimateStatistics <- function(x, deviations=FALSE, axisnames=NULL, nnull=10, cex.axis=0.7, cex.nullpoints=0.3, cex.obspoints=0.5, cex.names=1.5, ...)
{	
	if (deviations)
	{
		axes <- row.names(x@deviations$deviation_mismatchDirections)
		numaxes <- 2+length(axes)
		par(mfrow=c(ceiling(numaxes/2),2))
		
		plotDeviation(x,"Volume")
		plotDeviation(x,"Mismatch")
		
		for (i in 1:length(axes))
		{
			plotDeviation(x, axes[i])
		}
		par(mfcol=c(1,1))
	}
	else
	{
		obsStats <- x@obsStats
		nullStats <- x@nullStats
		
		nvar <- ncol(obsStats$meanNiches)
		
		obsStats$meanNiches <- obsStats$meanNiches[complete.cases(obsStats$meanNiches),]
		
		climateVariables <- names(obsStats$observedClimate)
	
		if (is.null(axisnames))
		{
			axisnames = climateVariables
		}
		
		allpoints <- NULL
		allpoints <- rbind(allpoints, obsStats$meanNiches, obsStats$inferredClimate, obsStats$observedClimate)
		for (i in 1:length(nullStats))
		{
			allpoints <- rbind(allpoints, nullStats[[i]]$meanNiches, nullStats[[i]]$inferredClimate, nullStats[[i]]$observedClimate)
		}
		allpoints <- as.numeric(allpoints)
		

		xlim = c(min(allpoints, na.rm=TRUE), max(allpoints, na.rm=TRUE))
		ylim = c(min(allpoints, na.rm=TRUE), max(allpoints, na.rm=TRUE))

		if (length(nullStats) < nnull)
		{
			nnull = length(nullStats)
			warning('Number of nulls plotted was reduced to match data')
		}
	
		
		opar = par(no.readonly=TRUE)
		par(oma=c(1,1,1,1))	
		par(mfrow=c(nvar, nvar))
		par(mar=c(0,0,0,0))
		par(oma=c(2,2,2,2))
		par(mgp=c(2,0.8,0))
		for (i in 1:length(climateVariables))
		{
			for (j in 1:length(climateVariables))
			{
				whichx <- climateVariables[j]
				whichy <- climateVariables[i]
				
				if (i < j)
				{				
					# observed points
					avg_xv_present <- obsStats$meanNiches[,whichx]
					avg_yv_present <- obsStats$meanNiches[,whichy]
					avg_xv_inferred <- obsStats$inferredClimate[whichx]
					avg_yv_inferred <- obsStats$inferredClimate[whichy]
					avg_xv_observed <- obsStats$observedClimate[whichx]
					avg_yv_observed <- obsStats$observedClimate[whichy]
					
					xlmin <- min(avg_xv_present, na.rm=TRUE)
					xlmax <- max(avg_xv_present, na.rm=TRUE)
					ylmin <- min(avg_yv_present, na.rm=TRUE)
					ylmax <- max(avg_yv_present, na.rm=TRUE)
					
					sf <- 0.1
	
					plot(0, 0,type='n',axes=FALSE,xlim=xlim,ylim=ylim)
					
					if (abs(i-j)==1)
					{
						axis(side=1,labels=TRUE, cex.axis=cex.axis, tck=-.05)
						axis(side=2,labels=TRUE, cex.axis=cex.axis, tck=-.05)
					}				
	
					# null points
					for (k in 1:nnull)
					{
						xv_present <- nullStats[[k]]$meanNiches[,whichx]
						yv_present <- nullStats[[k]]$meanNiches[,whichy]
						xv_inferred <- nullStats[[k]]$inferredClimate[whichx]
						yv_inferred <- nullStats[[k]]$inferredClimate[whichy]
						xv_observed <- nullStats[[k]]$observedClimate[whichx]
						yv_observed <- nullStats[[k]]$observedClimate[whichy]
						
						points(xv_present, yv_present,col=rgb(0.5,0.5,0.5,0.75),pch=16,cex=cex.nullpoints)				
						segments(xv_inferred, yv_inferred, xv_observed, yv_observed,col=rgb(0.0,0.0,0.0,0.75),lwd= 2*cex.nullpoints)
						points(xv_inferred, yv_inferred,col=rgb(0.0,0.0,0.0,0.5),cex=cex.nullpoints*2,lwd=0.5,pch=16)
											
						plotcircle(nullStats[[k]]$inferredClimate[i], nullStats[[k]]$inferredClimate[j],nullStats[[k]]$volumeMagnitude,col=rgb(0.5,0.5,0.5,0.75),lwd=2*cex.nullpoints)
			
						
					}
	
					# observed points
					points(avg_xv_present, avg_yv_present,pch=16,cex=cex.obspoints,col=rgb(1,0,0,0.25))
					segments(avg_xv_inferred, avg_yv_inferred, avg_xv_observed, avg_yv_observed,col=rgb(1,0,0,0.75),lwd= 2*cex.obspoints)
					points(avg_xv_inferred, avg_yv_inferred,pch=16,col=rgb(1,0,0,0.75),cex=2*cex.obspoints,lwd=0.5)					
					plotcircle(obsStats$inferredClimate[i], obsStats$inferredClimate[j],obsStats$volumeMagnitude,col=rgb(1,0,0,0.75),lwd= 2*cex.obspoints)
					
					
					# observed climate
	        		points(avg_xv_observed, avg_yv_observed,pch=1,cex=2*cex.obspoints,col='black',lwd=2*cex.obspoints)
					
					box()
					
					
				}
				else if (i==j)
				{
					plot(0,0,type='n',axes=FALSE,xlab='',ylab='')
					text(0,0,axisnames[i],cex=cex.names)
				}
				else
				{
					
					plot(0,0,type='n',axes=FALSE,xlab='',ylab='')
				}
			}
		}
		
		st <- data.frame(delta=x@deviations$deviation_volumeMagnitude, lambda=x@deviations$deviation_mismatchMagnitude)
		
		title(main=substitute(paste(delta == dval, ", ",~ lambda == lval, sep=""), list(dval=round(st$delta,digits=3), lval=round(st$lambda,digits=3))),outer=TRUE)
		
		par(opar)
	}
}

plot.CommunityClimateInput <- function(x, climateaxes=NULL, axisnames=NULL, cex.community=0.5, cex.pool=0.25, pch.community=16, pch.pool=16, colors="rainbow", ...)
{
	climateniches = x@climate_niches_tinf
	localcommunity = x@species_list_tinf
	regionalpool = x@regional_pool_tinf
	
	if (is.null(climateaxes))
	{
		climateaxes <- setdiff(names(climateniches),"taxon")
	}
	if (is.null(axisnames))
	{
		axisnames = climateaxes
	}
	
	colorfun <- match.fun(colors)
	
	pairs(climateniches[,climateaxes],col=colorfun(1.5*length(regionalpool))[climateniches$taxon],cex=ifelse(climateniches$taxon %in% localcommunity, cex.community, cex.pool),pch=ifelse(climateniches$taxon %in% localcommunity, pch.community, pch.pool),labels=axisnames,lower.panel=NULL)	
}


generatedemodata <- function(num_regionalpool = 50, num_community = 5, num_occurrences = 40, num_climateaxes = 3, observed=0)
{	
	climateniches <- NULL
	for (i in 1:num_regionalpool)
	{
		randdata = NULL
		for (j in 1:num_climateaxes)
		{
			meanpos = runif(num_climateaxes,min=2,max=4)
			tcol = rnorm(num_occurrences, mean=meanpos[j] + runif(n=1,min=-2,max=2), sd=runif(1, 0.2,0.4))
			randdata <- cbind(randdata, tcol)
		}
		
		randdata <- as.data.frame(randdata)
		
		names(randdata) <- paste("ClimateAxis", 1:num_climateaxes, sep='')
		randdata$taxon = paste("Species", i, collapse='')
		
		climateniches <- rbind(climateniches, randdata)
	}
	climateniches $taxon <- factor(climateniches$taxon)
	
	nichepos <- do.call("rbind",by(climateniches[,1:num_climateaxes], climateniches$taxon, function(x) {cm <- colMeans(x); return(data.frame(pos=sqrt(sum(cm^2))))}))
	
	# select for species on the lower edge of the climate space
	localcommunity <- row.names(nichepos)[order(nichepos,decreasing=FALSE)[1:num_community]]
	regionalpool <- row.names(nichepos)
	
	# suppose the observed climate is at a certain position
	observedclimate <- rep(observed, num_climateaxes); names(observedclimate) <- paste("ClimateAxis", 1:num_climateaxes, sep='')
	
	demodata <- inputcommunitydata(
		localcommunity = localcommunity, 
		regionalpool = regionalpool, 
		climateniches = climateniches, 
		observedclimate = observedclimate)
	
	return(demodata)
}

inputcommunitydata <- function(localcommunity, regionalpool, regionalpoolweights=numeric(), climateniches, observedclimate)
{
	object <- new("CommunityClimateInput", 
		species_list_tinf = localcommunity, 
		regional_pool_tinf = regionalpool, 
		regional_pool_weights_tinf = regionalpoolweights,
		climate_niches_tinf = climateniches, 
		observed_climate_tobs = observedclimate)
	
	return(object)	
}