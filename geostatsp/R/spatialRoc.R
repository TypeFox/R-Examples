
spatialRoc = function(fit, 
		rr=c(1,1.2, 1.5,2), truth, 
		border=NULL, random=FALSE){
	
	prob= 1-exp(seq(0,-12,len=200))
	
	if(any(names(fit)=='inla')){
		fit = list(fit)
	}
	
	
	breaks = c(-Inf, rr, Inf)
	
	if('raster'%in% names(truth))
		truth = truth$raster
	
	if(random) {
		truthVariable = 'random'
		truth = truth[[
				intersect(
						paste(truthVariable, c('',1:length(fit)), sep=''),
						names(truth)
				)
		]]
		truth = exp(truth)	
	} else { 
		truthVariable = 'relativeIntensity'
		
		truth = truth[[
				intersect(
						paste(truthVariable, c('',1:length(fit)), sep=''),
						names(truth)
				)
		]]
	}

	if(!is.null(border))
		truth = mask(truth, border)
	
	
	if(any(names(fit[[1]])=='raster')){
		
		template = fit[[1]]$raster[['space']]
		if(!is.null(border))
			template = mask(template, border)
		
		Srow= rowFromY(template, seq(ymax(truth), ymin(truth),
						len=nrow(truth)))
		Scol= colFromX(template, seq(xmin(truth), xmax(truth), 
						len=ncol(truth)))
		Scell = cellFromRowColCombine(template, rownr=Srow, colnr=Scol)
		Scell= values(template)[Scell]
		
		toKeep = which(!is.na(Scell) & !is.na(values(truth[[1]])))
	
		Sregion = na.omit(values(template))
		
	} else { # bym model
		Sregion = 1:length(fit[[1]]$data)
		regionRaster = rasterize(fit[[1]]$data, truth, 
				field=Sregion)
		if(!is.null(border))
			regionRaster = mask(regionRaster, border)
		Scell = values(regionRaster)
		toKeep = which(!is.na(Scell))
	
	}

	Nlevels = length(breaks)-1
	truthOver = cut(truth, breaks=breaks)
	names(truthOver) = gsub(
			paste("^",truthVariable, "$", sep=''),
			paste(truthVariable, 1, sep=''),
			names(truth))

	
	truthOver = as.data.frame(truthOver)[toKeep,,drop=FALSE]	
	truthOver$fitId = Scell[toKeep]
	truthOver = na.omit(truthOver)	
	
	SlevelsC = as.character(1:Nlevels)
	
	result = NULL
	for(Dsim in 1:length(fit)) {
		
		if(random) {
			marginals = fit[[Dsim]]$inla$marginals.bym 
			if(!length(marginals)) { # lgcp
				marginals = fit[[Dsim]]$inla$marginals.random$space 
			} else { # bym
				names(marginals) = paste("bym.", Sregion, sep="")
			}
		} else { 
			marginals = fit[[Dsim]]$inla$marginals.fitted.bym 
			if(!length(marginals)) {
				marginals = fit[[Dsim]]$inla$marginals.predict
			} else { # bym
				names(marginals) = paste("bym.", Sregion, sep="")
			}
			
		}
		
		
		truthFreqList = tapply(
				truthOver[,paste(truthVariable,Dsim,sep='')], 
				truthOver[,'fitId'],
				function(qq) 
					table(qq)[SlevelsC]
		)
		truthFreq = matrix(unlist(truthFreqList), ncol=Nlevels, byrow=TRUE,
				dimnames = list(names(truthFreqList), 
						paste("level", 1:Nlevels, sep="")))
		truthFreq[is.na(truthFreq)]=0
		truthCusum = t(apply(truthFreq, 1, cumsum))
		truthCusum = cbind(truthCusum, 
				fitId=as.numeric(rownames(truthFreq))
		)
		colnames(truthCusum) = gsub(paste("^level", Nlevels	, sep=''), "n",
				colnames(truthCusum))
		colnames(truthCusum) = gsub("^level", "under", colnames(truthCusum))
		
		underCols = grep("^under", colnames(truthCusum), value=TRUE)
		overCusum = NULL
		for(D in underCols){
			overCusum = cbind(overCusum, truthCusum[,'n'] - truthCusum[,D])
		}
		colnames(overCusum) = gsub("^under", "over", underCols)
		truthCusum = cbind(truthCusum, overCusum)
		
		x=NULL	

		for(Drr in rr ) {
			x = cbind(x,excProb(marginals, log(Drr)))
		}

		excCols = colnames(x) = paste('exc', 1:(Nlevels-1), sep='')		
		x = cbind(x, fitId = as.numeric(
						gsub("^[[:alpha:]]+\\.?", "", 
								names(marginals)))
		)
		x = merge(truthCusum, x, by='fitId')	
		colnames(x) = gsub("^level", "below", colnames(x))
		
		belowCols = grep("^under", colnames(x),value=TRUE)
		aboveCols = grep("^over", colnames(x),value=TRUE)
		
		freqMat = NULL
		for(Dprob in rev(prob)) {
			
			aboveP = x[,excCols]>=Dprob
			naboveP = !aboveP
			
			freqMat = abind::abind(cbind(
							fp = apply(aboveP * x[,belowCols,drop=FALSE], 2, sum,na.rm=TRUE), 
							tp = apply(aboveP * x[,aboveCols,drop=FALSE],2,sum,na.rm=TRUE),
							fn = apply(naboveP * x[,aboveCols,drop=FALSE],2,sum,na.rm=TRUE),	
							tn = apply(naboveP * x[,belowCols,drop=FALSE],2,sum,na.rm=TRUE)	
					), freqMat, along=3)
		}

    
		resD=
				abind::abind(
						sens = (freqMat[, 'tp', ,drop=FALSE] / 
                    (freqMat[,'tp',,drop=FALSE] + freqMat[,'fn',,drop=FALSE])),
						onemspec = (1-freqMat[, 'tn', ,drop=FALSE] / 
                    (freqMat[,'fp',,drop=FALSE] + freqMat[,'tn',,drop=FALSE])),
						along=4)

    dimnames(resD)[[3]]=as.character(prob)
		dimnames(resD)[[1]] = paste("exc", rr, sep='')
		
		result = abind::abind(result, resD, along=length(dim(resD))+1)
		
	} # end loop through fits
		dimnames(result)[[length(dim(result))]] = 
				paste('sim', 1:length(fit),sep='') 
    if(length(fit)>1) {
      result = abind::abind(result,
				mean=apply(result, seq(1,length(dim(result))-1), 
						mean),
				along=length(dim(result))
		)	
	}
	result = drop(result)

    result
}
