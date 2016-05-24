lcOneRow = function(thisrow, idxCol=NULL) {
	thisrow = thisrow[!is.na(thisrow)]
	if(length(thisrow)) {
		thisrow = sapply(thisrow, function(qq) list(list(weight=qq)))
		for(D  in idxCol)
			thisrow[[D]] = list(
				weight=1, 
				idx=thisrow[[D]]$weight
		)
		for(D in names(thisrow))
			thisrow[[D]] = thisrow[D]
		names(thisrow) = paste("v", 1:length(thisrow), sep="")
	}
	thisrow
}

setGeneric('glgm', 
		function(
				formula, data, grid, 
				covariates=NULL, 
				...) {
			standardGeneric("glgm")
		}
)

 # sort out formula
# null formula
setMethod("glgm", 
		signature("NULL"), 
    function(formula=NULL, data, grid, 
        covariates=NULL, ...) {
      formula =  1 
      callGeneric(formula, data, grid, covariates, ...)
    }
)


setMethod("glgm", 
		signature("numeric"),  
    function(formula, data, grid, 
        covariates=NULL, ...) {

      formula = names(data)[formula]
      callGeneric(formula, data, grid, covariates, ...)
    }
)

# change character to formula
setMethod("glgm", 
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

setMethod("glgm", 
		signature("formula", "ANY", "numeric", "ANY"),
    function(formula, data, grid, covariates=NULL, ...) {
      grid = squareRaster(data, grid)
      callGeneric(formula, data, grid, covariates, ...)
    }
)



# extrat covariates for data, convert covariates to a stack
setMethod("glgm", 
		signature("formula", "Raster", "Raster", "ANY"),
    function(
        formula, 
        data,  
        grid,
        covariates=NULL,
        buffer=0,
        ...) {
      
      dataCov = gm.dataRaster(
          formula, data,
          grid,
          covariates,
          buffer)
      
      callGeneric(
			  formula = dataCov$formula, 
          data = dataCov$data, 
		  grid = dataCov$grid, 
          covariates = dataCov$covariates, ...)
    }
)


setMethod("glgm", 
				signature("formula", "Spatial", "Raster", "ANY"),
        function(formula, 
            data, grid, 
            covariates=NULL, 
            buffer=0,...) {
          
          dataCov = gm.dataSpatial(
              formula, data, 
              grid, covariates, buffer)
          callGeneric(formula, 
              data=dataCov$data@data, 
							grid=dataCov$grid, 
              covariates=dataCov$covariates, ...)
        }
    )

#################
#### the real work
##################

setMethod("glgm", 
		signature("formula", "data.frame", "Raster", "data.frame"), 
		function(formula, data,  grid, 
				covariates=NULL, 
				shape=1, priorCI=NULL, 
				mesh=FALSE,...) {

    if(!any(names(grid)=='space'))
		warning("grid must have a layer called space with inla cell ID's")

	if(!all(all.vars(formula)%in% names(data)))
		warning("some covariates seem to be missing: formula ", paste(all.vars(formula), collapse=" "), ", data: ", paste(names(data), collapse=" "))

    cells = trim(grid[['space']])
	firstCell = values(cells)[1]
	cellDim = dim(cells)[1:2]
	# first cell = 2 * buffer^2 + ncolSmall * buffer + buffer
	# buffer = -(nrowSmall+1) + sqrt (  (nrowSmall+1)^2 + 8 firstCell / 4
	buffer = (-(cellDim[1]+1) + sqrt(  (cellDim[1]+1)^2 + 8* (firstCell-1) ))/4
	# data, cells, and covariates must have varilable called 'space'		
	# values of cells must be index numbers, and cells shouldnt include the buffer		
	thedots = list(...)
			
	# priors for spatial standard deviation and nugget std dev.
	sdNames = unique(c("sd",grep("^sd", names(priorCI), value=TRUE)))
	# if model is Gaussian, look for prior for sdNugget
	if(!any(names(thedots)=="family")) {
		thedots$family =  "gaussian"
	}
	if(thedots$family=="gaussian") {
		sdNames = unique(c(sdNames, "sdNugget"))
	}
	if(thedots$family=="gamma") {
		sdNames = unique(c(sdNames, "gammaShape"))
	}
	
  # list of prior distributions
      if(any(names(priorCI)=='distributions')){
        priorDistributions = priorCI$distributions
      } else {
        priorDistributions = list()
      }
  
	# priors for sd's (and precisions) 
	
	precPrior=list()
  
	for(Dsd in sdNames) {
    Dprec = gsub("^sd","precision",Dsd)
    if(any(names(priorDistributions)==Dprec)) {
    # distribution supplied
      
    if('scale' %in% names(priorDistributions[[Dprec]])){
      priorDistributions[[Dprec]]['rate'] = 1/ priorDistributions[[Dprec]]['scale']
    }
    
      if(all(c('shape','rate') %in% names(priorDistributions[[Dprec]]))) {

        precPrior[[Dsd]] = list(
						param = c(
            shape=as.numeric(priorDistributions[[Dprec]]['shape']), 
            rate=as.numeric(priorDistributions[[Dprec]]['rate'])),
				prior = 'loggamma')

      } else {
        precPrior[[Dsd]] = list(param=c(
            shape=priorDistributions[[Dprec]][1],
            rate=priorDistributions[[Dprec]][2]),
				prior = 'loggamma')
      }
    } else if(any(names(priorCI)==Dsd)) {
      # find distribution from interval supplied

		
		# if of length 1, it's pc prior u with alpha = 0.05
		if(length(priorCI[[Dsd]])==1){
			priorCI[[Dsd]] = c(
					u=as.numeric(priorCI[[Dsd]]),
					alpha = 0.05
					)
		}

		if(!length(names(priorCI[[Dsd]])))
			names(priorCI[[Dsd]]) = c('lower','upper')

		if(all(c('u','alpha') %in% names(priorCI[[Dsd]]))) {
			# pc priors
			precPrior[[Dsd]] = list(
					param=priorCI[[Dsd]],
					prior = 'pc.prec')
		} else {
			# gamma prior
		
    obj1 = sort(priorCI[[Dsd]]^-2)
		cifun = function(pars) {
				theci = 	pgamma(obj1, shape=pars[1], 
						rate=pars[2],log.p=T)
				
				(log(0.025) - theci[1])^2 +
				(2*(log(0.975) - theci[2]))^2		

			}
	
			
			
		precPrior2=optim(c(.5,.5/mean(obj1)), cifun, 
				lower=c(0.000001,0.0000001),method="L-BFGS-B")
		names(precPrior2$par) = c("shape","rate")

		precPrior[[Dsd]] = list(
				param = precPrior2$par,
				prior = 'loggamma')
		
				
 		#pgamma(obj1, shape= precPrior["shape"], rate=precPrior["rate"],log.p=F)
		#pgamma(obj1, shape= precPrior["shape"], rate=precPrior["rate"],log.p=T)
		#log(c(0.025, 0.975))
		#precPrior2
		#pgamma(obj1, shape=precPrior["shape"], rate=precPrior["rate"],log.p=T)
		#log(c(0.025, 0.975)) 
 		#1/sqrt(qgamma(c(0.975,0.025), shape=precPrior["shape"], rate=precPrior["rate"]))
		#priorCI$sd
		
		} # end gamma prior
		} else { # no prior supplied
			# default prior
			precPrior[[Dsd]] = list(
					param = c(shape=0.01, rate=0.01),
					prior = 'loggamma')
		}
	}
		
  if(any(names(priorDistributions)=='range')) {
    # distribution supplied
  
  if('scale' %in% names(priorDistributions[['range']])){
    priorDistributions[['range']]['rate'] = 1/ priorDistributions[['range']]['scale']
  }
  
  
    if(all(c('shape','rate') %in% names(priorDistributions$range))) {
      
      ratePrior = c(
          shape=as.numeric(priorDistributions$range['shape']), 
          rate=as.numeric(priorDistributions$range['rate']*xres(cells))
      )
      
    } else {
      ratePrior = c(
          shape=priorDistributions$range[1],
          rate = priorDistributions$range[2]*xres(cells)
      )
    }
  } else if("range" %in% names(priorCI)) {
		if(priorCI$range[1] < xres(cells)/4) {
			priorCI$range[1] = xres(cells)/4
			warning("lower bound of range CI too small, setting it to 1/4 cell size")
			
		}
		
		# rang parameter, in terms of cells, not km.
		obj1=sort(priorCI$range/xres(cells))
		
		cifun = function(pars) {
			theci = 		pgamma(obj1, shape=pars[1], rate=pars[2], log.p=T)
			
			(theci[1] - log(0.025))^2 +
					(theci[2] - log(0.925))^2 
		}
			
		ratePrior2=optim(c(2,2/mean(obj1)), cifun, 
				lower=c(0.001,0.001),method="L-BFGS-B")
		ratePrior = ratePrior2$par
		names(ratePrior ) = c("shape","rate")
	} else {
		ratePrior = c(shape=0.01, rate=0.01)
	}

	# prior for gamma shape
	# log-normal, priorCI is 4 standard deviations
	if("gammaShape" %in% names(priorCI)) {
		gammaShapePrior  = list(
				prior='gaussian',
				param=c(
						mean=as.numeric(mean(log(priorCI$gammaShape))),
						precision = as.numeric(abs(diff(log(priorCI$gammaShape)))[1]/4)^(-2)
						)
				)
	} else {
		gammaShapePrior = NULL
	}


  spaceFormula = paste(".~.+ f(space, model='matern2d', ",
				"nrow=", nrow(cells)+2*buffer, 
				", ncol=", ncol(cells)+2*buffer,
				", nu=", shape, 
				", hyper = list(",
				 "range=list( param=c(",
				      paste(ratePrior, collapse=","),
				"), prior='loggamma'),",
				"prec=list( param=c(",
				paste(precPrior$sd$param, collapse=","),
				"), prior='",precPrior$sd$prior,"')",
				" ) )", sep=""
			)
	
	formula = update.formula(formula,	as.formula(spaceFormula))

  # sort out factors
	thevars = rownames(attributes(terms(formula))$factors)
	thevars = grep("^factor\\(", thevars, value=TRUE)
	varsInData = apply(data, 2, is.factor)
	varsInData = names(data)[varsInData]
	thevars = c(varsInData, thevars)
	
	if(length(thevars)){
		thevars = gsub("^factor\\(|\\)", "", thevars)
		# loop through factors
		for(D in thevars){
			# biggest category is baseline
			thetable = table(data[,D])
			thebase = names(sort(thetable,decreasing=TRUE))[1]
			newLevels = unique(c(thebase, levels(factor(data[,D]))))
			data[,D] = factor(data[,D], levels=newLevels)
			covariates[,D] = factor(covariates[,D],
					levels=levels(data[,D]))
			
		}
	}
	theFactors = thevars

	
	# create linear combinations object for prediction.
	# create formula, strip out left variable and f(...) terms
 	formulaForLincombs = unlist(strsplit(as.character(formula), "~"))
	formulaForLincombs = formulaForLincombs[length(formulaForLincombs)]
	formulaForLincombs =
			gsub("\\+?[[:space:]]*f\\([[:print:]]*\\)[[:space:]]?($|\\+)", "+", formulaForLincombs)
	# strip out offsets
formulaForLincombs =
		gsub("\\+?[[:space:]]*offset\\([[:print:]]*\\)[[:space:]]?($|\\+)", "+", formulaForLincombs)

	# convert multiple + to a single +
	formulaForLincombs = gsub(
			"\\+[[:space:]]?\\+([[:space:]]?\\+)?", "+",
			formulaForLincombs)
	# strip out trailing +
formulaForLincombs = gsub("\\+[[:space:]]?$", "", formulaForLincombs)

	# if we have covariates in the formula and in the data
	if(nchar(formulaForLincombs) & nrow(covariates) ) {

		formulaForLincombs=as.formula(paste("~", formulaForLincombs))
	
	
		# variables in the model but not in prediction rasters
		thevars = rownames(attributes(terms(formulaForLincombs))$factors)
		thevars = gsub("^factor\\(|\\)", "", thevars)
		varsInPredict = thevars[thevars %in% names(covariates)]
		cantPredict = thevars[! thevars %in% names(covariates)]
		theFactors2 = grep("factor\\([[:print:]]*\\)", cantPredict)
		if(length(theFactors2)) {
			temp = cantPredict
			cantPredict = cantPredict[-theFactors2]
			theFactorsInFormula = temp[theFactors2]
		}
		if(length(cantPredict)){
				covariates[,cantPredict]= 0
		}
		covariates = covariates[,c("space", thevars),drop=FALSE]
		lincombMat = model.matrix(update.formula(
						formulaForLincombs, ~.+space),
					covariates, na.action=NULL)
		lincombMat[lincombMat==0] = NA
		
		spaceCol = grep("^space$", colnames(lincombMat), value=TRUE, ignore.case=TRUE)
	
		thelincombs <- apply(lincombMat, 1, lcOneRow, idxCol=spaceCol)
		names(thelincombs) = paste("c", lincombMat[,spaceCol],sep="")
		
	} else { # no covariates or no INLA
		thelincombs=list()	
		for(D in 1:ncell(cells)) {
			thelincombs[[D]] = list(
					list(
							"(Intercept)"=list(weight=1)
					),
					list(
							space=list(weight=1, idx=values(cells)[D])
					)
			)
		}
		names(thelincombs) = paste("c", values(cells),sep="")
		}

	# get rid of observations with NA's in covariates
	allVars = all.vars(formula)
	theNA = apply(data[,allVars], 1, function(qq) any(is.na(qq)))
 
	data = data[!theNA,]
	if(any(names(thedots)=='Ntrials'))
		thedots$Ntrials = thedots$Ntrials[!theNA]

	forInla = thedots
	forInla$lincomb = c(thelincombs, forInla$lincomb)
	forInla$data = data
	forInla$formula = formula
	
	

		
	# if model is gaussian, add prior for nugget
	if(!is.null(precPrior$sdNugget)) {
		forInla$control.family$hyper$prec =
				list(prior=precPrior$sdNugget$prior,
						param=precPrior$sdNugget$params
				) 
	}
	if(!is.null(gammaShapePrior)) {
		forInla$control.family$hyper$prec =
				gammaShapePrior 
	}
	
	
	# get rid of some elements of forInla that aren't required
	forInla = forInla[grep("^buffer$", names(forInla), invert=TRUE)]

#	return(forInla)
 

	if(requireNamespace("INLA", quietly=TRUE)) {
		inlaResult = do.call(INLA::inla, forInla) 
	} else {
		inlaResult = 
			list(logfile="INLA is not installed. \n install splines, numDeriv, Rgraphviz, graph,\n fields, rgl, mvtnorm, multicore, pixmap,\n splancs, orthopolynom \n then see www.r-inla.org")
	}

	
	if(all(names(inlaResult)=="logfile"))
		return(c(forInla, inlares=inlaResult))
	

 	# parameter priors for result


	params = list(
			range = list(userPriorCI=priorCI$range, 
					priorCI = 
							xres(cells)*
							qgamma(c(0.025,0.975), 
									shape=ratePrior["shape"], 
									rate=ratePrior["rate"]),
					priorCIcells = 
							qgamma(c(0.975,0.025), 
									shape=ratePrior["shape"], 
									rate=ratePrior["rate"]),
					params.intern = ratePrior))
	rangeLim = 	qgamma(c(0.001,0.999), 
			shape=ratePrior["shape"], 
			rate=ratePrior["rate"])
	rangeLim = xres(cells) *rangeLim
	rangeSeq = seq(min(rangeLim), max(rangeLim), len=1000)
	rangeSeqCells = rangeSeq/xres(cells)
	params$range$prior=cbind(
			x=rangeSeq,
			y=dgamma(rangeSeqCells, shape=ratePrior["shape"], 
					rate=ratePrior["rate"])  / xres(cells)
	)

  for(Dsd in names(precPrior)) {
		if(precPrior[[Dsd]]$prior == 'loggamma'){
		params[[Dsd]] = list(userPriorCI=priorCI[[Dsd]], 
			priorCI = 1/sqrt(
				qgamma(c(0.975,0.025), 
						shape=precPrior[[Dsd]]$param["shape"], 
						rate=precPrior[[Dsd]]$param["rate"])),
					params.intern=precPrior[[Dsd]])
	
	precLim = 	qgamma(c(0.999,0.001), 
			shape=precPrior[[Dsd]]$param["shape"], 
			rate=precPrior[[Dsd]]$param["rate"])
	sdLim = 1/sqrt(precLim)
	sdSeq = seq(min(sdLim), max(sdLim), len=1000)
	precSeq = sdSeq^(-2)
	params[[Dsd]]$prior=cbind(
			x=sdSeq,
			y=dgamma(precSeq, shape=precPrior[[Dsd]]$param["shape"], 
					rate=precPrior[[Dsd]]$param["rate"]) *2* (precSeq)^(3/2) 
	)
	} else { # pc prior
	
		params[[Dsd]] = list(userPriorCI=priorCI[[Dsd]], 
				priorCI = 1/sqrt(
						INLA::inla.pc.qprec(c(0.975,0.025),  
								u = precPrior[[Dsd]]$param['u'], 
								alpha = precPrior[[Dsd]]$param['alpha'])
				),
		params.intern=precPrior[[Dsd]]$param)
		
		precLim = INLA::inla.pc.qprec(c(0.999,0.001),  
						u = precPrior[[Dsd]]$param['u'], 
						alpha = precPrior[[Dsd]]$param['alpha'])
		sdLim = 1/sqrt(precLim)
		sdSeq = seq(min(sdLim), max(sdLim), len=1000)
		precSeq = sdSeq^(-2)
		params[[Dsd]]$prior=cbind(
				x=sdSeq,
				y=INLA::inla.pc.dprec(precSeq, 
						u = precPrior[[Dsd]]$param['u'], 
						alpha = precPrior[[Dsd]]$param['alpha']
					) * 2 * (precSeq)^(3/2) 
		)
		
	}
	}

	if(!is.null(gammaShapePrior)) {
		
		paramsGammaShape = 	c(
				gammaShapePrior$param["mean"], 
				sd=as.numeric(1/sqrt(gammaShapePrior$param["precision"]))
		)
		
		xLim = sort(exp(-qnorm(
								c(0.999,0.001), 
								mean=paramsGammaShape["mean"], 
								sd=paramsGammaShape["sd"])
						))
		xSeq  = seq(xLim[1], xLim[2], len=1000)
		
		
		params[['gammaShape']] = list(
				userPriorCI = priorCI[['gammaShape']],
				priorCI = sort(exp(-qnorm(c(0.975,0.025), 
								mean=paramsGammaShape["mean"], 
								sd=paramsGammaShape["sd"]))
				),
				params.intern=gammaShapePrior$param,
				params = paramsGammaShape,
				distribution = 'lognormal',
				prior = cbind(
						x=xSeq, 
						y = stats::dlnorm(xSeq, meanlog = paramsGammaShape['mean'],sdlog = paramsGammaShape['sd'])
						)
				)
	}
	
	
	# random into raster
# E exp(random)

if("summary.random" %in% names(inlaResult)) {

temp=unlist(
		lapply(inlaResult$marginals.random$space, function(qq) {
					sum(
							exp(qq[,"x"])*c(0,diff(qq[,"x"]))*qq[,"y"]	
					)
				})
)
inlaResult$summary.random[['space']][,"exp"] = temp



		forRast = 	as.matrix(inlaResult$summary.random[["space"]][values(cells),])
		resRasterRandom = 
				brick(extent(cells), nrows=nrow(cells),
						ncols=ncol(cells), crs=projection(cells),
						nl=dim(forRast)[2])
		names(resRasterRandom) = 
				paste("random.", colnames(forRast),sep="")
		
		values(resRasterRandom) = as.vector(forRast)
		
	} else {
		return(list(inla=inlaResult, parameters=params))
	}

	inlaResult$marginals.random$space = inlaResult$marginals.random$space[values(cells)]
	
	# E exp(lincombs)
	temp=unlist(
			lapply(inlaResult$marginals.lincomb.derived, function(qq) {
						sum(
								exp(qq[,"x"])*c(0,diff(qq[,"x"]))*qq[,"y"]	
						)
					})
	)
	inlaResult$summary.lincomb.derived[,"exp"] = temp
	
	# E inv logit(lincombs)
	if(length(grep("logit",inlaResult$misc$linkfunctions$names))) {
		temp=unlist(
				lapply(inlaResult$marginals.lincomb.derived, function(qq) {
							eqqx = exp(qq[,"x"])
							sum(
								eqqx/(1+eqqx)*c(0,diff(qq[,"x"]))*qq[,"y"]	
							)
						})
		)
		inlaResult$summary.lincomb.derived[,"invlogit"] = temp		
	}
 

	# lincombs into raster
theSpaceName = grep("^c[[:digit:]]+$", names(inlaResult$marginals.lincomb.derived), value=TRUE)
theSpace = as.integer(gsub("^c", "", theSpaceName))

	linc = inlaResult$summary.lincomb.derived[theSpaceName,]
	linc$space = theSpace
	inlaResult$marginals.predict = 
			inlaResult$marginals.lincomb.derived

	missingCells = values(cells)[! values(cells) %in% theSpace]

	if(length(missingCells)) {
		toadd = matrix(NA, length(missingCells), dim(linc)[2], 
					dimnames=list(
							paste("c", missingCells, sep=""), 
							colnames(linc)
			)
			)
		toadd[,"space"] = missingCells
	
		linc = rbind(linc, toadd)

		# Add in empty lists for the marginals of missing cells
		
		missingMarginals = vector("list", length(missingCells))
		names(missingMarginals) = rownames(toadd)
		
		inlaResult$marginals.predict = c(	
				inlaResult$marginals.predict,
				missingMarginals)
		
	}
	linc = as.matrix(linc[match( values(cells), linc$space),])
	inlaResult$marginals.predict = 
			inlaResult$marginals.predict[
					paste("c", values(cells), sep="")
			]
	
	
	resRasterFitted = 
			brick(extent(cells), nrows=nrow(cells),
					ncols=ncol(cells), crs=projection(cells),
					nl=ncol(linc))
	names(resRasterFitted) = 
			paste("predict.", colnames(linc),sep="")
	
	values(resRasterFitted) = as.vector(linc)
	
	
	
			
	
	
	# posterior distributions




params$range$posterior=inlaResult$marginals.hyperpar[["Range for space"]]
params$range$posterior[,"x"] =  xres(cells) * params$range$posterior[,"x"]
params$range$posterior[,"y"] = params$range$posterior[,"y"] / xres(cells)


# sum(c(0,diff(params$range$posterior[,"x"])) * params$range$posterior[,"y"])
# sum(c(0,diff(params$range$prior[,"x"])) * params$range$prior[,"y"])


params$summary = inlaResult$summary.fixed

params$summary = cbind(params$summary, 
		meanExp = unlist(
				lapply(inlaResult$marginals.fixed,
						function(qq) {
							sum(
									exp(qq[,"x"])*c(0,diff(qq[,"x"]))*qq[,"y"]	
							)
						}
				))
)

if(length(grep("logit",inlaResult$misc$linkfunctions$names))) {
	params$summary = cbind(params$summary, 
			meanInvLogit = unlist(
					lapply(inlaResult$marginals.fixed, function(qq) {
								eqqx = exp(qq[,"x"])
								sum(
										eqqx/(1+eqqx)*c(0,diff(qq[,"x"]))*qq[,"y"]	
								)
							}
							)
	))
}


thecols = paste(c("0.975", "0.5","0.025"), "quant", sep="")

thesd = c(
		sdNugget= grep("^Precision[[:print:]]*Gaussian observations$", 
				names(inlaResult$marginals.hyperpar), value=TRUE),
		gammaShape = grep("^Precision[[:print:]]*Gamma observations$", 
				names(inlaResult$marginals.hyperpar), value=TRUE),
		sd = grep("^Precision[[:print:]]*space$", 
				names(inlaResult$marginals.hyperpar), value=TRUE)
)

params$summary = rbind(params$summary,
		matrix(NA, nrow=length(thesd)+1, ncol=ncol(params$summary),
				dimnames = list(c("range", names(thesd)), 
						colnames(params$summary)))
)


# convert precisions to standard deviations
for(Dsd in grep("gammaShape", names(thesd), invert=TRUE, value=TRUE)) {
	
	params[[Dsd]]$posterior=
			inlaResult$marginals.hyperpar[[thesd[Dsd]]]
	params[[Dsd]]$posterior[,"y"] = params[[Dsd]]$posterior[,"y"] * 2*  
			params[[Dsd]]$posterior[,"x"]^(3/2) 
	params[[Dsd]]$posterior[,"x"] = 1/sqrt(params[[Dsd]]$posterior[,"x"])  
	params[[Dsd]]$posterior = params[[Dsd]]$posterior[
			seq(dim(params[[Dsd]]$posterior)[1],1),]		

	params$summary[Dsd, thecols] = 
				1/sqrt(inlaResult$summary.hyperpar[
								thesd[Dsd],rev(thecols)])
    params$summary[Dsd,"mode"] = 
        1/sqrt(inlaResult$summary.hyperpar[
            thesd[Dsd],'mode'])
    
		
	params$summary[Dsd,"mean"] =sum(
		1/sqrt(inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"x"])*
			c(0,diff(inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"x"]))*
				inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"y"]
	)
}

if(length(grep("gammaShape", names(thesd))) ) {
	Dsd = 'gammaShape'
params[[Dsd]]$posterior=
		inlaResult$marginals.hyperpar[[thesd[Dsd]]]
params[[Dsd]]$posterior[,"y"] = params[[Dsd]]$posterior[,"y"] *  
		params[[Dsd]]$posterior[,"x"]^2 
params[[Dsd]]$posterior[,"x"] = 1/params[[Dsd]]$posterior[,"x"]  
params[[Dsd]]$posterior = params[[Dsd]]$posterior[
		seq(dim(params[[Dsd]]$posterior)[1],1),]		

params$summary[Dsd, thecols] = 
		1/(inlaResult$summary.hyperpar[
						thesd[Dsd],rev(thecols)])
params$summary[Dsd,"mode"] = 
    1/(inlaResult$summary.hyperpar[
            thesd[Dsd],'mode'])


params$summary[Dsd,"mean"] =sum(
		1/(inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"x"])*
				c(0,diff(inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"x"]))*
				inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"y"]
)
}


# put range in summary, in units of distance, not numbers of cells
thecolsFull =c("mean","sd",thecols,"mode") 
params$summary["range",thecolsFull]=				
		xres(cells)*
		inlaResult$summary.hyperpar[
				"Range for space",
				thecolsFull
		]
dimnames(params$summary) = lapply(dimnames(params$summary),
		function(qq) {
			qq=gsub("_", "\\\\textunderscore~", qq)
			qq=gsub("\\$", "\\\\textdollar~", qq)
			qq=gsub("<", "\\\\textless~", qq)
			qq=gsub(">", "\\\\textgreater~", qq)
			qq
		}
)
params$summary = as.data.frame(params$summary)

for(Dvar in names(covariates)) {
	theLevels =levels(covariates[[Dvar]])[[1]]
	if(!is.null(nrow(theLevels))){
	for(D in 1:nrow(theLevels)) {
		rownames(params$summary) = gsub(
			paste("(factor)?(\\()?", Dvar, "(\\))?:?", 
					theLevels[D,1],"$",sep=""),
			paste(Dvar, ":",theLevels[D,2],sep=""), 
					rownames(params$summary))
	}
}
}



resRaster=stack(resRasterRandom, resRasterFitted, cells)

#	if(!is.null(smallBbox))
#		resRaster = crop(resRaster, extent(smallBbox))
		
	result=list(inla=inlaResult,
					raster=resRaster,
					parameters=params
			)

	result
	
}
)