### Function to implement Prais-Winsten correction for autocorrelation
### Author: Konstantin Kashin
### August 1, 2013

prais.correct <- function(method,env.base){	
	### pull objects from base environment
	lm.out <- get("lm.out",envir=env.base)
	yX <- get("yX",envir=env.base)
	units <- get("units",envir=env.base)
	panel.vec <- get("panel.vec",envir=env.base)
	e.mat <- get("e.mat",envir=env.base)
	panel.weight <- get("panel.weight",envir=env.base)
	N.times <- get("N.times",envir=env.base)
	N.units <- get("N.units",envir=env.base)
	bound.rho <- get("bound.rho",envir=env.base)
	rho.na.rm <- get("rho.na.rm",envir=env.base)
	rhotype <- get("rhotype",envir=env.base) 
	obs.mat <- get("obs.mat",envir=env.base)
	rank <- get("rank",envir=env.base)
	singular.ok <- get("singular.ok",envir=env.base)
	p <- ncol(yX)
	
	if(method=="none"){
		pw.output <- list(pw.lm = lm.out, pw.rho = NULL)
	} else{			
		### calculate rho by panel
		rhos <- apply(e.mat, MARGIN=2, function(e) est.rho(e, k=rank,rhotype=rhotype))
		
		### deal with missing rhos
		if(any(is.na(rhos))){
			if(!rho.na.rm){
				stop("Cannot estimate at least one panel-specific autocorrelation. Consider setting rho.na.rm to 'TRUE'.",call.=FALSE)
			} else {
				if(method=="ar1"){
					rhos <- na.omit(rhos)
				} else{
					rhos[is.na(rhos)] <- 0
					message("Setting panel-specific correlation to 0 for at least one panel because unable to estimate autocorrelation.")
				}
			}
		}
		
		### set bounds
		if(bound.rho & any(abs(rhos)>1)){
			rhos[rhos > 1] <- 1
			rhos[rhos < -1] <- -1
			message("Panel-specific correlations bounded to [-1,1]")
		}
		
		### create average rho if ar1
		if(method=="ar1"){
			# calculate total length of runs, not counting the first observation in each run
			weights <- apply(obs.mat,MARGIN=1,function(x) sum(rle(x)$length[rle(x)$values==TRUE]-1))
			# remove weights for rhos that were NA (and so that weight vectors matches rho vector)
			if(!is.null(attr(rhos,"na.action"))){
				weights <- weights[-attr(rhos,"na.action")]
			}
			
			if(panel.weight=="t"){
				# if panel.weight is "t", adjust by 1
				weights <- weights + 1
			}
			
			# calculate average rho
			rho.avg <- sum(rhos*weights)/sum(weights)
			rhos <- rep(rho.avg,N.units)
		}
		
		### do Prais-Winsten transformation here
		transformed.data <- lapply(1:N.units,function(el) prais.transform(el,rhos=rhos,p=p,N.times=N.times,units=units,panel.vec=panel.vec,yX=yX,obs.mat=obs.mat))
		transformed.df <- do.call(rbind,transformed.data)
		transformed.X <- transformed.df[,2:p]
		transformed.y <- transformed.df[,1]
		
		# 2nd stage regression
		pw.lm <- lm(transformed.y ~ transformed.X - 1, singular.ok= singular.ok)
		
		if(method=="ar1"){
			rho.out <- rho.avg
		} else{
			rho.out <- rhos
			}
			
		pw.output <- list(pw.lm = pw.lm, pw.rho = rho.out)
		} # end of autocorrelation else statement
		
	return(pw.output)
	}