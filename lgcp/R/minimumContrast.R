############## IMPORTANT NOTE ##################
# All functions below are merely considered 'helper' functions for the main functions 'minimum.contrast' and 'minimum.contrast.spatiotemporal' at the end of this file.
################################################

##' g.diff.single function
##'
##' A function to find the minimum contrast (squared discrepancy) value based on the pair correlation function, for one specific value of phi
##' (spatial scale) and one specific value of sigma^2 (spatial variance) for the LGCP.
##'
##' @param ps A numeric vector of length 2 giving the values of phi and sigma^2, in that order.
##' @param ghat A numeric vector giving the nonparametric estimate of the PCF at all distances specified in useq (see below)
##' @param useq An increasing, equally spaced numeric vector giving the spatial distances at which the contrast criterion is to be evaluated.
##' @param model A character string specifying the form of the theoretical spatial correlation function (matches 'model' argument for CovarianceFct in the RandomFields packages).
##' @param transform A scalar-valued function which performs a numerical transformation of its argument. Used for calibration of the contrast criterion, by transforming both parametric and nonparametric forms of the PCF.
##' @param power A scalar used for calibration of the contrast criterion: the power which to raise the parametric and nonparametric forms of the PCF to.
##' @param ... Additional arguments if required for definition of the correlation function as per 'model'. See ?CovarianceFct (RandomFields).
##' @return A single numeric value providing the minimum contrast value for the specified value of the ps argument.
##' @export

g.diff.single <- function(ps,ghat,useq,model,transform,power,...){
	if(any(ps<=0)) return(NA)
  	g.parametric <- exp(CovarianceFct(useq,model=model,param=c(mean=NA,variance=ps[2],nugget=0,scale=ps[1],...)))
	return(sum((transform(ghat)^power-transform(g.parametric)^power)^2))
}


##' K.val function
##'
##' An internal function used in computing the theoretical K function for the LGCP. See \link{K.u} for the theoretical K.
##'
##' @param val Spatial lag
##' @param phi Spatial scale parameter value
##' @param sig2 Spatial variance parameter value
##' @param model A character string specifying the form of the theoretical spatial correlation function (matches 'model' argument for CovarianceFct in the RandomFields packages)
##' @param ... Additional arguments if required for definition of the correlation function as per 'model'. See ?CovarianceFct (RandomFields)
##' @return A single numeric value representing a component of the theoretical K function 
##' @export

K.val <- function(val,phi,sig2,model,...){
    return(val*exp(CovarianceFct(val,model=model,param=c(mean=NA,variance=sig2,nugget=0,scale=phi,...))))
}


##' K.u function
##'
##' A function to compute the theoretical K function for the LGCP.
##'
##' @param u Spatial lag at which we wish to find the theoretical K function
##' @param phi Spatial scale parameter value
##' @param sig2 Spatial variance parameter value
##' @param model A character string specifying the form of the theoretical spatial correlation function (matches 'model' argument for CovarianceFct in the RandomFields packages) 
##' @param ... Additional arguments if required for definition of the correlation function as per 'model'. See ?CovarianceFct (RandomFields)
##' @return A single numeric value representing the theoretical K function evaluated at u.
##' @export

K.u <- function(u,phi,sig2,model,...){
    kseq <- c()
    for(i in 1:length(u)){
        intseq <- seq(0,u[i],length=500)
        kseq[i] <- 2*pi*sum(intseq[2]*K.val(intseq,phi,sig2,model,...))
    }
    return(kseq)
}



##' K.diff.single function
##'
##' A function to find the minimum contrast (squared discrepancy) value based on the K function, for one specific value of phi
##' (spatial scale) and one specific value of sigma^2 (spatial variance) for the LGCP.
##'
##' @param ps A numeric vector of length 2 giving the values of phi and sigma^2, in that order. 
##' @param khat A numeric vector giving the nonparametric estimate of the K function at all distances specified in useq (see below) 
##' @param useq An increasing, equally spaced numeric vector giving the spatial distances at which the contrast criterion is to be evaluated.
##' @param model A character string specifying the form of the theoretical spatial correlation function (matches 'model' argument for CovarianceFct in the RandomFields packages).
##' @param transform A scalar-valued function which performs a numerical transformation of its argument. Used for calibration of the contrast criterion, by transforming both parametric and nonparametric forms of the K function.
##' @param power A scalar used for calibration of the contrast criterion: the power which to raise the parametric and nonparametric forms of the K function to.
##' @param ... Additional arguments if required for definition of the correlation function as per 'model'. See ?CovarianceFct (RandomFields).
##' @return A single numeric value providing the minimum contrast value for the specified value of the ps argument.
##' @export

K.diff.single <- function(ps,khat,useq,model,transform,power,...){
    if(any(ps<=0)) return(NA)
    k.parametric <- K.u(useq,ps[1],ps[2],model,...)
    return(sum((transform(khat)^power-transform(k.parametric)^power)^2))
}


##' C.diff.single.im function
##'
##' A function to find the minimum contrast (squared discrepancy) value based on the the temporal autocorrelation function, for one specific value of theta (temporal scale) for the spatiotemporal LGCP.
##' Only the exponential form is considered for the theoretical temporal correlation function. This also depends upon a static pair of values for the spatial scale and spatial variance of the latent Gaussian process (usually estimated first).
##'
##' @param theta Single numeric value for the parameter controlling the scale of temporal dependence in the frequency of observations.
##' @param data Object of class stppp, giving the observed spatiotemporal data set.
##' @param ps A numeric vector of length 2 giving fixed values of phi and sigma^2, in that order. 
##' @param Chat A numeric vector giving the nonparametric estimate of the temporal autocorrelation function at all temporal lags specified by vseq.
##' @param vseq An increasing, equally spaced numeric vector giving the temporal distances at which the contrast criterion is to be evaluated.
##' @param spat A density estimate of the fixed, possibly inhomogeneous, density of the underlying spatial trend. An object of class 'im' (spatstat). May be unnormalised; in which case it will be scaled to integrate to 1 over the spatial study region. 
##' @param model A character string specifying the form of the theoretical spatial correlation function (matches 'model' argument for CovarianceFct in the RandomFields packages).
##' @return A single numeric value providing the minimum contrast value for the specified value of the theta argument.
##' @export

C.diff.single.im <- function(theta,data,ps,Chat,vseq,spat,model){
    myC <- Cvb(data,spat,100,model,covpars=NULL)
    theo <- sapply(vseq,myC,sigma=sqrt(ps[2]),phi=ps[1],theta=theta)
    return(sum((theo*(Chat[1]/theo[1])-Chat)^2))
}



##' my.ginhomAverage function
##'
##' A carbon-copy of ginhomAverage from package 'lgcp', with extra control over the printing of progress bars and other output to the console during execution. 
##' Computes the time-averaged version of the nonparametric PCF (for use with spatiotemporal data).
##'
##' @param xyt an object of class stppp.
##' @param spatial.intensity A spatialAtRisk object giving the possibly inhomogeneous underlying fixed spatial density of the data. 
##' @param temporal.intensity A temporalAtRisk object giving the possibly inhomogeneous underlying fixed temporal intensity of the data. 
##' @param time.window Time interval contained in the interval xyt$tlim over which to compute average. Useful if there is a lot of data over a lot of time points.
##' @param rvals Vector of values for the argument r at which g(r) should be evaluated (see ?pcfinhom). There is a sensible default. 
##' @param correction Choice of edge correction to use, see ?pcfinhom, default is Ripley isotropic correction. 
##' @param suppresswarnings Whether or not to suppress warnings generated by pcfinhom. 
##' @param verbose Whether or not to print function comments and progress to the console during execution. Defaults to TRUE.
##' @param ... Other parameters to be passed to pcfinhom, see ?pcfinhom.
##' @return A vector corresponding to the time-averaged PCF for spatiotemporal data, evaluated at spatial lags defined by 'rvals'.
##' @export

my.ginhomAverage <- function (xyt, spatial.intensity, temporal.intensity, time.window = xyt$tlim, rvals = NULL, correction = "iso", suppresswarnings = FALSE, verbose = TRUE, ...) 
{
    time.window <- sort(as.integer(time.window))
    verifyclass(spatial.intensity, "spatialAtRisk")
    verifyclass(temporal.intensity, "temporalAtRisk")
    density <- as.im(spatial.intensity)
    approxscale <- mean(sapply(xyt$t[xyt$t >= time.window[1] & 
        xyt$t <= time.window[2]], temporal.intensity))
    den <- density
    den$v <- den$v * approxscale
    xyt$t <- as.integer(xyt$t)
    numsamp <- min(c(200, xyt$n))
    if (!suppresswarnings) {
        if (is.null(rvals)) {
            gin <- pcfinhom(xyt[sample(1:xyt$n, numsamp, replace = FALSE)], 
                lambda = den, ...)
        }
        else {
            gin <- pcfinhom(xyt[sample(1:xyt$n, numsamp, replace = FALSE)], 
                lambda = den, r = rvals, ...)
        }
    }
    else {
        if (is.null(rvals)) {
            suppressWarnings(gin <- pcfinhom(xyt[sample(1:xyt$n, 
                numsamp, replace = FALSE)], lambda = den, ...))
        }
        else {
            suppressWarnings(gin <- pcfinhom(xyt[sample(1:xyt$n, 
                numsamp, replace = FALSE)], lambda = den, r = rvals, 
                ...))
        }
    }
    r <- gin$r
    nam <- names(gin)
    nam <- nam[nam != "r"]
    tls <- sapply(time.window[1]:time.window[2], function(x) {
        sum(xyt$t == x)
    })
    tls <- (time.window[1]:time.window[2])[tls != 0]
    ntls <- length(tls)
    if(verbose) pb <- txtProgressBar(min = tls[1], max = rev(tls)[1], style = 3)
    if (!suppresswarnings) {
        pcf <- lapply(tls, function(t) {
            if(verbose) setTxtProgressBar(pb, t)
            den <- density
            den$v <- den$v * temporal.intensity(t)
            try(pcfinhom(xyt[xyt$t == t], lambda = den, r = r, 
                ...))
        })
    }
    else {
        suppressWarnings(pcf <- lapply(tls, function(t) {
            if(verbose) setTxtProgressBar(pb, t)
            den <- density
            den$v <- den$v * temporal.intensity(t)
            try(pcfinhom(xyt[xyt$t == t], lambda = den, r = r, 
                ...))
        }))
    }
    if(verbose) close(pb)
    li <- as.list(pcf[[1]])
    ct <- 1
    if (ntls > 1) {
        for (i in 2:ntls) {
            if (!class(pcf[[i]])[1] == "try-error") {
                li <- add.list(li, as.list(pcf[[i]]))
                ct <- ct + 1
            }
        }
    }
    li <- smultiply.list(li, 1/ct)
    ginhom <- gin
    for (n in nam) {
        idx <- which(names(ginhom) == n)
        ginhom[[idx]] <- li[[idx]]
    }
    attr(ginhom, "correction") <- correction
    #cat("Returning an average of", ct, "curves\n")
    return(ginhom)
}



##' my.KinhomAverage function
##'
##' A carbon-copy of KinhomAverage from package 'lgcp', with extra control over the printing of progress bars and other output to the console during execution. 
##' Computes the time-averaged version of the nonparametric K function (for use with spatiotemporal data).
##'
##' @param xyt an object of class stppp.
##' @param spatial.intensity A spatialAtRisk object giving the possibly inhomogeneous underlying fixed spatial density of the data. 
##' @param temporal.intensity A temporalAtRisk object giving the possibly inhomogeneous underlying fixed temporal intensity of the data. 
##' @param time.window Time interval contained in the interval xyt$tlim over which to compute average. Useful if there is a lot of data over a lot of time points.
##' @param rvals Vector of values for the argument r at which g(r) should be evaluated (see ?Kinhom). There is a sensible default. 
##' @param correction Choice of edge correction to use, see ?Kinhom, default is Ripley isotropic correction. 
##' @param suppresswarnings Whether or not to suppress warnings generated by Kinhom. 
##' @param verbose Whether or not to print function comments and progress to the console during execution. Defaults to TRUE.
##' @return A vector corresponding to the time-averaged K function for spatiotemporal data, evaluated at spatial lags defined by 'rvals'.
##' @export

my.KinhomAverage <- function (xyt, spatial.intensity, temporal.intensity, time.window = xyt$tlim, rvals = NULL, correction = "iso", suppresswarnings = FALSE, verbose = TRUE) 
{
    time.window <- sort(as.integer(time.window))
    verifyclass(spatial.intensity, "spatialAtRisk")
    verifyclass(temporal.intensity, "temporalAtRisk")
    density <- as.im(spatial.intensity)
    approxscale <- mean(sapply(xyt$t[xyt$t >= time.window[1] & 
        xyt$t <= time.window[2]], temporal.intensity))
    den <- density
    den$v <- den$v * approxscale
    xyt$t <- as.integer(xyt$t)
    numsamp <- min(c(200, xyt$n))
    if (!suppresswarnings) {
        if (is.null(rvals)) {
            Kin <- Kinhom(xyt[sample(1:xyt$n, numsamp, replace = FALSE)], 
                lambda = den)
        }
        else {
            Kin <- Kinhom(xyt[sample(1:xyt$n, numsamp, replace = FALSE)], 
                lambda = den, r = rvals)
        }
    }
    else {
        if (is.null(rvals)) {
            suppressWarnings(Kin <- Kinhom(xyt[sample(1:xyt$n, 
                numsamp, replace = FALSE)], lambda = den))
        }
        else {
            suppressWarnings(Kin <- Kinhom(xyt[sample(1:xyt$n, 
                numsamp, replace = FALSE)], lambda = den, r = rvals))
        }
    }
    r <- Kin$r
    nam <- names(Kin)
    nam <- nam[nam != "r"]
    tls <- sapply(time.window[1]:time.window[2], function(x) {
        sum(xyt$t == x)
    })
    tls <- (time.window[1]:time.window[2])[tls != 0]
    ntls <- length(tls)
    if(verbose) pb <- txtProgressBar(min = tls[1], max = rev(tls)[1], style = 3)
    if (!suppresswarnings) {
        pcf <- lapply(tls, function(t) {
            if(verbose) setTxtProgressBar(pb, t)
            den <- density
            den$v <- den$v * temporal.intensity(t)
            try(Kinhom(xyt[xyt$t == t], lambda = den, r = r))
        })
    }
    else {
        suppressWarnings(pcf <- lapply(tls, function(t) {
            if(verbose) setTxtProgressBar(pb, t)
            den <- density
            den$v <- den$v * temporal.intensity(t)
            try(Kinhom(xyt[xyt$t == t], lambda = den, r = r))
        }))
    }
    if(verbose) close(pb)
    li <- as.list(pcf[[1]])
    ct <- 1
    if (ntls > 1) {
        for (i in 2:ntls) {
            if (!class(pcf[[i]])[1] == "try-error") {
                li <- add.list(li, as.list(pcf[[i]]))
                ct <- ct + 1
            }
        }
    }
    li <- smultiply.list(li, 1/ct)
    Kinhom <- Kin
    for (n in nam) {
        idx <- which(names(Kinhom) == n)
        Kinhom[[idx]] <- li[[idx]]
    }
    attr(Kinhom, "correction") <- correction
    #cat("Returning an average of", ct, "curves\n")
    return(Kinhom)
}
##################################################



##################### SPATIAL ####################


##' minimum.contrast function
##'
##' A function to provide minimum contrast (aka least squares) estimates of the spatial scale (phi) and spatial variance (sigma^2) assuming an LGCP modelling framework
##' for spatial data.
##'
##' @param data An object of class 'ppp' (package spatstat) with a polygonal window. May be univariate or multitype. 
##' @param model Assumed theoretical form of the spatial correlation function. Matches 'model' argument for 'CovarianceFct' in package RandomFields.
##' @param method Character string indicating which version of spatial minimum contrast to use: either "K" or "g".
##' @param intens Underlying deterministic spatial intensity. A single function f(x,y) or a single pixel image if univariate, a list of these objects if point pattern is multitype (order must correspond to order of ppp marks).
##' @param power Power to raise the functions to in the contrast criterion. Default 1.
##' @param transform Transformation to apply to the functions in the contrast criterion. Default no transformation. 
##' @param startvals Starting values for 'optim' in minimising the contrast criterion in the order c(phi,sigma2). A list of these if multitype. If NULL, the function automatically attempts to find suitable starting values, though no guarantee of 'optim' convergence can be given! 
##' @param verbose Boolean. Whether or not to print function progress. 
##' @param ... Additional arguments to be passed to 'param' in evaluation of 'CovarianceFct' (need dependent upon 'model').
##' @return Returned values are the minimum contrast estimates of phi and sigma^2, as well as the overall squared discrepancy between the parametric and nonparametric forms of the function used corresponding to these estimates. (This can be useful in deciding between several different theoretical forms of the correlation specified by 'model'). If the point pattern is multitype, each pair of parameters is estimated independently for each marginal (type-specific) data set.
##' @seealso \link{minimum.contrast.spatiotemporal}, link{chooseCellWidth}, \link{getpolyol}, \link{guessinterp}, \link{getZmat},
##' \link{addTemporalCovariates}, \link{lgcpPrior}, \link{lgcpInits}, \link{CovFunction}
##' \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}, \link{lgcpPredictSpatioTemporalPlusPars}, 
##' \link{lgcpPredictMultitypeSpatialPlusPars}
##' @export

minimum.contrast <- function(data,model,method="g",intens=NULL,power=1,transform=NULL,startvals=NULL,verbose=TRUE,...){
    
	if(!is.ppp(data)) stop("'data' must be an object of class 'ppp'")
	W <- data$window
	ntypes <- 1
	rec <- as.rectangle(W)
	
	
	if(!is.polygonal(W)&&!is.rectangle(W)) stop("'data' must have a polygon or rectangle window")
	
	if(is.null(transform)) transform <- function(x) return(x)
	if(!is.function(transform)) stop("'transform' must be a scalar-valued function")	
		
	if(!is.multitype(data)){
		if(is.null(startvals)){
			startvals <- c(min(diff(rec$xrange),diff(rec$yrange))/100,log(data$n)/2)
			#print(startvals)
		}
				
		data.m <- NULL
		if(!is.null(intens)&&!is.function(intens)&&!is.im(intens)) stop("'intens' must be a function f(x,y) OR a pixel image of class 'im'")
		if(is.im(intens)){
			intens.im <- intens
			intens <- function(x,y){
				rstr <- nearest.raster.point(x,y,w=as.mask(W,dimyx=dim(intens.im$v)))
				if(length(x)>1)	return(diag(intens.im$v[rstr$row,rstr$col]))
				else return(intens.im$v[rstr$row,rstr$col])
			}
		}
		
		if(verbose) cat("[Univariate spatial minimum contrast]\n")
		if(method=="g"){
			if(is.null(intens)){
				if(verbose) cat("Nonparametric homogeneous PCF estimation...")
				nonpar <- pcf(data,correction="Ripley",stoyan=0.15,nlarge=Inf)
				if(verbose) cat("done.\n")
			} else {
				if(verbose) cat("Nonparametric heterogenous PCF estimation...")
				nonpar <- try(pcfinhom(data,lambda=intens,correction="Ripley",stoyan=0.15,nlarge=Inf))
				if(inherits(nonpar,"try-error")){
				    if(!is.im(intens.im)){
				        nonpar <- try(pcfinhom(data,lambda=intens,correction="Ripley",stoyan=0.15,nlarge=Inf)) # will print the fail message from the previous line
				    }
				    nonpar <- try(pcfinhom(data,lambda=intens.im,correction="Ripley",stoyan=0.15,nlarge=Inf))
				}
				if(verbose) cat("done.\n")
			}
			if(verbose) cat(paste("Starting values are (",round(startvals[1],2),", ",round(startvals[2],2),"); optimising ",model," correlation function...",sep=""))
			MINCON <- optim(par=startvals,fn=g.diff.single,ghat=nonpar$iso[-1],useq=nonpar$r[-1],model=model,transform=transform,power=power,...)$par
			if(verbose) cat("done.\n")
			
			parametric <- exp(CovarianceFct(nonpar$r[-1],model=model,param=c(mean=NA,variance=MINCON[2],nugget=0,scale=MINCON[1],...)))
			disc.vec <- sum((transform(nonpar$iso[-1])^power-transform(parametric)^power)^2)*nonpar$r[2]
		} else if(method=="K"){
			if(is.null(intens)){
				if(verbose) cat("Nonparametric homogeneous K-function estimation...")
				nonpar <- Kest(data,correction="Ripley",nlarge=Inf,r=seq(0,rmax.rule("K",W=W,data$n/area.owin(W)),length=65))
				if(verbose) cat("done.\n")
			} else {
				if(verbose) cat("Nonparametric heterogeneous K-function estimation...")
				nonpar <- Kinhom(data,lambda=intens,correction="Ripley",nlarge=Inf,r=seq(0,rmax.rule("K",W=W,data$n/area.owin(W)),length=65))
				if(verbose) cat("done.\n")
			}
			if(verbose) cat(paste("Starting values are (",round(startvals[1],2),", ",round(startvals[2],2),"); optimising ",model," correlation function...",sep=""))
			MINCON <- optim(par=startvals,fn=K.diff.single,khat=nonpar$iso[-1],useq=nonpar$r[-1],model=model,transform=transform,power=power,...)$par
			if(verbose) cat("done.\n")
			
			parametric <- K.u(nonpar$r[-1],MINCON[1],MINCON[2],model=model,...)
			disc.vec <- sum((transform(nonpar$iso[-1])^power-transform(parametric)^power)^2)*nonpar$r[2]	
		} else {
			stop("'method' argument must be either \"g\" or \"K\"")
		}
	} else {
		data.types <- split(data)
		data.m <- as.character(unique(marks(data)))
		ntypes <- length(data.types)
		
		if(is.null(startvals)) startvals <- vector("list",ntypes)
		if(!is.list(startvals)) stop("'startvals' must be a list or NULL for multi-type point patterns")
		else if(length(startvals)!=ntypes) stop("'startvals' list length does not correspond to the number of observation types")

		
		if(is.null(intens)) intens <- vector("list",ntypes)
		if(!is.list(intens)) stop("'intens' must be a list or NULL for multi-type point patterns")
		else if(length(intens)!=ntypes) stop("'intens' list length does not correspond to the number of observation types")

		if(verbose) cat("[Multivariate (independent) spatial minimum contrast]\n")
		
		MINCON <- disc.vec <- c()
		intens.temp <- list()
		for(i in 1:ntypes){
			if(is.null(startvals[[i]])){
				startvals[[i]] <- c(min(diff(rec$xrange),diff(rec$yrange))/100,log(data.types[[i]]$n)/2)
				#print(startvals[[i]])
			}
			
			if(!is.null(intens[[i]])&&!is.function(intens[[i]])&&!is.im(intens[[i]])) stop("'intens' entries must be functions f(x,y) OR a pixel images of class 'im'")
			if(is.im(intens[[i]])){
				intens.im <- intens[[i]]
				intens.temp[[i]] <- intens[[i]]
				intens[[i]] <- function(x,y){
					rstr <- nearest.raster.point(x,y,w=as.mask(W,dimyx=dim(intens.im$v)))
					if(length(x)>1)	return(diag(intens.im$v[rstr$row,rstr$col]))
					else return(intens.im$v[rstr$row,rstr$col])
				}
			}
			
			if(method=="g"){
				if(is.null(intens[[i]])){
					if(verbose){
					    cat(paste("TYPE ",data.m[i],": Nonparametric homogeneous PCF estimation...",sep=""))
					}
					nonpar <- pcf(data.types[[i]],correction="Ripley",stoyan=0.15,nlarge=Inf)
					if(verbose){
					    cat("done.\n")
					}
				} 
				else {
					if(verbose){ 
					    cat(paste("TYPE ",data.m[i],": Nonparametric heterogeneous PCF estimation...",sep=""))
    					nonpar <- try(pcfinhom(data.types[[i]],lambda=intens[[i]],correction="Ripley",stoyan=0.15,nlarge=Inf))
    					if(inherits(nonpar,"try-error")){
        				    if(!is.im(intens.temp[[i]])){
        				        nonpar <- try(pcfinhom(data.types[[i]],lambda=intens[[i]],correction="Ripley",stoyan=0.15,nlarge=Inf)) # will print the fail message from the previous line
        				    }
        				    nonpar <- try(pcfinhom(data.types[[i]],lambda=intens.temp[[i]],correction="Ripley",stoyan=0.15,nlarge=Inf))
    				    }
				    }
				    if(verbose){ 
				        cat("done.\n")
				    }
				}
					
				
				if(verbose) cat(paste("TYPE ",data.m[i],": Starting values are (",round(startvals[[i]][1],2),", ",round(startvals[[i]][2],2),"); optimising ",model," correlation function...",sep=""))
				MINCON <- append(MINCON,optim(par=startvals[[i]],fn=g.diff.single,ghat=nonpar$iso[-1],useq=nonpar$r[-1],model=model,transform=transform,power=power,...)$par)
				if(verbose) cat("done.\n")

				parametric <- exp(CovarianceFct(nonpar$r[-1],model=model,param=c(mean=NA,variance=tail(MINCON,2)[2],nugget=0,scale=tail(MINCON,2)[1],...)))
				disc.vec[i] <- sum((transform(nonpar$iso[-1])^power-transform(parametric)^power)^2)*nonpar$r[2]

			} else if(method=="K"){
				if(is.null(intens[[i]])){
					if(verbose) cat(paste("TYPE ",data.m[i],": Nonparametric homogeneous K-function estimation...",sep=""))
					nonpar <- Kest(data.types[[i]],correction="Ripley",nlarge=Inf,r=seq(0,rmax.rule("K",W=W,data$n/area.owin(W)),length=65))
					if(verbose) cat("done.\n")
				} else {
					if(verbose) cat(paste("TYPE ",data.m[i],": Nonparametric heterogeneous K-function estimation...",sep=""))
					nonpar <- Kinhom(data.types[[i]],lambda=intens[[i]],correction="Ripley",nlarge=Inf,r=seq(0,rmax.rule("K",W=W,data$n/area.owin(W)),length=65))
					if(verbose) cat("done.\n")
				}
				if(verbose) cat(paste("TYPE ",data.m[i],": Starting values are (",round(startvals[[i]][1],2),", ",round(startvals[[i]][2],2),"); optimising ",model," correlation function...",sep=""))
				MINCON <- append(MINCON,optim(par=startvals[[i]],fn=K.diff.single,khat=nonpar$iso[-1],useq=nonpar$r[-1],model=model,transform=transform,power=power,...)$par)
				if(verbose) cat("done.\n")
				
				parametric <- K.u(nonpar$r[-1],tail(MINCON,2)[1],tail(MINCON,2)[2],model=model,...)
				disc.vec[i] <- sum((transform(nonpar$iso[-1])^power-transform(parametric)^power)^2)*nonpar$r[2]	
			} else {
				stop("'method' entries must be either \"g\" or \"K\"")
			}
		}
	}
		
	result <- matrix(MINCON,ntypes,2,byrow=T)
	dimnames(result) <- list(data.m,c("scale","variance"))
	discrepancy <- matrix(disc.vec)
	dimnames(discrepancy) <- list(data.m,"Squared discrepancy")
	return(list(estimates=result,discrepancy=discrepancy))
}



##################### SPATIOTEMPORAL ####################

##' minimum.contrast.spatiotemporal function
##'
##' A function to provide minimum contrast (aka least squares) estimates of the spatial scale (phi), spatial variance (sigma^2) and temporal scale (theta) assuming an LGCP modelling framework
##' for spatiotemporal data. Currently only implemented for univariate (i.e. unmarked) spatiotemporal point patterns
##'
##' @param data An object of class 'stppp' from package 'lgcp'. Must be univariate i.e. have ' data$markformat=="none" ' 
##' @param model Assumed theoretical form of the spatial correlation function. Matches 'model' argument for 'CovarianceFct' in package RandomFields.
##' @param method Character string indicating which version of spatial minimum contrast to use: either "K" or "g". 
##' @param spatial.dens An object of class 'spatialAtRisk', or a (possibly unnormalised) pixel image of class 'im', giving the underlying deterministic spatial density.   
##' @param temporal.intens An object of class 'temporalAtRisk', giving the deterministic, possibly inhomogeneous, temporal intensity.
##' @param power Power to raise the functions to in the spatial contrast criterion. Default 1.
##' @param transform Transformation to apply to the spatial functions in the contrast criterion. Default no transformation. 
##' @param spatial.startvals Starting values for 'optim' in minimising the contrast criterion in the order c(phi,sigma2). If NULL, the function automatically attempts to find suitable starting values, though no guarantee of 'optim' convergence can be given! 
##' @param temporal.interval Defaults to c(0.1,10) if NULL. An interval of the form 'c(lowerlimit,upperlimit)' to be passed to 'optimise'. This is the interval in which the function will search for an optimal value for theta (the scale parameter for temporal dependence). Note that only the exponential covariance model is implemented for temporal dependence. 
##' @param verbose Boolean. Whether or not to print function progress.  
##' @param ... Additional arguments to be passed to 'param' in evaluation of 'CovarianceFct' (need dependent upon 'model'). 
##' @return Returned values are the minimum contrast estimates of phi, sigma^2 and theta, as well as the overall squared discrepancy between the parametric and nonparametric forms of the spatial function used corresponding to these estimates. (This can be useful in deciding between several different theoretical forms of the spatial correlation specified by 'model').
##' @seealso \link{minimum.contrast}, link{chooseCellWidth}, \link{getpolyol}, \link{guessinterp}, \link{getZmat},
##' \link{addTemporalCovariates}, \link{lgcpPrior}, \link{lgcpInits}, \link{CovFunction}
##' \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}, \link{lgcpPredictSpatioTemporalPlusPars}, 
##' \link{lgcpPredictMultitypeSpatialPlusPars}
##' @export


minimum.contrast.spatiotemporal <- function(data,model,method="g",spatial.dens=NULL,temporal.intens=NULL,power=1,transform=NULL,spatial.startvals=NULL,temporal.interval=NULL,verbose=TRUE,...){
	if(class(data)[1]!="stppp") stop("'data' must be an object of class 'stppp'")
	W <- data$window
	uni.t <- unique(data$t)
	rec <- as.rectangle(W)
	xs <- seq(W$xrange[1],W$xrange[2],length=129)
	xc <- xs[1]+0.5*(xs[2]-xs[1])+0:127*(xs[2]-xs[1])
	ys <- seq(W$yrange[1],W$yrange[2],length=129)
	yc <- ys[1]+0.5*(ys[2]-ys[1])+0:127*(ys[2]-ys[1])

	if(data$markformat=="none"){
		if(verbose) cat("[Univariate spatio-temporal minimum contrast]\n")
		
		if(is.null(spatial.startvals)){
			spatial.startvals <- c(min(diff(rec$xrange),diff(rec$yrange))/100,log(data$n)/2)
		}
		
		if(is.null(temporal.intens)) temporal.intens <- temporalAtRisk(function(x) 1,data$tlim,data)
		if(is.null(spatial.dens)) spatial.dens <- im(matrix(1,128,128),xcol=xc,yrow=yc)
		if(!any(class(spatial.dens)=="spatialAtRisk")){
			if(is.im(spatial.dens)) spatial.dens <- spatialAtRisk(spatial.dens)
			else stop("'spatial.dens' must be an object of class 'spatialAtRisk' or 'im'")
		}					
		if(!any(class(temporal.intens)=="temporalAtRisk")) stop("'temporal.intens' must be an intensity-scaled object of class 'temporalAtRisk'")	

		if(method=="g"){
			if(verbose) cat("Spatial: Time-averaged PCF estimation...\n")
    		nonpar <- my.ginhomAverage(data,spatial.dens,temporal.intens,suppresswarnings=T,verbose=verbose)
    		if(verbose) cat(paste("Spatial: Starting values are (",round(spatial.startvals[1],2),", ",round(spatial.startvals[2],2),"); optimising ",model," correlation function...",sep=""))
			MINCON.SPATIAL <- optim(par=spatial.startvals,fn=g.diff.single,ghat=nonpar$iso[-1],useq=nonpar$r[-1],model=model,transform=transform,power=power,...)$par
			parametric <- exp(CovarianceFct(nonpar$r[-1],model=model,param=c(mean=NA,variance=MINCON.SPATIAL[2],nugget=0,scale=MINCON.SPATIAL[1],...)))
			disc.vec <- sum((transform(nonpar$iso[-1])^power-transform(parametric)^power)^2)*nonpar$r[2]
			if(verbose) cat("done.\n")
		} else if(method=="K"){
			if(verbose) cat("Spatial: Time-averaged K-function estimation...\n")   		
    		nonpar <- my.KinhomAverage(data,spatial.dens,temporal.intens,suppresswarnings=T,verbose=verbose)
    		if(verbose) cat(paste("Spatial: Starting values are (",round(spatial.startvals[1],2),", ",round(spatial.startvals[2],2),"); optimising ",model," correlation function...",sep=""))
			MINCON.SPATIAL <- optim(par=spatial.startvals,fn=K.diff.single,khat=nonpar$iso[-1],useq=nonpar$r[-1],model=model,transform=transform,power=power,...)$par
			parametric <- K.u(nonpar$r[-1],MINCON.SPATIAL[1],MINCON.SPATIAL[2],model=model,...)
			disc.vec <- sum((transform(nonpar$iso[-1])^power-transform(parametric)^power)^2)*nonpar$r[2]
			if(verbose) cat("done.\n")
		}
		
		if(verbose) cat("Temporal: Covariance function estimation...")
		uqt <- as.numeric(names(table(as.integer(data$t))))
		tvals <- c()
		for(i in 1:length(uni.t)) tvals[i] <- temporal.intens(uni.t[i])
    	autocov <- acf(table(as.integer(data$t))-tvals,type="covariance",plot=F)
    	v <- 0:5
    	blen <- 10
    	index <- 1
    	vseq <- Chat <- rep(NA,(blen-1)*max(v)+1)
    	for(i in 1:(length(v)-1)){
        	xs <- seq(v[i],v[i+1],length=blen)
        	if(i>1) xs <- xs[-1]
        	vseq[index:(index+length(xs)-1)] <- xs
        	Chat[index:(index+length(xs)-1)] <- approx(v[c(i,i+1)],autocov$acf[c(i,i+1)],xout=xs)$y
        	index <- index+length(xs)
    	}
		if(is.null(temporal.interval)) temporal.interval <- c(0.1,10)
		if(verbose) cat("done.\n")
		if(verbose) cat(paste("Temporal: Starting interval is (",round(temporal.interval[1],2),", ",round(temporal.interval[2],2),"); optimising for exponential dependence...",sep=""))
    	MINCON.TEMPORAL <- optimise(f=C.diff.single.im,interval=temporal.interval,data=data,ps=MINCON.SPATIAL,Chat=Chat,vseq=vseq,spat=im(t(spatial.dens$Zm),xcol=xc,yrow=yc),model=model)$minimum
		if(verbose) cat("done.\n")
		result <- matrix(c(MINCON.SPATIAL,MINCON.TEMPORAL),1,3,byrow=T,dimnames=list(NULL,c("scale (spatial)","variance (spatial)","scale (temporal)")))
		return(list(estimates=result,discrepancy=matrix(disc.vec,1,1,dimnames=list(NULL,"Squared discrepancy (spatial)"))))
	} else {
		if(verbose) cat("[Multivariate spatio-temporal minimum contrast]\n")
		
		stop("Code for multitype space time processes not yet implemented.")
	}	
}

