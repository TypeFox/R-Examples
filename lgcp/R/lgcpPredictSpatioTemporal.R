##' lgcpPredict function
##'
##' The function \code{lgcpPredict} performs spatiotemporal prediction for log-Gaussian Cox Processes
##'
##' The following is a mathematical description of a log-Gaussian Cox Process, it is best viewed in the pdf version of the manual.
##'
##' Let \eqn{\mathcal Y(s,t)}{\mathcal Y(s,t)} be a spatiotemporal Gaussian process, \eqn{W\subset R^2}{W\subset R^2} be an 
##' observation window in space and \eqn{T\subset R_{\geq 0}}{T\subset R_{\geq 0}} be an interval of time of interest. 
##' Cases occur at spatio-temporal positions \eqn{(x,t) \in W \times T}{(x,t) \in W \times T} 
##'  according to an inhomogeneous spatio-temporal Cox process,
##' i.e. a Poisson process with a stochastic intensity \eqn{R(x,t)}{R(x,t)},
##'   The number of cases, \eqn{X_{S,[t_1,t_2]}}{X_{S,[t_1,t_2]}}, arising in 
##'   any \eqn{S \subseteq W}{S \subseteq W} during the interval \eqn{[t_1,t_2]\subseteq T}{[t_1,t_2]\subseteq T} is 
##'   then Poisson distributed conditional on \eqn{R(\cdot)}{R(\cdot)},
##' \deqn{X_{S,[t_1,t_2]} \sim \mbox{Poisson}\left\{\int_S\int_{t_1}^{t_2} R(s,t)d sd t\right\}}{X_{S,[t_1,t_2]} \sim \mbox{Poisson}\left\{\int_S\int_{t_1}^{t_2} R(s,t)d sd t\right\}.}
##' Following Brix and Diggle (2001) and Diggle et al (2005), the intensity is decomposed multiplicatively as
##' \deqn{R(s,t) = \lambda(s)\mu(t)\exp\{\mathcal Y(s,t)\}.}{R(s,t) = \lambda(s)\mu(t)Exp\{\mathcal Y(s,t)\}.}
##' In the above, the fixed spatial component, \eqn{\lambda:R^2\mapsto R_{\geq 0}}{\lambda:R^2\mapsto R_{\geq 0}}, 
##' is a known function, proportional to the population at risk at each point in space and scaled so that
##' \deqn{\int_W\lambda(s)d s=1,}{\int_W\lambda(s)d s=1,}
##' whilst the fixed temporal component, 
##'  \eqn{\mu:R_{\geq 0}\mapsto R_{\geq 0}}{\mu:R_{\geq 0}\mapsto R_{\geq 0}}, is also a known function with
##' \deqn{\mu(t) \delta t = E[X_{W,\delta t}],}{\mu(t) \delta t = E[X_{W,\delta t}],}
##' for \eqn{t}{t} in a small interval of time, \eqn{\delta t}{\delta t}, over which the rate of the process over \eqn{W}{W} can be considered constant.
##'
##' \bold{
##'     NOTE: the xyt stppp object can be recorded in continuous time, but for the purposes of prediciton,    
##'     discretisation must take place. For the time dimension, this is achieved invisibly by \code{as.integer(xyt$t)} and
##'     \code{as.integer(xyt$tlim)}. Therefore, before running an analysis please make sure that this is commensurate
##'     with the physical inerpretation and requirements of your output. The spatial discretisation is
##'     chosen with the argument cellwidth (or gridsize). If the chosen discretisation in time and space is too coarse for a
##'     given set of parameters (sigma, phi and theta) then the proper correlation structures implied by the model will not
##'     be captured in the output.
##' }
##'
##' Before calling this function, the user must decide on the time point of interest, the
##' number of intervals of data to use, the parameters, spatial covariance model, spatial discretisation,
##' fixed spatial (\eqn{\lambda(s)}{\lambda(s)}) and temporal (\eqn{\mu(t)}{\mu(t)}) components, mcmc parameters, and whether or not any output is
##' required.
##'
##' @param xyt a spatio-temporal point pattern object, see ?stppp
##' @param T time point of interest
##' @param laglength specifies lag window, so that data from and including  time (T-laglength) to time T is used in the MALA algorithm
##' @param model.parameters values for parameters, see ?lgcppars
##' @param spatial.covmodel correlation type see ?CovarianceFct 
##' @param covpars vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct
##' @param cellwidth width of grid cells on which to do MALA (grid cells are square) in same units as observation window. Note EITHER gridsize OR cellwidth must be specified.
##' @param gridsize size of output grid required. Note EITHER gridsize OR cellwidthe must be specified.
##' @param spatial.intensity the fixed spatial component: an object of that can be coerced to one of class spatialAtRisk
##' @param temporal.intensity the fixed temporal component: either a numeric vector, or a function that can be coerced into an object of class temporalAtRisk
##' @param mcmc.control MCMC paramters, see ?mcmcpars
##' @param output.control output choice, see ?setoutput
##' @param missing.data.areas a list of owin objects (of length laglength+1) which has xyt$window as a base window, but with polygonal holes specifying spatial areas where there is missing data.
##' @param autorotate logical: whether or not to automatically do MCMC on optimised, rotated grid.   
##' @param gradtrunc truncation for gradient vector equal to H parameter Moller et al 1998 pp 473. Default is Inf, which means no gradient truncation. Set to NULL to estimate this automatically (though note that this may not necessarily be a good choice). The default seems to work in most settings.
##' @param ext integer multiple by which grid should be extended, default is 2. Generally this will not need to be altered, but if the spatial correlation decays very slowly (compared withe the size of hte observation window), increasing 'ext' may be necessary.
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window.
##' further notes on autorotate argument: If set to TRUE, and the argument spatial is not NULL, then the argument spatial must be computed in the original frame of reference (ie NOT in the rotated frame). 
##' Autorotate performs bilinear interpolation (via interp.im) on an inverse transformed grid; if there is no computational advantage in doing this, a warning message will be issued. Note that best accuracy 
##' is achieved by manually rotating xyt and then computing spatial on the transformed xyt and finally feeding these in as arguments to the function lgcpPredict. By default autorotate is set to FALSE.
##' @return the results of fitting the model in an object of class \code{lgcpPredict}
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor, Tilman M. Davies, Barry S. Rowlingson, Peter J. Diggle (2013). Journal of Statistical Software, 52(4), 1-40. URL http://www.jstatsoft.org/v52/i04/        
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##'     \item Wood ATA, Chan G (1994). Simulation of Stationary Gaussian Processes in [0,1]d. Journal of Computational and Graphical Statistics, 3(4), 409-432.
##'     \item Moller J, Syversveen AR, Waagepetersen RP (1998). Log Gaussian Cox Processes. Scandinavian Journal of Statistics, 25(3), 451-482.
##' }
##' @seealso \link{KinhomAverage}, \link{ginhomAverage}, \link{lambdaEst}, \link{muEst}, \link{spatialparsEst}, \link{thetaEst},  
##' \link{spatialAtRisk}, \link{temporalAtRisk}, \link{lgcppars}, \link{CovarianceFct}, \link{mcmcpars}, \link{setoutput} 
##' \link{print.lgcpPredict}, \link{xvals.lgcpPredict}, \link{yvals.lgcpPredict}, \link{plot.lgcpPredict}, \link{meanfield.lgcpPredict},
##' \link{rr.lgcpPredict}, \link{serr.lgcpPredict}, \link{intens.lgcpPredict},   
##' \link{varfield.lgcpPredict}, \link{gridfun.lgcpPredict}, \link{gridav.lgcpPredict}, \link{hvals.lgcpPredict}, \link{window.lgcpPredict},
##' \link{mcmctrace.lgcpPredict}, \link{plotExceed.lgcpPredict}, \link{quantile.lgcpPredict}, \link{identify.lgcpPredict}, \link{expectation.lgcpPredict},
##' \link{extract.lgcpPredict}, \link{showGrid.lgcpPredict}
##' @export 
    
lgcpPredict <- function(xyt,
					    T,
					    laglength,
					    model.parameters=lgcppars(),
					    spatial.covmodel="exponential",
					    covpars=c(),
					    cellwidth=NULL,
					    gridsize=NULL,
					    spatial.intensity,
					    temporal.intensity,					
					    mcmc.control,
					    output.control=setoutput(),
					    missing.data.areas=NULL,
					    autorotate=FALSE,
					    gradtrunc=Inf,
					    ext=2,
					    inclusion="touching"){

    
    starttime <- Sys.time()
    
    ###
    # Convert times into integer-valued vectors
    ###
    
    if (!inherits(T,"integer")){
	    warning("Converting T into integer value, see ?as.integer",immediate.=TRUE)
	    T <- as.integer(T) 
	}
	if (!inherits(laglength,"integer")){
	    warning("Converting laglength into integer values, see ?as.integer",immediate.=TRUE)
	    laglength <- as.integer(laglength) 
	}
	if (!inherits(xyt$tlim,"integer")){
	    warning("Converting xyt$tlim into integer values, see ?as.integer",immediate.=TRUE)
	    xyt$tlim <- as.integer(xyt$tlim) # convert times into integer values: they should already be in this form.
	}	
	if (!inherits(xyt$t,"integer")){
	    warning("Converting xyt$t into integer values, see ?as.integer",immediate.=TRUE)
	    xyt$t <- as.integer(xyt$t)
	}
	if(xyt$window$type=="rectangle"){
	    xyt$window <- as.polygonal(xyt$window)
	}
	
	###
	# select cellwidth if gridsize specified
	###
	
	if(is.null(cellwidth) & is.null(gridsize)){
	    stop("Either cell width OR grid size must be specified")
	}
	if(!is.null(cellwidth) & !is.null(gridsize)){
	    stop("Either cell width OR grid size must be specified")
	}
	if (!all(sapply(gridsize,is.pow2))){
	    stop("All elements of gridsize must be a power of 2")
	}
	if(!is.null(gridsize) & autorotate==TRUE){
	    warning("In order to use autorotate, you must specify a cell width instead ... SETTING autorotate=FALSE.",immediate.=TRUE)
	    autorotate <- FALSE
	}
	if(!is.null(gridsize)){
	    approxcw <- diff(xyt$window$xrange)/gridsize[1] # approx cell width
	    cwseq <- seq(approxcw/2,2*approxcw,length.out=500)
	    cwfun <- function(cw){
	        ow <- selectObsWindow(xyt,cw)
	        return(c(ow$M,ow$N))
	    }
	    gsmat <- t(sapply(cwseq,cwfun))
	    tf <- apply(gsmat,1,function(x){return(all(x==gridsize))})
	    if(sum(tf)==0){
	        stop("No sensible observation window found: either change gridsize, or specify cellwidth instead")
	    }
	    else{
	        cellwidth <- cwseq[min(which(tf))]
	    }
	}
    
    ###
    # Perform basic checks 
    ###    					

    if (!is.null(gradtrunc)){
        if(gradtrunc<0){
            stop("gradtrunc must be non-negative")
        }
    }
		
	if (!inherits(temporal.intensity,"temporalAtRisk")){
	    temporal.intensity <- temporalAtRisk(temporal.intensity,tlim=xyt$tlim,xyt=xyt)
	}
	else{
	    if(!all(as.integer(xyt$tlim)==attr(temporal.intensity,"tlim"))){
	        stop("Incompatible temporal.intensity, integer time limits (xyt$tlim and temporal.intensity$tlim) do not match")
	    }
	}
	
	if(laglength==0){
	    stop("laglength must be >= 1")
	}
	
	if(mcmc.control$burnin>mcmc.control$mala.length){
		stop("Number of burnin iterations must be less than the total number of iterations")
	}
	
	aggtimes <- T - laglength:0
	nobser <- 0
	for (i in 1:(laglength+1)){
	    nobser <- nobser + sum(xyt$t==aggtimes[i])
	}
	if(nobser==0){
	    cat("NOTE: time data should be integer-valued.\n")
		stop("No data in chosen time interval")
	}
	
	temporalfit <- sapply(aggtimes,temporal.intensity)
	if (any(is.na(temporalfit))){
	    stop("Missing temporal fitted values")
	}
	
	###
	# Initialise
	###
	
	sigma <- model.parameters$sigma
	phi <- model.parameters$phi
	mu <- model.parameters$mu
	theta <- model.parameters$theta
	
	tdiff <- c(Inf,diff(aggtimes)) # require Inf for target and gradient function to work
	numt <- length(tdiff)
    bt <- exp(-theta*tdiff)
	
	###
    # compute whether there is any efficiency gain in rotating window
    ###  
    
    if (!autorotate){
        test <- roteffgain(xyt,cellwidth)            
    }
    else{
        test <- roteffgain(xyt,cellwidth)
        if (!test){
            warning("There is no gain in efficiency by rotating, see ?roteffgain",immediate.=TRUE)
            cat("Not rotating window.\n")
        }
        else{
            rotmat <- getRotation(xyt)$rotation
            xyt <- affine(xyt,mat=rotmat)
        }
    }
    ow <- selectObsWindow(xyt,cellwidth)
	xyt <- ow$xyt
	M <- ow$M
	N <- ow$N
	if(!is.null(missing.data.areas)){
	    if(autorotate){
	        if(test){ # only bother if it was worthwhile rotating
        	    missing.data.areas <- lapply(missing.data.areas,affine,mat=rotmat)
        	    lapply(1:numt,function(i){missing.data.areas[[i]]$xrange<<-xyt$window$xrange;missing.data.areas[[i]]$yrange<<-xyt$window$yrange})
    	    }
	    }
	}
	
	tst <- mget("lgcpPredictstapptriggertestvalue",envir=parent.frame(),ifnotfound=FALSE)$lgcpPredictstapptriggertestvalue
	if (tst){
	    del1 <- (xyt$window$xrange[2]-xyt$window$xrange[1])/M
	    del2 <- (xyt$window$yrange[2]-xyt$window$yrange[1])/N 
	    mcens <- xyt$window$xrange[1]+.5*del1+(0:(M-1))*del1
	    ncens <- xyt$window$yrange[1]+.5*del2+(0:(N-1))*del2
        xls <- rep(mcens,N)
        yls <- rep(ncens,each=M)	    
	    spdf <- get("app",envir=parent.frame())$spdf
	    #EJP: olay <- overlay(SpatialPoints(cbind(xls,yls)),spdf)
	    olay <- over(SpatialPoints(cbind(xls,yls)), geometry(spdf))
	    if(length(table(olay))!=length(spdf)){
	        cat("\n")
	        warning("With chosen cell width, will not be able to generate aggregated inference for all regions.",.immediate=TRUE)
	        cat("\n")
	    }
	    olay[is.na(olay)] <- 0
	    olay <- matrix(olay,M,N)
	}
	
	if (M*N>=(256^2)){
	    Sys.sleep(1)
	    cat("\n")
	    warning("USING LARGE FFT GRID: COMPUTATION MAY BE SLOW ON SOME MACHINES ...",.immediate=TRUE)
	    cat("\n")
	}  
	
	cat(paste("FFT Grid size: [",ext*M," , ",ext*N,"]\n",sep=""))
	Sys.sleep(1)
    rm(ow)
    
    ###
    # Deal with spatial component and rotate, if necessary
    ###
					
	if(!any(class(spatial.intensity)=="spatialAtRisk")){			
        spatial <- spatialAtRisk(spatial.intensity)
    }
    else{
        spatial <- spatial.intensity
    }
    
    if(autorotate){
        if(test){
            spatial <- affine(spatial,mat=rotmat)
        }
    }
    
    ################################################################
    # Create grid and FFT objects
    ################################################################   
    
    study.region <- xyt$window
	
	## DEFINE LATTICE & CENTROIDS ##
		
	del1 <- (study.region$xrange[2]-study.region$xrange[1])/M
	del2 <- (study.region$yrange[2]-study.region$yrange[1])/N 
	
	Mext <- ext*M	
	Next <- ext*N
	
	mcens <- study.region$xrange[1]+.5*del1+(0:(Mext-1))*del1
	ncens <- study.region$yrange[1]+.5*del2+(0:(Next-1))*del2	
	
	cellarea <- del1*del2
	
	if(is.null(missing.data.areas)){
    	if(inclusion=="centroid"){
            cellInside <- inside.owin(x=rep(mcens,Next),y=rep(ncens,each=Mext),w=study.region)
        }        
        else if(inclusion=="touching"){
            cellInside <- touchingowin(x=mcens,y=ncens,w=study.region)
        }
        else{
            stop("Invlaid choice for argument 'inclusion'.")
        }
    	cellInside <- matrix(as.numeric(cellInside),Mext,Next)
    	cellInside <- rep(list(cellInside),numt)
	}
	else{
	    cellInside <- list()
	    for (i in 1:numt){
	        #
	        
	        if(inclusion=="centroid"){
                cellInside[[i]] <- inside.owin(x=rep(mcens,Next),y=rep(ncens,each=Mext),w=missing.data.areas[[i]])
            }            
            else if(inclusion=="touching"){
                cellInside[[i]] <- touchingowin(x=mcens,y=ncens,w=missing.data.areas[[i]])
            }
            else{
                stop("Invlaid choice for argument 'inclusion'.")
            }
    	    cellInside[[i]] <- matrix(as.numeric(cellInside[[i]]),Mext,Next)
	    }
	}
	
	## OBTAIN SPATIAL VALS ON LATTICE (LINEAR INTERPOLATION) ##
	
	if(is.null(missing.data.areas)){
    	spatialvals <- fftinterpolate(spatial,mcens,ncens,ext=ext)
    	spatialvals <- spatialvals*cellInside[[1]]
    	spatialvals <- spatialvals / (cellarea*sum(spatialvals))
    	spatialvals <- rep(list(spatialvals),numt)
	}
	else{
        if(inclusion=="centroid"){
            cellIns <- inside.owin(x=rep(mcens,Next),y=rep(ncens,each=Mext),w=study.region) # for purposes of computing normalising constant of lambda
        }        
        else if(inclusion=="touching"){
            cellIns <- touchingowin(x=mcens,y=ncens,w=study.region)
        }
        else{
            stop("Invlaid choice for argument 'inclusion'.")
        }
    	cellIns <- matrix(as.numeric(cellIns),Mext,Next)
	    spatialinterp <- fftinterpolate(spatial,mcens,ncens,ext=ext)
	    tempinterp <- spatialinterp*cellIns
	    NC <- cellarea*sum(tempinterp) # gives normalising constant for lambda over the whole observation window, with no missing areas (compare with the version where missing.data.area is null above)
	    spatialvals <- list()
	    for (i in 1:numt){ 
        	spatialvals[[i]] <- spatialinterp*cellInside[[i]]
        	spatialvals[[i]] <- spatialvals[[i]] / NC
	    }
	}
	
	# compute the base matrix of the covariance matrix
    bcb <- blockcircbase(x=mcens,y=ncens,sigma=sigma,phi=phi,model=spatial.covmodel,additionalparameters=covpars)
    
    Qeigs <- eigenfrombase(inversebase(bcb)) # eigenvalues of Q (the precision matrix)
    rqe <- sqrt(Qeigs) # square root of the eigenvalues (used in computation)
    irqe <- 1/rqe # reciprocal root (commputation)
    	
	################################################################
		
	
		
    ###
    # Set up MCMC loop, required to compute nsamp, below
    ###
    mLoop = mcmcLoop(N=mcmc.control$mala.length,burnin=mcmc.control$burnin,thin=mcmc.control$retain,progressor=mcmcProgressTextBar)	
	
	# issue warning if dumping information to disc
	nsamp <- floor((mLoop$N-mLoop$burnin)/mLoop$thin)
	if (!is.null(output.control$gridfunction) & class(output.control$gridfunction)[1]=="dump2dir"){
    	cat("WARNING: disk space required for saving is approximately ",round(nsamp*object.size(array(runif(M*N),dim=c(M,N)))/1024^2,2)," Mb, ",sep="")
        if (!output.control$gridfunction$forceSave){
            m <- menu(c("yes","no"),title="continue?")
            if(m==1){
                cat("Note: to bypass this menu, set forceSave=TRUE in dump2dir\n")
                Sys.sleep(2)
            }
            else{
                stop("Stopped")
            }
        }
    }
	
	nis <- list()
	for(i in 1:numt){
	    if (sum(xyt$t==aggtimes[i])>0){
		    nis[[i]] <- getCounts(xyt=xyt,subset=(xyt$t==aggtimes[i]),M=M,N=N,ext=ext)
		}
		else{
		    nis[[i]] <- matrix(0,ext*M,ext*N)
		}
		ct1 <- sum(nis[[i]])
		nis[[i]] <- nis[[i]] * (spatialvals[[i]]>0)
		ct2 <- sum(nis[[i]])
		if(ct2<ct1){
		    warning(paste("Time ",aggtimes[i],": ",ct1-ct2," data points lost due to discretisation.",sep=""),immediate.=TRUE)
		}
	}	
	
	###
	# Compute gradient truncation, if necessary
	###
    
    if(is.null(gradtrunc)){
        gradtrunc <- computeGradtruncSpatioTemporal(nsims=100,
                                                    scale=1,
                                                    nis=nis,
                                                    mu=mu,
                                                    rootQeigs=rqe,
                                                    invrootQeigs=irqe,
                                                    spatial=spatialvals,
                                                    temporal=temporalfit,
                                                    bt=bt,
                                                    cellarea=cellarea)
    }
	
	###
	# Run MALA
	##
	
	gridfun <- output.control$gridfunction
	if (is.null(gridfun)){
	    gridfun <- nullFunction()
	}
    gridav <- output.control$gridmeans
	if (is.null(gridav)){
	    gridav <- nullAverage()
	}    
    
    lg <- MALAlgcp( mcmcloop=mLoop,
                    inits=mcmc.control$inits,
                    adaptivescheme=mcmc.control$adaptivescheme,
                    M=M,
                    N=N,
                    Mext=Mext,
                    Next=Next,
                    sigma=sigma,
                    phi=phi,
                    theta=theta,
                    mu=mu,
                    nis=nis,
                    cellarea=cellarea,
                    spatialvals=spatialvals,
                    temporal.fitted=temporalfit,
                    tdiff=tdiff,
                    rootQeigs=rqe,
                    invrootQeigs=irqe,
                    cellInside=cellInside,
                    MCMCdiag=mcmc.control$MCMCdiag,
                    gradtrunc=gradtrunc,
                    gridfun=gridfun,
                    gridav=gridav,
                    mcens=mcens,
                    ncens=ncens,
                    aggtimes=aggtimes)
                                                     
	
	endtime <- Sys.time()
	timetaken <- endtime-starttime
	
	lg$xyt <- xyt
	lg$M <- M
	lg$N <- N
	lg$aggtimes <- aggtimes
	lg$tdiffs <- tdiff
	lg$vars <- bt
	lg$spatial <- spatial
	lg$temporal <- temporalfit
	lg$grid <- spatialvals
	lg$nis <- lgcpgrid(nis,xvals=mcens[1:M],yvals=ncens[1:N],zvals=aggtimes)
	lg$mcens <- mcens[1:M]
	lg$ncens <- ncens[1:N]
	lg$cellarea <- diff(mcens[1:2]) * diff(ncens[1:2]) 
	lg$sigma <- sigma
	lg$phi <- phi
	lg$mu <- mu
	lg$theta <- theta
	lg$mcmcpars <- mcmc.control
	lg$timetaken <- timetaken
	lg$ext <- ext
	lg$cellInside <- lapply(cellInside,function(x){x[1:M,1:N]})
	lg$spatialonly <- FALSE
	lg$inclusion <- inclusion	
	if (tst){
	    lg$overlay <- olay
	}
	
	class(lg) <- c("lgcpPredict","lgcpobject")	
	
	return(lg)															
					
}	
	

##' MALAlgcp function
##'
##' ADVANCED USE ONLY A function to perform MALA for the spatial only case
##'
##' @param mcmcloop an mcmcLoop object
##' @param inits initial values from mcmc.control
##' @param adaptivescheme adaptive scheme from mcmc.control
##' @param M number of cells in x direction on output grid
##' @param N number of cells in y direction on output grid
##' @param Mext number of cells in x direction on extended output grid
##' @param Next number of cells in y direction on extended output grid
##' @param sigma spatial covariance parameter sigma
##' @param phi spatial covariance parameter phi
##' @param theta temporal correlation parameter theta
##' @param mu spatial covariance parameter mu
##' @param nis cell counts matrix
##' @param cellarea area of cells
##' @param spatialvals spatial at risk, function lambda, interpolated onto the requisite grid
##' @param temporal.fitted temporal fitted values representing mu(t)
##' @param tdiff vecto of time differences with convention that the first element is Inf
##' @param scaleconst expected number of observations
##' @param rootQeigs square root of eigenvalues of precision matrix
##' @param invrootQeigs inverse square root of eigenvalues of precision matrix
##' @param cellInside logical matrix dictating whether cells are inside the observation window
##' @param MCMCdiag defunct
##' @param gradtrunc gradient truncation parameter
##' @param gridfun grid functions
##' @param gridav grid average functions
##' @param mcens x-coordinates of cell centroids 
##' @param ncens y-coordinates of cell centroids
##' @param aggtimes z-coordinates of cell centroids (ie time)
##' @return object passed back to lgcpPredictSpatial
##' @export
MALAlgcp <- function(mcmcloop,
                            inits,
                            adaptivescheme,
                            M,
                            N,
                            Mext,
                            Next,
                            sigma,
                            phi,
                            theta,
                            mu,
                            nis,
                            cellarea,
                            spatialvals,
                            temporal.fitted,
                            tdiff,
                            scaleconst,
                            rootQeigs,
                            invrootQeigs,
                            cellInside,
                            MCMCdiag,
                            gradtrunc,
                            gridfun,
                            gridav,
                            mcens,
                            ncens,
                            aggtimes){
                            
    SpatialOnlyMode <- FALSE
    ##ImprovedAlgorithm <- TRUE 
    SpatialPlusParameters <- FALSE
    SpatioTemporalPlusParameters <- FALSE   
    MultiTypeMode <- FALSE
    
    nlevs <- NULL # note this line is here for gridFunction and gridAverage methods and is not used otherwise
    
    n <- length(temporal.fitted)
    gt <- 1-exp(-2*theta*tdiff) # note tdiff[1] defined as Inf, for psimplifying code
    bt <- exp(-theta*tdiff)
    bt[1] <- NA # again, a computational assumption for simplifying code
      

    cellOutside <- lapply(cellInside,function(x){!as.logical(x)})
    logspatial <- lapply(1:n,function(i){log(temporal.fitted[i]*spatialvals[[i]])})
    lapply(1:n,function(i){logspatial[[i]][cellOutside[[i]] | spatialvals[[i]]==0] <<- 0}) # NOTE THIS IS FOR SIMPLIFYING THE COMPUTATION OF THE TARGET!!                               
                                               
    GFinitialise(gridfun) # note these two lines must come after M and N have been computed or defined
	GAinitialise(gridav) # note these two lines must come after M and N have been computed or defined
    
    h <- initialiseAMCMC(adaptivescheme)
    hrec <- h
    nsamp <- 0
    icount <- 0
    MCMCacc <- 0
    y.mean <- rep(list(matrix(0,M,N)),n)
    y.var <- rep(list(matrix(0,M,N)),n)
    EY.mean <- rep(list(matrix(0,M,N)),n)
    EY.var <- rep(list(matrix(0,M,N)),n)
     	    
    Gamma <- rep(list(matrix(0,Mext,Next)),n)  # initialise with mean                         
    oldtags <- target.and.grad.spatiotemporal(Gamma=Gamma,nis=nis,cellarea=cellarea,rootQeigs=rootQeigs,invrootQeigs=invrootQeigs,mu=mu,spatial=spatialvals,logspat=logspatial,temporal=temporal.fitted,bt=bt,gt=gt,gradtrunc=gradtrunc)    
    
    while(nextStep(mcmcloop)){
    
        propmeans <- lapply(1:n,function(i){Gamma[[i]] + (h/2)*oldtags$grad[[i]]})
        
        propGamma <- lapply(1:n,function(i){propmeans[[i]] + sqrt(h)*rnorm(Mext*Next)})
        proptags <- target.and.grad.spatiotemporal(Gamma=propGamma,nis=nis,cellarea=cellarea,rootQeigs=rootQeigs,invrootQeigs=invrootQeigs,mu=mu,spatial=spatialvals,logspat=logspatial,temporal=temporal.fitted,bt=bt,gt=gt,gradtrunc=gradtrunc)
        revpropmeans <- lapply(1:n,function(i){propGamma[[i]] + (h/2)*proptags$grad[[i]]})
        
        ac <- exp(proptags$logtarget-oldtags$logtarget-sum(sapply(1:n,function(i){sum((Gamma[[i]]-revpropmeans[[i]])^2)}))/(2*h) + sum(sapply(1:n,function(i){sum((propGamma[[i]]-propmeans[[i]])^2)}))/(2*h))
        ac <- min(ac,1)        
        
        icount <- icount + 1        
        MCMCacc <- ((icount-1)/icount) * MCMCacc + ac/icount
        
        ###if(iteration(mcmcloop)==5000){browser()}
        ###cat("\n",h,"\n")
        ###browser() 
             
        
        trigger <- FALSE
		if (proptags$logtarget==-Inf | is.na(ac) | is.nan(ac)){ # gradient truncation insufficient, so reduce
	        gradtrunc <- gradtrunc/2
	        cat("Reducing gradient truncation to:",gradtrunc,"\n")
	        oldtags <- target.and.grad.spatiotemporal(Gamma=Gamma,nis=nis,cellarea=cellarea,rootQeigs=rootQeigs,invrootQeigs=invrootQeigs,mu=mu,spatial=spatialvals,logspat=logspatial,temporal=temporal.fitted,bt=bt,gt=gt,gradtrunc=gradtrunc) 	        
	        if (!is.burnin(mcmcloop)){
	            cat("Gradient truncation currently",gradtrunc,"\n")
	            cat("Suggest reducing this further and setting gradtrunc manually in lgcpPredictSpatial.\n")
	            cat("Alternatively, a consistently failing algorithm can also suggest problems with the compatibility of the observed data with the spatialAtRisk, temporalAtRisk and observation window components\n")
	            stop(paste("Problem with gradient truncation at iteration",iteration(mcmcloop),"acceptance probability =",ac))
	        }
	        ac <- 0 # don't accept the move if a suitable gradient truncation has not been found.
	        trigger <- TRUE # this is set to true so that the adaptive scheme is halted for one iteration if a suitable gradient trunctation has not been found
	    } 
        
        if (ac>runif(1)){
            Gamma <- propGamma
            oldtags <- proptags
        }
    
        if (iteration(mcmcloop)>1){
	        hrec <- c(hrec,h)
	    }
	    
	    trigger <- FALSE
	    if (!trigger){ # ie if there was no problem with gradient truncation
	        h <- updateAMCMC(adaptivescheme)
        }
        ##cat("\n")
        ##cat(h,"\n")
        ##cat(oldtags$logtarget)
        ##cat("\n")
        
        if (is.retain(mcmcloop)){
	        nsamp <- nsamp + 1
        	y.mean <- lapply(1:n,function(i){((nsamp-1)/nsamp) * y.mean[[i]] + oldtags$Y[[i]][1:M,1:N]/nsamp})
        	EY.mean <- lapply(1:n,function(i){((nsamp-1)/nsamp) * EY.mean[[i]] + oldtags$expY[[i]][1:M,1:N]/nsamp})
        	if (nsamp>1){
    			y.var <- lapply(1:n,function(i){((nsamp-2)/(nsamp-1))*y.var[[i]] + (nsamp/(nsamp-1)^2)*(y.mean[[i]]-oldtags$Y[[i]][1:M,1:N])^2})
    			EY.var <- lapply(1:n,function(i){((nsamp-2)/(nsamp-1))*EY.var[[i]] + (nsamp/(nsamp-1)^2)*(EY.mean[[i]]-oldtags$expY[[i]][1:M,1:N])^2})
    		}                 	               
			GFupdate(gridfun)
			GAupdate(gridav)		    
		}
    } 
	
	retlist <- list(lasth=h,lastGAM=lgcpgrid(Gamma,xvals=mcens[1:M],yvals=ncens[1:N],zvals=aggtimes))
	
	GFfinalise(gridfun) # these two lines must appear after retlist has been initialised
	GAfinalise(gridav)  #	
	
	retlist$mcmcacc <- MCMCacc
	retlist$hrec <- hrec
    retlist$y.mean <- lgcpgrid(y.mean,xvals=mcens[1:M],yvals=ncens[1:N],zvals=aggtimes)
    retlist$y.var <- lgcpgrid(y.var,xvals=mcens[1:M],yvals=ncens[1:N],zvals=aggtimes)
    retlist$EY.mean <- lgcpgrid(EY.mean,xvals=mcens[1:M],yvals=ncens[1:N],zvals=aggtimes)
    retlist$EY.var <- lgcpgrid(EY.var,xvals=mcens[1:M],yvals=ncens[1:N],zvals=aggtimes)
    retlist$gridfunction <- GFreturnvalue(gridfun)
    retlist$gridaverage <- GAreturnvalue(gridav)
    retlist$mcmcinfo <- mcmcloop
    retlist$gradtrunc <- gradtrunc	
	
	return(retlist)                          
}
                        

##' target.and.grad.spatiotemporal function
##'
##' A function to compute the target and gradient for 'spatial only' MALA
##'
##' @param Gamma current state of the chain, Gamma
##' @param nis matrix of cell counts
##' @param cellarea area of cells, a positive number
##' @param rootQeigs square root of the eigenvectors of the precision matrix
##' @param invrootQeigs inverse square root of the eigenvectors of the precision matrix
##' @param mu parameter of the latent Gaussian field
##' @param spatial spatial at risk function, lambda, interpolated onto correct grid
##' @param logspat log of spatial at risk function, lambda*scaleconst, interpolated onto correct grid
##' @param temporal fitted temoporal values
##' @param bt in Brix and Diggle vector b(delta t)
##' @param gt in Brix and Diggle vector g(delta t) (ie the coefficient of R in G(t)), with convention that (deltat[1])=Inf
##' @param gradtrunc gradient truncation parameter
##' @return the back-transformed Y, its exponential, the log-target and gradient for use in MALAlgcp 
##' @export
target.and.grad.spatiotemporal <- function(Gamma,nis,cellarea,rootQeigs,invrootQeigs,mu,spatial,logspat,temporal,bt,gt,gradtrunc){
    
    ncells <- length(Gamma[[1]])
    tcons <- 1/gt
    n <- length(tcons)
    Y <- list() 
    expY <- list() 
    gradmult <- list()
    gradcomp <- list()
    grad <- list()
    ml <- get("mcmcloop",envir=parent.frame())
    iter <- iteration(ml)
    gradfun <- function(i){
        Y[[i]] <<- YfromGamma(Gamma[[i]],invrootQeigs=invrootQeigs,mu=mu)
        expY[[i]] <<- exp(Y[[i]])
        expYtrunc <- expY[[i]]
        expYtrunc[expYtrunc>gradtrunc] <- gradtrunc # gradient truncation
        grad[[i]] <<- (-1)*tcons[i]*Gamma[[i]]
        gradmult[[i]] <<- as.list(sqrt(gt[i])*c(rep(0,i-1),c(1,cumprod(bt[-c(1,1:i)])))) # computes the product of sqrt(gt(i)) * beta(delta t) terms in the summation
        gradcomp[[i]] <<- nis[[i]]-temporal[i]*spatial[[i]]*expYtrunc*cellarea
    }
    sapply(1:n,gradfun)

    for(i in 1:n){
        gradcum <- 0
        for(j in i:n){
            gradcum <- gradcum + gradmult[[i]][[j]] * gradcomp[[j]]
        }
        grad[[i]] <- grad[[i]] + (1/ncells)*Re(fft(invrootQeigs*fft(gradcum,inverse=TRUE)))
    }
    
    ###logtarget <- -(1/2)*sum(tcons*sapply(Gamma,function(x){sum(x*x)})) + sum(sapply(1:n,function(i){sum((Y[[i]]+log(temporal[i])+logspat[[i]])*nis[[i]] - temporal[i]*spatial[[i]]*expY[[i]]*cellarea)})) # note that both nis=0, logspat=0 and spatial=0 outside of observation window, so this effectively limits summation to the observation window only
    logtarget <- -(1/2)*sum(tcons*sapply(Gamma,function(x){sum(x*x)})) + sum(sapply(1:n,function(i){sum(Y[[i]]*nis[[i]] - temporal[i]*spatial[[i]]*expY[[i]]*cellarea)})) # ... a more computationally efficient way of doing this # note that both nis=0, logspat=0 and spatial=0 outside of observation window, so this effectively limits summation to the observation window only
    
    return(list(Y=Y,expY=expY,logtarget=logtarget,grad=grad))
}

##' lgcpPredictAggregated function
##'
##' The function \code{lgcpPredict} performs spatiotemporal prediction for log-Gaussian Cox Processes for point process data where counts
##' have been aggregated to the regional level. This is achieved by imputation of the regional counts onto a spatial continuum; if something
##' is known about the underlying spatial density of cases, then this information can be added to improve the quality of the imputation, 
##' without this, the counts are distributed uniformly within regions.
##'
##' The following is a mathematical description of a log-Gaussian Cox Process, it is best viewed in the pdf version of the manual.
##'
##' Let \eqn{\mathcal Y(s,t)}{\mathcal Y(s,t)} be a spatiotemporal Gaussian process, \eqn{W\subset R^2}{W\subset R^2} be an 
##' observation window in space and \eqn{T\subset R_{\geq 0}}{T\subset R_{\geq 0}} be an interval of time of interest. 
##' Cases occur at spatio-temporal positions \eqn{(x,t) \in W \times T}{(x,t) \in W \times T} 
##'  according to an inhomogeneous spatio-temporal Cox process,
##' i.e. a Poisson process with a stochastic intensity \eqn{R(x,t)}{R(x,t)},
##'   The number of cases, \eqn{X_{S,[t_1,t_2]}}{X_{S,[t_1,t_2]}}, arising in 
##'   any \eqn{S \subseteq W}{S \subseteq W} during the interval \eqn{[t_1,t_2]\subseteq T}{[t_1,t_2]\subseteq T} is 
##'   then Poisson distributed conditional on \eqn{R(\cdot)}{R(\cdot)},
##' \deqn{X_{S,[t_1,t_2]} \sim \mbox{Poisson}\left\{\int_S\int_{t_1}^{t_2} R(s,t)d sd t\right\}}{X_{S,[t_1,t_2]} \sim \mbox{Poisson}\left\{\int_S\int_{t_1}^{t_2} R(s,t)d sd t\right\}.}
##' Following Brix and Diggle (2001) and Diggle et al (2005), the intensity is decomposed multiplicatively as
##' \deqn{R(s,t) = \lambda(s)\mu(t)\exp\{\mathcal Y(s,t)\}.}{R(s,t) = \lambda(s)\mu(t)Exp\{\mathcal Y(s,t)\}.}
##' In the above, the fixed spatial component, \eqn{\lambda:R^2\mapsto R_{\geq 0}}{\lambda:R^2\mapsto R_{\geq 0}}, 
##' is a known function, proportional to the population at risk at each point in space and scaled so that
##' \deqn{\int_W\lambda(s)d s=1,}{\int_W\lambda(s)d s=1,}
##' whilst the fixed temporal component, 
##'  \eqn{\mu:R_{\geq 0}\mapsto R_{\geq 0}}{\mu:R_{\geq 0}\mapsto R_{\geq 0}}, is also a known function with
##' \deqn{\mu(t) \delta t = E[X_{W,\delta t}],}{\mu(t) \delta t = E[X_{W,\delta t}],}
##' for \eqn{t}{t} in a small interval of time, \eqn{\delta t}{\delta t}, over which the rate of the process over \eqn{W}{W} can be considered constant.
##'
##' \bold{
##'     NOTE: the xyt stppp object can be recorded in continuous time, but for the purposes of prediciton,    
##'     discretisation must take place. For the time dimension, this is achieved invisibly by \code{as.integer(xyt$t)} and
##'     \code{as.integer(xyt$tlim)}. Therefore, before running an analysis please make sure that this is commensurate
##'     with the physical inerpretation and requirements of your output. The spatial discretisation is
##'     chosen with the argument cellwidth (or gridsize). If the chosen discretisation in time and space is too coarse for a
##'     given set of parameters (sigma, phi and theta) then the proper correlation structures implied by the model will not
##'     be captured in the output.
##' }
##'
##' Before calling this function, the user must decide on the time point of interest, the
##' number of intervals of data to use, the parameters, spatial covariance model, spatial discretisation,
##' fixed spatial (\eqn{\lambda(s)}{\lambda(s)}) and temporal (\eqn{\mu(t)}{\mu(t)}) components, mcmc parameters, and whether or not any output is
##' required.
##'
##' @param app a spatio-temporal aggregated point pattern object, see ?stapp
##' @param popden a spatialAtRisk object of class 'fromFunction' describing the population density, if known. Default is NULL, which gives a uniform density on each region.
##' @param T time point of interest
##' @param laglength specifies lag window, so that data from and including  time (T-laglength) to time T is used in the MALA algorithm
##' @param model.parameters values for parameters, see ?lgcppars
##' @param spatial.covmodel correlation type see ?CovarianceFct 
##' @param covpars vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct
##' @param cellwidth width of grid cells on which to do MALA (grid cells are square). Note EITHER gridsize OR cellwidthe must be specified.
##' @param gridsize size of output grid required. Note EITHER gridsize OR cellwidthe must be specified.
##' @param spatial.intensity the fixed spatial component: an object of that can be coerced to one of class spatialAtRisk
##' @param temporal.intensity the fixed temporal component: either a numeric vector, or a function that can be coerced into an object of class temporalAtRisk
##' @param mcmc.control MCMC paramters, see ?mcmcpars
##' @param output.control output choice, see ?setoutput
##' @param autorotate logical: whether or not to automatically do MCMC on optimised, rotated grid.   
##' @param gradtrunc truncation for gradient vector equal to H parameter Moller et al 1998 pp 473. Set to NULL to estimate this automatically (default). Set to zero for no gradient truncation.
##' @param n parameter for as.stppp. If popden is NULL, then this parameter controls the resolution of the uniform. Otherwise if popden is of class 'fromFunction', it controls the size of the imputation grid used for sampling. Default is 100.
##' @param dmin parameter for as.stppp. If any reginal counts are missing, then a set of polygonal 'holes' in the observation window will be computed for each. dmin is the parameter used to control the simplification of these holes (see ?simplify.owin). default is zero.
##' @param check logical parameter for as.stppp. If any reginal counts are missing, then roughly speaking, check specifies whether to check the 'holes'. 
##' further notes on autorotate argument: If set to TRUE, and the argument spatial is not NULL, then the argument spatial must be computed in the original frame of reference (ie NOT in the rotated frame). 
##' Autorotate performs bilinear interpolation (via interp.im) on an inverse transformed grid; if there is no computational advantage in doing this, a warning message will be issued. Note that best accuracy 
##' is achieved by manually rotating xyt and then computing spatial on the transformed xyt and finally feeding these in as arguments to the function lgcpPredict. By default autorotate is set to FALSE.
##' @return the results of fitting the model in an object of class \code{lgcpPredict}
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor, Tilman M. Davies, Barry S. Rowlingson, Peter J. Diggle (2013). Journal of Statistical Software, 52(4), 1-40. URL http://www.jstatsoft.org/v52/i04/
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##'     \item Wood ATA, Chan G (1994). Simulation of Stationary Gaussian Processes in [0,1]d. Journal of Computational and Graphical Statistics, 3(4), 409-432.
##'     \item Moller J, Syversveen AR, Waagepetersen RP (1998). Log Gaussian Cox Processes. Scandinavian Journal of Statistics, 25(3), 451-482.
##' }
##' @seealso \link{KinhomAverage}, \link{ginhomAverage}, \link{lambdaEst}, \link{muEst}, \link{spatialparsEst}, \link{thetaEst},  
##' \link{spatialAtRisk}, \link{temporalAtRisk}, \link{lgcppars}, \link{CovarianceFct}, \link{mcmcpars}, \link{setoutput} 
##' \link{print.lgcpPredict}, \link{xvals.lgcpPredict}, \link{yvals.lgcpPredict}, \link{plot.lgcpPredict}, \link{meanfield.lgcpPredict},
##' \link{rr.lgcpPredict}, \link{serr.lgcpPredict}, \link{intens.lgcpPredict},   
##' \link{varfield.lgcpPredict}, \link{gridfun.lgcpPredict}, \link{gridav.lgcpPredict}, \link{hvals.lgcpPredict}, \link{window.lgcpPredict},
##' \link{mcmctrace.lgcpPredict}, \link{plotExceed.lgcpPredict}, \link{quantile.lgcpPredict}, \link{identify.lgcpPredict}, \link{expectation.lgcpPredict},
##' \link{extract.lgcpPredict}, \link{showGrid.lgcpPredict}
##' @export 
    
lgcpPredictAggregated <- function(  app,
                                    popden=NULL,
            					    T,
            					    laglength,
            					    model.parameters=lgcppars(),
            					    spatial.covmodel="exponential",
            					    covpars=c(),
            					    cellwidth=NULL,
            					    gridsize=NULL,
            					    spatial.intensity,
            					    temporal.intensity,					
            					    mcmc.control,
            					    output.control=setoutput(),
            					    autorotate=FALSE,
        					        gradtrunc=NULL,
        					        n=100,
        					        dmin=0,
        					        check=TRUE){

    cat("Imputing data ...\n")
    xyt <- as.stppp(app,popden=popden,n=n,dmin=dmin,check=check) # generate imputed xyt object
    olay <- attr(xyt,"overlay")
    
    if(mget("lgcpPredictstapptriggertestvalue",envir=parent.frame(),ifnotfound="no-object")$lgcpPredictstapptriggertestvalue!="no-object"){
        stop("Please rename object 'lgcpPredictstapptriggertestvalue' before using lgcpPredict.stapp.")
    }
    lgcpPredictstapptriggertestvalue <- TRUE # used in lgcpPredict.stppp
    
    mda <- NULL
    if(!is.null(attr(xyt,"owinlist"))){
        mda <- attr(xyt,"owinlist")
    }
        					    
    lg <- lgcpPredict(  xyt=xyt,
					    T=T,
					    laglength=laglength,
					    model.parameters=model.parameters,
					    spatial.covmodel=spatial.covmodel,
					    covpars=covpars,
					    cellwidth=cellwidth,
					    gridsize=gridsize,
					    spatial.intensity=spatial.intensity,
					    temporal.intensity=temporal.intensity,					
					    mcmc.control=mcmc.control,
					    output.control=output.control,
					    missing.data.areas=mda[(T-laglength):T],
					    autorotate=autorotate,
					    gradtrunc=gradtrunc)
    
    lg$app <- app # the original aggregated data
    nreg <- length(lg$app$spdf) # number of regions
    lg$RegCounts <- sapply(1:nreg,function(x){sum(olay==x)}) # regional counts, some may be zero reflects accuracy of aggregated operations   					    
    attr(lg$xyt,"overlay") <- olay # since this attribute is stripped off inside lgcpPredict.stppp
    attr(lg$xyt,"owinlist") <- mda
    					           					    
    class(lg) <- c("aggregatedPredict","lgcpPredict","lgcpobject")
    return(lg)	        					    
}
