##' lgcpPredictSpatial function
##'
##' The function \code{lgcpPredictSpatial} performs spatial prediction for log-Gaussian Cox Processes
##'
##' The following is a mathematical description of a log-Gaussian Cox Process, it is best viewed in the pdf version of the manual.
##'
##' Let \eqn{\mathcal Y(s)}{\mathcal Y(s)} be a spatial Gaussian process and \eqn{W\subset R^2}{W\subset R^2} be an 
##' observation window in space. 
##' Cases occur at spatial positions \eqn{x \in W}{x \in W} 
##'  according to an inhomogeneous spatial Cox process,
##' i.e. a Poisson process with a stochastic intensity \eqn{R(x)}{R(x)},
##'   The number of cases, \eqn{X_{S}}{X_{S}}, arising in 
##'   any \eqn{S \subseteq W}{S \subseteq W} is 
##'   then Poisson distributed conditional on \eqn{R(\cdot)}{R(\cdot)},
##' \deqn{X_{S} \sim \mbox{Poisson}\left\{\int_S R(s)ds\right\}}{X_{S} \sim \mbox{Poisson}\left\{\int_S R(s)ds\right\}.}
##' Following Brix and Diggle (2001) and Diggle et al (2005) (but ignoring temporal variation), the intensity is decomposed multiplicatively as
##' \deqn{R(s,t) = \lambda(s)\exp\{\mathcal Y(s,t)\}.}{R(s,t) = \lambda(s)Exp\{\mathcal Y(s,t)\}.}
##' In the above, the fixed spatial component, \eqn{\lambda:R^2\mapsto R_{\geq 0}}{\lambda:R^2\mapsto R_{\geq 0}}, 
##' is a known function, proportional to the population at risk at each point in space and scaled so that
##' \deqn{\int_W\lambda(s)d s=1.}{\int_W\lambda(s)d s=1.}
##'
##'
##' Before calling this function, the user must decide on the parameters, spatial covariance model, spatial discretisation,
##' fixed spatial (\eqn{\lambda(s)}{\lambda(s)}) component, mcmc parameters, and whether or not any output is
##' required. Note there is no autorotate option for this function.
##'
##' @param sd a spatial point pattern object, see ?ppp
##' @param model.parameters values for parameters, see ?lgcppars
##' @param spatial.covmodel correlation type see ?CovarianceFct 
##' @param covpars vector of additional parameters for certain classes of covariance function (eg Matern), these must be supplied in the order given in ?CovarianceFct
##' @param cellwidth width of grid cells on which to do MALA (grid cells are square) in same units as observation window. Note EITHER gridsize OR cellwidthe must be specified.
##' @param gridsize size of output grid required. Note EITHER gridsize OR cellwidthe must be specified.
##' @param spatial.intensity the fixed spatial component: an object of that can be coerced to one of class spatialAtRisk
##' @param spatial.offset Numeric of length 1. Optional offset parameter, corresponding to the expected number of cases. NULL by default, in which case, this is estimateed from teh data.
##' @param mcmc.control MCMC paramters, see ?mcmcpars
##' @param output.control output choice, see ?setoutput   
##' @param gradtrunc truncation for gradient vector equal to H parameter Moller et al 1998 pp 473. Default is Inf, which means no gradient truncation. Set to NULL to estimate this automatically (though note that this may not necessarily be a good choice). The default seems to work in most settings.
##' @param ext integer multiple by which grid should be extended, default is 2. Generally this will not need to be altered, but if the spatial correlation decays slowly, increasing 'ext' may be necessary.
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former, the default, includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window. 
##' @return the results of fitting the model in an object of class \code{lgcpPredict}
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor, Tilman M. Davies, Barry S. Rowlingson, Peter J. Diggle (2013). Journal of Statistical Software, 52(4), 1-40. URL http://www.jstatsoft.org/v52/i04/
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##'     \item Wood ATA, Chan G (1994). Simulation of Stationary Gaussian Processes in [0,1]d. Journal of Computational and Graphical Statistics, 3(4), 409-432.
##'     \item Moller J, Syversveen AR, Waagepetersen RP (1998). Log Gaussian Cox Processes. Scandinavian Journal of Statistics, 25(3), 451-482.
##' }
##' @seealso \link{lgcpPredict} \link{KinhomAverage}, \link{ginhomAverage}, \link{lambdaEst}, \link{muEst}, \link{spatialparsEst}, \link{thetaEst},  
##' \link{spatialAtRisk}, \link{temporalAtRisk}, \link{lgcppars}, \link{CovarianceFct}, \link{mcmcpars}, \link{setoutput} 
##' \link{print.lgcpPredict}, \link{xvals.lgcpPredict}, \link{yvals.lgcpPredict}, \link{plot.lgcpPredict}, \link{meanfield.lgcpPredict},
##' \link{rr.lgcpPredict}, \link{serr.lgcpPredict}, \link{intens.lgcpPredict},   
##' \link{varfield.lgcpPredict}, \link{gridfun.lgcpPredict}, \link{gridav.lgcpPredict}, \link{hvals.lgcpPredict}, \link{window.lgcpPredict},
##' \link{mcmctrace.lgcpPredict}, \link{plotExceed.lgcpPredict}, \link{quantile.lgcpPredict}, \link{identify.lgcpPredict}, \link{expectation.lgcpPredict},
##' \link{extract.lgcpPredict}, \link{showGrid.lgcpPredict}
##' @export 
    
lgcpPredictSpatial <- function( sd,
        					    model.parameters=lgcppars(),
        					    spatial.covmodel="exponential",
        					    covpars=c(),
        					    cellwidth=NULL,
        					    gridsize=NULL,
        					    spatial.intensity,
        					    spatial.offset=NULL,				
        					    mcmc.control,
        					    output.control=setoutput(),
        					    gradtrunc=Inf,
        					    ext=2,
        					    inclusion="touching"){
    
    starttime <- Sys.time()
    
    GMRF <- FALSE
    if(!is.null(attr(sd,"covbase"))){
        GMRF <- TRUE
    }
    
    if(sd$window$type=="rectangle"){
	    sd$window <- as.polygonal(sd$window)
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
	if(!is.null(gridsize)){
	    approxcw <- diff(sd$window$xrange)/gridsize[1] # approx cell width
	    cwseq <- seq(approxcw/2,2*approxcw,length.out=500)
	    cwfun <- function(cw){
	        ow <- selectObsWindow(sd,cw)
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
	
	if(mcmc.control$burnin>mcmc.control$mala.length){
		stop("Number of burnin iterations must be less than the total number of iterations")
	}
	
	
	###
	# Initialise
	###		

	sigma <- model.parameters$sigma
	phi <- model.parameters$phi
	mu <- model.parameters$mu
	theta <- model.parameters$theta
	if(!is.null(spatial.offset)){
	    scaleconst <- spatial.offset 
	}
	else{
	    scaleconst <- sd$n/exp(mu+sigma^2/2) # ML estimate of scaling constant
    }
	
    ow <- selectObsWindow(sd,cellwidth) 
	sd <- ow$xyt
	M <- ow$M
	N <- ow$N
	
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
    
    ################################################################
    # Create grid and FFT objects
    ################################################################
    
    study.region <- sd$window
	
	## DEFINE LATTICE & CENTROIDS ##
		
	del1 <- (study.region$xrange[2]-study.region$xrange[1])/M
	del2 <- (study.region$yrange[2]-study.region$yrange[1])/N 
	
	Mext <- ext*M	
	Next <- ext*N
	
	mcens <- study.region$xrange[1]+.5*del1+(0:(Mext-1))*del1
	ncens <- study.region$yrange[1]+.5*del2+(0:(Next-1))*del2	
	
	cellarea <- del1*del2
	
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
	
	## OBTAIN SPATIAL VALS ON LATTICE (LINEAR INTERPOLATION) ##
	
	spatialvals <- fftinterpolate(spatial,mcens,ncens,ext=ext)
	spatialvals <- spatialvals*cellInside
	spatialvals <- spatialvals / (cellarea*sum(spatialvals))
	#### NOT NECESSARY spatialvals[cellInside & spatialvals==0] <- 1e-200 # impute a very small number into cells inside the observation window with apparently zero risk
	
	# compute the base matrix of the covariance matrix
	if(GMRF){ # this is invoked if data have been simulated from lgcpSimSpatialGMRF
        bcb <- attr(sd,"covbase")
        Qeigs <- eigenfrombase(attr(bcb,"precbase"))
    }
    else{
        bcb <- blockcircbase(x=mcens,y=ncens,sigma=sigma,phi=phi,model=spatial.covmodel,additionalparameters=covpars)    
        Qeigs <- eigenfrombase(inversebase(bcb)) # eigenvalues of Q (the precision matrix)
    }
    
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
    	cat("WARNING: disk space required for saving is approximately ",round(nsamp*object.size(array(runif(N*N),dim=c(M,N)))/1024^2,2)," Mb, ",sep="")
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
	
	nis  <- getCounts(xyt=sd,M=M,N=N,ext=ext) 
    ct1 <- sum(nis)
	nis <- nis * cellInside
	ct2 <- sum(nis)
	if(ct2<ct1){
	    warning(paste(ct1-ct2," data points lost due to discretisation.",sep=""),immediate.=TRUE)
	}	
	
	###
	# Compute gradient truncation, if necessary
	###
    
    if(is.null(gradtrunc)){
        gradtrunc <- computeGradtruncSpatial(   nsims=100,
                                                scale=1,
                                                nis=nis,
                                                mu=mu,
                                                rootQeigs=rqe,
                                                invrootQeigs=irqe,
                                                scaleconst=scaleconst,
                                                spatial=spatialvals,
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
    
    lg <- MALAlgcpSpatial(  mcmcloop=mLoop,
                            inits=mcmc.control$inits,
                            adaptivescheme=mcmc.control$adaptivescheme,
                            M=M,
                            N=N,
                            Mext=Mext,
                            Next=Next,
                            sigma=sigma,
                            phi=phi,
                            mu=mu,
                            nis=nis,
                            cellarea=cellarea,
                            spatialvals=spatialvals,
                            scaleconst=scaleconst,
                            rootQeigs=rqe,
                            invrootQeigs=irqe,
                            cellInside=cellInside,
                            MCMCdiag=mcmc.control$MCMCdiag,
                            gradtrunc=gradtrunc,
                            gridfun=gridfun,
                            gridav=gridav,
                            mcens=mcens,
                            ncens=ncens)
                                                  
	
	endtime <- Sys.time()
	timetaken <- difftime(endtime,starttime,units="mins")
	
	lg$xyt <- sd
	lg$M <- M
	lg$N <- N
	lg$aggtimes <- NA
	lg$tdiffs <- NA
	lg$vars <- NA
	lg$spatial <- spatial
	lg$temporal <- scaleconst
	lg$grid <- list(spatialvals)
	lg$nis <- lgcpgrid(nis,xvals=mcens[1:M],yvals=ncens[1:N])
	lg$mcens <- mcens[1:M]
	lg$ncens <- ncens[1:N]
	lg$sigma <- sigma
	lg$phi <- phi
	lg$mu <- mu
	lg$theta <- theta
	lg$mcmcpars <- mcmc.control
	lg$timetaken <- timetaken
	lg$spatialonly <- TRUE
	lg$ext <- ext
	lg$cellInside <- cellInside[1:M,1:N]
	lg$inclusion <- inclusion
	
	class(lg) <- c("lgcpPredict","lgcpobject")	
	
	return(lg)															
					
}	
	

##' MALAlgcpSpatial function
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
##' @param mu spatial covariance parameter mu
##' @param nis cell counts matrix
##' @param cellarea area of cells
##' @param spatialvals spatial at risk, function lambda, interpolated onto the requisite grid
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
##' @return object passed back to lgcpPredictSpatial
##' @export
MALAlgcpSpatial <- function(mcmcloop,
                            inits,
                            adaptivescheme,
                            M,
                            N,
                            Mext,
                            Next,
                            sigma,
                            phi,
                            mu,
                            nis,
                            cellarea,
                            spatialvals,
                            scaleconst,
                            rootQeigs,
                            invrootQeigs,
                            cellInside,
                            MCMCdiag,
                            gradtrunc,
                            gridfun,
                            gridav,
                            mcens,
                            ncens){
                            
    SpatialOnlyMode <- TRUE
    ##ImprovedAlgorithm <- TRUE
    SpatialPlusParameters <- FALSE
    SpatioTemporalPlusParameters <- FALSE
    MultiTypeMode <- FALSE


    cellOutside <- !as.logical(cellInside)
    logspatial <- log(scaleconst*spatialvals)
    logspatial[cellOutside|spatialvals==0] <- 0 # NOTE THIS IS FOR SIMPLIFYING THE COMPUTATION OF THE TARGET!!
                                  
    temporal.fitted <- Inf # note this line is here for gridFunction and gridAverage methods and is not used otherwise
    nlevs <- NULL # note this line is here for gridFunction and gridAverage methods and is not used otherwise                            
    GFinitialise(gridfun) # note these two lines must come after M and N have been computed or defined
	GAinitialise(gridav) # note these two lines must come after M and N have been computed or defined
    
    h <- initialiseAMCMC(adaptivescheme)
    hrec <- h
    nsamp <- 0
    icount <- 0
    MCMCacc <- 0
    y.mean <- matrix(0,M,N)
    y.var <- matrix(0,M,N)
    EY.mean <- matrix(0,M,N)
    EY.var <- matrix(0,M,N)	                 
             
    Gamma <- matrix(0,Mext,Next) # initialise with mean     
    oldtags <- target.and.grad.spatial(Gamma=Gamma,nis=nis,cellarea=cellarea,rootQeigs=rootQeigs,invrootQeigs=invrootQeigs,mu=mu,spatial=spatialvals,logspat=logspatial,scaleconst=scaleconst,gradtrunc=gradtrunc)
    
    logtarget <- c()
    
    while(nextStep(mcmcloop)){  
    
        propmeans <- Gamma + (h/2)*oldtags$grad
        
        propGamma <- propmeans + sqrt(h)*rnorm(Mext*Next)
        proptags <- target.and.grad.spatial(Gamma=propGamma,nis=nis,cellarea=cellarea,rootQeigs=rootQeigs,invrootQeigs=invrootQeigs,mu=mu,spatial=spatialvals,logspat=logspatial,scaleconst=scaleconst,gradtrunc=gradtrunc)
        revpropmeans <- propGamma + (h/2)*proptags$grad
        
        ac <- exp(proptags$logtarget-oldtags$logtarget-sum((Gamma-revpropmeans)^2)/(2*h) + sum((propGamma-propmeans)^2)/(2*h))
        ac <- min(ac,1)
        
        icount <- icount + 1        
        MCMCacc <- ((icount-1)/icount) * MCMCacc + ac/icount
        
        ###browser()
        
        trigger <- FALSE
		if (proptags$logtarget==-Inf | is.na(ac) | is.nan(ac)){ # gradient truncation insufficient, so reduce
	        gradtrunc <- gradtrunc/2
	        cat("Reducing gradient truncation to:",gradtrunc,"\n")
	        oldtags <- target.and.grad.spatial(Gamma=Gamma,nis=nis,cellarea=cellarea,rootQeigs=rootQeigs,invrootQeigs=invrootQeigs,mu=mu,spatial=spatialvals,logspat=logspatial,scaleconst=scaleconst,gradtrunc=gradtrunc) 	        
	        if (!is.burnin(mcmcloop)){
	            cat("Gradient truncation currently",gradtrunc,"\n")
	            cat("Suggest reducing this further and setting gradtrunc manually in lgcpPredictSpatial.\n")
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
        
        if (is.retain(mcmcloop)){
	        nsamp <- nsamp + 1
        	y.mean <- ((nsamp-1)/nsamp) * y.mean + oldtags$Y[1:M,1:N]/nsamp
        	EY.mean <- ((nsamp-1)/nsamp) * EY.mean + oldtags$expY[1:M,1:N]/nsamp
        	if (nsamp>1){
    			y.var <- ((nsamp-2)/(nsamp-1))*y.var + (nsamp/(nsamp-1)^2)*(y.mean-oldtags$Y[1:M,1:N])^2
    			EY.var <- ((nsamp-2)/(nsamp-1))*EY.var + (nsamp/(nsamp-1)^2)*(EY.mean-oldtags$expY[1:M,1:N])^2
    		}                 	               
			GFupdate(gridfun)
			GAupdate(gridav)
			logtarget <- c(logtarget,oldtags$logtarget)		    
		}
    } 
	
	retlist <- list(lasth=h,lastGAM=oldtags$Gamma)

    GFfinalise(gridfun) # these two lines must appear after retlist has been initialised
	GAfinalise(gridav)  #	
	
	retlist$mcmcacc <- MCMCacc
	retlist$hrec <- hrec
    retlist$y.mean <- lgcpgrid(list(y.mean),xvals=mcens[1:M],yvals=ncens[1:N])
    retlist$y.var <- lgcpgrid(list(y.var),xvals=mcens[1:M],yvals=ncens[1:N])
    retlist$EY.mean <- lgcpgrid(list(EY.mean),xvals=mcens[1:M],yvals=ncens[1:N])
    retlist$EY.var <- lgcpgrid(list(EY.var),xvals=mcens[1:M],yvals=ncens[1:N])
    retlist$gridfunction <- GFreturnvalue(gridfun)
    retlist$gridaverage <- GAreturnvalue(gridav)
    retlist$mcmcinfo <- mcmcloop
    retlist$gradtrunc <- gradtrunc
    retlist$logtarget <- logtarget	
	
	return(retlist)                          
}
                        

