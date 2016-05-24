##' lgcpPredictSpatioTemporalPlusPars function
##'
##' A function to deliver fully Bayesian inference for the spatiotemporal log-Gaussian Cox process.\cr
##'
##' See the vignette "Bayesian_lgcp" for examples of this code in use.\cr
##'
##' The model for the data is as follows:\cr
##' \cr
##' X(s) ~ Poisson[R(s,t)]\cr
##' \cr
##' R(s) = C_A lambda(s,t) exp[Z(s,t)beta+Y(s,t)]\cr
##' \cr
##'
##' Here X(s,t) is the number of events in the cell of the computational grid containing s, R(s,t) is the Poisson rate, 
##' C_A is the cell area, lambda(s,t) is a known offset, Z(s,t) is a vector of measured covariates and Y(s,t) is the 
##' latent Gaussian process on the computational grid. The other parameters in the model are beta, the covariate effects; 
##' and eta=[log(sigma),log(phi),log(theta)], the parameters of the process Y on an appropriately transformed (in this case log) scale.\cr
##'
##' We recommend the user takes the following steps before running this method: 
##'
##' \enumerate{
##'    \item Compute approximate values of the parameters, eta, of the process Y using the function minimum.contrast. 
##'        These approximate values are used for two main reasons: (1) to help inform the size of the computational grid, since we 
##'        will need to use a cell width that enables us to capture the dependence properties of Y and (2) to help inform the 
##'        proposal kernel for the MCMC algorithm.
##'    \item Choose an appropriate grid on which to perform inference using the function chooseCellwidth; this will partly be determined 
##'        by the results of the first stage and partly by the available computational resource available to perform inference. 
##'    \item Using the function getpolyol, construct the computational grid and polygon overlays, as required. As this can be an expensive step, 
##'        we recommend that the user saves this object after it has been 
##'        constructed and in future reference to the data, reloads this object, rather than having to re-compute it (provided the 
##'        computational grid has not changed).
##'    \item Decide on which covariates are to play a part in the analysis and use the lgcp function getZmat to interpolate these 
##'        onto the computational grid. Note that having saved the results from the previous step, this is a relatively quick operation, 
##'        and allows the user to quickly construct different design matrices, Z, from different candidate models for the data
##'    \item If required, set up the population offset using SpatialAtRisk functions (see the vignette "Bayesian_lgcp"); specify the priors 
##'        using lgcpPrior; and if desired, the initial values for the MCMC, using the function lgcpInits.
##'    \item Run the MCMC algorithm and save the output to disk. We recommend dumping information to disk using the dump2dir function 
##'        in the output.control argument because it offers much greater flexibility in terms of MCMC diagnosis and post-processing.
##'    \item Perform post-processing analyses including MCMC diagnostic checks and produce summaries of the posterior expectations 
##'        we require for presentation. (see the vignette "Bayesian_lgcp" for further details). Functions of use in this step include
##'        traceplots, autocorr, parautocorr, ltar, parsummary, priorpost, postcov, textsummary, expectation, exceedProbs and lgcp:::expectation.lgcpPredict
##' }
##'
##' The user must provide a list of design matrices to use this function. In the interpolation step above, there are three cases to consider 
##'
##' \enumerate{
##'    \item where Z(s,t) cannot be decomposed, i.e., Z are true spatiotemporal covariates. In this case, each element of the list must 
##'        be constructed separately using the function getZmat on the covariates for each time point.
##'    \item Z(s,t)beta = Z_1(s)beta_1 + Z_2(t)beta_2: the spatial and temporal effects are separable; 
##'        in this case use the function addTemporalCovariates, to aid in the construction of the list.
##'    \item Z(s,t)beta = Z(s)beta, in which case the user only needs to perform the interpolation using getZmat 
##'        once, each of the elements of the  list will then be identical.
##'    \item Z(s,t)beta =  Z(t)beta in this case we follow the procedure for the separable case above.
##'        For example, if dotw is a temporal covariate we would use formula <- X ~ dotw for the main algorithm, formula.spatial <- X ~ 1 to 
##'        interpolate the spatial covariates using getZmat, followed by temporal.formula <- t ~ dotw - 1 using addTemporalCovariates
##'        to construct the list of design matrices, Zmat.
##' }
##'
##'
##' @param formula a formula object of the form X ~ var1 + var2 etc. The name of the dependent variable must be "X". Only accepts 'simple' formulae, such as the example given.
##' @param xyt An object of class stppp
##' @param T the time point of interest
##' @param laglength the number of previous time points to include in the analysis
##' @param ZmatList A list of design matrices Z constructed with getZmat and possibly addTemporalCovariates see the details below and Bayesian_lgcp vignette for details on how to construct this.
##' @param model.priors model priors, set using lgcpPrior
##' @param model.inits model initial values. The default is NULL, in which case lgcp will use the prior mean to initialise eta and beta will be initialised from an oversispersed glm fit to the data. Otherwise use lgcpInits to specify.
##' @param spatial.covmodel choice of spatial covariance function. See ?CovFunction
##' @param cellwidth the width of computational cells
##' @param poisson.offset A list of SpatialAtRisk objects (of length the number of types) defining lambda_k (see below)
##' @param mcmc.control MCMC paramters, see ?mcmcpars
##' @param output.control output choice, see ?setoutput   
##' @param gradtrunc truncation for gradient vector equal to H parameter Moller et al 1998 pp 473. Default is Inf, which means no gradient truncation, which seems to work in most settings.
##' @param ext integer multiple by which grid should be extended, default is 2. Generally this will not need to be altered, but if the spatial correlation decays slowly, increasing 'ext' may be necessary.
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former, the default, includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window.
##' @return an object of class lgcpPredictSpatioTemporalPlusParameters
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor, Tilman M. Davies, Barry S. Rowlingson, Peter J. Diggle. Bayesian Inference and Data Augmentation Schemes for Spatial, Spatiotemporal and Multivariate Log-Gaussian Cox Processes in R. Submitted. 
##'     \item Benjamin M. Taylor, Tilman M. Davies, Barry S. Rowlingson, Peter J. Diggle (2013). Journal of Statistical Software, 52(4), 1-40. URL http://www.jstatsoft.org/v52/i04/
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##'     \item Wood ATA, Chan G (1994). Simulation of Stationary Gaussian Processes in [0,1]d. Journal of Computational and Graphical Statistics, 3(4), 409-432.
##'     \item Moller J, Syversveen AR, Waagepetersen RP (1998). Log Gaussian Cox Processes. Scandinavian Journal of Statistics, 25(3), 451-482.
##' }
##' @seealso \link{minimum.contrast}, \link{minimum.contrast.spatiotemporal}, link{chooseCellWidth}, \link{getpolyol}, \link{guessinterp}, \link{getZmat},
##' \link{addTemporalCovariates}, \link{lgcpPrior}, \link{lgcpInits}, \link{CovFunction}
##' \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}, 
##' \link{lgcpPredictMultitypeSpatialPlusPars},
##' \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export

lgcpPredictSpatioTemporalPlusPars <- function( formula,
                                        xyt,
                                        T,
					                    laglength,
                                        ZmatList=NULL,       
                					    model.priors,
                					    model.inits=lgcpInits(),
                					    spatial.covmodel,
                					    cellwidth=NULL,
                					    poisson.offset=NULL,				
                					    mcmc.control,
                					    output.control=setoutput(),
                					    gradtrunc=Inf,
                					    ext=2,
                					    inclusion="touching"){
    
    regionalcovariates <- NULL
    pixelcovariates <- NULL     
    nprobe <- 100000    
    gridsize <- NULL
    starttime <- Sys.time()
    
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
	
	if(class(model.priors)!="lgcpPrior"){
	    stop("Argument model.priors must be of class lgcpPrior, see ?lgcpPrior")
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
	
	tdiff <- c(Inf,diff(aggtimes)) # require Inf for target and gradient function to work
	numt <- length(tdiff)
	
    ow <- selectObsWindow(xyt,cellwidth) 
	xyt <- ow$xyt
	M <- ow$M # note for this function, M and N are powers of 2 
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
    
    if(is.null(poisson.offset)){
        poisson.offset <- list()
        for (i in 1:numt){
            poisson.offset[[i]] <- list(X=seq(xyt$window$xrange[1],xyt$window$xrange[2],length.out=100),
                                        Y=seq(xyt$window$yrange[1],xyt$window$yrange[2],length.out=100),
                                        Zm=matrix(1,100,100))
        }
    }
    else{
        if(!is.list(poisson.offset)){ # assume the offset is only spatially varying
            po <- list()
            for (i in 1:numt){
                po[[i]] <- poisson.offset
            }
            poisson.offset <- po
            rm(po)
            gc()
        }
        else{
            if(length(poisson.offset)!=numt){
                stop(paste("Poisson offset should have length",numt))
            }
        }
    } 
    
    spatial <- list()
    for(i in 1:numt){
        if(!inherits(poisson.offset[[i]],"spatialAtRisk")){			
            spatial[[i]] <- spatialAtRisk(poisson.offset[[i]])
        }
        else{
            spatial[[i]] <- poisson.offset[[i]]
        }
        
        if(inherits(spatial[[i]],"fromXYZ")){
            spatial[[i]]$Zm <- spatial[[i]]$Zm*attr(spatial[[i]],"NC") # put back in 'normalising constant' so that spatialAtRisk acts as an offset (ie it no longer integrates to 1 over the study region.)
        }
        if(inherits(spatial[[i]],"fromSPDF")){
            spatial[[i]]$atrisk <- spatial[[i]]$atrisk*attr(spatial[[i]],"NC")
            spatial[[i]]$spdf$atrisk <- spatial[[i]]$atrisk
        }
    }   	
					
	
	
	## no longer required as here this is an offset: spatialvals <- spatialvals / (cellarea*sum(spatialvals))
    
    ################################################################
    # Create grid and FFT objects
    ################################################################
    
    study.region <- xyt$window
	
	## DEFINE LATTICE & CENTROIDS ##
	
	if(!is.null(attr(ZmatList[[1]],"gridobj"))){
	    gridobj <- attr(ZmatList[[1]],"gridobj")
	}
	else{
	    gridobj <- genFFTgrid(study.region=study.region,M=M,N=N,ext=ext,inclusion=inclusion)
	}
	del1 <- gridobj$del1
	del2 <- gridobj$del2
	Mext <- gridobj$Mext
	Next <- gridobj$Next
	mcens <- gridobj$mcens
	ncens <- gridobj$ncens
	cellarea <- gridobj$cellarea
	cellInside <- gridobj$cellInside
	
	## COMPUTE GRID DISTANCES ##
	
	x <- gridobj$mcens
    y <- gridobj$ncens    
    xidx <- rep(1:Mext,Next)
    yidx <- rep(1:Next,each=Mext)
    dxidx <- pmin(abs(xidx-xidx[1]),Mext-abs(xidx-xidx[1]))
    dyidx <- pmin(abs(yidx-yidx[1]),Next-abs(yidx-yidx[1]))
    d <- sqrt(((x[2]-x[1])*dxidx)^2+((y[2]-y[1])*dyidx)^2)
    
    ## OBTAIN POISSON OFFSET ON LATTICE (LINEAR INTERPOLATION) ##
	
	spatial.offset <- list()
	for(i in 1:numt){
    	spatial.offset[[i]] <- fftinterpolate(spatial[[i]],mcens,ncens,ext=ext)
    	spatial.offset[[i]] <- spatial.offset[[i]]*cellInside
	}
	
	## SET UP SPATIAL COVARIATES ON GRID
	
	spatialOnlyCovariates <- FALSE
	if(is.null(ZmatList)){
	    ZmatList <- list()
	    if(!inherits(regionalcovariates,"list")&!inherits(pixelcovariates,"list")){ # ie if covariates do not change over time	        
	        ZmatList <- cov.interp.fft(formula=formula,W=study.region,regionalcovariates=regionalcovariates,pixelcovariates=pixelcovariates,mcens=mcens[1:M],ncens=ncens[1:N],cellInside=cellInside[1:M,1:N])
	        spatialOnlyCovariates <- TRUE
	    }
	    else{
            if((!inherits(regionalcovariates,"list")|!inherits(pixelcovariates,"list"))){
                stop("regionalcovariates and pixelcovariates must EITHER both be list objects OR SpatialPolygonsDataFrame and SpatialPixelsDataFrame objects respectively.")
            }
            else{
                ZmatList[[i]] <- cov.interp.fft(formula=formula,W=study.region,regionalcovariates=regionalcovariates[[i]],pixelcovariates=pixelcovariates[[i]],mcens=mcens[1:M],ncens=ncens[1:N],cellInside=cellInside[1:M,1:N])
            }	    
	    }
	}
	else{
	    for(i in 1:numt){
	        if(inherits(ZmatList,"matrix")){
	            if(!isTRUE(all.equal(mcens[1:M],attr(ZmatList,"mcens")))|!isTRUE(all.equal(ncens[1:N],attr(ZmatList,"ncens")))){
                    stop(paste("FFT grid and ZmatList[[",i,"]] are on different grids. Please recompute ZmatList using 'getZmat'.",sep=""))
                }
                spatialOnlyCovariates <- TRUE
	        }
	        else{
                if(!isTRUE(all.equal(mcens[1:M],attr(ZmatList[[i]],"mcens")))|!isTRUE(all.equal(ncens[1:N],attr(ZmatList[[i]],"ncens")))){
                    stop(paste("FFT grid and ZmatList[[",i,"]] are on different grids. Please recompute ZmatList using 'getZmat'.",sep=""))
                }
            }
        }    	
	}
		
	
	
	
    	
	################################################################
		
	
		
    ###
    # Set up MCMC loop, required to compute nsamp, below
    ###
    mLoop = mcmcLoop(N=mcmc.control$mala.length,burnin=mcmc.control$burnin,thin=mcmc.control$retain,progressor=mcmcProgressTextBar)	
	
	# issue warning if dumping information to disc
	nsamp <- floor((mLoop$N-mLoop$burnin)/mLoop$thin)
	if (!is.null(output.control$gridfunction) & class(output.control$gridfunction)[1]=="dump2dir"){
    	cat("WARNING: disk space required for saving is approximately ",round(nsamp*object.size(array(runif((M)*(N)*(length(aggtimes))),dim=c((M),(N),(length(aggtimes)))))/1024^2,2)," Mb, ",sep="")
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
		nis[[i]] <- nis[[i]] * (spatial.offset[[i]]>0)
		ct2 <- sum(nis[[i]])
		if(ct2<ct1){
		    warning(paste("Time ",aggtimes[i],": ",ct1-ct2," data points lost due to discretisation.",sep=""),immediate.=TRUE)
		}
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
    
    lg <- MALAlgcpSpatioTemporal.PlusPars( mcmcloop=mLoop,
                                    inits=mcmc.control$inits,
                                    adaptivescheme=mcmc.control$adaptivescheme,
                                    M=M,
                                    N=N,
                                    Mext=Mext,
                                    Next=Next,
                                    mcens=mcens,
                                    ncens=ncens,
                                    formula=formula,
                                    ZmatList=ZmatList,
                                    model.priors=model.priors,
                                    model.inits=model.inits,
                                    fftgrid=gridobj,
                                    spatial.covmodel=spatial.covmodel,
                                    tdiff=tdiff,
                                    nis=nis,
                                    cellarea=cellarea,
                                    spatialvals=spatial.offset,
                                    cellInside=cellInside,
                                    MCMCdiag=mcmc.control$MCMCdiag,
                                    gradtrunc=gradtrunc,
                                    gridfun=gridfun,
                                    gridav=gridav,
                                    d=d,
                                    aggtimes=aggtimes,
                                    spatialOnlyCovariates=spatialOnlyCovariates)
                                                  
	
	endtime <- Sys.time()
	timetaken <- endtime-starttime
	
	lg$xyt <- xyt
	lg$M <- M
	lg$N <- N
	lg$aggtimes <- aggtimes
	lg$tdiffs <- tdiff
	lg$vars <- NA
	lg$spatial <- spatial
	lg$temporal <- NA
	lg$grid <- gridobj
	lg$nis <- nis
	lg$mcens <- mcens[1:M]
	lg$ncens <- ncens[1:N]
	lg$mcmcpars <- mcmc.control
	lg$timetaken <- timetaken
	lg$spatialonly <- TRUE
	lg$spatialonlyplusparameters <- TRUE
	lg$ext <- ext
	lg$cellInside <- cellInside[1:M,1:N]
	lg$inclusion <- inclusion
    lg$poisson.offset <- spatial.offset
	lg$priors <- model.priors
	lg$covFct <- spatial.covmodel
	
	class(lg) <- c("lgcpPredictSpatioTemporalPlusParameters","lgcpPredict","lgcpobject")	
	
	return(lg)															
					
}	
	

##' MALAlgcpSpatioTemporal.PlusPars function
##'
##' A function to run the MCMC algorithm for spatiotemporal point process data. Not for general purpose use.
##'
##' @param mcmcloop details of the mcmc loop
##' @param inits initial values
##' @param adaptivescheme the adaptive MCMC scheme 
##' @param M number of grid cells in x direction
##' @param N number of grid cells in y direction
##' @param Mext number of extended grid cells in x direction 
##' @param Next number of extended grid cells in y direction
##' @param mcens centroids in x direction
##' @param ncens  centroids in y direction
##' @param formula a formula object of the form X ~ var1 + var2 etc. 
##' @param ZmatList list of design matrices constructed using getZmat
##' @param model.priors model priors, constructed using lgcpPrior 
##' @param model.inits initial values for the MCMC
##' @param fftgrid an objects of class FFTgrid, see genFFTgrid
##' @param spatial.covmodel spatial covariance model, consructed with CovFunction
##' @param nis cell counts on the etended grid
##' @param tdiff vector of time differences
##' @param cellarea the cell area
##' @param spatialvals inerpolated poisson offset on fft grid 
##' @param cellInside 0-1 matrix indicating inclusion in the observation window 
##' @param MCMCdiag not used
##' @param gradtrunc gradient truncation parameter
##' @param gridfun used to specify other actions to be taken, e.g. dumping MCMC output to disk.
##' @param gridav used for computing Monte Carlo expectations online
##' @param d matrix of toral distances
##' @param aggtimes the aggregate times
##' @param spatialOnlyCovariates whether this is a 'spatial' only problem
##' @return output from the MCMC run
##' @export

MALAlgcpSpatioTemporal.PlusPars <- function(   mcmcloop,
                                        inits,
                                        adaptivescheme,
                                        M,
                                        N,
                                        Mext,
                                        Next,
                                        mcens,
                                        ncens,
                                        formula,
                                        ZmatList,
                                        model.priors,
                                        model.inits,
                                        fftgrid,
                                        spatial.covmodel,
                                        nis,
                                        tdiff,
                                        cellarea,
                                        spatialvals,
                                        cellInside,
                                        MCMCdiag,
                                        gradtrunc,
                                        gridfun,
                                        gridav,
                                        d,
                                        aggtimes,
                                        spatialOnlyCovariates){
                            
    SpatialOnlyMode <- FALSE
    SpatialPlusParameters <- FALSE
    SpatioTemporalPlusParameters <- TRUE
    MultiTypeMode <- FALSE
    
    numt <- length(tdiff)
    
    cellInsideLogical <- as.logical(cellInside)
      
    M <- M
    N <- N
    temporal.fitted <- rep(Inf,length(aggtimes)) # note this line is here for gridFunction and gridAverage methods and is not used otherwise
    nlevs <- NULL # note this line is here for gridFunction and gridAverage methods and is not used otherwise
    GFinitialise(gridfun) # note these two lines must come after M and N have been computed or defined
	GAinitialise(gridav) # note these two lines must come after M and N have been computed or defined
    
    nsamp <- 0
    icount <- 0
    MCMCacc <- 0
    y.mean <- as.list(rep(0,numt)) #matrix(0,M,N)
    y.var <- as.list(rep(0,numt)) #matrix(0,M,N)
    EY.mean <- as.list(rep(0,numt)) #matrix(0,M,N)
    EY.var <- as.list(rep(0,numt)) #matrix(0,M,N)
    
    ###
	# Initialise
	###	
	
	etainvtrans <- model.priors$etaprior$inverse_transform
	
	if(is.null(model.inits$etainit)){
	    etaval <- model.priors$etaprior$mean
	}
	else{
	    etaval <- model.inits$etainit
	}
	
	etainv <- etainvtrans(etaval)
	cp <- CovParameters(list(sigma=etainv[1],phi=etainv[2],theta=etainv[3]))
	ETA0 <- cp
	
	off <- c() # required so the the call to glm passes CRAN check
	rm(off)
    if(!is.null(model.inits$betainit)){
        betaval <- model.inits$betainit	    
    }
    else{
    	betaval <- 0
    	if(spatialOnlyCovariates){
    	    for(i in 1:numt){    	    
    	        dfr <- attr(ZmatList,"data.frame")
    	        dfr$X <- nis[[i]][cellInsideLogical]
            	dfr$off <- log(cellarea*spatialvals[[i]][cellInsideLogical])
            	dfr$off[is.infinite(dfr$off)] <- NA # this excludes cells where the rate is zero
            	mod <- glm(formula,data=dfr,family=quasipoisson,offset=off)
            	betaval <- betaval+coefficients(mod)/numt # gives mean beta over each time point 
    	    }
    	} 
    	else{
    	    dfr <- c()
    	    for(i in 1:numt){
    	        dfrtemp <- attr(ZmatList[[i]],"data.frame") 
    	        dfrtemp$X <- nis[[i]][cellInsideLogical]
            	dfrtemp$off <- log(cellarea*spatialvals[[i]][cellInsideLogical])
            	dfr$off[is.infinite(dfr$off)] <- NA # this excludes cells where the rate is zero
            	dfr <- rbind(dfr,dfrtemp)
        	}
        	mod <- glm(formula,data=dfr,family=quasipoisson,offset=off)
        	betaval <- coefficients(mod) # gives mean beta over each time point 
	    }
	     
	    if(any(is.na(betaval))){
            stop("Initial value of beta, as computed using glm contains NA values.")
        } 	
    }
    betaval <- BetaParameters(betaval)
    

    tm <- matrix(FALSE,Mext,Next)
    tm[1:M,1:N] <- TRUE
    if(spatialOnlyCovariates){
        Z <- matrix(0,Next*Mext,ncol=ncol(ZmatList))
        if(!all((1:length(betaval))==match(names(betaval),colnames(ZmatList)))){
            warning("Re-ordering variables in Z matrix to match order in spatiotemporal formula.",immediate.=TRUE)
        }
        ZmatList <- ZmatList[,match(names(betaval),colnames(ZmatList))]        
        Z[as.vector(tm),] <- ZmatList
        Zt <- t(Z)
    } 
    else{
        Z <- list()
        Zt <- list()
        if(!all((1:length(betaval))==match(names(betaval),colnames(ZmatList[[1]])))){
            warning("Re-ordering variables in Z matrix to match order in spatiotemporal formula.",immediate.=TRUE)
        }
        for(i in 1:numt){
            Z[[i]] <- matrix(0,Next*Mext,ncol=ncol(ZmatList[[i]]))
            ZmatList[[i]] <- ZmatList[[i]][,match(names(betaval),colnames(ZmatList[[i]]))]                        
            Z[[i]][as.vector(tm),] <- ZmatList[[i]] 
            Zt[[i]] <- t(Z[[i]])
        }
    } 
    
    
    glmfitted <- list()
    for(i in 1:numt){
        if(inherits(Z,"matrix")){
            Zbeta <- matrix(as.vector(Z%*%betaval),Mext,Next)        
            glmfitted[[i]] <- cellarea*spatialvals[[i]]*exp(Zbeta)
        }
        else{ 
            Zbeta <- matrix(as.vector(Z[[i]]%*%betaval),Mext,Next)
            glmfitted[[i]] <- cellarea*spatialvals[[i]]*exp(Zbeta)   
        }
        glmfitted[[i]] <- glmfitted[[i]][1:M,1:N]
     }
    #browser()
    ###
	# Now loop
	###
	
    neta <- length(etaval)
	nbeta <- length(betaval)	
	
	# compute approximate variance matrix of Y
	y <- list()
	gammatemp <- list()
    GPdummy <- GPrealisation(gamma=matrix(rnorm(Mext*Next),Mext,Next),fftgrid=fftgrid,covFunction=spatial.covmodel,covParameters=cp,d=d)
	if(spatialOnlyCovariates){	
    	for(i in 1:numt){
        	y[[i]] <- log(nis[[i]]/(cellarea*spatialvals[[i]]* exp(matrix(as.vector(Z%*%betaval),Mext,Next))))
            y[[i]][is.na(y[[i]]) | is.infinite(y[[i]])] <- -cp$sigma^2/2
            gammainit <- GammafromY(Y=y[[i]],rootQeigs=GPdummy$rootQeigs,mu=GPdummy$CovParameters$mu) 
            gammatemp[[i]] <- gammainit  
        }
    }
    else{
        #stop()
        for(i in 1:numt){
        	y[[i]] <- log(nis[[i]]/(cellarea*spatialvals[[i]]* exp(matrix(as.vector(Z[[i]]%*%betaval),Mext,Next))))
            y[[i]][is.na(y[[i]]) | is.infinite(y[[i]])] <- -cp$sigma^2/2
            gammainit <- GammafromY(Y=y[[i]],rootQeigs=GPdummy$rootQeigs,mu=GPdummy$CovParameters$mu) 
            gammatemp[[i]] <- gammainit  
        }
    }   
    #  
    
    GP <- stGPrealisation(gammatemp,fftgrid=fftgrid,covFunction=spatial.covmodel,covParameters=cp,d=d,tdiff=tdiff)   
    
    gammaVar <- list()
    rootgammaVar <- list()
    CB <- list()
    CB <- matrix(0,nbeta,nbeta)
    
    for(i in 1:numt){                   
        if(spatialOnlyCovariates){
            eMat <- cellarea*spatialvals[[i]]*exp(matrix(as.vector(Z%*%betaval),nrow(GP$gamma[[i]]),ncol(GP$gamma[[i]])))*GP$expY[[i]]
        }
        else{
            eMat <- cellarea*spatialvals[[i]]*exp(matrix(as.vector(Z[[i]]%*%betaval),nrow(GP$gamma[[i]]),ncol(GP$gamma[[i]])))*GP$expY[[i]]
        }
        fi <- (1/length(GP$Y[[i]]))*Re(fft(GP$invrootQeigs*fft((1/length(GP$Y[[i]]))*Re(fft(GP$invrootQeigs*fft(eMat,inverse=TRUE))),inverse=TRUE)))
        g_i <- 1-exp(-2*tdiff[i]*cp$theta)
        gammaVar[[i]] <- (1.65^2/((Mext*Next)^(1/3)))*(-1)/(-(1/g_i)-fi)
        rootgammaVar[[i]] <- sqrt(gammaVar[[i]])
        # compute variance matrix of beta
        if(spatialOnlyCovariates){
            EZ <- as.vector(cellarea*spatialvals[[i]]*exp(matrix(as.vector(Z%*%betaval),nrow(GP$gamma[[i]]),ncol(GP$gamma[[i]]))+GP$Y[[i]]))
        }
        else{
            EZ <- as.vector(cellarea*spatialvals[[i]]*exp(matrix(as.vector(Z[[i]]%*%betaval),nrow(GP$gamma[[i]]),ncol(GP$gamma[[i]]))+GP$Y[[i]]))
        }
        
        for(j in 1:nbeta){
            for(k in 1:nbeta){
                if(spatialOnlyCovariates){
                    CB[j,k] <- CB[j,k] + sum(Z[,j]*Z[,k]*EZ)  # observed information
                }
                else{
                    CB[j,k] <- CB[j,k] + sum(Z[[i]][,j]*Z[[i]][,k]*EZ)  # observed information
                }
            }
        } 
    }
    CB <- CB + model.priors$betaprior$precision


    additional_scaleconst <- (0.234/0.574)     

    cat("Computing posterior variance w.r.to eta ...\n")

    lensq <- 10
    sqsigma <- seq(model.priors$etaprior$mean[1]-2*sqrt(model.priors$etaprior$variance[1,1]),model.priors$etaprior$mean[1]+2*sqrt(model.priors$etaprior$variance[1,1]),length.out=lensq)
    sqphi <- seq(model.priors$etaprior$mean[2]-2*sqrt(model.priors$etaprior$variance[2,2]),model.priors$etaprior$mean[2]+2*sqrt(model.priors$etaprior$variance[2,2]),length.out=lensq)
    sqtheta <- seq(model.priors$etaprior$mean[3]-2*sqrt(model.priors$etaprior$variance[3,3]),model.priors$etaprior$mean[3]+2*sqrt(model.priors$etaprior$variance[3,3]),length.out=lensq)
    #sqsigma <- sort(rnorm(lensq,model.priors$etaprior$mean[1],sqrt(model.priors$etaprior$variance[1,1])))
    #sqphi <- sort(rnorm(lensq,model.priors$etaprior$mean[2],sqrt(model.priors$etaprior$variance[2,2])))
    #sqtheta <- sort(rnorm(lensq,model.priors$etaprior$mean[3],sqrt(model.priors$etaprior$variance[3,3])))
    
    ltarmat <- array(NA,dim=c(lensq,lensq,lensq))
    pb <- txtProgressBar(min=0,max=lensq,style=3)
    for (i in 1:lensq){
        setTxtProgressBar(pb,i)
        for(j in 1:lensq){
            for(k in 1:lensq){                        
                cpA <- CovParameters(list(sigma=exp(sqsigma[i]),phi=exp(sqphi[j]),theta=exp(sqtheta[k])))      
                gpA <- stGPrealisation(gammatemp,fftgrid=fftgrid,covFunction=spatial.covmodel,covParameters=cpA,d=d,tdiff=tdiff)  
                matent <- try(target.and.grad.SpatioTemporalPlusPars(GP=gpA,prior=model.priors,Z=Z,Zt=Zt,eta=c(sqsigma[i],sqphi[j],sqtheta[k]),beta=betaval,nis=nis,cellarea=cellarea,spatial=spatialvals,gradtrunc=Inf,ETA0=ETA0,tdiff=tdiff)$logtarget)
                if(class(matent)!="try-error"){
                    ltarmat[i,j,k] <- matent 
                }
            }
        }        
    }
    close(pb)
    ltarmat[is.infinite(ltarmat)] <- NA
    dffit <- data.frame(ltar=as.vector(ltarmat))
    exgr <- expand.grid(sqsigma,sqphi,sqtheta)
    dffit$sigma <- exgr[,1]
    dffit$sigma2 <- exgr[,1]^2
    dffit$phi <- exgr[,2]
    dffit$phi2 <- exgr[,2]^2
    dffit$theta <- exgr[,3]
    dffit$theta2 <- exgr[,3]^2
    dffit$sigmaphi <- exgr[,1]*exgr[,2]
    dffit$sigmatheta <- exgr[,1]*exgr[,3]
    dffit$phitheta <- exgr[,2]*exgr[,3]
    
    
    try(tarmod <- lm(ltar~sigma2+sigma+phi2+phi+theta2+theta+sigmaphi+sigmatheta+phitheta,data=dffit))
    try(coef <- coefficients(tarmod))
    mmm <- matrix(c(2*coef[2],coef[8],coef[9],  
                    coef[8],2*coef[4],coef[10],
                    coef[9],coef[10],2*coef[6]),3,3)
    mmm[3,1:2] <- 0 
    mmm[1:2,3] <- 0                        
    try(etaCovMat <- additional_scaleconst*solve((-1)*mmm))    
    etaCovMattest2 <- FALSE    
    if(class(etaCovMat)!="try-error"){
        etaCovMattest2 <- !all(eigen(etaCovMat)$values>0)
    }
    if(class(etaCovMat)=="try-error" | etaCovMattest2){
        warning("Cannot find good approximation of posterior variance w.r.to eta: using the following variance instead:",immediate.=TRUE)
        etaCovMat <- (1/100)*model.priors$etaprior$variance
    }
    else{
       cat("Approximate of posterior variance w.r.to eta found:\n")
    }         

	sigma_eta <- (2.38^2/length(etaval))*etaCovMat 

	SIGMA_ETA <- sigma_eta
	Q_eta <- solve(sigma_eta)
	sigma_beta <- (1.65^2/(length(betaval)^(1/3)))*solve(CB) 
	Q_beta <- solve(sigma_beta)
	
	sigma_eta_chol <- t(chol(sigma_eta))
	sigma_beta_chol <- t(chol(sigma_beta))
	
    print(sigma_eta)
    print(sigma_beta)	
	
	betarec <- c()
	etarec <- c()
	
	sigmaetarec <- list(sigma_eta)
	sigmabetarec <- list(sigma_beta)
    
    adapt_h <- adaptivescheme
    h <- initialiseAMCMC(adapt_h)
    
    adapt_heta <- andrieuthomsh(inith=1,alpha=0.5,C=1,targetacceptance=0.574)
    heta <- initialiseAMCMC(adapt_heta)


	betaetarec <- c()
     
	
    GP <- stGPrealisation(gamma=lapply(1:numt,function(x){matrix(0,Mext,Next)}),fftgrid=fftgrid,covFunction=spatial.covmodel,covParameters=cp,d=d,tdiff=tdiff)
    oldtags <- target.and.grad.SpatioTemporalPlusPars(GP=GP,prior=model.priors,Z=Z,Zt=Zt,eta=etaval,beta=betaval,nis=nis,cellarea=fftgrid$cellarea,spatial=spatialvals,gradtrunc=gradtrunc,ETA0=ETA0,tdiff=tdiff)
    
    tarrec <- oldtags$logtarget 
    hallrec <- h  
    
    reject_its <- c()   
       
    while(nextStep(mcmcloop)){
    
        propmeans_beta <- betaval + (h/2)*sigma_beta%*%oldtags$gradbeta     
        propbeta <- BetaParameters(as.vector(propmeans_beta+sqrt(h)*sigma_beta_chol%*%rnorm(nbeta))) 
    
        propmeans_gamma <- list()
        for (i in 1:numt){
            propmeans_gamma[[i]] <- GP$gamma[[i]] + (h/2)*gammaVar[[i]]*oldtags$gradgamma[[i]]
        }
        propmeans_eta <- etaval 
        propeta <- as.vector(propmeans_eta + sqrt(h)*sigma_eta_chol%*%rnorm(neta))
        

        
        propetainv <- etainvtrans(propeta) 
        propcp <- CovParameters(list(sigma=propetainv[1],phi=propetainv[2],theta=propetainv[3]))
        propgamma <- list()
        for(i in 1:numt){
            propgamma[[i]] <- propmeans_gamma[[i]] + sqrt(h)*rootgammaVar[[i]]*rnorm(Mext*Next) 
        }
        propGP <- stGPrealisation(gamma=propgamma,fftgrid=fftgrid,covFunction=spatial.covmodel,covParameters=propcp,d=d,tdiff=tdiff)
                
                
          
        proptags <- target.and.grad.SpatioTemporalPlusPars(GP=propGP,prior=model.priors,Z=Z,Zt=Zt,eta=propeta,beta=propbeta,nis=nis,cellarea=fftgrid$cellarea,spatial=spatialvals,gradtrunc=gradtrunc,ETA0=ETA0,tdiff=tdiff)    
        
        revpropmeans_gamma <- list()
        for(i in 1:numt){
            revpropmeans_gamma[[i]] <- propGP$gamma[[i]] + (h/2)*gammaVar[[i]]*proptags$gradgamma[[i]]
        }
        revpropmeans_eta <- propeta 
        revpropmeans_beta <- propbeta + (h/2)*sigma_beta%*%proptags$gradbeta

        gammacontrib <- 0
        for(i in 1:numt){
            gammacontrib <- gammacontrib - sum((GP$gamma[[i]]-revpropmeans_gamma[[i]])^2/(2*h*gammaVar[[i]])) + sum((propGP$gamma[[i]]-propmeans_gamma[[i]])^2/(2*h*gammaVar[[i]]))
        }
   
        
        ac <- exp(proptags$logtarget-oldtags$logtarget + gammacontrib +
                            (-(0.5/h)*t(betaval-revpropmeans_beta)%*%Q_beta%*%(betaval-revpropmeans_beta)) -
                            (-(0.5/h)*t(propbeta-propmeans_beta)%*%Q_beta%*%(propbeta-propmeans_beta)))
                            
       
        ac <- min(ac,1)

        icount <- icount + 1        
        MCMCacc <- ((icount-1)/icount) * MCMCacc + ac/icount
        

		if (proptags$logtarget==-Inf | is.na(ac) | is.nan(ac)){ # gradient truncation insufficient, so reduce

	        warning("One possible cause of this warning is that there may be evidence in the data for quite large values of the spatial correlation parameter. If this is the case, then this warning means that the MCMC chain has wandered into a region of the phi-space that causes the variance matrix of Y (computed by the DFT) to become non positive-definite. One possible way of rectifying this issue is to restart the chain using a larger value of 'ext' in the call to lgcpPredictSpatioTemporalPlusPars. You should do this if the warning message is repeated many times. If this warning message appears at all then you should be warned that inference from this run may not be reliable: the proposed move at this iteration was rejected.",immediate.=TRUE)
	        ac <- 0
	        reject_its <- c(reject_its,iteration(mcmcloop))
	    } 
	    
        
        if (ac>runif(1)){
            GP <- propGP
            oldtags <- proptags
            etaval <- propeta
            betaval <- propbeta
        }
        
        h <- updateAMCMC(adapt_h)
        if(iteration(mcmcloop)%%100==0){print(h)}
        

       
         
        
        if (is.retain(mcmcloop)){

            hallrec <- c(hallrec,h)
            tarrec <- c(tarrec,oldtags$logtarget)
            betarec <- rbind(betarec,betaval)
	        etarec <- rbind(etarec,etaval)
      
	        nsamp <- nsamp + 1
	        
	        for(i in 1:numt){
            	y.mean[[i]] <- ((nsamp-1)/nsamp) * y.mean[[i]] + GP$Y[[i]][1:(M),1:(N)]/nsamp
            	EY.mean[[i]] <- ((nsamp-1)/nsamp) * EY.mean[[i]] + GP$expY[[i]][1:(M),1:(N)]/nsamp
            	if (nsamp>1){
        			y.var[[i]] <- ((nsamp-2)/(nsamp-1))*y.var[[i]] + (nsamp/(nsamp-1)^2)*(y.mean[[i]]-GP$Y[[i]][1:(M),1:(N)])^2
        			EY.var[[i]] <- ((nsamp-2)/(nsamp-1))*EY.var[[i]] + (nsamp/(nsamp-1)^2)*(EY.mean[[i]]-GP$expY[[i]][1:(M),1:(N)])^2
        		}
    		}                 	               
			GFupdate(gridfun)
			GAupdate(gridav)		    
		}
    } 
	
	retlist <- list(lasth=rev(hallrec)[1],lastGAM=oldtags$Gamma)
	retlist$mcmcacc <- MCMCacc
    retlist$y.mean <- lgcpgrid(y.mean)
    retlist$y.var <- lgcpgrid(y.var)
    retlist$EY.mean <- lgcpgrid(EY.mean)
    retlist$EY.var <- lgcpgrid(EY.var)
    retlist$gridfunction <- GFreturnvalue(gridfun)
    retlist$gridaverage <- GAreturnvalue(gridav)
    retlist$mcmcinfo <- mcmcloop
    retlist$gradtrunc <- gradtrunc
    retlist$etarec <- etarec
    retlist$betarec <- betarec  
    retlist$glmfit <- mod
    #retlist$dfr <- dfr
    retlist$glmfittedvals <- glmfitted
    retlist$sigmaetarec <- sigmaetarec
    retlist$sigmabetarec <- sigmabetarec
    retlist$hallrec <- hallrec
    retlist$tarrec <- tarrec
    retlist$Z <- Z
    retlist$reject_its <- reject_its
	
	return(retlist)                          
}                        

##' target.and.grad.SpatioTemporalPlusPars function
##'
##' A function to compute the target and gradient for the Bayesian spatiotemporal LGCP.
##'
##' @param GP an object created using the stGPrealisation function
##' @param prior the priors for hte model, created using lgcpPrior
##' @param Z the design matrix on the FFT grid
##' @param Zt the transpose of the design matrix
##' @param eta the paramers eta
##' @param beta the parameters beta
##' @param nis the cell counts on the FFT grid
##' @param cellarea the cell area
##' @param spatial the poisson offset
##' @param gradtrunc the gradient truncation parameter
##' @param ETA0 the initial value of eta
##' @param tdiff vector of time differences between time points
##' @return the target and gradient for the spatiotemporal model.
##' @export

target.and.grad.SpatioTemporalPlusPars <- function(GP,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc,ETA0,tdiff){

    if(inherits(Z,"matrix")){
        Zbeta <- matrix(as.vector(Z%*%beta),nrow(GP$gamma[[1]]),ncol(GP$gamma[[1]]))
        expZbeta <- exp(Zbeta)
    }
    else{ 
        Zbeta <- lapply(Z,function(x){matrix(as.vector(x%*%beta),nrow(GP$gamma[[1]]),ncol(GP$gamma[[1]]))})
        expZbeta <- lapply(Zbeta,exp)   
    }
    
    pri <- EvaluatePrior(etaParameters=eta,betaParameters=beta,prior=prior)    
    
    ncells <- length(GP$gamma[[1]])
    tcons <- 1/(1-exp(-2*ETA0$theta*tdiff))
    n <- length(tcons)
    
    e <- list()
    NminusE <- list()
    if(inherits(Z,"matrix")){
        lapply(1:n,function(i){e[[i]] <<- spatial[[i]]*expZbeta*GP$expY[[i]]*cellarea; NminusE[[i]] <<- nis[[i]]-e[[i]]})
    }
    else{
        lapply(1:n,function(i){e[[i]] <<- spatial[[i]]*expZbeta[[i]]*GP$expY[[i]]*cellarea; NminusE[[i]] <<- nis[[i]]-e[[i]]} )   
    }    
    

    gradmult <- list()
    gradcomp <- list()
    grad <- list()

    expYtrunc <- GP$expY
    if(!is.infinite(gradtrunc)){
        expYtrunc <- lapply(GP$expY,function(x){x[x>gradtrunc]<-gradtrunc; return(x)})
    }    
    
    if(is.infinite(gradtrunc)){
        gradfun <- function(i){
            grad[[i]] <<- (-1)*tcons[i]*GP$gamma[[i]]
            gradmult[[i]] <<- as.list(sqrt(GP$gt[i])*c(rep(0,i-1),c(1,cumprod(GP$bt[-c(1,1:i)])))) # computes the product of sqrt(gt(i)) * beta(delta t) terms in the summation
            gradcomp[[i]] <<- NminusE[[i]]
        }
    }
    else{    
        if(inherits(Z,"matrix")){
            gradfun <- function(i){
                grad[[i]] <<- (-1)*tcons[i]*GP$gamma[[i]]
                gradmult[[i]] <<- as.list(sqrt(GP$gt[i])*c(rep(0,i-1),c(1,cumprod(GP$bt[-c(1,1:i)])))) # computes the product of sqrt(gt(i)) * beta(delta t) terms in the summation
                gradcomp[[i]] <<- nis[[i]]-spatial[[i]]*expZbeta*expYtrunc[[i]]*cellarea
            }
        }
        else if(inherits(Z,"list")){
            gradfun <- function(i){
                grad[[i]] <<- (-1)*tcons[i]*GP$gamma[[i]]
                gradmult[[i]] <<- as.list(sqrt(GP$gt[i])*c(rep(0,i-1),c(1,cumprod(GP$bt[-c(1,1:i)])))) # computes the product of sqrt(gt(i)) * beta(delta t) terms in the summation
                gradcomp[[i]] <<- nis[[i]]-spatial[[i]]*expZbeta[[i]]*expYtrunc[[i]]*cellarea
            }
        }
        else{
            stop("Error in target.and.grad.SpatioTemporalPlusPars, print out the function to see where this occurred")
        }
    }
    
    sapply(1:n,gradfun)


    gradbeta <- pri$betacontrib$gradcontrib
    for(i in 1:n){
        gradcum <- 0
        for(j in i:n){
            gradcum <- gradcum + gradmult[[i]][[j]] * gradcomp[[j]]
        }        
        grad[[i]] <- grad[[i]] + (1/ncells)*Re(fft(GP$invrootQeigs*fft(gradcum,inverse=TRUE)))        
        if(inherits(Z,"matrix")){
            gradbeta <- gradbeta + as.vector(Zt%*%as.vector(NminusE[[i]]))
        }
        else{
            gradbeta <- gradbeta + as.vector(Zt[[i]]%*%as.vector(NminusE[[i]]))
        }
    }
      
    gradeta <- NULL
    
    if(inherits(Z,"matrix")){
        logtarget <- -(1/2)*sum(tcons*sapply(GP$gamma,function(x){sum(x*x)})) + sum(sapply(1:n,function(i){sum((Zbeta+GP$Y[[i]])*nis[[i]] - e[[i]])}))
    }
    else{
        logtarget <- -(1/2)*sum(tcons*sapply(GP$gamma,function(x){sum(x*x)})) + sum(sapply(1:n,function(i){sum((Zbeta[[i]]+GP$Y[[i]])*nis[[i]] - e[[i]])})) 
    }
    logtarget <- logtarget + pri$etacontrib$loglik + pri$betacontrib$loglik 
    
    return(list(logtarget=logtarget,gradgamma=grad,gradbeta=gradbeta,gradeta=gradeta,e=e))
}




