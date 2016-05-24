##' lgcpPredictMultitypeSpatialPlusPars function
##'
##' A function to deliver fully Bayesian inference for a multitype spatial log-Gaussian Cox process.\cr
##'
##' See the vignette "Bayesian_lgcp" for examples of this code in use.\cr
##'
##' We suppose there are K point types of interest. The model for point-type k is as follows:\cr
##' \cr
##' X_k(s) ~ Poisson[R_k(s)]\cr
##' \cr
##' R_k(s) = C_A lambda_k(s) exp[Z_k(s)beta_k+Y_k(s)]\cr
##' \cr
##' 
##' Here X_k(s) is the number of events of type k in the computational grid cell containing the
##' point s, R_k(s) is the Poisson rate, C_A is the cell area, lambda_k(s) is a known offset, Z_k(s) is a vector
##' of measured covariates and Y_i(s) where i = 1,...,K+1 are latent Gaussian processes on the
##' computational grid. The other parameters in the model are beta_k , the covariate effects for the
##' kth type; and eta_i = [log(sigma_i),log(phi_i)], the parameters of the process Y_i for i = 1,...,K+1 on
##' an appropriately transformed (again, in this case log) scale.
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
##' @param formulaList an object of class formulaList, see ?formulaList. A list of formulae of the form t1 ~ var1 + var2 etc. The name of the dependent variable must correspond to the name of the point type. Only accepts 'simple' formulae, such as the example given.
##' @param sd a marked ppp object, the mark of interest must be able to be coerced to a factor variable
##' @param typemark if there are multiple marks, thrun the MCMC algorithm for spatial point process data. Not for general purpose use.is sets the name of the mark by which 
##' @param Zmat design matrix including all covariate effects from each point type, constructed with getZmat
##' @param model.priorsList model priors, a list object of length the number of types, each element set using lgcpPrior
##' @param model.initsList list of model initial values (of length the number of types). The default is NULL, in which case lgcp will use the prior mean to initialise eta and beta will be initialised from an oversispersed glm fit to the data. Otherwise use lgcpInits to specify.
##' @param spatial.covmodelList list of spatial covariance functions (of length the number of types). See ?CovFunction
##' @param cellwidth the width of computational cells
##' @param poisson.offset A list of SpatialAtRisk objects (of length the number of types) defining lambda_k (see below)
##' @param mcmc.control MCMC paramters, see ?mcmcpars
##' @param output.control output choice, see ?setoutput   
##' @param gradtrunc truncation for gradient vector equal to H parameter Moller et al 1998 pp 473. Default is Inf, which means no gradient truncation, which seems to work in most settings.
##' @param ext integer multiple by which grid should be extended, default is 2. Generally this will not need to be altered, but if the spatial correlation decays slowly, increasing 'ext' may be necessary.
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former, the default, includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window.
##' @return an object of class lgcpPredictMultitypeSpatialPlusParameters
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
##' \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}, \link{lgcpPredictSpatioTemporalPlusPars},
##' \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export
  
lgcpPredictMultitypeSpatialPlusPars <- function(formulaList,
                                                sd,
                                                typemark=NULL,
                                                Zmat=NULL,       
                        					    model.priorsList, # gives normal prior with mean betahat and variance |betahat|*betaerror where betahat is estimated using glm in the first instance
                        					    model.initsList=NULL,
                        					    spatial.covmodelList,
                        					    cellwidth=NULL,
                        					    poisson.offset=NULL,				
                        					    mcmc.control,
                        					    output.control=setoutput(),
                        					    gradtrunc=Inf,
                        					    ext=2,
                        					    inclusion="touching"){
    
    spatial.offsetList <- poisson.offset
    regionalcovariates <- NULL
    pixelcovariates <- NULL  
    gridsize <- NULL                    					    
    formula <- aggregateformulaList(formulaList)                         					    
    
    starttime <- Sys.time()
    
    if(is.null(typemark)){
        if(!is.null(dim(sd$marks))){
            stop("Since sd$marks is not a vector, you must specify which of the marks defines the point type by setting the argument typemark.")
        }
    }
    else{
        if(!any(names(sd$marks)==typemark)){
            stop(paste("None of the marks is named",typemark))
        }
    }
    
    if(is.null(dim(sd$marks))){
        if(!is.factor(sd$marks)){
            warning(paste("Converting",typemark,"to a factor"),.immediate=TRUE)            
        }
        marks <- as.factor(sd$marks)
    }
    else{
        if(!is.factor(sd$marks$typemark)){
            warning(paste("Converting",typemark,"to a factor"),.immediate=TRUE)            
        }
        marks <- as.factor(sd$marks$typemark)        
    }
    ntypes <- length(formulaList)
    
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
	
	
	
	
    ow <- selectObsWindow(sd,cellwidth) 
	sd <- ow$xyt
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
    
    spatial <- list()	
	
	if(is.null(spatial.offsetList)){
	    for(i in 1:ntypes){
            spatial[[i]] <- spatialAtRisk(list(X=seq(sd$window$xrange[1],sd$window$xrange[2],length.out=M),Y=seq(sd$window$yrange[1],sd$window$yrange[2],length.out=N),Zm=matrix(1,M,N)))
        }
    }
    else if (class(spatial.offsetList)=="SpatialAtRisk"){
        for(i in 1:ntypes){
            spatial[[i]] <- spatial.offsetList
        }
    }
    else if (class(spatial.offsetList)=="list"){
        for(i in 1:ntypes){
            if(!any(class(spatial.offsetList[[i]])=="spatialAtRisk")){		
                spatial[[i]] <- spatialAtRisk(spatial.offsetList[[i]])
            }
            else{
                spatial[[i]] <- spatial.offsetList[[i]]
            }
        }
    }
    else{
        stop("Invalid spatial.offsetList")    
    }
    
    
    funk <- function(spatial){
        if(any(class(spatial)=="fromXYZ")){
            spatial$Zm <- spatial$Zm*attr(spatial,"NC") # put back in 'normalising constant' so that spatialAtRisk acts as an offset (ie it no longer integrates to 1 over the study region.)
        }
        if(any(class(spatial)=="fromSPDF")){
            spatial$atrisk <- spatial$atrisk*attr(spatial,"NC")
            spatial$spdf$atrisk <- spatial$atrisk
        }
        return(spatial)
    }
    
    spatial <- lapply(spatial,funk)
    
    ################################################################
    # Create grid and FFT objects
    ################################################################
    
    study.region <- sd$window
	
	## DEFINE LATTICE & CENTROIDS ##
	
	if(!is.null(attr(Zmat,"gridobj"))){
	    gridobj <- attr(Zmat,"gridobj")
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
	
	## SET UP SPATIAL COVARIATES ON GRID

	
	# NOT NEEDED ANYMORE Zmat <- cov.interp.fft(formula=formula,W=study.region,regionalcovariates=regionalcovariates,pixelcovariates=pixelcovariates,mcens=mcens[1:M],ncens=ncens[1:N],nprobe=nprobe,cellInside=cellInside[1:M,1:N])
    if(is.null(Zmat)){
	    Zmat <- cov.interp.fft(formula=formula,W=study.region,regionalcovariates=regionalcovariates,pixelcovariates=pixelcovariates,mcens=mcens[1:M],ncens=ncens[1:N],cellInside=cellInside[1:M,1:N])
	}
	else{
        if(!isTRUE(all.equal(mcens[1:M],attr(Zmat,"mcens")))|!isTRUE(all.equal(ncens[1:N],attr(Zmat,"ncens")))){
            stop("FFT grid and Zmat are on different grids. Please recompute Zmat using 'getZmat'.")
        }    	
	}	
	zml <- getZmats(Zmat=Zmat,formulaList=formulaList)	
	
	## OBTAIN POISSON OFFSET ON LATTICE (LINEAR INTERPOLATION) ##
	
	spatialvals <- list()	
	for(i in 1:ntypes){
    	spatialvals[[i]] <- fftinterpolate(spatial[[i]],mcens,ncens,ext=ext)
    	spatialvals[[i]] <- spatialvals[[i]]*cellInside ## note that this is not scaled, as we wish for lambda to act as an offset / (cellarea*sum(spatialvals))
	}
    	
	################################################################
		
    ###
    # Set up MCMC loop, required to compute nsamp, below
    ###
    mLoop = mcmcLoop(N=mcmc.control$mala.length,burnin=mcmc.control$burnin,thin=mcmc.control$retain,progressor=mcmcProgressTextBar)	
	
	# issue warning if dumping information to disc
	nsamp <- floor((mLoop$N-mLoop$burnin)/mLoop$thin)
	if (!is.null(output.control$gridfunction) & class(output.control$gridfunction)[1]=="dump2dir"){
    	cat("WARNING: disk space required for saving is approximately ",round(nsamp*object.size(array(runif((M)*(N)*(ntypes+1)),dim=c((M),(N),ntypes+1)))/1024^2,2)," Mb, ",sep="")
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
	
	mlev <- sapply(formulaList,function(x){variablesinformula(x)[1]})
	nis <- list()
	ct1 <- 0
	ct2 <- 0
	for (i in 1:ntypes){
    	nis[[i]]  <- getCounts(xyt=sd[marks==mlev[i],],M=M,N=N,ext=ext) 
        ct1 <- ct1+ sum(nis[[i]])
    	nis[[i]] <- nis[[i]] * cellInside
    	ct2 <- ct2 + sum(nis[[i]])
    }	
    if(ct2<ct1){
	    warning(paste(ct1-ct2," data points lost due to discretisation.",sep=""),immediate.=TRUE)
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
    
    lg <- MALAlgcpMultitypeSpatial.PlusPars( mcmcloop=mLoop,
                                    inits=mcmc.control$inits,
                                    adaptivescheme=mcmc.control$adaptivescheme,
                                    M=M,
                                    N=N,
                                    Mext=Mext,
                                    Next=Next,
                                    mcens=mcens,
                                    ncens=ncens,
                                    formulaList=formulaList,
                                    zml=zml,
                                    Zmat=Zmat,
                                    model.priorsList=model.priorsList,
                                    model.initsList=model.initsList,
                                    fftgrid=gridobj,
                                    spatial.covmodelList=spatial.covmodelList,
                                    nis=nis,
                                    cellarea=cellarea,
                                    spatialvals=spatialvals,
                                    cellInside=cellInside,
                                    MCMCdiag=mcmc.control$MCMCdiag,
                                    gradtrunc=gradtrunc,
                                    gridfun=gridfun,
                                    gridav=gridav,
                                    marks=marks,
                                    ntypes=ntypes,
                                    d=d)
                                                  
	
	endtime <- Sys.time()
	timetaken <- endtime-starttime
	
	lg$xyt <- sd
	lg$M <- M
	lg$N <- N
	lg$aggtimes <- NA
	lg$tdiffs <- NA
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
	lg$poisson.offset <- spatialvals
	lg$priors <- model.priorsList
	lg$covFct <- spatial.covmodelList
	
	class(lg) <- c("lgcpPredictMultitypeSpatialPlusParameters","lgcpPredict","lgcpobject")	
	
	return(lg)															
					
}	
	

##' MALAlgcpMultitypeSpatial.PlusPars function
##'
##' A function to run the MCMC algorithm for multivariate spatial point process data. Not for general purpose use.
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
##' @param formulaList a list of formula objects of the form X ~ var1 + var2 etc. 
##' @param zml list of design matrices
##' @param Zmat a design matrix constructed using getZmat
##' @param model.priorsList list of model priors, see lgcpPriors
##' @param model.initsList list of model initial values, see lgcpInits
##' @param fftgrid an objects of class FFTgrid, see genFFTgrid
##' @param spatial.covmodelList list of spatial covariance models constructed using CovFunction
##' @param nis cell counts on the etended grid
##' @param cellarea the cell area
##' @param spatialvals inerpolated poisson offset on fft grid 
##' @param cellInside 0-1 matrix indicating inclusion in the observation window 
##' @param MCMCdiag not used
##' @param gradtrunc gradient truncation parameter
##' @param gridfun used to specify other actions to be taken, e.g. dumping MCMC output to disk.
##' @param gridav used for computing Monte Carlo expectations online
##' @param marks the marks from the marked ppp object
##' @param ntypes the number of types being analysed
##' @param d matrix of toral distances
##' @return output from the MCMC run
##' @export

MALAlgcpMultitypeSpatial.PlusPars <- function(  mcmcloop,
                                                inits,
                                                adaptivescheme,
                                                M,
                                                N,
                                                Mext,
                                                Next,
                                                mcens,
                                                ncens,
                                                formulaList,
                                                zml,
                                                Zmat,
                                                model.priorsList,
                                                model.initsList,
                                                fftgrid,
                                                spatial.covmodelList,
                                                nis,
                                                cellarea,
                                                spatialvals,
                                                cellInside,
                                                MCMCdiag,
                                                gradtrunc,
                                                gridfun,
                                                gridav,
                                                marks,
                                                ntypes,
                                                d){
                            
    SpatialOnlyMode <- TRUE
    ##ImprovedAlgorithm <- TRUE
    SpatialPlusParameters <- TRUE
    SpatioTemporalPlusParameters <- FALSE
    MultiTypeMode <- TRUE
    
    cellInsideLogical <- as.logical(cellInside)      
    
    cidx <- ntypes + 1
    nlevs <- length(formulaList)  
    M <- M
    N <- N
    temporal.fitted <- Inf # note this line is here for gridFunction and gridAverage methods and is not used otherwise                            
    GFinitialise(gridfun) # note these two lines must come after M and N have been computed or defined
	GAinitialise(gridav) # note these two lines must come after M and N have been computed or defined
    
    nsamp <- 0
    icount <- 0
    MCMCacc <- 0
    y.mean <- array(0,dim=c(M,N,ntypes+1))
    y.var <- array(0,dim=c(M,N,ntypes+1))
    EY.mean <- array(0,dim=c(M,N,ntypes+1))
    EY.var <- array(0,dim=c(M,N,ntypes+1))
    
    ###
	# Initialise
	###	
	
	etainvtrans <- list()
	etaval <- list()
    etainv <- list()
    cp <- list()
	for(i in 1:(ntypes+1)){
    	etainvtrans[[i]] <- model.priorsList[[i]]$etaprior$inverse_transform
    	
    	if(is.null(model.initsList)){
    	    etaval[[i]] <- model.priorsList[[i]]$etaprior$mean
    	}
    	else{
    	    etaval[[i]] <- model.initsList[[i]]$etainit
    	}
    	etainv[[i]] <- etainvtrans[[i]](etaval[[i]])
    	cp[[i]] <- CovParameters(list(sigma=etainv[[i]][1],phi=etainv[[i]][2]))
	}
	
	mlev <- sapply(formulaList,function(x){variablesinformula(x)[1]})
	SigmaBeta <- list()
	betavals <- list()

	dfr <- attr(Zmat,"data.frame")
	modls <- list()
	off <- c() # required so the the call to glm passes CRAN check
	rm(off)
	for(i in 1:ntypes){    	
    	dfr[[mlev[i]]] <- nis[[i]][cellInsideLogical] 
    	dfr$off <- log(cellarea*spatialvals[[i]][cellInsideLogical])
    	dfr$off[is.infinite(dfr$off)] <- NA # this excludes cells where the rate is zero
    	mod <- glm(formulaList[[i]],data=dfr,family=quasipoisson,offset=off)
    	modls[[i]] <- mod    	
        betavals[[i]] <- BetaParameters(coefficients(mod))
        SigmaBeta[[i]] <- vcov(mod)
    }
   
    Zlist <- list()
    Ztlist <- list()
    tm <- matrix(FALSE,Mext,Next)
    tm[1:M,1:N] <- TRUE
    for(i in 1:ntypes){
        Zlist[[i]] <- matrix(0,Next*Mext,ncol=ncol(zml[[i]]))        
        Zlist[[i]][as.vector(tm),] <- zml[[i]]
        Ztlist[[i]] <- t(Zlist[[i]])
    }  
    
    neta <- list()
    nbeta <- list()
    for(i in 1:(ntypes+1)){
        neta[[i]] <- length(etaval[[i]])
    }
    for(i in 1:ntypes){
	    nbeta[[i]] <- length(betavals[[i]])        
    }
    
    gammaVar <- list()
    CB <- list()
    etaCovMat <- list()  
      
    
    GPls <- list()
    gmm <- 0
    for (k in 1:ntypes){
        # compute approximate variance matrix of Y        
    	y <- log(nis[[k]]/(cellarea*spatialvals[[k]]* exp(matrix(as.vector(Zlist[[k]]%*%betavals[[k]]),Mext,Next))))
        y[is.na(y) | is.infinite(y)] <- -cp[[k]]$sigma^2/2
        
        GPdummy <- GPrealisation(gamma=matrix(0,Mext,Next),fftgrid=fftgrid,covFunction=spatial.covmodelList[[k]],covParameters=cp[[k]],d=d)
        gammainit <- GammafromY(Y=y,rootQeigs=GPdummy$rootQeigs,mu=cp[[k]]$mu)        
        GP <- GPrealisation(gamma=gammainit,fftgrid=fftgrid,covFunction=spatial.covmodelList[[k]],covParameters=cp[[k]],d=d)
        gmm <- gmm + gammainit
        GPls[[k]] <- GP
        eMat <- cellarea*spatialvals[[k]]*exp(matrix(as.vector(Zlist[[k]]%*%betavals[[k]]),nrow(GP$gamma),ncol(GP$gamma)))*GP$expY
        fi <- (1/length(GP$Y))*Re(fft(GP$invrootQeigs*fft((1/length(GP$Y))*Re(fft(GP$invrootQeigs*fft(eMat,inverse=TRUE))),inverse=TRUE)))
        gammaVar[[k]] <- 2 * (1.65^2/((Mext*Next)^(1/3)))*(-1)/(-1-fi)

        # compute approximate variance matrix of beta
        EZ <- as.vector(cellarea*spatialvals[[k]]*exp(matrix(as.vector(Zlist[[k]]%*%betavals[[k]]),nrow(GP$gamma),ncol(GP$gamma))+GP$Y))    
        CB[[k]] <- matrix(NA,nbeta[[k]],nbeta[[k]]) 
           
        for(i in 1:nbeta[[k]]){
            for(j in 1:nbeta[[k]]){
                CB[[k]][i,j] <- sum(Zlist[[k]][,i]*Zlist[[k]][,j]*EZ) + model.priorsList[[k]]$betaprior$precision[i,j]# observed information
            }
        }
        
    }
    GPls[[ntypes+1]] <- GPrealisation(gamma=matrix(0,Mext,Next),fftgrid=fftgrid,covFunction=spatial.covmodelList[[ntypes+1]],covParameters=cp[[ntypes+1]],d=d)
    #for (k in 1:ntypes){    # compute approximate variance matrix of eta              
    #    etaCovMat[[k]] <- (0.234/0.574)*solve((-1)*GPdrv2_Multitype(GPlist=GPls,priorlist=model.priorsList,Zlist=Zlist,Ztlist=Ztlist,etalist=etaval,betalist=betavals,nis=nis,cellarea=cellarea,spatial=spatialvals,gradtrunc=Inf,fftgrid=fftgrid,covfunction=spatial.covmodelList,d=d,eps=1e-6,k=k)$hess)       
    #}
    
    if(TRUE){
        tmpls <- list()
        additional_scaleconst <- (0.234/0.574) 
        for(k in 1:(ntypes+1)){
            lensq <- 10
            sqsigma <- seq(model.priorsList[[k]]$etaprior$mean[1]-2*sqrt(model.priorsList[[k]]$etaprior$variance[1,1]),model.priorsList[[k]]$etaprior$mean[1]+2*sqrt(model.priorsList[[k]]$etaprior$variance[1,1]),length.out=lensq)
            sqphi <- seq(model.priorsList[[k]]$etaprior$mean[2]-2*sqrt(model.priorsList[[k]]$etaprior$variance[2,2]),model.priorsList[[k]]$etaprior$mean[2]+2*sqrt(model.priorsList[[k]]$etaprior$variance[2,2]),length.out=lensq)
            ltarmat <- matrix(NA,lensq,lensq)
            for (i in 1:lensq){
                for(j in 1:lensq){                        
                    cpA <- CovParameters(list(sigma=exp(sqsigma[i]),phi=exp(sqphi[j])))                    
                    gpA <- GPls
                    etatmp <- etaval
                    etatmp[[k]] <- c(sqsigma[i],sqphi[j])      
                    gpA[[k]] <- GPrealisation(gamma=GPls[[k]]$gamma,fftgrid=fftgrid,covFunction=spatial.covmodelList[[k]],covParameters=cpA,d=d)
                    #matent <- try(target.and.grad.spatialPlusPars(GP=gpA,prior=model.priors,Z=Z,Zt=Zt,eta=c(sqsigma[i],sqphi[j]),beta=betaval,nis=nis,cellarea=cellarea,spatial=spatialvals,gradtrunc=Inf)$logtarget)
                    matent <- try(target.and.grad.MultitypespatialPlusPars(GPlist=gpA,priorlist=model.priorsList,Zlist=Zlist,Ztlist=Ztlist,eta=etatmp,beta=betavals,nis=nis,cellarea=fftgrid$cellarea,spatial=spatialvals,gradtrunc=gradtrunc)$logtarget)
                    if(class(matent)!="try-error"){
                        ltarmat[i,j] <- matent 
                    }
                }        
            }
            ltarmat[is.infinite(ltarmat)] <- NA
            dffit <- data.frame(ltar=as.vector(ltarmat))
            exgr <- expand.grid(sqsigma,sqphi)
            dffit$sigma <- exgr[,1]
            dffit$sigma2 <- exgr[,1]^2
            dffit$phi <- exgr[,2]
            dffit$phi2 <- exgr[,2]^2
            dffit$sigmaphi <- exgr[,1]*exgr[,2]
            try(tarmod <- lm(ltar~sigma2+sigma+phi2+phi+sigmaphi,data=dffit))
            try(coef <- coefficients(tarmod))
            try(etaCovMat[[k]] <- additional_scaleconst*solve((-1)*matrix(c(2*coef[2],coef[6],coef[6],2*coef[4]),2,2)))
        }
    }
    
    gammaVar[[cidx]] <- 0
    #etaCovMat[[cidx]] <- 0
    for(i in 1:ntypes){
        gammaVar[[cidx]] <- gammaVar[[cidx]] + (1/ntypes)*gammaVar[[i]]
    #    etaCovMat[[cidx]] <- etaCovMat[[cidx]] + (1/ntypes)*etaCovMat[[i]]
    }
    gammaVar[[cidx]] <- gammaVar[[cidx]]/ntypes
    #etaCovMat[[cidx]] <- etaCovMat[[cidx]]/ntypes

    rootgammaVar <- lapply(gammaVar,sqrt)


    ###
	# Now loop
	###
	
	sigma_eta <- list()
	SIGMA_ETA <- list()
	Q_eta <- list()
	sigma_beta <- list()
	Q_beta <- list()
	sigma_eta_chol <- list()
	sigma_beta_chol <- list()

	
	for(i in 1:(ntypes+1)){
    	sigma_eta[[i]] <- (2.38^2/length(etaval[[i]]))*etaCovMat[[i]]
    	SIGMA_ETA[[i]] <- sigma_eta[[i]] 
    	Q_eta[[i]] <- solve(sigma_eta[[i]])    	
    	sigma_eta_chol[[i]] <- t(chol(sigma_eta[[i]])) 
    	print(sigma_eta[[i]])   	
	}
	for(i in 1:ntypes){
    	sigma_beta[[i]] <- (1.65^2/(length(betavals[[i]])^(1/3)))*solve(CB[[i]])
    	Q_beta[[i]] <- as.matrix(solve(sigma_beta[[i]]))
	    sigma_beta_chol[[i]] <- as.matrix(t(chol(sigma_beta[[i]])))
	    print(sigma_beta[[i]])
	}

	
	
	
	betarec <- c()
	etarec <- c()

    adapt_h <- adaptivescheme
    h <- initialiseAMCMC(adapt_h)  

    gammacount <- etacount <- betacount <- 1
    
    GPlist <- list()
    for(i in 1:(ntypes+1)){
        GPlist[[i]] <- GPrealisation(gamma=matrix(0,Mext,Next),fftgrid=fftgrid,covFunction=spatial.covmodelList[[i]],covParameters=cp[[i]],d=d)
    }
    oldtags <- target.and.grad.MultitypespatialPlusPars(GPlist=GPlist,priorlist=model.priorsList,Zlist=Zlist,Ztlist=Ztlist,eta=etaval,beta=betavals,nis=nis,cellarea=fftgrid$cellarea,spatial=spatialvals,gradtrunc=gradtrunc)
	
    trigger1 <- TRUE

    tarrec <- oldtags$logtarget 
    hallrec <- h   
    
    reject_its <- c()   
    
    while(nextStep(mcmcloop)){

        propmeans_gamma <- list()
        propmeans_eta <- list()
        propeta <- list()
        propetainv <- list()
        propcp <- list()
        propGPlist <- list()
        propmeans_beta <- list()
        propbeta <- list()
        for (i in (ntypes+1):1){ #note this loop is backwards (for a reason!)
            propmeans_gamma[[i]] <- GPlist[[i]]$gamma + (h/2)*gammaVar[[i]]*oldtags$gradgamma[[i]]
            propmeans_eta[[i]] <- etaval[[i]] 
            propeta[[i]] <- as.vector(propmeans_eta[[i]] + sqrt(h)*sigma_eta_chol[[i]]%*%rnorm(neta[[i]])) 
            propetainv[[i]] <- etainvtrans[[i]](propeta[[i]])
            propcp[[i]] <- CovParameters(list(sigma=propetainv[[i]][1],phi=propetainv[[i]][2]))
            propGPlist[[i]] <- GPrealisation(gamma=propmeans_gamma[[i]] + sqrt(h)*rootgammaVar[[i]]*rnorm(Mext*Next),fftgrid=fftgrid,covFunction=spatial.covmodelList[[i]],covParameters=propcp[[i]],d=d)        
        }
        for (i in 1:ntypes){            
            propmeans_beta[[i]] <- betavals[[i]] + (h/2)*sigma_beta[[i]]%*%oldtags$gradbeta[[i]]                    
            propbeta[[i]] <- BetaParameters(as.vector(propmeans_beta[[i]]+sqrt(h)*sigma_beta_chol[[i]]%*%rnorm(nbeta[[i]])))         
        }
        
          
        proptags <- target.and.grad.MultitypespatialPlusPars(GPlist=propGPlist,priorlist=model.priorsList,Zlist=Zlist,Ztlist=Ztlist,eta=propeta,beta=propbeta,nis=nis,cellarea=fftgrid$cellarea,spatial=spatialvals,gradtrunc=gradtrunc)    
        
        revpropmeans_gamma <- list()
        revpropmeans_eta <- list()
        revpropmeans_beta <- list()
        gamma_acceptance_contrib <- 0
        eta_acceptance_contrib <- 0
        beta_acceptance_contrib <- 0
        for (i in (ntypes+1):1){ # again, doen in reverse order
            revpropmeans_gamma[[i]] <- propGPlist[[i]]$gamma + (h/2)*gammaVar[[i]]*proptags$gradgamma[[i]]
            revpropmeans_eta[[i]] <- propeta[[i]] 
            
            gamma_acceptance_contrib <- gamma_acceptance_contrib -sum((GPlist[[i]]$gamma-revpropmeans_gamma[[i]])^2/(2*h*gammaVar[[i]])) + sum((propGPlist[[i]]$gamma-propmeans_gamma[[i]])^2/(2*h*gammaVar[[i]]))
            
        }
        for (i in 1:ntypes){
            revpropmeans_beta[[i]] <- propbeta[[i]] + (h/2)*sigma_beta[[i]]%*%proptags$gradbeta[[i]]           
            beta_acceptance_contrib <- beta_acceptance_contrib + (-(0.5/h)*t(betavals[[i]]-revpropmeans_beta[[i]])%*%Q_beta[[i]]%*%(betavals[[i]]-revpropmeans_beta[[i]])) - (-(0.5/h)*t(propbeta[[i]]-propmeans_beta[[i]])%*%Q_beta[[i]]%*%(propbeta[[i]]-propmeans_beta[[i]]))
        }        
        
        ac <- exp(proptags$logtarget-oldtags$logtarget + gamma_acceptance_contrib + eta_acceptance_contrib + beta_acceptance_contrib)

        ac <- min(ac,1)

        icount <- icount + 1        
        MCMCacc <- ((icount-1)/icount) * MCMCacc + ac/icount
        

		if (proptags$logtarget==-Inf | is.na(ac) | is.nan(ac)){ # gradient truncation insufficient, so reduce

	        warning("One possible cause of this warning is that there may be evidence in the data for quite large values of the spatial correlation parameter. If this is the case, then this warning means that the MCMC chain has wandered into a region of the phi-space that causes the variance matrix of Y (computed by the DFT) to become non positive-definite. One possible way of rectifying this issue is to restart the chain using a larger value of 'ext' in the call to lgcpPredictMultitypeSpatialPlusPars. You should do this if the warning message is repeated many times. If this warning message appears at all then you should be warned that inference from this run may not be reliable: the proposed move at this iteration was rejected.",immediate.=TRUE)
	        ac <- 0
	        reject_its <- c(reject_its,iteration(mcmcloop))
	    } 
        
        if (ac>runif(1)){
            GPlist <- propGPlist
            oldtags <- proptags
            etaval <- propeta
            betavals <- propbeta
        }
        
       
        if(iteration(mcmcloop)%%100==0){print(h)}
        h <- updateAMCMC(adapt_h)
             
        
        if (is.retain(mcmcloop)){

            hallrec <- c(hallrec,h)
            tarrec <- c(tarrec,oldtags$logtarget)
            betarec <- rbind(betarec,unlist(betavals))
	        etarec <- rbind(etarec,unlist(etaval))
	        
	        Y <- GPlist2array(GPlist=GPlist,element="Y")
	        expY <- GPlist2array(GPlist=GPlist,element="expY")
      
	        nsamp <- nsamp + 1
        	y.mean <- ((nsamp-1)/nsamp) * y.mean + Y[1:M,1:N,]/nsamp
        	EY.mean <- ((nsamp-1)/nsamp) * EY.mean + expY[1:M,1:N,]/nsamp
        	if (nsamp>1){
    			y.var <- ((nsamp-2)/(nsamp-1))*y.var + (nsamp/(nsamp-1)^2)*(y.mean-Y[1:M,1:N,])^2
    			EY.var <- ((nsamp-2)/(nsamp-1))*EY.var + (nsamp/(nsamp-1)^2)*(EY.mean-expY[1:M,1:N,])^2
    		}                 	               
			GFupdate(gridfun)
			GAupdate(gridav)		    
		}
    } 
	
	retlist <- list(lasth=rev(hallrec)[1],lastGAM=oldtags$Gamma)
	retlist$mcmcacc <- MCMCacc
    retlist$y.mean <- lgcpgrid(y.mean,xvals=mcens[1:M],yvals=ncens[1:N])
    retlist$y.var <- lgcpgrid(y.var,xvals=mcens[1:M],yvals=ncens[1:N])
    retlist$EY.mean <- lgcpgrid(EY.mean,xvals=mcens[1:M],yvals=ncens[1:N])
    retlist$EY.var <- lgcpgrid(EY.var,xvals=mcens[1:M],yvals=ncens[1:N])
    retlist$gridfunction <- GFreturnvalue(gridfun)
    retlist$gridaverage <- GAreturnvalue(gridav)
    retlist$mcmcinfo <- mcmcloop
    retlist$gradtrunc <- gradtrunc
    retlist$etarec <- etarec
    retlist$betarec <- betarec  
    retlist$glmfit <- modls
    retlist$Zlist <- Zlist
    retlist$hallrec <- hallrec
    retlist$tarrec <- tarrec
    retlist$cinslogical <- cellInsideLogical
    retlist$reject_its <- reject_its
	
	return(retlist)                          
}                        

##' target.and.grad.MultitypespatialPlusPars function
##'
##' A function to compute the taget an gradient for the Bayesian multivariate lgcp
##'
##' @param GPlist list of Gaussian processes
##' @param priorlist list of priors
##' @param Zlist list of design matrices on the FFT gridd
##' @param Ztlist list of transposed design matrices
##' @param eta LGCP model parameter eta
##' @param beta LGCP model parameter beta
##' @param nis matrix of cell counts on the extended grid
##' @param cellarea the cell area
##' @param spatial the poisson offset interpolated onto the correcy grid
##' @param gradtrunc gradient truncation paramter
##' @return the target and gradient
##' @export

target.and.grad.MultitypespatialPlusPars <- function(GPlist,priorlist,Zlist,Ztlist,eta,beta,nis,cellarea,spatial,gradtrunc){

    cidx <- length(GPlist)
    ntypes <- cidx-1
    Zbeta <- list()
    expZbeta <- list()
    pri <- list()
    pri[[cidx]] <- EvaluatePrior(etaParameters=eta[[cidx]],betaParameters=NULL,prior=priorlist[[cidx]]) 
    e <- list()
    NminusE <- list()
    logtarget <- -(1/2)*sum(GPlist[[cidx]]$gamma^2) + pri[[cidx]]$etacontrib$loglik
    gradgamma <- list()
    gradbeta <- list()
    gradeta <- list()
    gradgamma[[cidx]] <- (-1)*GPlist[[cidx]]$gamma
    for(i in 1:ntypes){
        Zbeta[[i]] <- matrix(as.vector(Zlist[[i]]%*%beta[[i]]),nrow(GPlist[[i]]$gamma),ncol(GPlist[[i]]$gamma))
        expZbeta[[i]] <- exp(Zbeta[[i]])
        pri[[i]] <- EvaluatePrior(etaParameters=eta[[i]],betaParameters=beta[[i]],prior=priorlist[[i]])    
        e[[i]] <- spatial[[i]]*expZbeta[[i]]*GPlist[[cidx]]$expY*GPlist[[i]]$expY*cellarea    
        NminusE[[i]] <- nis[[i]] - e[[i]]

        logtarget <- logtarget -(1/2)*sum(GPlist[[i]]$gamma^2) + sum((Zbeta[[i]]+GPlist[[cidx]]$Y+GPlist[[i]]$Y)*nis[[i]] - e[[i]]) + pri[[i]]$etacontrib$loglik + pri[[i]]$betacontrib$loglik # note that both nis=0, logspat=0 and spatial=0 outside of observation window, so this effectively limits summation to the observation window only
    
        if(!is.infinite(gradtrunc)){
            expYtrunc <- GPlist[[cidx]]$expY
            expStrunc <- GPlist[[i]]$expY
            expYtrunc[expYtrunc>gradtrunc] <- gradtrunc
            expStrunc[expStrunc>gradtrunc] <- gradtrunc  
            gradgamma[[i]] <- (-1)*GPlist[[i]]$gamma + (1/length(GPlist[[i]]$Y))*Re(fft(GPlist[[i]]$invrootQeigs*fft(nis[[i]]-spatial[[i]]*expZbeta[[i]]*expYtrunc*expStrunc*cellarea,inverse=TRUE)))
            gradgamma[[cidx]] <- gradgamma[[cidx]] + (1/length(GPlist[[cidx]]$Y))*Re(fft(GPlist[[cidx]]$invrootQeigs*fft(nis[[i]]-spatial[[i]]*expZbeta[[i]]*expYtrunc*expStrunc*cellarea,inverse=TRUE)))
        }
        else{
            gradgamma[[i]] <- (-1)*GPlist[[i]]$gamma + (1/length(GPlist[[i]]$Y))*Re(fft(GPlist[[i]]$invrootQeigs*fft(NminusE[[i]],inverse=TRUE)))
            gradgamma[[cidx]] <- gradgamma[[cidx]] + (1/length(GPlist[[cidx]]$Y))*Re(fft(GPlist[[cidx]]$invrootQeigs*fft(NminusE[[i]],inverse=TRUE)))
        }
        
        gradbeta[[i]] <- as.vector(Ztlist[[i]]%*%as.vector(NminusE[[i]])) + as.vector(pri[[i]]$betacontrib$gradcontrib)
     }
              
    logtarget <- as.vector(logtarget)

    return(list(logtarget=logtarget,gradgamma=gradgamma,gradbeta=gradbeta,gradeta=gradeta))
}



##' lgcpSimMultitypeSpatialCovariates function
##'
##' A function to Simulate multivariate point process models
##'
##' @param formulaList a list of formulae objetcs 
##' @param owin a spatstat owin object on which to simulate the data    
##' @param regionalcovariates a SpatialPolygonsDataFrame object
##' @param pixelcovariates a SpatialPixelsDataFrame object
##' @param betaList list of beta parameters
##' @param spatial.offsetList list of poisson offsets
##' @param cellwidth cellwidth
##' @param model.parameters model parameters, a list eg list(sigma=1,phi=0.2)
##' @param spatial.covmodel the choice of spatial covariance model, can be anything from the RandomFields covariance function, CovariacenFct.
##' @param covpars additional covariance parameters, for the chosen model, optional.
##' @param ext number of times to extend the simulation window
##' @param plot whether to plot the results automatically
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former, the default, includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window.
##' @return a marked ppp object, the simulated data
##' @export

lgcpSimMultitypeSpatialCovariates <- function(  formulaList,
                                                owin,
                                                regionalcovariates,
                                                pixelcovariates,
                                                betaList,
                                                spatial.offsetList=NULL,
                                                cellwidth,
                                                model.parameters,
                                                spatial.covmodel="exponential",
                                                covpars=c(),
                                                ext=2,
                                                plot=FALSE,
                                                inclusion="touching"){
                                                
    ntypes <- length(formulaList)     
    vartypes <- getLHSformulaList(formulaList)
    formula <- aggregateformulaList(formulaList)                                           
                                                                                      
    sigma <- model.parameters$sigma
	phi <- model.parameters$phi
	mu <- model.parameters$mu

    if(is.null(owin)){
        owin <- owin()
    }    
    
    const0 <- 0.05 # level below which correlation function must drop to be considered not important
    
    # check space discretisation
    c2 <- -phi*log(const0) # the number at which exp(-theta*x) = 0.05, use this to choose time discretisation 
    if (cellwidth>c2/2){
        warning(paste("cellwidth should be at least",c2/2,"to get accurate results."))
    }
     
    xyt <- ppp(window=owin)
      
    ow <- selectObsWindow(xyt,cellwidth)
	xyt <- ow$xyt
	M <- ow$M
	N <- ow$N
	cat(paste("FFT Grid size: [",ext*M," , ",ext*N,"]\n",sep=""))

    spatial <- list()	
	
	if(is.null(spatial.offsetList)){
	    for(i in 1:ntypes){
            spatial[[i]] <- spatialAtRisk(list(X=seq(xyt$window$xrange[1],xyt$window$xrange[2],length.out=M),Y=seq(xyt$window$yrange[1],xyt$window$yrange[2],length.out=N),Zm=matrix(1,M,N)))
        }
    }
    else if (class(spatial.offsetList)=="SpatialAtRisk"){
        for(i in 1:ntypes){
            spatial[[i]] <- spatial.offsetList
        }
    }
    else if (class(spatial.offsetList)=="list"){
        for(i in 1:ntypes){
            if(!any(class(spatial.offsetList[[i]])=="spatialAtRisk")){		
                spatial[[i]] <- spatialAtRisk(spatial.offsetList[[i]])
            }
            else{
                spatial[[i]] <- spatial.offsetList[[i]]
            }
        }
    }
    else{
        stop("Invalid spatial.offsetList")    
    }
    
    
    funk <- function(spatial){
        if(any(class(spatial)=="fromXYZ")){
            spatial$Zm <- spatial$Zm*attr(spatial,"NC") # put back in 'normalising constant' so that spatialAtRisk acts as an offset (ie it no longer integrates to 1 over the study region.)
        }
        if(any(class(spatial)=="fromSPDF")){
            spatial$atrisk <- spatial$atrisk*attr(spatial,"NC")
            spatial$spdf$atrisk <- spatial$atrisk
        }
        return(spatial)
    }
    
    spatial <- lapply(spatial,funk)
    
    ################################################################
    # Create grid and FFT objects
    ################################################################
    
    study.region <- xyt$window
	
	## DEFINE LATTICE & CENTROIDS ##
		
	gridobj <- genFFTgrid(study.region=study.region,M=M,N=N,ext=ext,inclusion=inclusion)
	del1 <- gridobj$del1
	del2 <- gridobj$del2
	Mext <- gridobj$Mext
	Next <- gridobj$Next
	mcens <- gridobj$mcens
	ncens <- gridobj$ncens
	cellarea <- gridobj$cellarea
	cellInside <- gridobj$cellInside[1:M,1:N]
	
	xg <- mcens[1:M]
	yg <- ncens[1:N]
	
	## OBTAIN SPATIAL VALS ON LATTICE (LINEAR INTERPOLATION) ##

   
    spatialvals <- list()	
	for(i in 1:ntypes){
    	spatialvals[[i]] <- fftinterpolate(spatial[[i]],mcens,ncens,ext=ext)
    	spatialvals[[i]] <- spatialvals[[i]][1:M,1:N]
    	spatialvals[[i]] <- spatialvals[[i]]*cellInside ## note that this is not scaled, as we wish for lambda to act as an offset / (cellarea*sum(spatialvals))
	}
	
	# compute the base matrix of the covariance matrix
    bcb <- blockcircbase(x=mcens,y=ncens,sigma=sigma,phi=phi,model=spatial.covmodel,additionalparameters=covpars)    
    Qeigs <- eigenfrombase(inversebase(bcb)) # eigenvalues of Q (the precision matrix)
    rqe <- sqrt(Qeigs) # square root of the eigenvalues (used in computation)
    irqe <- 1/rqe # reciprocal root (commputation)
    
    # Interpolate covariate info onto FFT grid
    
    Zmat <- cov.interp.fft(formula=formula,W=study.region,regionalcovariates=regionalcovariates,pixelcovariates=pixelcovariates,mcens=mcens[1:M],ncens=ncens[1:N],cellInside=cellInside[1:M,1:N])
    
    zml <- getZmats(Zmat=Zmat,formulaList=formulaList)
    

    	
	################################################################
    # Simulate the data	
	################################################################	
	
	s <- Sys.time()
	
    commonfield <- YfromGamma(matrix(rnorm(Mext*Next),Mext,Next),invrootQeigs=irqe,mu=mu)[1:M,1:N]
    
    otherfields <- list()
    ratelist <- list() 
    ncases <- 0
    types <- c()
    cmatlist <- list()
    for(i in 1:ntypes){
        otherfields[[i]] <- YfromGamma(matrix(rnorm(Mext*Next),Mext,Next),invrootQeigs=irqe,mu=mu)[1:M,1:N]
        ratelist[[i]] <- as.vector(spatialvals[[i]]*cellarea*exp(matrix(zml[[i]]%*%betaList[[i]],M,N) + commonfield + otherfields[[i]]))
        
        if(any(is.na(ratelist[[i]]))){
            stop("Choice of beta gives Inf Poisson rates.")
        }
        cmatlist[[i]] <- matrix(rpois(M*N,ratelist[[i]]),M,N)
        nc <- sum(cmatlist[[i]])
        ncases <- ncases + nc
        
        if(nc==0){
            stop(paste("Choice of beta gives expected number of cases for type ",vartypes[i]," over region as",sum(ratelist[[i]])))
    	    
    	}
    	
        types <- c(types,rep(vartypes[i],nc))
    }
    
	if(ncases>100000){
	    warning("Number of cases for chosen parameters is large, ncases > 100000")
	}
	
	cases <- c()
	for(i in 1:ntypes){
        caseidx <- which(cmatlist[[i]]>0)
        caseidx <- unlist(sapply(caseidx,function(x){rep(x,cmatlist[[i]][x])}))
        cases <- rbind(cases,cbind(rep(xg,length(yg)),rep(yg,each=length(xg)))[caseidx,] + cbind(runif(sum(cmatlist[[i]]),-del1/2,del1/2),runif(sum(cmatlist[[i]]),-del2/2,del2/2)))
    } 
	
	xyt <- ppp(x=cases[,1],y=cases[,2],marks=types,window=owin)
	attr(xyt,"rejects") <- NULL # get rid of rejects: these are due to discrete approximation
	attr(xyt,"spatialoffset") <- spatial
	attr(xyt,"spatialoffsetGRID") <- spatialvals
	attr(xyt,"cellwidth") <- cellwidth
	attr(xyt,"sigma") <- sigma
	attr(xyt,"phi") <- phi
	attr(xyt,"spatial.covmodel") <- spatial.covmodel
	attr(xyt,"covpars") <- covpars
	attr(xyt,"ext") <- ext
    attr(xyt,"xvals") <- xg
    attr(xyt,"yvals") <- yg
    attr(xyt,"ratelist") <- lapply(ratelist,matrix,nrow=M,ncol=N)
    attr(xyt,"commonfield") <- commonfield
    attr(xyt,"otherfields") <- otherfields
    attr(xyt,"formula") <- formula
    attr(xyt,"Zmatlist") <- zml
    attr(xyt,"betaList") <- betaList 
    
    class(xyt) <- c("lgcpSimMultitypeSpatialPlusParameters","ppp")
    return(xyt)
}

##' intens.lgcpSimMultitypeSpatialPlusParameters function
##'
##' A function to return the cellwise Poisson intensity used during in constructing the simulated data.
##'
##' @method intens lgcpSimMultitypeSpatialPlusParameters
##' @param obj an object of class lgcpSimMultitypeSpatialPlusParameters
##' @param ... other parameters
##' @return the Poisson intensity
##' @usage "intens(obj, ...)"
##' @export
intens.lgcpSimMultitypeSpatialPlusParameters <- function(obj,...){
    M <- dim(attr(obj,"spatialoffsetGRID"))[1]
    N <- dim(attr(obj,"spatialoffsetGRID"))[2]
    spatialvals <- attr(obj,"spatialoffsetGRID")[1:M,1:N]
    cellarea <- diff(attr(obj,"xvals")[1:2])*diff(attr(obj,"yvals")[1:2])
    Zmat <- attr(obj,"Zmat")
    beta <- attr(obj,"beta")
    truefield <- attr(obj,"truefield")
    return(spatialvals*cellarea*exp(matrix(Zmat%*%beta,M,N) + truefield))
}

##' autocorrMultitype function
##'
##' A function to compute cell-wise autocorrelation in the latent field at specifiec lags 
##'
##' @param x an object of class lgcpPredictMultitypeSpatialPlusParameters 
##' @param lags the lags at which to compute the autocorrelation
##' @param fieldno the index of the lateyt field, the i in Y_i, see the help file for lgcpPredictMultitypeSpatialPlusParameters. IN diagnostic checking ,this command should be called for each field in the model.
##' @param inWindow an observation owin window on which to compute the autocorrelations, can speed up calculation. Default is x$xyt$window, set to NULL for full grid.
##' @param crop2parentwindow logical: whether to only compute autocorrelations for cells inside x$xyt$window (the 'parent window')
##' @param ... other arguments
##' @return an array, the [,,i]th slice being the grid of cell-wise autocorrelations.
##' @export

autocorrMultitype <- function(x,lags,fieldno,inWindow=x$xyt$window,crop2parentwindow=TRUE,...){
    tidx <- NULL
    if (is.null(x$gridfunction$dirname)){
        stop("dump2dir not specified, MCMC output must have be dumped to disk to use this function. See ?dump2dir.")
    }
    if (length(tidx)>1){
        stop("tidx should either be NULL, or a vector of length 1")
    }
    fn <- paste(x$gridfunction$dirname,"simout.nc",sep="")
    ncdata <- nc_open(fn)
    datadim <- ncdata$var$simrun$varsize
    if (is.null(tidx)){
        tidx <- datadim[3]
    }
    if(!is.null(inWindow)){
        if (crop2parentwindow){
            grinw <- matrix(as.logical(gridInWindow(x$mcens,x$ncens,inWindow) * gridInWindow(x$mcens,x$ncens,x$xyt$window)),length(xvals(x)),length(yvals(x)))
        }
        else{
            grinw <- gridInWindow(x$mcens,x$ncens,inWindow)
        }
    }
    result <- array(dim=c(datadim[1],datadim[2],length(lags)))
    sampcount <- datadim[4]
    pb <- txtProgressBar(min=1,max=datadim[1]*datadim[2],style=3)
    trigger <- TRUE
    for (i in 1:datadim[1]){
        for(j in 1:datadim[2]){
            setTxtProgressBar(pb,j+(i-1)*datadim[2])
            if(!is.null(inWindow)){
                if (isTRUE(grinw[i,j])){
                    tr <- ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(i,j,tidx,fieldno,1), count=c(1,1,1,1,-1))
                    result[i,j,] <- acf(tr,plot=FALSE)$acf[lags+1]
                }
                else{
                    result[i,j,] <- NA
                }
            }
            else{
                tr <- ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(i,j,tidx,fieldno,1), count=c(1,1,1,1,-1))
                result[i,j,] <- acf(tr,plot=FALSE)$acf[lags+1]
            }
            
            if(trigger){
                ltst <- length(acf(ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(1,1,tidx,fieldno,1), count=c(1,1,1,1,-1)),plot=FALSE)$acf)
                tst <- (lags+1)>ltst
                if(any(tst)){
                    stop(paste("Cannot return lag",ltst,"or above."))
                }
                trigger <- FALSE
            }
        }
    }
    close(pb)
    nc_close(ncdata)
    attr(result,"lags") <- lags
    attr(result,"xcoords") <- x$mcens
    attr(result,"ycoords") <- x$ncens
    attr(result,"window") <- NULL
    if(!is.null(inWindow)){
        attr(result,"window") <- inWindow
    }
    attr(result,"ParentWindow") <- x$xyt$window
    class(result) <- c("lgcpAutocorr","array")
    return(result)   
}
