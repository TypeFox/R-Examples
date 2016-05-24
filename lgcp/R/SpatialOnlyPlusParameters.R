##' lgcpPredictSpatialPlusPars function
##'
##' A function to deliver fully Bayesian inference for the spatial log-Gaussian Cox process.\cr
##'
##' See the vignette "Bayesian_lgcp" for examples of this code in use.\cr
##'
##' The model for the data is as follows:\cr
##' \cr
##' X(s) ~ Poisson[R(s)]\cr
##' \cr
##' R(s) = C_A lambda(s) exp[Z(s)beta+Y(s)]\cr
##' \cr
##'
##' Here X(s) is the number of events in the cell of the computational grid containing s, R(s) is the Poisson rate, 
##' C_A is the cell area, lambda(s) is a known offset, Z(s) is a vector of measured covariates and Y(s) is the 
##' latent Gaussian process on the computational grid. The other parameters in the model are beta, the covariate effects; 
##' and eta=[log(sigma),log(phi)], the parameters of the process Y on an appropriately transformed (in this case log) scale.
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
##' @param formula a formula object of the form X ~ var1 + var2 etc. The name of the dependent variable must be "X". Only accepts 'simple' formulae, such as the example given.
##' @param sd a spatstat ppp object
##' @param Zmat design matrix Z (see below) constructed with getZmat
##' @param model.priors model priors, set using lgcpPrior
##' @param model.inits model initial values. The default is NULL, in which case lgcp will use the prior mean to initialise eta and beta will be initialised from an oversispersed glm fit to the data. Otherwise use lgcpInits to specify.
##' @param spatial.covmodel choice of spatial covariance function. See ?CovFunction
##' @param cellwidth the width of computational cells
##' @param poisson.offset A SpatialAtRisk object defining lambda (see below)
##' @param mcmc.control MCMC paramters, see ?mcmcpars
##' @param output.control output choice, see ?setoutput   
##' @param gradtrunc truncation for gradient vector equal to H parameter Moller et al 1998 pp 473. Default is Inf, which means no gradient truncation, which seems to work in most settings.
##' @param ext integer multiple by which grid should be extended, default is 2. Generally this will not need to be altered, but if the spatial correlation decays slowly, increasing 'ext' may be necessary.
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former, the default, includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window.
##' @return an object of class lgcpPredictSpatialOnlyPlusParameters
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
##' \link{lgcpPredictAggregateSpatialPlusPars}, \link{lgcpPredictSpatioTemporalPlusPars}, 
##' \link{lgcpPredictMultitypeSpatialPlusPars},
##' \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export 
   
lgcpPredictSpatialPlusPars <- function( formula,
                                        sd,
                                        Zmat=NULL,       
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
    
    if(sd$window$type=="rectangle"){
	    sd$window <- as.polygonal(sd$window)
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
    
    if(is.null(poisson.offset)){
        poisson.offset <- list( X=seq(sd$window$xrange[1],sd$window$xrange[2],length.out=100),
                                Y=seq(sd$window$yrange[1],sd$window$yrange[2],length.out=100),
                                Zm=matrix(1,100,100))
    }	
					
	if(!any(class(poisson.offset)=="spatialAtRisk")){			
        spatial <- spatialAtRisk(poisson.offset)
    }
    else{
        spatial <- poisson.offset
    }
    
    if(any(class(spatial)=="fromXYZ")){
        spatial$Zm <- spatial$Zm*attr(spatial,"NC") # put back in 'normalising constant' so that spatialAtRisk acts as an offset (ie it no longer integrates to 1 over the study region.)
    }
    if(any(class(spatial)=="fromSPDF")){
        spatial$atrisk <- spatial$atrisk*attr(spatial,"NC")
        spatial$spdf$atrisk <- spatial$atrisk
    }
    
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
	
	spts <- SpatialPoints(cbind(rep(mcens[1:M],N),rep(ncens[1:N],each=M)))
	spatialcovs <- data.frame()
	
	if(is.null(Zmat)){
	    Zmat <- cov.interp.fft(formula=formula,W=study.region,regionalcovariates=regionalcovariates,pixelcovariates=pixelcovariates,mcens=mcens[1:M],ncens=ncens[1:N],cellInside=cellInside[1:M,1:N])
	}
	else{
        if(!isTRUE(all.equal(mcens[1:M],attr(Zmat,"mcens")))|!isTRUE(all.equal(ncens[1:N],attr(Zmat,"ncens")))){
            stop("FFT grid and Zmat are on different grids. Please recompute Zmat using 'getZmat'.")
        }    	
	}
		
	
	## OBTAIN POISSON OFFSET ON LATTICE (LINEAR INTERPOLATION) ##
	
	spatial.offset <- fftinterpolate(spatial,mcens,ncens,ext=ext)
	spatial.offset <- spatial.offset*cellInside
	## no longer required as here this is an offset: spatialvals <- spatialvals / (cellarea*sum(spatialvals))
	
    	
	################################################################
		
	
		
    ###
    # Set up MCMC loop, required to compute nsamp, below
    ###
    mLoop = mcmcLoop(N=mcmc.control$mala.length,burnin=mcmc.control$burnin,thin=mcmc.control$retain,progressor=mcmcProgressTextBar)	
	
	# issue warning if dumping information to disc
	nsamp <- floor((mLoop$N-mLoop$burnin)/mLoop$thin)
	if (!is.null(output.control$gridfunction) & class(output.control$gridfunction)[1]=="dump2dir"){
    	cat("WARNING: disk space required for saving is approximately ",round(nsamp*object.size(array(runif((M)*(N)),dim=c((M),(N))))/1024^2,2)," Mb, ",sep="")
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
	
    
    lg <- MALAlgcpSpatial.PlusPars( mcmcloop=mLoop,
                                    inits=mcmc.control$inits,
                                    adaptivescheme=mcmc.control$adaptivescheme,
                                    M=M,
                                    N=N,
                                    Mext=Mext,
                                    Next=Next,
                                    mcens=mcens,
                                    ncens=ncens,
                                    formula=formula,
                                    Zmat=Zmat,
                                    model.priors=model.priors,
                                    model.inits=model.inits,
                                    fftgrid=gridobj,
                                    spatial.covmodel=spatial.covmodel,
                                    nis=nis,
                                    cellarea=cellarea,
                                    spatialvals=spatial.offset,
                                    cellInside=cellInside,
                                    MCMCdiag=mcmc.control$MCMCdiag,
                                    gradtrunc=gradtrunc,
                                    gridfun=gridfun,
                                    gridav=gridav,
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
    lg$poisson.offset <- spatial.offset
    lg$priors <- model.priors
    lg$covFct <- spatial.covmodel
	
	class(lg) <- c("lgcpPredictSpatialOnlyPlusParameters","lgcpPredict","lgcpobject")	
	
	return(lg)															
					
}	
	

##' MALAlgcpSpatial.PlusPars function
##'
##' A function to run the MCMC algorithm for spatial point process data. Not for general purpose use.
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
##' @param Zmat design matrix constructed using getZmat
##' @param model.priors model priors, constructed using lgcpPrior 
##' @param model.inits initial values for the MCMC
##' @param fftgrid an objects of class FFTgrid, see genFFTgrid
##' @param spatial.covmodel spatial covariance model, consructed with CovFunction
##' @param nis cell counts on the etended grid
##' @param cellarea the cell area
##' @param spatialvals inerpolated poisson offset on fft grid 
##' @param cellInside 0-1 matrix indicating inclusion in the observation window 
##' @param MCMCdiag not used
##' @param gradtrunc gradient truncation parameter
##' @param gridfun used to specify other actions to be taken, e.g. dumping MCMC output to disk.
##' @param gridav used for computing Monte Carlo expectations online
##' @param d matrix of toral distances
##' @return output from the MCMC run
##' @export

MALAlgcpSpatial.PlusPars <- function(   mcmcloop,
                                        inits,
                                        adaptivescheme,
                                        M,
                                        N,
                                        Mext,
                                        Next,
                                        mcens,
                                        ncens,
                                        formula,
                                        Zmat,
                                        model.priors,
                                        model.inits,
                                        fftgrid,
                                        spatial.covmodel,
                                        nis,
                                        cellarea,
                                        spatialvals,
                                        cellInside,
                                        MCMCdiag,
                                        gradtrunc,
                                        gridfun,
                                        gridav,
                                        d){
                            
    SpatialOnlyMode <- TRUE
    SpatialPlusParameters <- TRUE
    SpatioTemporalPlusParameters <- FALSE
    MultiTypeMode <- FALSE
    
    

    
    cellInsideLogical <- as.logical(cellInside)
      
    M <- M
    N <- N
    temporal.fitted <- Inf # note this line is here for gridFunction and gridAverage methods and is not used otherwise
    nlevs <- NULL # note this line is here for gridFunction and gridAverage methods and is not used otherwise
    GFinitialise(gridfun) # note these two lines must come after M and N have been computed or defined
	GAinitialise(gridav) # note these two lines must come after M and N have been computed or defined
    
    nsamp <- 0
    icount <- 0
    MCMCacc <- 0
    y.mean <- matrix(0,M,N)
    y.var <- matrix(0,M,N)
    EY.mean <- matrix(0,M,N)
    EY.var <- matrix(0,M,N)
    
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
	cp <- CovParameters(list(sigma=etainv[1],phi=etainv[2]))
	
	off <- c() # required so the the call to glm passes CRAN check
	rm(off)
	dfr <- attr(Zmat,"data.frame")
	dfr$X <- nis[cellInsideLogical]
	dfr$off <- log(cellarea*spatialvals[cellInsideLogical])
	dfr$off[is.infinite(dfr$off)] <- NA # this excludes cells where the rate is zero
	mod <- glm(formula,data=dfr,family=quasipoisson,offset=off)    	
    if(is.null(model.inits$betainit)){
    	betaval <- BetaParameters(coefficients(mod))
    }
    else{
        betaval <- BetaParameters(model.inits$betainit)
    }
    
    Z <- matrix(0,Next*Mext,ncol=ncol(Zmat))
    tm <- matrix(FALSE,Mext,Next)
    tm[1:M,1:N] <- TRUE
    Z[as.vector(tm),] <- Zmat
    Zt <- t(Z)  
    

    ###
	# Now loop
	###
	
    neta <- length(etaval)
	nbeta <- length(betaval)	
	
	# compute approximate variance matrix of Y
	y <- log(nis/(cellarea*spatialvals* exp(matrix(as.vector(Z%*%betaval),Mext,Next))))
    y[is.na(y) | is.infinite(y)] <- -cp$sigma^2/2   
    GPdummy <- GPrealisation(gamma=matrix(rnorm(Mext*Next),Mext,Next),fftgrid=fftgrid,covFunction=spatial.covmodel,covParameters=cp,d=d)
   
    gammainit <- GammafromY(Y=y,rootQeigs=GPdummy$rootQeigs,mu=GPdummy$CovParameters$mu)       
    GP <- GPrealisation(gamma=gammainit,fftgrid=fftgrid,covFunction=spatial.covmodel,covParameters=cp,d=d)
    eMat <- cellarea*spatialvals*exp(matrix(as.vector(Z%*%betaval),nrow(GP$gamma),ncol(GP$gamma)))*GP$expY
    fi <- (1/length(GP$Y))*Re(fft(GP$invrootQeigs*fft((1/length(GP$Y))*Re(fft(GP$invrootQeigs*fft(eMat,inverse=TRUE))),inverse=TRUE)))
    gammaVar <- (1.65^2/((Mext*Next)^(1/3)))*(-1)/(-1-fi)
    GAMMAVAR <- gammaVar
    rootgammaVar <- sqrt(gammaVar)
    # compute variance matrix of beta
    EZ <- as.vector(cellarea*spatialvals*exp(matrix(as.vector(Z%*%betaval),nrow(GP$gamma),ncol(GP$gamma))+GP$Y))    
    CB <- matrix(NA,nbeta,nbeta)    
    for(i in 1:nbeta){
        for(j in 1:nbeta){
            CB[i,j] <- sum(Z[,i]*Z[,j]*EZ) + model.priors$betaprior$precision[i,j] # observed information
        }
    } 

    
    # compute variance matrix of eta

 

    additional_scaleconst <- (0.234/0.574)     
    
    etaCovMat <- try(additional_scaleconst*solve((-1)*GPdrv2(GP=GP,prior=model.priors,Z=Z,Zt=Zt,eta=etaval,beta=betaval,nis=nis,cellarea=cellarea,spatial=spatialvals,gradtrunc=Inf,fftgrid=fftgrid,covfunction=spatial.covmodel,d=d,eps=1e-6)$hess)) #mean(diff(mcens[1:2]),diff(ncens[1:2]))/100)$hess) # eps=mean(diff(mcens[1:2]),diff(ncens[1:2]))/100 is approx 1/100 cellwidth
    etaCovMattest <- FALSE 
    #browser()   
    if(class(etaCovMat)!="try-error"){
        etaCovMattest <- !all(eigen(etaCovMat)$values>0)
    }
    if(class(etaCovMat)=="try-error" | etaCovMattest){

        cat("Computing posterior variance w.r.to eta via finite differencing failed, trying global alternative ...\n")
        lensq <- 10
        sqsigma <- seq(model.priors$etaprior$mean[1]-2*sqrt(model.priors$etaprior$variance[1,1]),model.priors$etaprior$mean[1]+2*sqrt(model.priors$etaprior$variance[1,1]),length.out=lensq)
        sqphi <- seq(model.priors$etaprior$mean[2]-2*sqrt(model.priors$etaprior$variance[2,2]),model.priors$etaprior$mean[2]+2*sqrt(model.priors$etaprior$variance[2,2]),length.out=lensq)
        ltarmat <- matrix(NA,lensq,lensq)
        for (i in 1:lensq){
            for(j in 1:lensq){                        
                cpA <- CovParameters(list(sigma=exp(sqsigma[i]),phi=exp(sqphi[j])))      
                gpA <- GPrealisation(gamma=GP$gamma,fftgrid=fftgrid,covFunction=spatial.covmodel,covParameters=cpA,d=d)
                matent <- try(target.and.grad.spatialPlusPars(GP=gpA,prior=model.priors,Z=Z,Zt=Zt,eta=c(sqsigma[i],sqphi[j]),beta=betaval,nis=nis,cellarea=cellarea,spatial=spatialvals,gradtrunc=Inf)$logtarget)
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
        try(etaCovMat <- additional_scaleconst*solve((-1)*matrix(c(2*coef[2],coef[6],coef[6],2*coef[4]),2,2)))
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
    }
    else{
        cat("Approximate of posterior variance w.r.to eta found:\n")
    }


	sigma_eta <- (2.38^2/length(etaval))*etaCovMat #model.priors$etaprior$variance #diag(c(var(log(vecrec)),model.priors$etaprior$variance[2,2]))
	SIGMA_ETA <- sigma_eta
	Q_eta <- solve(sigma_eta)

	sigma_beta <- (1.65^2/(length(betaval)^(1/3)))*solve(CB) #vcov(mod)	
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
    
    cons <- 1
   
    consrec <- rep(NA,1000)

	betaetarec <- c()
      
	
    GP <- GPrealisation(gamma=matrix(0,Mext,Next),fftgrid=fftgrid,covFunction=spatial.covmodel,covParameters=cp,d=d)
    oldtags <- target.and.grad.spatialPlusPars(GP=GP,prior=model.priors,Z=Z,Zt=Zt,eta=etaval,beta=betaval,nis=nis,cellarea=fftgrid$cellarea,spatial=spatialvals,gradtrunc=gradtrunc)
    
    tarrec <- oldtags$logtarget
    hallrec <- h 
    
    reject_its <- c()   
       
    while(nextStep(mcmcloop)){
    
        propmeans_beta <- betaval + (h/2)*sigma_beta%*%oldtags$gradbeta     
        propbeta <- BetaParameters(as.vector(propmeans_beta+sqrt(h)*sigma_beta_chol%*%rnorm(nbeta)))
    
        propmeans_gamma <- GP$gamma + (h/2)*gammaVar*oldtags$gradgamma
        propmeans_eta <- etaval 
        propeta <- as.vector(propmeans_eta + sqrt(h)*sigma_eta_chol%*%rnorm(neta)) 

        
        propetainv <- etainvtrans(propeta) 
        propcp <- CovParameters(list(sigma=propetainv[1],phi=propetainv[2]))
        propGP <- GPrealisation(gamma=propmeans_gamma+ sqrt(h)*rootgammaVar*rnorm(Mext*Next),fftgrid=fftgrid,covFunction=spatial.covmodel,covParameters=propcp,d=d)
                
                
          
        proptags <- target.and.grad.spatialPlusPars(GP=propGP,prior=model.priors,Z=Z,Zt=Zt,eta=propeta,beta=propbeta,nis=nis,cellarea=fftgrid$cellarea,spatial=spatialvals,gradtrunc=gradtrunc)    
        
        revpropmeans_gamma <- propGP$gamma + (h/2)*gammaVar*proptags$gradgamma
        revpropmeans_eta <- propeta 
        revpropmeans_beta <- propbeta + (h/2)*sigma_beta%*%proptags$gradbeta

        
        ac <- exp(proptags$logtarget-oldtags$logtarget-sum((GP$gamma-revpropmeans_gamma)^2/(2*h*gammaVar)) + 
                            sum((propGP$gamma-propmeans_gamma)^2/(2*h*gammaVar)) +
                            (-(0.5/h)*t(betaval-revpropmeans_beta)%*%Q_beta%*%(betaval-revpropmeans_beta)) -
                            (-(0.5/h)*t(propbeta-propmeans_beta)%*%Q_beta%*%(propbeta-propmeans_beta)))
                            
 
                            
        ac <- min(ac,1)

        icount <- icount + 1        
        MCMCacc <- ((icount-1)/icount) * MCMCacc + ac/icount
        

		if (proptags$logtarget==-Inf | is.na(ac) | is.nan(ac)){ # gradient truncation insufficient, so reduce

	        warning("One possible cause of this warning is that there may be evidence in the data for quite large values of the spatial correlation parameter. If this is the case, then this warning means that the MCMC chain has wandered into a region of the phi-space that causes the variance matrix of Y (computed by the DFT) to become non positive-definite. One possible way of rectifying this issue is to restart the chain using a larger value of 'ext' in the call to lgcpPredictSpatialPlusPars. You should do this if the warning message is repeated many times. If this warning message appears at all then you should be warned that inference from this run may not be reliable: the proposed move at this iteration was rejected.",immediate.=TRUE)
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
        	y.mean <- ((nsamp-1)/nsamp) * y.mean + GP$Y[1:(M),1:(N)]/nsamp
        	EY.mean <- ((nsamp-1)/nsamp) * EY.mean + GP$expY[1:(M),1:(N)]/nsamp
        	if (nsamp>1){
    			y.var <- ((nsamp-2)/(nsamp-1))*y.var + (nsamp/(nsamp-1)^2)*(y.mean-GP$Y[1:(M),1:(N)])^2
    			EY.var <- ((nsamp-2)/(nsamp-1))*EY.var + (nsamp/(nsamp-1)^2)*(EY.mean-GP$expY[1:(M),1:(N)])^2
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
    retlist$sigmaetarec <- sigmaetarec
    retlist$sigmabetarec <- sigmabetarec
    retlist$hallrec <- hallrec
    retlist$tarrec <- tarrec
    retlist$Z <- Z
    retlist$reject_its <- reject_its
	
	return(retlist)                          
}                        

##' target.and.grad.spatialPlusPars function
##'
##' A function to compute the target and gradient for the Bayesian spatial LGCP
##'
##' @param GP an object created using GPrealisation
##' @param prior the model priors, created using lgcpPrior
##' @param Z the design matrix on the FFT grid
##' @param Zt transpose of the design matrix
##' @param eta the paramters, eta
##' @param beta the parameters, beta
##' @param nis cell counts on the FFT grid
##' @param cellarea the cell area
##' @param spatial poisson offset
##' @param gradtrunc the gradient truncation parameter
##' @return the target and graient for this model
##' @export

target.and.grad.spatialPlusPars <- function(GP,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc){

    Zbeta <- matrix(as.vector(Z%*%beta),nrow(GP$gamma),ncol(GP$gamma))
    expZbeta <- exp(Zbeta)
    pri <- EvaluatePrior(etaParameters=eta,betaParameters=beta,prior=prior)    

    e <- spatial*expZbeta*GP$expY*cellarea
    
    NminusE <- nis - e
    
    logtarget <- -(1/2)*sum(GP$gamma^2) + sum((Zbeta+GP$Y)*nis - e) + pri$etacontrib$loglik + pri$betacontrib$loglik # note that both nis=0, logspat=0 and spatial=0 outside of observation window, so this effectively limits summation to the observation window only
    logtarget <- as.vector(logtarget)
    
    if(!is.infinite(gradtrunc)){
        expYtrunc <- GP$expY
        expYtrunc[expYtrunc>gradtrunc] <- gradtrunc  
        gradgamma <- (-1)*GP$gamma + (1/length(GP$Y))*Re(fft(GP$invrootQeigs*fft(nis-spatial*expZbeta*expYtrunc*cellarea,inverse=TRUE))) 
    }
    else{
        gradgamma <- (-1)*GP$gamma + (1/length(GP$Y))*Re(fft(GP$invrootQeigs*fft(NminusE,inverse=TRUE)))
    }
    
    gradbeta <- as.vector(Zt%*%as.vector(NminusE) + pri$betacontrib$gradcontrib)
     
    gradeta <- NULL    

    return(list(logtarget=logtarget,gradgamma=gradgamma,gradbeta=gradbeta,gradeta=gradeta,e=e))
}


###############################################################
#
# !! NOTE THIS IS NOT THE DOCUMENTATION OBJECT, IT IS BELOW !!
#
###############################################################
# lgcpSimSpatialCovariates function
#
# A function to simulate from a log gaussian process
#
# @param owin observation window
# @param regionalcovariates an optional SpatialPolygonsDataFrame object giving covariate information measured at the regional level. None of the data columns should have name 'X'.
# @param pixelcovariates an optional SpatialPixelsData frame giving covariate information measured on a grid of points over the region. None of the data columns should have name 'X'.
# @param beta vector of regression parameters. The length of the vector depends on the number of paramters, whether there is an intercept and wethere there are any factor variables.
# The first entry of the beta vector should be the intercept term (if applicable), followed by the parameters for each of the variables as
# they appear in the model formula. For factor-like variables (which should ideally be prepared as factors in advance) there should be a 
# parameter for EACH of the levels of the factor.  
# @param spatial.offset an object that can be coerced to one of class spatialAtRisk
# @param cellwidth width of cells
# @param model.parameters parameters of model, see ?lgcppars. Only set sigma and phi for spatial model.
# @param spatial.covmodel spatial covariance function, default is exponential, see ?CovarianceFct
# @param covpars vector of additional parameters for spatial covariance function, in order they appear in chosen model in ?CovarianceFct
# @param nprobe number of interpolation points over which to compute covariate information on FFT grid
# @param ext how much to extend the parameter space by. Default is 2.
# @param plot logical, wheter to plot the latent field.
# @param inclusion
# @return a ppp object containing the data
# @export

##' lgcpSimSpatialCovariates function
##'
##' A function to  simulate a spatial LGCP.
##'
##' @param formula a formula of the form X ~ var1 + var2 etc.
##' @param owin the observation window on which to do the simulation
##' @param regionalcovariates an optional object of class SpatialPolygonsDataFrame containing covariates
##' @param pixelcovariates an optional object of class SpatialPixelsDataFrame containing covariates
##' @param Zmat optional design matrix, if the polygon/polygon overlays have already been computed
##' @param beta the parameters, beta for the model
##' @param poisson.offset the poisson offet, created using a SpatialAtRisk.fromXYZ class of objects
##' @param cellwidth the with of cells on which to do the simulation
##' @param model.parameters the paramters of the model eg list(sigma=1,phi=0.2)
##' @param spatial.covmodel the choice of spatial covariance model, can be anything from the RandomFields covariance function, CovariacenFct.
##' @param covpars additional covariance parameters, for the chosen model, optional.
##' @param ext the amount by which to extend the observation grid in each direction, default is 2
##' @param plot whether to plot the resulting data
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former, the default, includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window.
##' @return a ppp onject containing the simulated data
##' @export

lgcpSimSpatialCovariates <- function(   formula,    
                                        owin,
                                        regionalcovariates=NULL,
                                        pixelcovariates=NULL,
                                        Zmat=NULL,
                                        beta,
                                        poisson.offset=NULL,
                                        cellwidth,
                                        model.parameters,
                                        spatial.covmodel="exponential",
                                        covpars=c(),
                                        ext=2,
                                        plot=FALSE,
                                        inclusion="touching"){
                                        
    nprobe <- 100000                                                                                      
    spatial.offset <- poisson.offset    
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
	if(is.null(spatial.offset)){
        spatial <- spatialAtRisk(list(X=seq(xyt$window$xrange[1],xyt$window$xrange[2],length.out=M),Y=seq(xyt$window$yrange[1],xyt$window$yrange[2],length.out=N),Zm=matrix(1,M,N)))
    }
    else{
        if(!any(class(spatial.offset)=="spatialAtRisk")){		
            spatial <- spatialAtRisk(spatial.offset)
        }
        else{
            spatial <- spatial.offset
        }
    }
    
    if(any(class(spatial)=="fromXYZ")){
        spatial$Zm <- spatial$Zm*attr(spatial,"NC") # put back in 'normalising constant' so that spatialAtRisk acts as an offset (ie it no longer integrates to 1 over the study region.)
    }
    if(any(class(spatial)=="fromSPDF")){
        spatial$atrisk <- spatial$atrisk*attr(spatial,"NC")
        spatial$spdf$atrisk <- spatial$atrisk
    }
    
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
	
	spatialvals <- fftinterpolate(spatial,mcens,ncens,ext=ext)
	spatialvals <- spatialvals[1:M,1:N]
	spatialvals <- spatialvals*cellInside ## note that this is not scaled, as we wish for lambda to act as an offset / (cellarea*sum(spatialvals))
	
	# compute the base matrix of the covariance matrix
    bcb <- blockcircbase(x=mcens,y=ncens,sigma=sigma,phi=phi,model=spatial.covmodel,additionalparameters=covpars)    
    Qeigs <- eigenfrombase(inversebase(bcb)) # eigenvalues of Q (the precision matrix)
    rqe <- sqrt(Qeigs) # square root of the eigenvalues (used in computation)
    irqe <- 1/rqe # reciprocal root (commputation)
    
    # Interpolate covariate info onto FFT grid
    
    #Zmat <- cov.interp.fft(formula=formula,W=study.region,regionalcovariates=regionalcovariates,pixelcovariates=pixelcovariates,mcens=mcens[1:M],ncens=ncens[1:N],nprobe=nprobe,cellInside=cellInside[1:M,1:N]) 
    if(is.null(Zmat)){
	    Zmat <- cov.interp.fft(formula=formula,W=study.region,regionalcovariates=regionalcovariates,pixelcovariates=pixelcovariates,mcens=mcens[1:M],ncens=ncens[1:N],cellInside=cellInside[1:M,1:N])
	}
	else{
        if(!isTRUE(all.equal(mcens[1:M],attr(Zmat,"mcens")))|!isTRUE(all.equal(ncens[1:N],attr(Zmat,"ncens")))){
            stop("FFT grid and Zmat are on different grids. Please recompute Zmat using 'getZmat'.")
        }    	
	}    
    	
	################################################################
    # Simulate the data	
	################################################################	

	s <- Sys.time()
    truefield <- YfromGamma(matrix(rnorm(Mext*Next),Mext,Next),invrootQeigs=irqe,mu=mu)[1:M,1:N]
    rate <- as.vector(spatialvals*cellarea*exp(matrix(Zmat%*%beta,M,N) + truefield))
    
    if(any(is.na(rate))){
        stop("Choice of beta gives Inf Poisson rates.")
    }
    cmat <- matrix(rpois(M*N,rate),M,N)
    ncases <- sum(cmat)
    if(ncases==0){
        cat("Choice of beta gives expected number of cases over region as",sum(rate),"\n")
	    stop("No data generated for chosen parameters")
	}
	if(ncases>100000){
	    warning("Number of cases for chosen parameters is large, ncases > 100000")
	}
    caseidx <- which(cmat>0)
    caseidx <- unlist(sapply(caseidx,function(x){rep(x,cmat[x])}))
    cases <- cbind(rep(xg,length(yg)),rep(yg,each=length(xg)))[caseidx,] + cbind(runif(ncases,-del1/2,del1/2),runif(ncases,-del2/2,del2/2))
    if(plot){
        rate[rate==0] <- NA
        image(xg,yg,matrix(rate,M,N)^0.25)
        points(cases,pch="+",cex=0.5)
    }    
	#browser()
	xyt <- ppp(x=cases[,1],y=cases[,2],window=owin)
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
    attr(xyt,"rate") <- matrix(rate,M,N)
    attr(xyt,"truefield") <- truefield
    attr(xyt,"formula") <- formula
    attr(xyt,"Zmat") <- Zmat
    names(beta) <- colnames(Zmat)
    attr(xyt,"beta") <- beta 
    
    class(xyt) <- "lgcpSimSpatialPlusParameters"
    return(xyt)
}

##' intens.lgcpSimSpatialPlusParameters function
##'
##' A function to return the cellwise Poisson intensity used during in constructing the simulated data.
##'
##' @method intens lgcpSimSpatialPlusParameters
##' @param obj an object of class lgcpSimSpatialPlusParameters
##' @param ... other parameters
##' @return the Poisson intensity
##' @export
intens.lgcpSimSpatialPlusParameters <- function(obj,...){
    M <- dim(attr(obj,"spatialoffsetGRID"))[1]
    N <- dim(attr(obj,"spatialoffsetGRID"))[2]
    spatialvals <- attr(obj,"spatialoffsetGRID")[1:M,1:N]
    cellarea <- diff(attr(obj,"xvals")[1:2])*diff(attr(obj,"yvals")[1:2])
    Zmat <- attr(obj,"Zmat")
    beta <- attr(obj,"beta")
    truefield <- attr(obj,"truefield")
    return(spatialvals*cellarea*exp(matrix(Zmat%*%beta,M,N) + truefield))
}


