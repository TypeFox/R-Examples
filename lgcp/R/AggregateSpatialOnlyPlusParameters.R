##' lgcpPredictAggregateSpatialPlusPars function
##'
##' A function to deliver fully Bayesian inference for the aggregated spatial log-Gaussian Cox process.\cr
##'
##' See the vignette "Bayesian_lgcp" for examples of this code in use.\cr
##'
##' In this case, we OBSERVE case counts in the regions of a SpatialPolygonsDataFrame; the counts are stored as a variable, X.
##' The model for the UNOBSERVED data, X(s), is as follows:\cr
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
##' @param spdf a SpatialPolygonsDataFrame object with variable "X", the event counts per region.
##' @param Zmat design matrix Z (see below) constructed with getZmat
##' @param overlayInZmat if the covariate information in Zmat also comes from spdf, set to TRUE to avoid replicating the overlay operations. Default is FALSE. 
##' @param model.priors model priors, set using lgcpPrior
##' @param model.inits model initial values. The default is NULL, in which case lgcp will use the prior mean to initialise eta and beta will be initialised from an oversispersed glm fit to the data. Otherwise use lgcpInits to specify.
##' @param spatial.covmodel choice of spatial covariance function. See ?CovFunction
##' @param cellwidth the width of computational cells
##' @param poisson.offset A SpatialAtRisk object defining lambda (see below)
##' @param mcmc.control MCMC paramters, see ?mcmcpars
##' @param output.control output choice, see ?setoutput   
##' @param gradtrunc truncation for gradient vector equal to H parameter Moller et al 1998 pp 473. Default is Inf, which means no gradient truncation, which seems to work in most settings.
##' @param ext integer multiple by which grid should be extended, default is 2. Generally this will not need to be altered, but if the spatial correlation decays slowly, increasing 'ext' may be necessary.
##' @param Nfreq the sampling frequency for the cell counts. Default is every 101 iterations.
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former, the default, includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window.
##' @return an object of class lgcpPredictAggregateSpatialPlusParameters
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
##' \link{lgcpPredictSpatialPlusPars},  \link{lgcpPredictSpatioTemporalPlusPars}, 
##' \link{lgcpPredictMultitypeSpatialPlusPars},
##' \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export
  
lgcpPredictAggregateSpatialPlusPars <- function( formula,
                                                spdf,
                                                Zmat=NULL,
                                                overlayInZmat=FALSE,       
                        					    model.priors,
                        					    model.inits=lgcpInits(),
                        					    spatial.covmodel,
                        					    cellwidth=NULL,
                        					    poisson.offset=NULL,				
                        					    mcmc.control,
                        					    output.control=setoutput(),
                        					    gradtrunc=Inf,
                        					    ext=2,
                        					    Nfreq = 101,
                        					    inclusion="touching"){

    nprobe <- 100000    
    gridsize <- NULL
    if(overlayInZmat){
        ol <- attr(Zmat,"polygonOverlay")
    }    
    else{
        fftpoly <- grid2spoly(attr(Zmat,"mcens"),attr(Zmat,"ncens"))
        cat("Overlaying fft grid onto spdf ...\n")   
        ol <- gOverlay(fftpoly,spdf)
    }
    
    cellcts <- rep(0,attr(Zmat,"M")*attr(Zmat,"N"))
    grinf <- expand.grid(attr(Zmat,"mcens"),attr(Zmat,"ncens"))
    spts <- c()
    #browser()
    for(i in 1:length(spdf)){
        if(spdf$X[i]==0){
            next
        }
        else{
            cellinfo <- ol$info[ol$info$polyidx==i,]
            if(nrow(cellinfo)>1){
                idx <- sample(cellinfo$grididx,spdf$X[i],replace=TRUE,prob=cellinfo$area)
            }
            else{
                idx <- rep(cellinfo$grididx,spdf$X[i])
            }
            idxunq <- unique(idx)
            sapply(idxunq,function(x){spts <<- rbind(spts,matrix(rep(grinf[x,],sum(idx==x)),ncol=2,byrow=TRUE))})
        }
    }

    spatstat.options(checkpolygons=FALSE)
    W <- as(gUnaryUnion(spdf),"owin")
    spatstat.options(checkpolygons=TRUE)

    sd <- suppressWarnings(ppp(x=unlist(spts[,1]),y=unlist(spts[,2]),window=W)) # NOTE THIS IS REFINED BELOW                      					    
    
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
	
	## OBTAIN POISSON OFFSET ON LATTICE (LINEAR INTERPOLATION) ##
	
	spatial.offset <- fftinterpolate(spatial,mcens,ncens,ext=ext)
	spatial.offset <- spatial.offset*cellInside
	## no longer required as here this is an offset: spatialvals <- spatialvals / (cellarea*sum(spatialvals))
	
	cts <- rep(0,N*M)
	soffvec <- as.vector(spatial.offset[1:M,1:N])
	for(i in 1:length(spdf)){
        if(spdf$X[i]==0){
            next
        }
        else{
            cellinfo <- ol$info[ol$info$polyidx==i,]
            if(nrow(cellinfo)>1){
                idx <- sample(cellinfo$grididx,spdf$X[i],replace=TRUE,prob=cellinfo$area*soffvec[cellinfo$grididx])
            }
            else{
                idx <- rep(cellinfo$grididx,spdf$X[i])
            }
            idxunq <- unique(idx)
            sapply(idxunq,function(x){cts[x] <<- cts[x] + sum(idx==x)})
        }
    }
    
    

    nis <- matrix(0,Mext,Next)
    nis[1:M,1:N] <- cts    
    	
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
    
    lg <- MALAlgcpAggregateSpatial.PlusPars( mcmcloop=mLoop,
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
                                    d=d,
                                    spdf=spdf,
                                    ol=ol,
                                    Nfreq=Nfreq)
                                                  
	
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
	
	class(lg) <- c("lgcpPredictAggregateSpatialPlusParameters","lgcpPredict","lgcpobject")	
	
	return(lg)															
					
}	
	

##' MALAlgcpAggregateSpatial.PlusPars function
##'
##' A function to run the MCMC algorithm for aggregated spatial point process data. Not for general purpose use.
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
##' @param spdf the SpatialPolygonsDataFrame containing the aggregate counts as a variable X 
##' @param ol overlay of fft grid onto spdf 
##' @param Nfreq frequency at which to resample nis
##' @return output from the MCMC run
##' @export

MALAlgcpAggregateSpatial.PlusPars <- function(   mcmcloop,
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
                                        d,
                                        spdf,
                                        ol,
                                        Nfreq){
                                        
    
    
    #browser()                                           
                            
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
    
    #h <- initialiseAMCMC(adaptivescheme)
    #hrec <- h
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
    
    if(length(betaval)!=ncol(Z)){
        stop(paste("Check prior for beta, the number of parameters should be:",ncol(Z)))
    }    
    
    if(length(betaval)!=ncol(Zmat)){
        stop("Prior for beta has incorrect number of variables.")    
    }  
    
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
    eMat <- spatialvals*exp(matrix(as.vector(Z%*%betaval),nrow(GP$gamma),ncol(GP$gamma)))*GP$expY*cellarea
    fi <- (1/length(GP$Y))*Re(fft(GP$invrootQeigs*fft((1/length(GP$Y))*Re(fft(GP$invrootQeigs*fft(eMat,inverse=TRUE))),inverse=TRUE)))
    gammaVar <- (1.65^2/((Mext*Next)^(1/3)))*(-1)/(-1-fi)
    GAMMAVAR <- gammaVar
    rootgammaVar <- sqrt(gammaVar)
    # compute variance matrix of beta
    EZ <- as.vector(cellarea*spatialvals*exp(matrix(as.vector(Z%*%betaval),nrow(GP$gamma),ncol(GP$gamma))+GP$Y))    
    CB <- matrix(NA,nbeta,nbeta)    

    for(i in 1:nbeta){
        for(j in 1:nbeta){
            CB[i,j] <- sum(Z[,i]*Z[,j]*EZ)  + model.priors$betaprior$precision[i,j] # observed information
        }
    }   
    
    
    etaCovMat <- try((0.234/0.574)*solve((-1)*GPdrv2(GP=GP,prior=model.priors,Z=Z,Zt=Zt,eta=etaval,beta=betaval,nis=nis,cellarea=cellarea,spatial=spatialvals,gradtrunc=Inf,fftgrid=fftgrid,covfunction=spatial.covmodel,d=d,eps=1e-6)$hess)) #mean(diff(mcens[1:2]),diff(ncens[1:2]))/100)$hess) # eps=mean(diff(mcens[1:2]),diff(ncens[1:2]))/100 is approx 1/100 cellwidth
    etaCovMattest <- FALSE    
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
        try(etaCovMat <- (0.234/0.574)*solve((-1)*matrix(c(2*coef[2],coef[6],coef[6],2*coef[4]),2,2)))
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
    print(etaCovMat)
	

	sigma_eta <- (2.38^2/length(etaval))*etaCovMat 
	SIGMA_ETA <- sigma_eta
	Q_eta <- solve(sigma_eta)
	sigma_beta <- (1.65^2/(length(betaval)^(1/3)))*solve(CB) #vcov(mod)	
	Q_beta <- solve(sigma_beta)
	
	sigma_eta_chol <- t(chol(sigma_eta))
	sigma_beta_chol <- t(chol(sigma_beta))
	
	betarec <- c()
	etarec <- c()
	
	sigmaetarec <- list(sigma_eta)
	sigmabetarec <- list(sigma_beta)
    
    adapt_h <- adaptivescheme
    h <- initialiseAMCMC(adapt_h)          
	
    GP <- GPrealisation(gamma=matrix(0,Mext,Next),fftgrid=fftgrid,covFunction=spatial.covmodel,covParameters=cp,d=d)    
    oldtags <- target.and.grad.AggregateSpatialPlusPars(GP=GP,prior=model.priors,Z=Z,Zt=Zt,eta=etaval,beta=betaval,nis=nis,cellarea=fftgrid$cellarea,spatial=spatialvals,gradtrunc=gradtrunc)

    trigger1 <- TRUE
    
    nsav <- (mcmcloop$N-mcmcloop$burnin)/mcmcloop$thin

    nisrec <- array(dim=c(M,N,nsav))
    niscount <- 0
    
    tarrec <- oldtags$logtarget
    hallrec <- h
    
    reject_its <- c()   
    
    while(nextStep(mcmcloop)){
    
        if(iteration(mcmcloop)%%Nfreq==1){
            oldtags <- target.and.grad.AggregateSpatialPlusPars(GP=GP,prior=model.priors,Z=Z,Zt=Zt,eta=etaval,beta=betaval,nis=nis,cellarea=fftgrid$cellarea,spatial=spatialvals,gradtrunc=gradtrunc)
        }
    
        propmeans_gamma <- GP$gamma + (h/2)*gammaVar*oldtags$gradgamma
        propmeans_eta <- etaval 
        propeta <- as.vector(propmeans_eta + sqrt(h)*sigma_eta_chol%*%rnorm(neta)) 
        
        propetainv <- etainvtrans(propeta) 
        propcp <- CovParameters(list(sigma=propetainv[1],phi=propetainv[2]))
        propGP <- GPrealisation(gamma=propmeans_gamma+ sqrt(h)*rootgammaVar*rnorm(Mext*Next),fftgrid=fftgrid,covFunction=spatial.covmodel,covParameters=propcp,d=d)
                
        propmeans_beta <- betaval + (h/2)*sigma_beta%*%oldtags$gradbeta     
        propbeta <- BetaParameters(as.vector(propmeans_beta+sqrt(h)*sigma_beta_chol%*%rnorm(nbeta)))         
          
        proptags <- target.and.grad.AggregateSpatialPlusPars(GP=propGP,prior=model.priors,Z=Z,Zt=Zt,eta=propeta,beta=propbeta,nis=nis,cellarea=fftgrid$cellarea,spatial=spatialvals,gradtrunc=gradtrunc)    
        
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
        
		if (proptags$logtarget==-Inf | is.na(ac) | is.nan(ac)){ 

	        warning("One possible cause of this warning is that there may be evidence in the data for quite large values of the spatial correlation parameter. If this is the case, then this warning means that the MCMC chain has wandered into a region of the phi-space that causes the variance matrix of Y (computed by the DFT) to become non positive-definite. One possible way of rectifying this issue is to restart the chain using a larger value of 'ext' in the call to lgcpPredictAggregateSpatialPlusPars. You should do this if the warning message is repeated many times. If this warning message appears at all then you should be warned that inference from this run may not be reliable: the proposed move at this iteration was rejected.",immediate.=TRUE)
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
        
        
        if(iteration(mcmcloop)%%Nfreq==0){ # update nis
            
            cts <- rep(0,N*M)
        	e <-  as.vector(oldtags$e[1:M,1:N])
        	for(i in 1:length(spdf)){
                if(spdf$X[i]==0){
                    next
                }
                else{
                    cellinfo <- ol$info[ol$info$polyidx==i,]
                    if(nrow(cellinfo)>1){
                        idx <- sample(cellinfo$grididx,spdf$X[i],replace=TRUE,prob=cellinfo$area*e[cellinfo$grididx])
                    }
                    else{
                        idx <- rep(cellinfo$grididx,spdf$X[i])
                    }
                    idxunq <- unique(idx)
                    sapply(idxunq,function(x){cts[x] <<- cts[x] + sum(idx==x)})
                }
            }
        
            nis <- matrix(0,Mext,Next)
            nis[1:M,1:N] <- cts 
        }

        
        
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
			
			niscount <- niscount + 1
            nisrec[,,niscount] <- nis[1:M,1:N]		    
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
    retlist$nisrec <- nisrec
    retlist$reject_its <- reject_its
    
	
	return(retlist)                          
}                        

##' target.and.grad.AggregateSpatialPlusPars function
##'
##' A function to compute the target and gradient for the Bayesian aggregated point process model. Not for general use. 
##'
##' @param GP an object constructed using GPrealisation
##' @param prior the prior, created using lgcpPrior
##' @param Z the design matrix on the full FFT grid
##' @param Zt the transpose of the design matrix
##' @param eta the model parameter, eta
##' @param beta the model parameters, beta 
##' @param nis cell counts on the FFT grid
##' @param cellarea the cell area 
##' @param spatial the poisson offset
##' @param gradtrunc the gradient truncation parameter
##' @return the target and gradient
##' @export

target.and.grad.AggregateSpatialPlusPars <- function(GP,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc){

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


