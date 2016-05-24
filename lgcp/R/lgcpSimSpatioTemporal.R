##' lgcpSim function
##'
##' Approximate simulation from a spatiotemoporal log-Gaussian Cox Process. Returns an stppp object.
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
##' @param owin polygonal observation window
##' @param tlim time interval on which to simulate data
##' @param spatial.intensity object that can be coerced into a spatialAtRisk object. if NULL then uniform spatial is chosen 
##' @param temporal.intensity the fixed temporal component: either a numeric vector, or a function that can be coerced into an object of class temporalAtRisk
##' @param cellwidth width of cells  in same units as observation window
##' @param model.parameters parameters of model, see ?lgcppars. 
##' @param spatial.covmodel spatial covariance function, default is exponential, see ?CovarianceFct
##' @param covpars vector of additional parameters for spatial covariance function, in order they appear in chosen model in ?CovarianceFct
##' @param returnintensities logigal, whether to return the spatial intensities and true field Y at each time. Default FALSE.
##' @param progressbar logical, whether to print a progress bar. Default TRUE.
##' @param ext how much to extend the parameter space by. Default is 2.
##' @param plot logical, whether to plot intensities.
##' @param ratepow power that intensity is raised to for plotting purposes (makes the plot more pleasign to the eye), defaul 0.25
##' @param sleeptime time in seconds to sleep between plots
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window.
##' @return an stppp object containing the data
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor, Tilman M. Davies, Barry S. Rowlingson, Peter J. Diggle (2013). Journal of Statistical Software, 52(4), 1-40. URL http://www.jstatsoft.org/v52/i04/
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##'     \item Wood ATA, Chan G (1994). Simulation of Stationary Gaussian Processes in [0,1]d. Journal of Computational and Graphical Statistics, 3(4), 409-432.
##'     \item Moller J, Syversveen AR, Waagepetersen RP (1998). Log Gaussian Cox Processes. Scandinavian Journal of Statistics, 25(3), 451-482.
##' }
##' @seealso \link{lgcpPredict}, \link{showGrid.stppp}, \link{stppp}
##' @examples xyt <- lgcpSim()
##' @export

lgcpSim <- function(owin=NULL,
                    tlim=as.integer(c(0,10)),
                    spatial.intensity=NULL,
                    temporal.intensity=NULL,
                    cellwidth = 0.05,
                    model.parameters=lgcppars(sigma=2,phi=0.2,theta=1),
                    spatial.covmodel="exponential",
                    covpars=c(),
                    returnintensities=FALSE,
                    progressbar=TRUE,
                    ext=2,
                    plot=FALSE,
                    ratepow=0.25,
                    sleeptime=0,
                    inclusion="touching"){
                                    
    if (!inherits(tlim,"integer")){
	    warning("Converting tlim into integer values, see ?as.integer")
	    tlim <- as.integer(tlim) # convert times into integer values: they should already be in this form.
	}
	tlim <- sort(tlim)
	if (tlim[1]==tlim[2]){
	    stop("Length of time interval given by as.integer(tlim) must be >= 1")
	}
	toffset <- tlim[1]
	maxt <- tlim[2] - toffset
                                                  
    sigma <- model.parameters$sigma
	phi <- model.parameters$phi
	mu <- model.parameters$mu
	theta <- model.parameters$theta

    if(is.null(owin)){
        owin <- owin()
    }
    
    if (is.null(temporal.intensity)){
        temporal.intensity <- constantInTime(100,tlim)
    }
    else{
        if (!inherits(temporal.intensity,"temporalAtRisk")){
            temporal.intensity <- temporalAtRisk(temporal.intensity,tlim)
        }
        if(!all(tlim==attr(temporal.intensity,"tlim"))){
	        stop("Incompatible temporal.intensity, integer time limits (tlim and temporal.intensity$tlim) do not match")
	    }
    } 
    
    ndivs <- diff(tlim) # just choose "integer" time divisions: gives greater consistency between data used in simulation and inference (compared with the above)
    tdiff = maxt/ndivs
    times <- tdiff/2 + tdiff*(0:(ndivs-1)) 
    mut <- sapply(times+toffset,temporal.intensity)
    
    if (progressbar){
        pb <- txtProgressBar(min=1,max=ndivs,style=3)
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

	if(is.null(spatial.intensity)){
        spatial <- spatialAtRisk(list(X=seq(xyt$window$xrange[1],xyt$window$xrange[2],length.out=M),Y=seq(xyt$window$yrange[1],xyt$window$yrange[2],length.out=N),Zm=matrix(1/(M*N),M,N)))
    }
    else{
        if(!any(class(spatial.intensity)=="spatialAtRisk")){		
            spatial <- spatialAtRisk(spatial.intensity)
        }
        else{
            spatial <- spatial.intensity
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
	
	xg <- mcens[1:M]
	yg <- ncens[1:N]
	
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
	cellInside <- as.numeric(matrix(as.logical(cellInside),Mext,Next)[1:M,1:N])
	
	## OBTAIN SPATIAL VALS ON LATTICE (LINEAR INTERPOLATION) ##
	
	spatialvals <- fftinterpolate(spatial,mcens,ncens,ext=ext)
	spatialvals <- spatialvals[1:M,1:N]
	spatialvals <- spatialvals*cellInside
	spatialvals <- spatialvals / (cellarea*sum(spatialvals))
	
	# compute the base matrix of the covariance matrix
    bcb <- blockcircbase(x=mcens,y=ncens,sigma=sigma,phi=phi,model=spatial.covmodel,additionalparameters=covpars)
    
    Qeigs <- eigenfrombase(inversebase(bcb)) # eigenvalues of Q (the precision matrix)
    rqe <- sqrt(Qeigs) # square root of the eigenvalues (used in computation)
    irqe <- 1/rqe # reciprocal root (commputation)    
    	
	################################################################
    # Simulate the data	
	################################################################

    if(returnintensities){
	    intensities <- array(NA,c(M,N,ndivs))
	    truefield <- array(NA,c(M,N,ndivs))
	}
	else{
	    intensities <- NULL
	    truefield <- NULL
	}
	
    
    
    cases <- NULL
    t <- NULL
    Y <- YfromGamma(matrix(rnorm(Mext*Next),Mext,Next),invrootQeigs=irqe,mu=mu)[1:M,1:N]
    rate <- as.vector(mut[1]*spatialvals*cellarea*exp(Y))
    if(returnintensities){
	    intensities[,,1] <- rate
	    truefield[,,1] <- Y
	}
    cmat <- matrix(rpois(M*N,rate),M,N)
    ncases <- sum(cmat)
    if(ncases>0){
        caseidx <- which(cmat>0)
        caseidx <- unlist(sapply(caseidx,function(x){rep(x,cmat[x])}))
        cases <- cbind(rep(xg,length(yg)),rep(yg,each=length(xg)))[caseidx,] + cbind(runif(ncases,-del1/2,del1/2),runif(ncases,-del2/2,del2/2))
        t <- sort(runif(ncases,times[1]-tdiff/2,times[1]+tdiff/2))
    }
    if(plot){
        rate[rate==0] <- NA
        image.plot(xg,yg,matrix(rate,M,N)^ratepow)
        points(cases,pch="+",cex=0.5)
        Sys.sleep(sleeptime)
    }
    for(i in 2:ndivs){
        Y <- mu*(1-exp(-theta)) + exp(-theta)*Y + sqrt(1-exp(-2*theta))*YfromGamma(matrix(rnorm(Mext*Next),Mext,Next),invrootQeigs=irqe,mu=0)[1:M,1:N] # note delta t is integral and = 1 here
        rate <- as.vector(mut[i]*spatialvals*cellarea*exp(Y))
        cmat <- matrix(rpois(M*N,rate),M,N)
        ncases <- sum(cmat)
        if(ncases>0){
            caseidx <- which(cmat>0)
            caseidx <- unlist(sapply(caseidx,function(x){rep(x,cmat[x])}))
            newcases <- cbind(rep(xg,length(yg)),rep(yg,each=length(xg)))[caseidx,] + cbind(runif(ncases,-del1/2,del1/2),runif(ncases,-del2/2,del2/2))
            cases <- rbind(cases,newcases)
            t <- c(t,sort(runif(ncases,times[i]-tdiff/2,times[i]+tdiff/2)))
            if(plot){
                rate[rate==0] <- NA
                image.plot(xg,yg,matrix(rate,M,N)^ratepow)
                points(newcases,pch="+",cex=0.5)
                Sys.sleep(sleeptime)
            }    
        }
        
        if(returnintensities){
    	    intensities[,,i] <- rate
    	    truefield[,,i] <- Y
    	}
        if (progressbar){
            setTxtProgressBar(pb,i)
        }
    }
    if (progressbar){
	    close(pb)
	}
	if(is.null(t)){
	    stop("No data generated for chosen parameters")
	}    
	
	if(!all(inside.owin(cases[,1],cases[,2],owin))){
	    remidx <- which(!inside.owin(cases[,1],cases[,2],owin))
	    cases <- cases[-remidx,]
	    t <- t[-remidx]
	}
	
	xyt <- stppp(ppp(x=cases[,1],y=cases[,2],window=owin),t=(t+toffset),tlim=tlim)

	attr(xyt,"rejects") <- NULL # get rid of rejects: these are due to discrete approximation
	attr(xyt,"spatialatrisk") <- spatial
	attr(xyt,"temporalfitted") <- mut
	attr(xyt,"cellwidth") <- cellwidth
	attr(xyt,"sigma") <- sigma
	attr(xyt,"phi") <- phi
	attr(xyt,"theta") <- theta
    attr(xyt,"temporalintensity") <- temporal.intensity
    attr(xyt,"temporalfitted") <- mut
    attr(xyt,"spatialcovmodel") <- spatial.covmodel
    attr(xyt,"covpars") <- covpars
	attr(xyt,"ext") <- ext
    attr(xyt,"xvals") <- xg
    attr(xyt,"yvals") <- yg
    attr(xyt,"intensities") <- intensities
    attr(xyt,"truefield") <- truefield 
    attr(xyt,"inclusion") <- inclusion  
    return(xyt)
}

