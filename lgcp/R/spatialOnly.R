##' blockcircbase function
##'
##' Compute the base matrix of a continuous Gaussian field. Computed as a block circulant matrix on a torus where x and y is the 
##' x and y centroids (must be equally spaced)
##'
##' @param x x centroids, an equally spaced vector
##' @param y y centroids, an equally spaced vector
##' @param sigma spatial variance parameter
##' @param phi spatial decay parameter
##' @param model covariance model, see ?CovarianceFct
##' @param additionalparameters additional parameters for chosen covariance model. See ?CovarianceFct
##' @param inverse logical. Whether to return the base matrix of the inverse covariance matrix (ie the base matrix for the precision matrix), default is FALSE
##' @return the base matrix of a block circulant matrix representing a stationary covariance function on a toral grid.
##' @export
blockcircbase <- function(x,y,sigma,phi,model,additionalparameters,inverse=FALSE){
    M <- length(x)
    N <- length(y)
    xidx <- rep(1:M,N)
    yidx <- rep(1:N,each=M)
    dxidx <- pmin(abs(xidx-xidx[1]),M-abs(xidx-xidx[1]))
    dyidx <- pmin(abs(yidx-yidx[1]),N-abs(yidx-yidx[1]))
    d <- sqrt(((x[2]-x[1])*dxidx)^2+((y[2]-y[1])*dyidx)^2)
    covbase <- matrix(gu(d,sigma=sigma,phi=phi,model=model,additionalparameters=additionalparameters),M,N)
    if(!inverse){ 
        return(covbase)
    }
    else{
        return(inversebase(covbase))
    }
}


##' inversebase function
##'
##' A function to compute the base of the inverse os a block circulant matrix, given the base of the matrix 
##'
##' @param x the base matrix of a block circulant matrix
##' @return the base matrix of the inverse of the circulant matrix
##' @export
inversebase <- function(x){    
    return(Re((1/prod(dim(x)))*fft(1/fft(x),inverse=TRUE)))
}

##' eigenfrombase function
##'
##' A function to compute the eigenvalues of an SPD block circulant matrix given the base matrix. 
##'
##' @param x the base matrix
##' @return the eigenvalues
##' @export
eigenfrombase <- function(x){
    return(Re(fft(x))) # do not need scaling as in Rue and Held (2005) pp64
}


##' is.SPD function
##'
##' A function to compute whether a block circulant matrix is symmetric positive definite (SPD), given its base matrix.
##'
##' @param base base matrix of a block circulant matrix
##' @return logical, whether the circulant matrix the base represents is SPD
##' @export
is.SPD <- function(base){
    if(prod(sign(eigenfrombase(base)))>0){
        return(TRUE)
    }
    return(FALSE)
}

##' GammafromY function
##'
##' A function to change Ys (spatially correlated noise) into Gammas (white noise). Used in the MALA algorithm.
##'
##' @param Y Y matrix
##' @param rootQeigs square root of the eigenvectors of the precision matrix
##' @param mu parameter of the latent Gaussian field
##' @return Gamma
##' @export
GammafromY <- function(Y,rootQeigs,mu){
    nc <- dim(rootQeigs)[2]
    nb <- length(Y)
    return((1/nb)*Re(fft(fft(Y-mu)*rootQeigs,inverse=TRUE)))    
}


##' YfromGamma function
##'
##' A function to change Gammas (white noise) into Ys (spatially correlated noise). Used in the MALA algorithm. 
##'
##' @param Gamma Gamma matrix
##' @param invrootQeigs inverse square root of the eigenvectors of the precision matrix
##' @param mu parameter of the latent Gaussian field
##' @return Y
##' @export
YfromGamma <- function(Gamma,invrootQeigs,mu){
    nc <- dim(invrootQeigs)[2]
    nb <- length(Gamma)
    return(mu + (1/nb)*Re(fft(invrootQeigs*fft(Gamma,inverse=TRUE))))
}

##' target.and.grad.spatial function
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
##' @param scaleconst the expected number of cases
##' @param gradtrunc gradient truncation parameter
##' @return the back-transformed Y, its exponential, the log-target and gradient for use in MALAlgcpSpatial 
##' @export
target.and.grad.spatial <- function(Gamma,nis,cellarea,rootQeigs,invrootQeigs,mu,spatial,logspat,scaleconst,gradtrunc){
    
    Y <- YfromGamma(Gamma=Gamma,invrootQeigs=invrootQeigs,mu=mu)
    expY <- exp(Y)    
    
    ###
    ### note logspat not necessary ...logtarget <- -(1/2)*sum(Gamma^2) + sum((Y+logspat)*nis - scaleconst*spatial*expY*cellarea) # note that both nis=0, logspat=0 and spatial=0 outside of observation window, so this effectively limits summation to the observation window only 
    ###    
     
    logtarget <- -(1/2)*sum(Gamma^2) + sum(Y*nis - scaleconst*spatial*expY*cellarea) # note that both nis=0, logspat=0 and spatial=0 outside of observation window, so this effectively limits summation to the observation window only

    expYtrunc <- expY
    expYtrunc[expYtrunc>gradtrunc] <- gradtrunc
    grad <- (-1)*Gamma + (1/length(Y))*Re(fft(invrootQeigs*fft(nis-scaleconst*spatial*expYtrunc*cellarea,inverse=TRUE)))         
    
    return(list(Y=Y,expY=expY,logtarget=logtarget,grad=grad))
}


##' circulant function
##'
##' generic function for constructing circulant matrices
##'
##' @param x an object
##' @param ... additional arguments
##' @return method circulant
##' @export
circulant <- function(x,...){
    UseMethod("circulant")
}



##' circulant.numeric function
##'
##' returns a circulant matrix with base x
##'
##' @method circulant numeric
##' @param x an numeric object
##' @param ... additional arguments
##' @return a circulant matrix with base x
##' @export
circulant.numeric <- function(x,...){
    n <- length(x)
    M <- matrix(NA,n,n)
    idx <- 1:n
    for(i in 1:n){
        M[i,] <- x[idx]
        idx <- c(idx[n],idx[-n])
    }
    return(M)
}



##' circulant.matrix function
##'
##' If x is a matrix whose columns are the bases of the sub-blocks of a block circulant matrix, then this function returns the 
##' block circulant matrix of interest.
##'
##' @method circulant matrix
##' @param x a matrix object
##' @param ... additional arguments
##' @return If x is a matrix whose columns are the bases of the sub-blocks of a block circulant matrix, then this function returns the block circulant matrix of interest.
##' @export
circulant.matrix <- function(x,...){
    M <- dim(x)[1]
    N <- dim(x)[2]
    submats <- t(apply(x,2,circulant))
    mat <- matrix(NA,M*N,M*N)
    idx <- 1:N
    ct <- 1
    for (i in 1:N){
        for (j in 1:N){ 
            xstart <- M*floor((ct-1)/N) + 1
            mult <- ct%%N
            if (mult==0){
                mult <- N
            }
            ystart <- M*(mult-1) + 1
            mat[xstart:(xstart+M-1),ystart:(ystart+M-1)] <- submats[idx[j],]
            ct <- ct + 1
        }
        idx <- c(idx[N],idx[-N])
    }
    return(mat)
}

##' toral.cov.mat function
##'
##' A function to compute the covariance matrix of a stationary process on a torus.
##'
##' @param xg x grid
##' @param yg y grid
##' @param sigma spatial variability parameter
##' @param phi spatial decay parameter
##' @param model model for covariance, see ?CovarianceFct
##' @param additionalparameters additional parameters for covariance structure
##' @return circulant covariacne matrix
##' @export
toral.cov.mat <- function(xg,yg,sigma,phi,model,additionalparameters){
    bcb <- blockcircbase(x=xg,y=yg,sigma=sigma,phi=phi,model=model,additionalparameters=additionalparameters)
    return(circulant(bcb))
}



##' lgcpSimSpatial function
##'
##' A function to simulate from a log gaussian process
##'
##' @param owin observation window
##' @param spatial.intensity an object that can be coerced to one of class spatialAtRisk
##' @param expectednumcases the expected number of cases
##' @param cellwidth width of cells in same units as observation window
##' @param model.parameters parameters of model, see ?lgcppars. Only set sigma and phi for spatial model.
##' @param spatial.covmodel spatial covariance function, default is exponential, see ?CovarianceFct
##' @param covpars vector of additional parameters for spatial covariance function, in order they appear in chosen model in ?CovarianceFct
##' @param ext how much to extend the parameter space by. Default is 2.
##' @param plot logical, whether to plot the latent field.
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window.
##' @return a ppp object containing the data
##' @export
lgcpSimSpatial <- function( owin=NULL,
                            spatial.intensity=NULL,
                            expectednumcases=100,
                            cellwidth = 0.05,
                            model.parameters=lgcppars(sigma=2,phi=0.2),
                            spatial.covmodel="exponential",
                            covpars=c(),
                            ext=2,
                            plot=FALSE,
                            inclusion="touching"){
                                                  
    sigma <- model.parameters$sigma
	phi <- model.parameters$phi
	mu <- model.parameters$mu
	theta <- model.parameters$theta

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
	
	s <- Sys.time()
    truefield <- YfromGamma(matrix(rnorm(Mext*Next),Mext,Next),invrootQeigs=irqe,mu=mu)[1:M,1:N]
    rate <- as.vector(expectednumcases*spatialvals*cellarea*exp(truefield))
    cmat <- matrix(rpois(M*N,rate),M,N)
    ncases <- sum(cmat)
    if(ncases==0){
	    stop("No data generated for chosen parameters")
	}
    caseidx <- which(cmat>0)
    caseidx <- unlist(sapply(caseidx,function(x){rep(x,cmat[x])}))
    cases <- cbind(rep(xg,length(yg)),rep(yg,each=length(xg)))[caseidx,] + cbind(runif(ncases,-del1/2,del1/2),runif(ncases,-del2/2,del2/2))
    if(plot){
        rate[rate==0] <- NA
        image(xg,yg,matrix(rate,M,N)^0.25)
        points(cases,pch="+",cex=0.5)
    }    
	
	xyt <- ppp(x=cases[,1],y=cases[,2],window=owin)
	attr(xyt,"rejects") <- NULL # get rid of rejects: these are due to discrete approximation
	attr(xyt,"spatialatrisk") <- spatial
	attr(xyt,"expectednumcases") <- expectednumcases
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
    attr(xyt,"inclusion") <- inclusion  
    return(xyt)
}

##' rgauss function
##'
##' A function to simulate a Gaussian field on a regular square lattice, the returned object is of class lgcpgrid.
##'
##' @param n the number of realisations to generate. Default is 1.
##' @param range a vector of length 2, defining the left-most and right most cell centroids in the x-direction. Note that the centroids in the y-direction are the same as those in the x-direction.
##' @param ncells the number of cells, typially a power of 2
##' @param spatial.covmodel spatial covariance function, default is exponential, see ?CovarianceFct 
##' @param model.parameters parameters of model, see ?lgcppars. Only set sigma and phi for spatial model.
##' @param covpars vector of additional parameters for spatial covariance function, in order they appear in chosen model in ?CovarianceFct
##' @param ext how much to extend the parameter space by. Default is 2.
##' @return an lgcp grid object containing the simulated field(s).
##' @export

rgauss <- function(n=1,range=c(0,1),ncells=128,spatial.covmodel="exponential",model.parameters=lgcppars(sigma=2,phi=0.1),covpars=c(),ext=2){    
    
    sigma <- model.parameters$sigma
	phi <- model.parameters$phi
	mu <- model.parameters$mu    
    
    mcens <- seq(range[1],range[1]+ext*diff(range)+(ext-1)*diff(range)/(ncells-1),length.out=ext*ncells)
    ncens <- seq(range[1],range[1]+ext*diff(range)+(ext-1)*diff(range)/(ncells-1),length.out=ext*ncells)        
    
    bcb <- blockcircbase(x=mcens,y=ncens,sigma=sigma,phi=phi,model=spatial.covmodel,additionalparameters=covpars)    
    Qeigs <- eigenfrombase(inversebase(bcb)) # eigenvalues of Q (the precision matrix)
    rqe <- sqrt(Qeigs) # square root of the eigenvalues (used in computation)
    irqe <- 1/rqe # reciprocal root (commputation)
    slist <- list()
    for(i in 1:n){
        slist[[i]] <- YfromGamma(matrix(rnorm((ncells*ext)^2),ncells*ext,ncells*ext),invrootQeigs=irqe,mu=mu)[1:ncells,1:ncells]
    }    
    return(lgcpgrid(slist,xvals=mcens[1:ncells],yvals=ncens[1:ncells]))        
}

