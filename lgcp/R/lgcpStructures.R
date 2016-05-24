##' genFFTgrid function
##'
##' A function to generate an FFT grid and associated quantities including cell dimensions,
##' size of extended grid, centroids, cell area, cellInside matrix (a 0/1 matrix: is the centroid of the cell inside the observation window?)
##'
##' @param study.region an owin object
##' @param M number of cells in x direction
##' @param N number of cells in y direction
##' @param ext multiplying constant: the size of the extended grid: ext*M by ext*N
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window.
##' @return a list
##' @export

genFFTgrid <- function(study.region,M,N,ext,inclusion="touching"){
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
    
    obj <- list(del1=del1,del2=del2,Mext=Mext,Next=Next,mcens=mcens,ncens=ncens,cellarea=cellarea,cellInside=cellInside,inclusion=inclusion)    
    class(obj) <- "FFTgrid"
    
    return(obj)   
}

##' blockcircbaseFunction function
##'
##' Compute the base matrix of a continuous Gaussian field. Computed as a block circulant matrix on a torus where x and y is the 
##' x and y centroids (must be equally spaced). This is an extension of the function blockcircbase to extend the range of covariance functions
##' that can be fitted to the model.
##'
##' @param x x centroids, an equally spaced vector
##' @param y y centroids, an equally spaced vector
##' @param CovFunction a function of distance, returning the covariance between points that distance apart
##' @param CovParameters an object of class CovParamters, see ?CovParameters 
##' @param inverse logical. Whether to return the base matrix of the inverse covariance matrix (ie the base matrix for the precision matrix), default is FALSE
##' @return the base matrix of a block circulant matrix representing a stationary covariance function on a toral grid.
##' @seealso \link{minimum.contrast}, \link{minimum.contrast.spatiotemporal}, \link{chooseCellwidth}, \link{getpolyol}, \link{guessinterp}, \link{getZmat},
##' \link{addTemporalCovariates}, \link{lgcpPrior}, \link{lgcpInits},
##' \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}, \link{lgcpPredictSpatioTemporalPlusPars}, 
##' \link{lgcpPredictMultitypeSpatialPlusPars}
##' @export
blockcircbaseFunction <- function(x,y,CovFunction,CovParameters,inverse=FALSE){
    M <- length(x)
    N <- length(y)
    xidx <- rep(1:M,N)
    yidx <- rep(1:N,each=M)
    dxidx <- pmin(abs(xidx-xidx[1]),M-abs(xidx-xidx[1]))
    dyidx <- pmin(abs(yidx-yidx[1]),N-abs(yidx-yidx[1]))
    d <- sqrt(((x[2]-x[1])*dxidx)^2+((y[2]-y[1])*dyidx)^2)
    covbase <- matrix(CovFunction(d,CovParameters=CovParameters),M,N)
    if(!inverse){ 
        return(covbase)
    }
    else{
        return(inversebase(covbase))
    }
}


##' GPrealisation function
##'
##' A function to store a realisation of a spatial gaussian process for use in MCMC algorithms that include Bayesian parameter estimation.
##' Stores not only the realisation, but also computational quantities.
##'
##' @param gamma the transformed (white noise) realisation of the process
##' @param fftgrid an object of class FFTgrid, see ?genFFTgrid
##' @param covFunction an object of class function returning the spatial covariance
##' @param covParameters an object of class CovParamaters, see ?CovParamaters
##' @param d matrix of grid distances
##' @return a realisation of a spatial Gaussian process on a regular grid
##' @export
GPrealisation <- function(gamma,fftgrid,covFunction,covParameters,d){
    if(!inherits(gamma,"matrix")){
        stop("argument 'gamma' must be a matrix")
    }
    if(class(fftgrid)!="FFTgrid"){
        stop("argument 'fftgrid' must be an object of class FFTgrid, see ?genFFTgrid")
    }
    if(!any(class(covFunction)=="CovFunction")){
        stop("argument 'covFunction' must be an object of class CovFunction")
    }
    if(class(covParameters)!="CovParamaters"){
        stop("argument 'covParamaters' must be an object of class CovParamaters, see ?CovParamaters")
    }
    
    Mext <- length(fftgrid$mcens)
    Next <- length(fftgrid$ncens)
    covbase <- covFunction(d=d,CovParameters=covParameters)
    covbase <- lapply(covbase,matrix,nrow=Mext,ncol=Next)
    ###############

    gp <- list()
    gp$gamma <- gamma
    gp$CovFunction <- covFunction 
    gp$CovParameters <- covParameters
    gp$covbase <- covbase 
    #gp$fftcovbase <- lapply(covbase,function(x){Re(fft(x))})
    gp$rootQeigs <- sqrt(1/Re(fft(covbase[[1]])))#sqrt(1/gp$fftcovbase$eval)
    gp$invrootQeigs <- 1/gp$rootQeigs
    gp$Y <- YfromGamma(gamma,invrootQeigs=gp$invrootQeigs,mu=covParameters$mu)
    gp$expY <- exp(gp$Y) # for computational purposes
    gp$mcens <- fftgrid$mcens
    gp$ncens <- fftgrid$ncens
    
    class(gp) <- "GPrealisation"    
    return(gp)    
    
}

##' stGPrealisation function
##'
##' A function to store a realisation of a spatiotemporal gaussian process for use in MCMC algorithms that include Bayesian parameter estimation.
##' Stores not only the realisation, but also computational quantities.
##'
##' @param gamma the transformed (white noise) realisation of the process
##' @param fftgrid an object of class FFTgrid, see ?genFFTgrid
##' @param covFunction an object of class function returning the spatial covariance
##' @param covParameters an object of class CovParamaters, see ?CovParamaters
##' @param d matrix of grid distances
##' @param tdiff vector of time differences 
##' @return a realisation of a spatiotemporal Gaussian process on a regular grid
##' @export



stGPrealisation <- function(gamma,fftgrid,covFunction,covParameters,d,tdiff){
    if(class(gamma)!="list"){
        stop("argument 'gamma' must be a list")
    }
    if(class(fftgrid)!="FFTgrid"){
        stop("argument 'fftgrid' must be an object of class FFTgrid, see ?genFFTgrid")
    }
    if(!any(class(covFunction)=="CovFunction")){
        stop("argument 'covFunction' must be an object of class CovFunction")
    }
    if(class(covParameters)!="CovParamaters"){
        stop("argument 'covParamaters' must be an object of class CovParamaters, see ?CovParamaters")
    }
    
    Mext <- length(fftgrid$mcens)
    Next <- length(fftgrid$ncens)
    covbase <- covFunction(d=d,CovParameters=covParameters)
    covbase <- lapply(covbase,matrix,nrow=Mext,ncol=Next)
    ###############
    

    gp <- list()
    gp$gamma <- gamma
    gp$CovFunction <- covFunction 
    gp$CovParameters <- covParameters
    gp$covbase <- covbase 
    gp$rootQeigs <- sqrt(1/Re(fft(covbase[[1]])))
    gp$invrootQeigs <- 1/gp$rootQeigs
    gp$at <- covParameters$mu*(1-exp(-covParameters$theta*tdiff))
    gp$bt <- exp(-covParameters$theta*tdiff)
    gp$gt <- 1-exp(-2*covParameters$theta*tdiff)
    gp$tdiff <- tdiff
    gp$Y <- list()
    gp$Y[[1]] <- YfromGamma(gamma[[1]],invrootQeigs=gp$invrootQeigs,mu=covParameters$mu)
    sapply(2:length(tdiff),function(i){gp$Y[[i]] <<- gp$at[[i]]+gp$bt[[i]]*gp$Y[[i-1]]+sqrt(gp$gt[[i]])*YfromGamma(gamma[[i]],invrootQeigs=gp$invrootQeigs,mu=0)})       
    gp$expY <- lapply(gp$Y,exp) # for computational purposes
    gp$mcens <- fftgrid$mcens
    gp$ncens <- fftgrid$ncens
    
    class(gp) <- "stGPrealisation"    
    return(gp)    
    
}

##' GPdrv2 function
##'
##' A function to compute the second derivative of the log target with respect to the paramters of the latent field. Not intended for general purpose use.
##'
##' @param GP an object of class GPrealisation 
##' @param prior priors for the model 
##' @param Z design matirix on the FFT grid
##' @param Zt transpose of the design matrix
##' @param eta vector of parameters, eta
##' @param beta vector of parameters, beta
##' @param nis cell counts on the extended grid
##' @param cellarea the cell area 
##' @param spatial the poisson offset
##' @param gradtrunc gradient truncation parameter 
##' @param fftgrid an object of class FFTgrid
##' @param covfunction the choice of covariance function, see ?CovFunction 
##' @param d matrix of toral distances
##' @param eps the finite difference step size 
##' @return first and second derivatives of the log target at the specified paramters Y, eta and beta
##' @export

GPdrv2 <- function(GP,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc,fftgrid,covfunction,d,eps=1e-6){

    tar <- target.and.grad.spatialPlusPars(GP,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc)$logtarget
              
    cp1 <- CovParameters(list(sigma=exp(eta[1]+eps),phi=exp(eta[2])))      
    gp1 <- GPrealisation(gamma=GP$gamma,fftgrid=fftgrid,covFunction=covfunction,covParameters=cp1,d=d)
    tar1 <- target.and.grad.spatialPlusPars(GP=gp1,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc)$logtarget
    cp11 <- CovParameters(list(sigma=exp(eta[1]-eps),phi=exp(eta[2])))      
    gp11 <- GPrealisation(gamma=GP$gamma,fftgrid=fftgrid,covFunction=covfunction,covParameters=cp11,d=d)
    tar11 <- target.and.grad.spatialPlusPars(GP=gp11,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc)$logtarget
    
    cp2 <- CovParameters(list(sigma=exp(eta[1]),phi=exp(eta[2]+eps)))      
    gp2 <- GPrealisation(gamma=GP$gamma,fftgrid=fftgrid,covFunction=covfunction,covParameters=cp2,d=d)
    tar2 <- target.and.grad.spatialPlusPars(GP=gp2,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc)$logtarget
    cp22 <- CovParameters(list(sigma=exp(eta[1]),phi=exp(eta[2]-eps)))      
    gp22 <- GPrealisation(gamma=GP$gamma,fftgrid=fftgrid,covFunction=covfunction,covParameters=cp22,d=d)
    tar22 <- target.and.grad.spatialPlusPars(GP=gp22,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc)$logtarget
    
    dlsigma <- (tar1-tar11)/(2*eps)
    dlphi <- (tar2-tar22)/(2*eps)    
    
    cpA <- CovParameters(list(sigma=exp(eta[1]+eps),phi=exp(eta[2]+eps)))      
    gpA <- GPrealisation(gamma=GP$gamma,fftgrid=fftgrid,covFunction=covfunction,covParameters=cpA,d=d)
    tarA <- target.and.grad.spatialPlusPars(GP=gpA,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc)$logtarget
    
    cpB <- CovParameters(list(sigma=exp(eta[1]+eps),phi=exp(eta[2]-eps)))      
    gpB <- GPrealisation(gamma=GP$gamma,fftgrid=fftgrid,covFunction=covfunction,covParameters=cpB,d=d)
    tarB <- target.and.grad.spatialPlusPars(GP=gpB,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc)$logtarget
    
    cpC <- CovParameters(list(sigma=exp(eta[1]-eps),phi=exp(eta[2]+eps)))      
    gpC <- GPrealisation(gamma=GP$gamma,fftgrid=fftgrid,covFunction=covfunction,covParameters=cpC,d=d)
    tarC <- target.and.grad.spatialPlusPars(GP=gpC,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc)$logtarget
    
    cpD <- CovParameters(list(sigma=exp(eta[1]-eps),phi=exp(eta[2]-eps)))      
    gpD <- GPrealisation(gamma=GP$gamma,fftgrid=fftgrid,covFunction=covfunction,covParameters=cpD,d=d)
    tarD <- target.and.grad.spatialPlusPars(GP=gpD,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc)$logtarget  
    
    d2lsigma2 <- (tar1-2*tar+tar11)/(eps^2)
    d2lphi2 <- (tar2-2*tar+tar22)/(eps^2)
    d2lsigmalphi <- (tarA-tarB-tarC+tarD)/(4*eps^2)
    hess <- matrix(c(d2lsigma2,d2lsigmalphi,d2lsigmalphi,d2lphi2),2,2)
    
    return(list(dlsigma=dlsigma,dlphi=dlphi,hess=hess))
}




##' GPdrv2_Multitype function
##'
##' A function to compute the second derivatives of the log target for the multivariate model with respect to the paramters of the latent field. Not intended for general use.
##'
##' @param GPlist a list of objects of class GPrealisation 
##' @param priorlist list of priors for the model
##' @param Zlist list of design matirices on the FFT grid
##' @param Ztlist list of transpose design matrices
##' @param etalist list of parameters, eta, for each realisation
##' @param betalist clist of parameters, beta, for each realisation
##' @param nis cell counts of each type the extended grid
##' @param cellarea  the cell area  
##' @param spatial list of poisson offsets for each type
##' @param gradtrunc  gradient truncation parameter
##' @param fftgrid an object of class FFTgrid
##' @param covfunction list giving the choice of covariance function for each type, see ?CovFunction  
##' @param d matrix of toral distances
##' @param eps the finite difference step size 
##' @param k index of type for which to compute the gradient and hessian
##' @return first and second derivatives of the log target  for tyupe k at the specified paramters Y, eta and beta
##' @export

GPdrv2_Multitype <- function(GPlist,priorlist,Zlist,Ztlist,etalist,betalist,nis,cellarea,spatial,gradtrunc,fftgrid,covfunction,d,eps=1e-6,k){ # k is list index

    tar <- target.and.grad.MultitypespatialPlusPars(GPlist=GPlist,priorlist=priorlist,Zlist=Zlist,Ztlist=Ztlist,eta=etalist,beta=betalist,nis=nis,cellarea=cellarea,spatial=spatial,gradtrunc=gradtrunc)$logtarget
    
    
    etatemp <- etalist
    gptemp <- GPlist
              
    cp1 <- CovParameters(list(sigma=exp(etalist[[k]][1]+eps),phi=exp(etalist[[k]][2])))      
    gp1 <- GPrealisation(gamma=GPlist[[k]]$gamma,fftgrid=fftgrid,covFunction=covfunction[[k]],covParameters=cp1,d=d)
    gplist1 <- gptemp
    gplist1[[k]] <- gp1
    tar1 <- target.and.grad.MultitypespatialPlusPars(GPlist=gplist1,priorlist=priorlist,Zlist=Zlist,Ztlist=Ztlist,eta=etalist,beta=betalist,nis=nis,cellarea=cellarea,spatial=spatial,gradtrunc=gradtrunc)$logtarget
    cp11 <- CovParameters(list(sigma=exp(etalist[[k]][1]-eps),phi=exp(etalist[[k]][2])))      
    gp11 <- GPrealisation(gamma=GPlist[[k]]$gamma,fftgrid=fftgrid,covFunction=covfunction[[k]],covParameters=cp11,d=d)
    gplist11 <- gptemp
    gplist11[[k]] <- gp11
    tar11 <- target.and.grad.MultitypespatialPlusPars(GPlist=gplist11,priorlist=priorlist,Zlist=Zlist,Ztlist=Ztlist,eta=etalist,beta=betalist,nis=nis,cellarea=cellarea,spatial=spatial,gradtrunc=gradtrunc)$logtarget
    
    cp2 <- CovParameters(list(sigma=exp(etalist[[k]][1]),phi=exp(etalist[[k]][2]+eps)))      
    gp2 <- GPrealisation(gamma=GPlist[[k]]$gamma,fftgrid=fftgrid,covFunction=covfunction[[k]],covParameters=cp2,d=d)
    gplist2 <- gptemp
    gplist2[[k]] <- gp2
    tar2 <- target.and.grad.MultitypespatialPlusPars(GPlist=gplist2,priorlist=priorlist,Zlist=Zlist,Ztlist=Ztlist,eta=etalist,beta=betalist,nis=nis,cellarea=cellarea,spatial=spatial,gradtrunc=gradtrunc)$logtarget
    cp22 <- CovParameters(list(sigma=exp(etalist[[k]][1]),phi=exp(etalist[[k]][2]-eps)))      
    gp22 <- GPrealisation(gamma=GPlist[[k]]$gamma,fftgrid=fftgrid,covFunction=covfunction[[k]],covParameters=cp22,d=d)
    gplist22 <- gptemp
    gplist22[[k]] <- gp22
    tar22 <- target.and.grad.MultitypespatialPlusPars(GPlist=gplist22,priorlist=priorlist,Zlist=Zlist,Ztlist=Ztlist,eta=etalist,beta=betalist,nis=nis,cellarea=cellarea,spatial=spatial,gradtrunc=gradtrunc)$logtarget
    
    dlsigma <- (tar1-tar11)/(2*eps)
    dlphi <- (tar2-tar22)/(2*eps)    
    
    cpA <- CovParameters(list(sigma=exp(etalist[[k]][1]+eps),phi=exp(etalist[[k]][2]+eps)))      
    gpA <- GPrealisation(gamma=GPlist[[k]]$gamma,fftgrid=fftgrid,covFunction=covfunction[[k]],covParameters=cpA,d=d)
    gplistA <- gptemp
    gplistA[[k]] <- gpA
    tarA <- target.and.grad.MultitypespatialPlusPars(GPlist=gplistA,priorlist=priorlist,Zlist=Zlist,Ztlist=Ztlist,eta=etalist,beta=betalist,nis=nis,cellarea=cellarea,spatial=spatial,gradtrunc=gradtrunc)$logtarget
    
    cpB <- CovParameters(list(sigma=exp(etalist[[k]][1]+eps),phi=exp(etalist[[k]][2]-eps)))      
    gpB <- GPrealisation(gamma=GPlist[[k]]$gamma,fftgrid=fftgrid,covFunction=covfunction[[k]],covParameters=cpB,d=d)
    gplistB <- gptemp
    gplistB[[k]] <- gpB
    tarB <- target.and.grad.MultitypespatialPlusPars(GPlist=gplistB,priorlist=priorlist,Zlist=Zlist,Ztlist=Ztlist,eta=etalist,beta=betalist,nis=nis,cellarea=cellarea,spatial=spatial,gradtrunc=gradtrunc)$logtarget
    
    cpC <- CovParameters(list(sigma=exp(etalist[[k]][1]-eps),phi=exp(etalist[[k]][2]+eps)))      
    gpC <- GPrealisation(gamma=GPlist[[k]]$gamma,fftgrid=fftgrid,covFunction=covfunction[[k]],covParameters=cpC,d=d)
    gplistC <- gptemp
    gplistC[[k]] <- gpC
    tarC <- target.and.grad.MultitypespatialPlusPars(GPlist=gplistC,priorlist=priorlist,Zlist=Zlist,Ztlist=Ztlist,eta=etalist,beta=betalist,nis=nis,cellarea=cellarea,spatial=spatial,gradtrunc=gradtrunc)$logtarget
    
    cpD <- CovParameters(list(sigma=exp(etalist[[k]][1]-eps),phi=exp(etalist[[k]][2]-eps)))      
    gpD <- GPrealisation(gamma=GPlist[[k]]$gamma,fftgrid=fftgrid,covFunction=covfunction[[k]],covParameters=cpD,d=d)
    gplistD <- gptemp
    gplistD[[k]] <- gpD
    tarD <- target.and.grad.MultitypespatialPlusPars(GPlist=gplistD,priorlist=priorlist,Zlist=Zlist,Ztlist=Ztlist,eta=etalist,beta=betalist,nis=nis,cellarea=cellarea,spatial=spatial,gradtrunc=gradtrunc)$logtarget
    
        
    
    d2lsigma2 <- (tar1-2*tar+tar11)/(eps^2)
    d2lphi2 <- (tar2-2*tar+tar22)/(eps^2)
    d2lsigmalphi <- (tarA-tarB-tarC+tarD)/(4*eps^2)
    hess <- matrix(c(d2lsigma2,d2lsigmalphi,d2lsigmalphi,d2lphi2),2,2)
    
    return(list(dlsigma=dlsigma,dlphi=dlphi,hess=hess))
}

##' GPdrv function
##'
##' A function to compute the first derivatives of the log target with respect to the paramters of the latent field. Not intended for general purpose use.
##'
##' @param GP an object of class GPrealisation 
##' @param prior priors for the model 
##' @param Z design matirix on the FFT grid
##' @param Zt transpose of the design matrix
##' @param eta vector of parameters, eta
##' @param beta vector of parameters, beta
##' @param nis cell counts on the extended grid
##' @param cellarea the cell area 
##' @param spatial the poisson offset
##' @param gradtrunc gradient truncation parameter 
##' @param fftgrid an object of class FFTgrid
##' @param covfunction the choice of covariance function, see ?CovFunction 
##' @param d matrix of toral distances
##' @param eps the finite difference step size 
##' @return first  derivatives of the log target at the specified paramters Y, eta and beta
##' @export

GPdrv <- function(GP,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc,fftgrid,covfunction,d,eps=1e-6){

    tar <- target.and.grad.spatialPlusPars(GP,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc)$logtarget
              
    cp1 <- CovParameters(list(sigma=exp(eta[1]+eps),phi=exp(eta[2])))      
    gp1 <- GPrealisation(gamma=GP$gamma,fftgrid=fftgrid,covFunction=covfunction,covParameters=cp1,d=d)
    tar1 <- target.and.grad.spatialPlusPars(GP=gp1,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc)$logtarget

    cp2 <- CovParameters(list(sigma=exp(eta[1]),phi=exp(eta[2]+eps)))      
    gp2 <- GPrealisation(gamma=GP$gamma,fftgrid=fftgrid,covFunction=covfunction,covParameters=cp2,d=d)
    tar2 <- target.and.grad.spatialPlusPars(GP=gp2,prior,Z,Zt,eta,beta,nis,cellarea,spatial,gradtrunc)$logtarget
      
    dlsigma <- (tar1-tar)/(eps)
    dlphi <- (tar2-tar)/(eps)
    
    return(list(dlsigma=dlsigma,dlphi=dlphi))
}





##' PriorSpec function
##'
##' Generic for declaring that an object is of valid type for use as as prior in lgcp. For further details and examples, see the vignette "Bayesian_lgcp". 
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method PriorSpec
##' @seealso \link{PriorSpec.list}
##' @export

PriorSpec <- function(obj,...){
    UseMethod("PriorSpec")
}



##' PriorSpec.list function
##'
##' Method for declaring a Bayesian prior density in lgcp. Checks to confirm that the object obj has the requisite components for functioning as a prior.
##'
##' @method PriorSpec list
##' @param obj a list object defining a prior , see ?GaussianPrior and ?LogGaussianPrior
##' @param ... additional arguments
##' @return an object suitable for use in a call to the MCMC routines
##' @seealso \link{GaussianPrior}, \link{LogGaussianPrior}
##' @examples 
##' \dontrun{PriorSpec(LogGaussianPrior(mean=log(c(1,500)),variance=diag(0.15,2)))}
##' \dontrun{PriorSpec(GaussianPrior(mean=rep(0,9),variance=diag(10^6,9)))}
##' @export

PriorSpec.list <- function(obj,...){
    if(is.null(obj$mean)){
        stop("obj must have an element $mean, the prior mean")
    }
    if(is.null(obj$variance)){
        stop("obj must have an element $variance, the prior variance")
    }
    if(is.null(obj$eval)){
        stop("obj must have an element $eval, a function returning the log likelihood")
    }
    if(is.null(obj$inverse_transform)){
        stop("obj must have an element $inverse_transform, a function returning the inverse of the parameter transform (which can be the identity function, see ?identity)")
    }    
    
    class(obj) <- c("PriorSpec",class(obj))
    return(obj)
}


##' LogGaussianPrior function
##'
##' A function to create a Gaussian prior on the log scale
##'
##' @param mean a vector of length 2 representing the mean (on the log scale)
##' @param variance a 2x2 matrix representing the variance (on the log scale)
##' @return an object of class LogGaussianPrior that can be passed to the function PriorSpec.
##' @seealso \link{GaussianPrior}, link{PriorSpec.list}
##' @examples 
##' \dontrun{LogGaussianPrior(mean=log(c(1,500)),variance=diag(0.15,2))}
##' @export

LogGaussianPrior <- function(mean,variance){
    prior <- list()
    prior$mean <- mean
    prior$variance <- as.matrix(variance)
    prior$determinant <- det(prior$variance)
    prior$precision <- solve(prior$variance)
    prior$dim <- length(mean)
    prior$const <- -(prior$dim/2)*log(2*pi) - (1/2)*log(prior$determinant)
    
    prior$eval <- function(x){
        vec <- prior$precision%*%(x-prior$mean)
        value <- prior$const -(1/2)*t(x-prior$mean)%*%vec
        return(list(loglik=value,gradcontrib=vec))
    }
    
    prior$transform <- log # ie the log function
    prior$inverse_transform <- exp # ie the exponential function
    prior$jacobian <- exp # ie (dsigma/dlogsigma,dphi/dlogphi), a function of (logsigma,logphi)
    class(prior) <- c("LogGaussianPrior","list")
    return(prior)
}

##' GaussianPrior function
##'
##' A function to create a Gaussian prior.
##'
##' @param mean a vector of length 2 representing the mean.
##' @param variance a 2x2 matrix representing the variance.
##' @return an object of class LogGaussianPrior that can be passed to the function PriorSpec.
##' @seealso \link{LogGaussianPrior}, link{PriorSpec.list}
##' @examples 
##' \dontrun{GaussianPrior(mean=rep(0,9),variance=diag(10^6,9))}
##' @export

GaussianPrior <- function(mean,variance){
    prior <- list()
    prior$mean <- mean
    prior$variance <- as.matrix(variance)
    prior$determinant <- det(prior$variance)
    prior$precision <- solve(prior$variance)
    prior$dim <- length(mean)
    prior$const <- -(prior$dim/2)*log(2*pi) - (1/2)*log(prior$determinant)
    
    prior$eval <- function(x){
        vec <- prior$precision%*%(x-prior$mean)
        value <- prior$const -(1/2)*t(x-prior$mean)%*%vec
        return(list(loglik=value,gradcontrib=vec))
    }
    
    prior$inverse_transform <- identity # ie the identity function
    class(prior) <- c("GaussianPrior","list")
    return(prior)
}


##' lgcpPrior function
##'
##' A function to create the prior for beta and eta ready for a run of the MCMC algorithm. 
##'
##' @param etaprior an object of class PriorSpec defining the prior for the parameters of the latent field, eta. See ?PriorSpec.list.
##' @param betaprior etaprior an object of class PriorSpec defining the prior for the parameters of main effects, beta. See ?PriorSpec.list.
##' @return an R structure representing the prior density ready for a run of the MCMC algorithm.
##' @seealso \link{GaussianPrior}, \link{LogGaussianPrior}, \link{PriorSpec.list}, \link{minimum.contrast}, \link{minimum.contrast.spatiotemporal}, \link{chooseCellwidth}, \link{getpolyol}, \link{guessinterp}, \link{getZmat},
##' \link{addTemporalCovariates}, \link{lgcpPrior}, \link{lgcpInits}, \link{CovFunction}
##' \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}, \link{lgcpPredictSpatioTemporalPlusPars}, 
##' \link{lgcpPredictMultitypeSpatialPlusPars}
##' @examples 
##' lgcpPrior(etaprior=PriorSpec(LogGaussianPrior(mean=log(c(1,500)),
##'     variance=diag(0.15,2))),betaprior=PriorSpec(GaussianPrior(mean=rep(0,9),
##'     variance=diag(10^6,9))))
##' @export

lgcpPrior <- function(etaprior=NULL,betaprior=NULL){
    # eta are parameters for the latent stochastic process(es). eta minimally consists of sigma and phi
    # beta are parameters for the regression part of the model.
    if(!is.null(etaprior)){
        if(!inherits(etaprior,"PriorSpec")){
            stop("etaprior must be an object of class PriorSpec")
        }
    }
    if(!is.null(betaprior)){
        if(!inherits(betaprior,"PriorSpec")){
            stop("betaprior must be an object of class PriorSpec")
        }
    }
    obj <- list()
    obj$etaprior <- etaprior 
    obj$betaprior <- betaprior
    class(obj) <- "lgcpPrior"
    return(obj)
}

##' CovParameters function
##'
##' A function to provide a structure for the parameters of the latent field. Not intended for general use.
##'
##' @param list a list
##' @return an object used in the MCMC routine.
##' @export

CovParameters <- function(list){ # spatial or spatiotemporal parameters
    if(is.null(list$sigma) | is.null(list$phi)){
        stop("argument 'list' must be a named list object containing elements $sigma and $phi as a minimum.")
    }
    
    if(!is.null(list$mu)){
        stop("argument 'list' must NOT contain an element $mu, as this is set by default")
    }
    else{
        list$mu <- -list$sigma^2/2
    }
    class(list) <- "CovParamaters"
    return(list)
}

##' CovFunction function
##'
##' A Generic method used to specify the choice of covariance function for use in the MCMC algorithm. For further details and examples, 
##' see the vignette "Bayesian_lgcp".
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method CovFunction
##' @seealso \link{CovFunction.function}, \link{exponentialCovFct}, \link{RandomFieldsCovFct}, \link{SpikedExponentialCovFct}
##' @export

CovFunction <- function(obj,...){
    UseMethod("CovFunction")
}



##' CovFunction.function function
##'
##' A function used to define the covariance function for the latent field prior to running the MCMC algorithm
##'
##' @method CovFunction function
##' @param obj a function object
##' @param ... additional arguments
##' @return the covariance function ready to run the MCMC routine.
##' @seealso \link{exponentialCovFct}, \link{RandomFieldsCovFct}, \link{SpikedExponentialCovFct}, \link{CovarianceFct}
##' @examples 
##' \dontrun{cf1 <- CovFunction(exponentialCovFct)}
##' \dontrun{cf2 <- CovFunction(RandomFieldsCovFct(model="matern",additionalparameters=1))}
##' @export

CovFunction.function <- function(obj,...){
    ans <- obj(0,list(sigma=1,phi=1))
    if(is.null(ans$eval)){
        stop("Invalid covariance function definition, function must return a list including an element $eval, the evaluation function")
    }
    class(obj) <- c("CovFunction","function")
    return(obj)
}

##' exponentialCovFct function
##'
##' A function to declare and also evaluate an exponential covariance function.
##'
##' @param d toral distance
##' @param CovParameters parameters of the latent field, an object of class "CovParamaters".
##' @return the exponential covariance function
##' @seealso \link{CovFunction.function}, \link{RandomFieldsCovFct}, \link{SpikedExponentialCovFct}
##' @export

exponentialCovFct <- function(d,CovParameters){
    x <- exp(-d/CovParameters$phi)
    ans <- list()
    ans$eval <- CovParameters$sigma^2*x
    return(ans)
}

##' SpikedExponentialCovFct function
##'
##' A function to declare and also evaluate a spiked exponential covariance function. Note that the present version of
##' lgcp only offers estimation for sigma and phi, the additional parameter 'spikevar' is treated as fixed.
##'
##' @param d toral distance
##' @param CovParameters parameters of the latent field, an object of class "CovParamaters".
##' @param spikevar the additional variance at distance 0
##' @return the spiked exponential covariance function; note that the spikevariance is currently not estimated as part of the MCMC routine, and is thus treated as a fixed parameter.
##' @seealso \link{CovFunction.function}, \link{exponentialCovFct}, \link{RandomFieldsCovFct}
##' @export

SpikedExponentialCovFct <- function(d,CovParameters,spikevar=1){
    x <- exp(-d/CovParameters$phi)
    ans <- list()
    ans$eval <- CovParameters$sigma^2*x
    ans$eval[d==0] <- ans$eval[d==0] + spikevar 
    return(ans)
}

##' RandomFieldsCovFct function
##'
##' A function to declare and also evaluate an covariance function from the RandomFields Package. See ?CovarianceFct. Note that the present version of
##' lgcp only offers estimation for sigma and phi, any additional paramters are treated as fixed.
##'
##' @param model the choice of model e.g.  "matern"
##' @param additionalparameters additional parameters for chosen covariance model. See ?CovarianceFct 
##' @return a covariance function from the RandomFields package
##' @seealso \link{CovFunction.function}, \link{exponentialCovFct}, \link{SpikedExponentialCovFct}, \link{CovarianceFct}
##' @examples 
##' \dontrun{RandomFieldsCovFct(model="matern",additionalparameters=1)}
##' @export

RandomFieldsCovFct <- function(model,additionalparameters=c()){
    force(model)
    force(additionalparameters)
    funk <- function(d,CovParameters){
        ans <- list()
        ans$eval <- gu(u=d,sigma=CovParameters$sigma,phi=CovParameters$phi,model=model,additionalparameters=additionalparameters)
        return(ans)
    }
    return(funk)
}


##' getCovParameters function
##'
##' Internal function for retrieving covariance parameters. not indended for general use.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method getCovParameters
##' @export

getCovParameters <- function(obj,...){
    UseMethod("getCovParameters")
}



##' getCovParameters.GPrealisation function
##'
##' Internal function for retrieving covariance parameters. not indended for general use.
##'
##' @method getCovParameters GPrealisation
##' @param obj an GPrealisation object
##' @param ... additional arguments
##' @return ...
##' @export

getCovParameters.GPrealisation <- function(obj,...){
    pars <- unlist(obj$CovParameters)
    pars <- pars[names(pars)!="mu"]
    return(pars)
}



##' getCovParameters.list function
##'
##' Internal function for retrieving covariance parameters. not indended for general use.
##'
##' @method getCovParameters list
##' @param obj an list object
##' @param ... additional arguments
##' @return ...
##' @export

getCovParameters.list <- function(obj,...){
    pars <- unlist(lapply(obj,function(x){return(x$CovParameters)}))
    pars <- pars[names(pars)!="mu"]
    return(pars)
}

##' BetaParameters function
##'
##' An internal function to declare a vector a parameter vector for the main effects.
##'
##' @param beta a vector
##' @return ...
##' @export

BetaParameters <- function(beta){
    obj <- beta
    class(obj) <- "BetaParameters"
    return(obj)
}

##' EvaluatePrior function
##'
##' An internal function used in the MCMC routine to evaluate the prior for a given set of parameters
##'
##' @param etaParameters the paramter eta
##' @param betaParameters the parameter beta
##' @param prior the prior
##' @return the prior evaluated at the given values.
##' @export

EvaluatePrior <- function(etaParameters,betaParameters,prior){
    if(!is.null(betaParameters)){
        if(class(betaParameters)!="BetaParameters"){
            stop("betaParameters must be an object of class BetaParameters")
        }
    }
    if(class(prior)!="lgcpPrior"){
        stop("prior must be an object of class lgcpPrior")
    }
 
    
    if(!is.null(betaParameters)){
        ans <- list(etacontrib=prior$etaprior$eval(etaParameters),betacontrib=prior$betaprior$eval(betaParameters))
        class(ans) <- "EvalPrior"
    }
    else{
        ans <- list(etacontrib=prior$etaprior$eval(etaParameters),betacontrib=NULL)
        class(ans) <- "EvalPrior"
    }
    return(ans)
}

##' lgcpInits function
##'
##' A function to declare initial values for a run of the MCMC routine. If specified, the MCMC algorithm will calibrate the proposal 
##' density using these as provisional estimates of the parameters.
##'
##' It is not necessary to supply intial values to the MCMC routine, by default the functions lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, 
##' lgcpPredictSpatioTemporalPlusPars and lgcpPredictMultitypeSpatialPlusPars will initialise the MCMC as follows. For eta, if no initial value is 
##' specified then the initial value of eta in the MCMC run will be the prior mean. For beta, if no initial value is specified then 
##' the initial value of beta in the MCMC run will be estimated from an overdispersed Poisson fit to the cell counts, ignoring spatial correlation. The user cannot
##' specify an initial value of Y (or equivalently Gamma), as a sensible value is chosen by the MCMC function.
##'
##' A secondary function of specifying initial values is to help design the MCMC proposal matrix, which is based on these initial estimates.
##'
##' @param etainit a vector, the initial value of eta to use
##' @param betainit a vector, the initial value of beta to use, this vector must have names the same as the variable names in the formula in use, and in the same order.
##' @return an object of class lgcpInits used in the MCMC routine.
##' @seealso \link{minimum.contrast}, \link{minimum.contrast.spatiotemporal}, \link{chooseCellwidth}, \link{getpolyol}, \link{guessinterp}, \link{getZmat},
##' \link{addTemporalCovariates}, \link{lgcpPrior}, \link{CovFunction},
##' \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}, \link{lgcpPredictSpatioTemporalPlusPars}, 
##' \link{lgcpPredictMultitypeSpatialPlusPars}
##' @examples 
##' \dontrun{INITS <- lgcpInits(etainit=log(c(sqrt(1.5),275)), betainit=NULL)}
##' @export

lgcpInits <- function(etainit=NULL,betainit=NULL){ # beta has a sensible initial value by default, via glm
    obj <- list()
    obj$etainit <- etainit
    obj$betainit <- betainit
    class(obj) <- "lgcpInits"
    return(obj)
}

##' GPlist2array function
##'
##' An internal function for turning a list of GPrealisation objects into an an array by a particular common element of the GPrealisation object
##'
##' @param GPlist an object of class GPrealisation
##' @param element the name of the element of GPlist[[1]] (for example) to extract, e.g. "Y"
##' @return an array
##' @export

GPlist2array <- function(GPlist,element){
    n <- length(GPlist)
    arr <- array(dim=c(dim(GPlist[[1]][[element]]),n))
    for(i in 1:n){
        arr[,,i] <- GPlist[[i]][[element]]
    }
    return(arr)
}

##' formulaList function
##'
##' A function to creat an object of class "formulaList" from a list of "formula" objects; use to define the model for the main effects
##' prior to running the multivariate MCMC algorithm.
##'
##' @param X a list object, each element of which is a formula
##' @return an object of class "formulaList"
##' @export

formulaList <- function(X){
    if(class(X)!="list"){
        stop("X must be a list object, each element bing of class 'formula'")
    }
    for(i in 1:length(X)){
        if(class(X[[i]])!="formula"){
            stop("Each element of X must be of class formula")
        }
    }
    class(X) <- c("formulaList","list")
    return(X)
}

##' aggregateformulaList function
##'
##' An internal function to collect terms from a formulalist. Not intended for general use. 
##'
##' @param x an object of class "formulaList"
##' @param ... other arguments
##' @return a formula of the form X ~ var1 + var2 tec.
##' @export
aggregateformulaList <- function(x,...){
    vrs <- unique(unlist(lapply(x,function(z){termsinformula(z)})))
    rhs <- paste(vrs,collapse="+")
    if(rhs==""){
        rhs <- "1"
    }
    f <- formula(paste("X~",rhs))
    attr(f, ".Environment")= .GlobalEnv
    return(f)
}  

##' getLHSformulaList function
##'
##' A function to retrieve the dependent variables from a formulaList object. Not intended for general use. 
##'
##' @param fl an object of class "formulaList"
##' @return the indepentdent variables
##' @export

getLHSformulaList <- function(fl){
    return(sapply(fl,function(x){as.character(x[2])}))
}  

##' getZmats function
##'
##' An internal function to create Z_k from an lgcpZmat object, for use in the multivariate MCMC algorithm. Not intended for general use.  
##'
##' @param Zmat an objecty of class "lgcpZmat"
##' @param formulaList an object of class "formulaList"
##' @return design matrices for each of the point types
##' @export

getZmats <- function(Zmat,formulaList){      
    dmat <- attr(Zmat,"data.frame")
    cellInside <- attr(Zmat,"cellInside")
    anymiss <- attr(Zmat,"anymiss")
    M <- attr(Zmat,"M")
    N <- attr(Zmat,"N")
    
    lhs <- getLHSformulaList(formulaList)
    
    Zmatlist <- list()
    for(i in 1:length(formulaList)){
        dmat[[lhs[i]]] <- dmat$X
        DM <- model.matrix(formulaList[[i]],data=dmat)
        Zmatlist[[i]] <- matrix(0,M*N,dim(DM)[2])
        colnames(Zmatlist[[i]]) <- colnames(DM)
        if (!is.null(anymiss)){
            Zmatlist[[i]][as.logical(cellInside),][!anymiss,] <-  DM # [!anymiss] because NAs are removed by model.matrix
        }
        else{
            Zmatlist[[i]][as.logical(cellInside),] <-  DM
        }
    }
    return(Zmatlist)  
}

## Redundant:
## we need to decide exactly what we're going to do here for the sharing of parameters
##SH <- function(expr,group){
##    f <- identity(expr) # a function for sharing parameters between models
##    attr(f,"group") <- group
##    return(f)
##}

##' fftmultiply function
##'
##' A function to pre-multiply a vector by a block cirulant matrix 
##'
##' @param efb eigenvalues of the matrix
##' @param vector the vector
##' @return a vector: the product of the matrix and the vector.
##' @export

fftmultiply <- function(efb,vector){ # efb == "eigenfrombase(base)" == Re(fft(base)) ... if constant, then compute only once
    return((1/length(vector))*Re(fft(efb*fft(vector,inverse=TRUE))))
}

##' samplePosterior function
##'
##' A function to draw a sample from the posterior of a spatial LGCP. Randomly selects an index i, and returns the ith value of eta, 
##' the ith value of beta and the ith value of Y as a named list.
##'
##' @param x an object of class lgcpPredictSpatialOnlyPlusParameters or lgcpPredictAggregateSpatialPlusParameters
##' @return a sample from the posterior named list object with names elements "eta", "beta" and "Y".
##' @export

samplePosterior <- function(x){
    n <- dim(x$betarec)[1]
    idx <- sample(1:n,1)
    ans <- list()
    ans$beta <- x$betarec[idx,]
    ans$eta <- x$etarec[idx,]
    
    fn <- paste(x$gridfunction$dirname,"simout.nc",sep="")
    ncdata <- nc_open(fn)    
        
    ans$Y <- ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(1,1,1,idx), count=c(-1,-1,-1,1))
    
    nc_close(ncdata)
    
    return(ans)
}

##' plotit function
##'
##' A function to plot various objects. A developmental tool: not intended for general use
##'
##' @param x an a list, matrix, or GPrealisation object.
##' @return plots the objects.
##' @export

plotit <- function(x){
    if(class(x)=="list"){
        lapply(x,function(x){image.plot(x);browser()})
        return()
    }
    if(class(x)=="matrix"){
        image.plot(x)
        return()
    }
    if(class(x)=="GPrealisation"){
        image.plot(x$mcens,x$ncens,x$Y)
    }    
}

##' getup function
##'
##' A function to get an object from a parent frame.
##'
##' @param n a character string, the name of the object 
##' @param lev how many levels up the hierarchy to go (see the argument "envir" from the function "get"), default is 1.
##' @return ...
##' @export

getup <- function(n,lev=1){
    return(get(n,envir=parent.frame(lev+1)))
}




##' checkObsWin function
##'
##' A function to run on an object generated by the "selectObsWindow" function. Plots the observation window with grid, use as a visual aid to check
##' the choice of cell width is correct.
##'
##' @param ow an object generated by selectObsWindow, see ?selectObsWindow
##' @return a plot of the observation window and grid
##' @seealso \link{chooseCellwidth}
##' @export

checkObsWin <- function(ow){
    xv <- ow$xvals
    yv <- ow$yvals
    gr <- expand.grid(xv,yv)
    plot(ow$xyt$window)
    points(gr,pch="+",cex=0.5,col="red")
    title(xlab=length(xv),ylab=length(yv))
}

##' tempRaster function
##'
##' A function to create a temporary raster object from an x-y regular grid of cell centroids. Useful for projection from one raster to another.
##'
##' @param mcens vector of equally-spaced coordinates of cell centroids in x-direction
##' @param ncens vector of equally-spaced coordinates of cell centroids in y-direction
##' @return an empty raster object
##' @export

tempRaster <- function(mcens,ncens){
    M <- length(mcens)
    N <- length(ncens)
    
    dx <- diff(mcens[1:2])
    dy <- diff(ncens[1:2])
    
    return(raster(nrows=M,ncols=N,xmn=mcens[1]-dx/2,xmx=mcens[M]+dx/2,ymn=ncens[1]-dy/2,ymx=ncens[N]+dy/2))
}

##' gOverlay function
##'
##' A function to overlay the FFT grid, a SpatialPolygons object, onto a SpatialPolygonsDataFrame object. 
##'
##' this code was adapted from Roger Bivand:\cr
##' https://stat.ethz.ch/pipermail/r-sig-geo/2011-June/012099.html
##'
##' @param grid the FFT grid, a SpatialPolygons object
##' @param spdf  a SpatialPolygonsDataFrame object
##' @return a matrix describing the features of the overlay: the originating indices of grid and spdf (all non-trivial intersections) and the area of each intersection.
##' @export


gOverlay <- function(grid,spdf){
    int <- gIntersects(spdf,grid,byid=TRUE)  
    cnt <- 0
    idx <- c()
    cellnum <- 0
    pb <- txtProgressBar(min=0,max=nrow(int)*ncol(int),style=3)
    for (i in 1:nrow(int)){
        for (j in 1:ncol(int)){
            cellnum <- cellnum + 1           
            if(int[i,j]){
                cnt <- cnt + 1
                vec <- try(gIntersection(grid[i,],spdf[j,], byid=TRUE),silent=TRUE)
                if(inherits(vec,"try-error")|inherits(vec,"SpatialPoints")|inherits(vec,"SpatialLines")){
                    cnt <- cnt - 1
                }
                else{
                    idx <- rbind(idx,c(i,j,sum(sapply(slot(vec,"polygons"),slot,"area"))))
                }
            }
            setTxtProgressBar(pb,cellnum)
        }
    }
    close(pb) 
    idx <- data.frame(idx)
    names(idx) <- c("grididx","polyidx","area")
    ans <- list()
    ans$info <- idx
    return(ans)
}

##' getZmat function
##'
##' A function to construct a design matrix for use with the Bayesian MCMC routines in lgcp. See the vignette "Bayesian_lgcp" for further details on
##' how to use this function.\cr
##' 
##' For example, a spatial LGCP model for the would have the form:\cr
##' \cr
##' X(s) ~ Poisson[R(s)]\cr
##' \cr
##' R(s) = C_A lambda(s) exp[Z(s)beta+Y(s)]\cr
##' \cr 
##'
##' The function getZmat helps create the matrix Z. The returned object is passed onto an MCMC function, for example lgcpPredictSpatialPlusPars or
##' lgcpPredictAggregateSpatialPlusPars. This function can also be used to help construct Z for use with lgcpPredictSpatioTemporalPlusPars and
##' lgcpPredictMultitypeSpatialPlusPars, but these functions require a list of such objects: see the vignette "Bayesian_lgcp" for examples. 
##'
##' @param formula a formula object of the form X ~ var1 + var2 etc. The name of the dependent variable must be "X". Only accepts 'simple' formulae, such as the example given.
##' @param data the data to be analysed (using, for example lgcpPredictSpatialPlusPars). Either an object of class ppp, or an object of class SpatialPolygonsDataFrame
##' @param regionalcovariates an optional SpatialPolygonsDataFrame object containing covariate information, if applicable
##' @param pixelcovariates an optional SpatialPixelsDataFrame object containing covariate information, if applicable
##' @param cellwidth the width of computational cells
##' @param ext integer multiple by which grid should be extended, default is 2. Generally this will not need to be altered, but if the spatial correlation decays slowly, increasing 'ext' may be necessary.
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former, the default, includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window.
##' @param overl an object of class "lgcppolyol", created by the function getpolyol. Such an object contains the FFT grid and a polygon/polygon overlay and speeds up computation massively.
##' @return a design matrix for passing on to the Bayesian MCMC functions
##' @seealso \link{minimum.contrast}, \link{minimum.contrast.spatiotemporal}, \link{chooseCellwidth}, \link{getpolyol}, \link{guessinterp},
##' \link{addTemporalCovariates}, \link{lgcpPrior}, \link{lgcpInits}, \link{CovFunction}
##' \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}, \link{lgcpPredictSpatioTemporalPlusPars}, 
##' \link{lgcpPredictMultitypeSpatialPlusPars}
##' @export

getZmat <- function(formula,data,regionalcovariates=NULL,pixelcovariates=NULL,cellwidth,ext=2,inclusion="touching",overl=NULL){
    if(!is.null(overl)){
        cat("Using 'cellwidth' and 'ext' from overl\n")
        cellwidth <- overl$cellwidth
        ext <- overl$ext
    }
    if(inherits(data,"SpatialPolygonsDataFrame")){
        spatstat.options(checkpolygons=FALSE)
        W <- as(gUnaryUnion(data),"owin")
        spatstat.options(checkpolygons=TRUE)
        sd <- ppp(window=W)         
    }
    else{    
        sd <- data
    }
    ow <- selectObsWindow(sd,cellwidth) 
	sd <- ow$xyt
	M <- ow$M # note for this function, M and N are powers of 2 
	N <- ow$N
	study.region <- sd$window
	if(!is.null(overl)){
	    gridobj <- overl$gridobj
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
	
	ans <- cov.interp.fft(formula=formula,W=study.region,regionalcovariates=regionalcovariates,pixelcovariates=pixelcovariates,mcens=mcens[1:M],ncens=ncens[1:N],cellInside=cellInside[1:M,1:N],overl=overl)
	attr(ans,"gridobj") <- gridobj
	attr(ans,"inclusion") <- inclusion
	attr(ans,"ext") <- ext
	attr(ans,"cellwidth") <- cellwidth 
	class(ans) <- c("lgcpZmat","matrix")
    return(ans)
}

##' chooseCellwidth function
##'
##' A function to help choose the cell width (the parameter "cellwidth" in lgcpPredictSpatialPlusPars, for example) prior to setting up the FFT grid, 
##' before an MCMC run. 
##' 
##' Ideally this function should be used after having made a preliminary guess at the parameters of the latent field.The idea is to run chooseCellwidth 
##' several times, adjusting the parameter "cwinit" so as to balance available computational resources with output grid size. 
##'
##' @param obj an object of class ppp, stppp, SpatialPolygonsDataFrame, or owin
##' @param cwinit the cell width
##' @return produces a plot of the observation window and computational grid. 
##' @seealso \link{minimum.contrast}, \link{minimum.contrast.spatiotemporal}, \link{getpolyol}, \link{guessinterp}, \link{getZmat},
##' \link{addTemporalCovariates}, \link{lgcpPrior}, \link{lgcpInits}, \link{CovFunction}
##' \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}, \link{lgcpPredictSpatioTemporalPlusPars}, 
##' \link{lgcpPredictMultitypeSpatialPlusPars}
##' @export

chooseCellwidth <- function(obj,cwinit){

    if(inherits(obj,"SpatialPolygons")){
        spatstat.options(checkpolygons=FALSE)
        W <- as(gUnaryUnion(obj),"owin")
        spatstat.options(checkpolygons=TRUE)
    }
    else if(inherits(obj,"ppp") | inherits(obj,"stppp")){
        W <- obj$window
    }
    else if(inherits(obj,"owin")){
        W <- obj
    }
    else{
        stop("Cannot extract observation window.")
    }
    
    OW <- selectObsWindow(ppp(window=W),cellwidth=cwinit)           
    plot(NULL,xlim=range(OW$xvals),ylim=range(OW$yvals),xlab="",ylab="",asp=1)
    plot(W,add=TRUE)
    points(expand.grid(OW$xvals,OW$yvals),pch="+",cex=0.5)
    title(sub=paste("Output Grid: M=",length(OW$xvals),"N=",length(OW$yvals)))  
}

##' priorpost function
##'
##' A function to plot the prior and posterior densities of the model parameters eta and beta. The prior appears as a red line
##' and the posterior appears as a histogram.
##'
##' @param obj an object produced by a call to lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars or lgcpPredictMultitypeSpatialPlusPars
##' @param breaks "breaks" paramter from the function "hist"
##' @param xlab optional label for x-axis, there is a sensible default.
##' @param ylab optional label for y-axis, there is a sensible default.
##' @param main optional title of the plot, there is a sensible default.
##' @param ask the paramter "ask", see ?par
##' @param ... other arguments passed to the function "hist"
##' @return plots of the prior and posterior of the model parameters eta and beta.
##' @seealso \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export

priorpost <- function(obj,breaks=30,xlab=NULL,ylab="Density",main="",ask=TRUE,...){

    XLAB <- xlab

    par(ask=ask)
    
    nms <- c("sigma","phi","theta")
    
    if(!inherits(obj$priors,"list")){
        if(inherits(obj$priors$etaprior,"LogGaussianPrior")){
            for (i in 1:length(obj$priors$etaprior$mean)){             
                rg <- diff(range(exp(obj$etarec[,i])))           
                denx <- seq(min(exp(obj$etarec[,i]))-rg/2,max(exp(obj$etarec[,i]))+rg/2,length.out=1000)
                if(is.null(XLAB)){
                    xlab <- parse(text=paste(nms[i],sep=""))
                }
                hist(exp(obj$etarec[,i]),breaks=breaks,freq=FALSE,xlim=range(denx),xlab=xlab,ylab=ylab,main=main,...)
                lines(denx,dlnorm(denx,obj$priors$etaprior$mean[i],sqrt(obj$priors$etaprior$variance[i,i])),col="red",lwd=2)
            }
        }
        else{
            stop("prior for eta currently unsupported by priorpost")
        }
        
        if(inherits(obj$priors$betaprior,"GaussianPrior")){
            for (i in 1:length(obj$priors$betaprior$mean)){
                rg <- diff(range(obj$betarec[,i]))            
                denx <- seq(min(obj$betarec[,i])-rg/2,max(obj$betarec[,i])+rg/2,length.out=1000)
                if(is.null(XLAB)){
                    xlab <- parse(text=paste("beta[",names(coefficients(obj$glmfit))[i],"]",sep=""))
                }
                hist(obj$betarec[,i],breaks=breaks,freq=FALSE,xlim=range(denx),xlab=xlab,ylab=ylab,main=main,...)
                lines(denx,dnorm(denx,obj$priors$betaprior$mean[i],sqrt(obj$priors$betaprior$variance[i,i])),col="red",lwd=2)
            }
        }
        else{
            stop("prior for beta currently unsupported by priorpost")
        }
    }
    else{ # multivariate
        etaidx <- 1        
        for(k in 1:length(obj$priors)){
            if(inherits(obj$priors[[k]]$etaprior,"LogGaussianPrior")){
                for (i in 1:length(obj$priors[[k]]$etaprior$mean)){             
                    rg <- diff(range(exp(obj$etarec[,etaidx])))           
                    denx <- seq(min(exp(obj$etarec[,etaidx]))-rg/2,max(exp(obj$etarec[,etaidx]))+rg/2,length.out=1000)
                    if(is.null(XLAB)){
                        xlab <- parse(text=paste(nms[i],"[",k,"]",sep=""))
                    }
                    hist(exp(obj$etarec[,etaidx]),breaks=breaks,freq=FALSE,xlim=range(denx),xlab=xlab,ylab=ylab,main=main,...)
                    lines(denx,dlnorm(denx,obj$priors[[k]]$etaprior$mean[i],sqrt(obj$priors[[k]]$etaprior$variance[i,i])),col="red",lwd=2)
                    etaidx <- etaidx + 1
                }
            }
            else{
                stop("prior for eta currently unsupported by priorpost")
            }
        }
        
        betaidx <- 1     
        for(k in 1:(length(obj$priors)-1)){    
            if(inherits(obj$priors[[k]]$betaprior,"GaussianPrior")){
                for (i in 1:length(obj$priors[[k]]$betaprior$mean)){
                    rg <- diff(range(obj$betarec[,betaidx]))            
                    denx <- seq(min(obj$betarec[,betaidx])-rg/2,max(obj$betarec[,betaidx])+rg/2,length.out=1000)
                    if(is.null(XLAB)){
                        xlab <- parse(text=paste("beta[(",k,")][",names(coefficients(obj$glmfit[[k]]))[i],"]",sep=""))
                    }
                    hist(obj$betarec[,betaidx],breaks=breaks,freq=FALSE,xlim=range(denx),xlab=xlab,ylab=ylab,main=main,...)
                    lines(denx,dnorm(denx,obj$priors[[k]]$betaprior$mean[i],sqrt(obj$priors[[k]]$betaprior$variance[i,i])),col="red",lwd=2)
                    betaidx <- betaidx + 1
                }
            }
            else{
                stop("prior for beta currently unsupported by priorpost")
            }
        }
    }   
    
}

##' transblue function
##'
##' A function to return a transparent blue colour.
##'
##' @param alpha transparency parameter, see ?rgb 
##' @return character string of colour
##' @export

transblue <- function(alpha=0.1){
    return(rgb(0,0,1,alpha=alpha))
}


##' transgreen function
##'
##' A function to return a transparent green colour.
##'
##' @param alpha transparency parameter, see ?rgb  
##' @return character string of colour
##' @export

transgreen <- function(alpha=0.1){
    return(rgb(0,1,0,alpha=alpha))
}


##' transred function
##'
##' A function to return a transparent red colour.
##'
##' @param alpha transparency parameter, see ?rgb  
##' @return character string of colour
##' @export

transred <- function(alpha=0.1){
    return(rgb(1,0,0,alpha=alpha))
}


##' transblack function
##'
##' A function to return a transparent black colour.
##'
##' @param alpha transparency parameter, see ?rgb 
##' @return character string of colour
##' @export

transblack <- function(alpha=0.1){
    return(rgb(0,0,0,alpha=alpha))
}


##' osppp2latlon function
##'
##' A function to transform a ppp object in the OSGB projection (epsg:27700) to a ppp object in the latitude/longitude (epsg:4326) projection. 
##'
##' @param obj a ppp object in OSGB
##' @return a pppobject in Lat/Lon
##' @export

osppp2latlon <- function(obj){
    win <- as(obj$window,"SpatialPolygons")
    proj4string(win) <- "+init=epsg:27700"
    win <- spTransform(win,CRS("+init=epsg:4326"))
    pts <- SpatialPoints(cbind(obj$x,obj$y),proj4string=CRS("+init=epsg:27700"))
    pts <- coordinates(spTransform(pts,CRS("+init=epsg:4326")))
    spatstat.options(checkpolygons=FALSE)
    win <- as(win,"owin")
    spatstat.options(checkpolygons=TRUE)
    return(ppp(x=pts[,1],y=pts[,2],window=win,marks=obj$marks))
}


##' osppp2merc function
##'
##' A function to transform a ppp object in the OS GB projection (epsg:27700) to a ppp object in the Mercator (epsg:3857) projection. 
##'
##' @param obj a ppp object in OSGB 
##' @return a ppp object in Mercator
##' @export

osppp2merc <- function(obj){
    win <- as(obj$window,"SpatialPolygons")
    proj4string(win) <- "+init=epsg:27700"
    win <- spTransform(win,CRS("+init=epsg:3857"))
    pts <- SpatialPoints(cbind(obj$x,obj$y),proj4string=CRS("+init=epsg:27700"))
    pts <- coordinates(spTransform(pts,CRS("+init=epsg:3857")))
    spatstat.options(checkpolygons=FALSE)
    win <- as(win,"owin")
    spatstat.options(checkpolygons=TRUE)
    return(ppp(x=pts[,1],y=pts[,2],window=win,marks=obj$marks))
}



##' numCases function
##'
##' A function used in conjunction with the function "expectation" to compute the expected number of cases in each computational grid cell. Currently
##' only implemented for spatial processes (lgcpPredictSpatialPlusPars and lgcpPredictAggregateSpatialPlusPars).
##'
##' @param Y the latent field 
##' @param beta the main effects 
##' @param eta the parameters of the latent field
##' @param Z the design matrix 
##' @param otherargs other arguments to the function (see vignette "Bayesian_lgcp" for an explanation)
##' @return the number of cases in each cell
##' @seealso \link{expectation},  \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}
##' @examples
##' \dontrun{ex <- expectation(lg,numCases)[[1]] # lg is output from spatial LGCP MCMC}
##' @export

numCases <- function(Y,beta,eta,Z,otherargs){
    ca <- diff(otherargs$mcens[1:2])*diff(otherargs$ncens[1:2])
    return(ca*otherargs$poisson.offset[1:otherargs$M,1:otherargs$N]*exp(matrix(Z%*%t(beta),otherargs$M*otherargs$ext,otherargs$N*otherargs$ext)[1:otherargs$M,1:otherargs$N] + Y))
}


##' covEffects function
##'
##' A function used in conjunction with the function "expectation" to compute the main covariate effects,\cr
##'     lambda(s) exp[Z(s)beta] \cr
##' in each computational grid cell. Currently
##' only implemented for spatial processes (lgcpPredictSpatialPlusPars and lgcpPredictAggregateSpatialPlusPars).
##'
##' @param Y the latent field 
##' @param beta the main effects 
##' @param eta the parameters of the latent field
##' @param Z the design matrix 
##' @param otherargs other arguments to the function (see vignette "Bayesian_lgcp" for an explanation)
##' @return the main effects
##' @seealso \link{expectation},  \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}
##' @examples
##' \dontrun{ex <- expectation(lg,covEffects)[[1]] # lg is output from spatial LGCP MCMC}
##' @export

covEffects <- function(Y,beta,eta,Z,otherargs){
    ca <- diff(otherargs$mcens[1:2])*diff(otherargs$ncens[1:2])
    return(ca*otherargs$poisson.offset[1:otherargs$M,1:otherargs$N]*exp(matrix(Z%*%t(beta),otherargs$M*otherargs$ext,otherargs$N*otherargs$ext)[1:otherargs$M,1:otherargs$N]))
}


##' traceplots function
##'
##' A function to produce trace plots for the paramerers beta and eta from a call to the function lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars or lgcpPredictMultitypeSpatialPlusPars
##'
##' @param obj an object produced by a call to lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars orlgcpPredictMultitypeSpatialPlusPars
##' @param xlab optional label for x-axis, there is a sensible default.
##' @param ylab optional label for y-axis, there is a sensible default.
##' @param main optional title of the plot, there is a sensible default.
##' @param ask the paramter "ask", see ?par
##' @param ... other arguments passed to the function "hist"
##' @return produces MCMC trace plots of the parameters beta and eta
##' @seealso \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export

traceplots <- function(obj,xlab="Sample No.",ylab=NULL,main="",ask=TRUE,...){

    YLAB <- ylab

    par(ask=ask)
    
    nms <- c("sigma","phi","theta")
    
    if(!inherits(obj$priors,"list")){
        if(inherits(obj$priors$etaprior,"LogGaussianPrior")){
            for (i in 1:length(obj$priors$etaprior$mean)){
                if(is.null(YLAB)){
                    ylab <- parse(text=paste(nms[i],sep=""))
                }
                plot(exp(obj$etarec[,i]),type="s",xlab=xlab,ylab=ylab,main=main,...)
            }
        }
        else{
            stop("prior for eta currently unsupported by traceplots")
        }
        
        if(inherits(obj$priors$betaprior,"GaussianPrior")){
            for (i in 1:length(obj$priors$betaprior$mean)){
                if(is.null(YLAB)){
                    ylab <- parse(text=paste("beta[",names(coefficients(obj$glmfit))[i],"]",sep=""))
                }
                plot(obj$betarec[,i],type="s",xlab=xlab,ylab=ylab,main=main,...)
            }
        }
        else{
            stop("prior for beta currently unsupported by traceplots")
        }
    }
    else{ # multivariate
        etaidx <- 1        
        for(k in 1:length(obj$priors)){
            if(inherits(obj$priors[[k]]$etaprior,"LogGaussianPrior")){
                for (i in 1:length(obj$priors[[k]]$etaprior$mean)){             
                    if(is.null(YLAB)){
                        ylab <- parse(text=paste(nms[i],"[",k,"]",sep=""))
                    }
                    plot(exp(obj$etarec[,etaidx]),type="s",xlab=xlab,ylab=ylab,main=main,...)
                    etaidx <- etaidx + 1
                }
            }
            else{
                stop("prior for eta currently unsupported by traceplots")
            }
        }
        
        betaidx <- 1     
        for(k in 1:(length(obj$priors)-1)){    
            if(inherits(obj$priors[[k]]$betaprior,"GaussianPrior")){
                for (i in 1:length(obj$priors[[k]]$betaprior$mean)){
                    if(is.null(YLAB)){
                        ylab <- parse(text=paste("beta[(",k,")][",names(coefficients(obj$glmfit[[k]]))[i],"]",sep=""))
                    }
                    plot(obj$betarec[,betaidx],type="s",xlab=xlab,ylab=ylab,main=main,...)
                    betaidx <- betaidx + 1
                }
            }
            else{
                stop("prior for beta currently unsupported by traceplots")
            }
        }
    }   
    
}


##' parautocorr function
##'
##' A function to produce autocorrelation plots for the paramerers beta and eta from a call to the function lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars or lgcpPredictMultitypeSpatialPlusPars
##'
##' @param obj an object produced by a call to lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars orlgcpPredictMultitypeSpatialPlusPars
##' @param xlab optional label for x-axis, there is a sensible default.
##' @param ylab optional label for y-axis, there is a sensible default.
##' @param main optional title of the plot, there is a sensible default.
##' @param ask the paramter "ask", see ?par
##' @param ... other arguments passed to the function "hist"
##' @return produces autocorrelation plots of the parameters beta and eta
##' @seealso \link{ltar}, \link{autocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export

parautocorr <- function(obj,xlab="Lag",ylab=NULL,main="",ask=TRUE,...){

    YLAB <- ylab

    par(ask=ask)
    
    nms <- c("log(sigma)","log(phi)","log(theta)")
    
    if(!inherits(obj$priors,"list")){
        if(inherits(obj$priors$etaprior,"LogGaussianPrior")){
            for (i in 1:length(obj$priors$etaprior$mean)){
                if(is.null(YLAB)){
                    ylab <- parse(text=paste(nms[i],sep=""))
                }
                acf(obj$etarec[,i],xlab=xlab,ylab=ylab,main=main,...)
            }
        }
        else{
            stop("prior for eta currently unsupported by parautocorr")
        }
        
        if(inherits(obj$priors$betaprior,"GaussianPrior")){
            for (i in 1:length(obj$priors$betaprior$mean)){
                if(is.null(YLAB)){
                    ylab <- parse(text=paste("beta[",names(coefficients(obj$glmfit))[i],"]",sep=""))
                }
                acf(obj$betarec[,i],xlab=xlab,ylab=ylab,main=main,...)
            }
        }
        else{
            stop("prior for beta currently unsupported by parautocorr")
        }
    }
    else{ # multivariate
        etaidx <- 1        
        for(k in 1:length(obj$priors)){
            if(inherits(obj$priors[[k]]$etaprior,"LogGaussianPrior")){
                for (i in 1:length(obj$priors[[k]]$etaprior$mean)){             
                    if(is.null(YLAB)){
                        ylab <- parse(text=paste(nms[i],"[",k,"]",sep=""))
                    }
                    acf(obj$etarec[,etaidx],xlab=xlab,ylab=ylab,main=main,...)
                    etaidx <- etaidx + 1
                }
            }
            else{
                stop("prior for eta currently unsupported by priorpost")
            }
        }
        
        betaidx <- 1     
        for(k in 1:(length(obj$priors)-1)){    
            if(inherits(obj$priors[[k]]$betaprior,"GaussianPrior")){
                for (i in 1:length(obj$priors[[k]]$betaprior$mean)){
                    if(is.null(YLAB)){
                        ylab <- parse(text=paste("beta[(",k,")][",names(coefficients(obj$glmfit[[k]]))[i],"]",sep=""))
                    }
                    acf(obj$betarec[,betaidx],xlab=xlab,ylab=ylab,main=main,...)
                    betaidx <- betaidx + 1
                }
            }
            else{
                stop("prior for beta currently unsupported by parautocorr")
            }
        }
    }   
    
}



##' parsummary function
##'
##' A function to produce a summary table for the parameters beta and eta from a call to the function lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars or lgcpPredictMultitypeSpatialPlusPars
##'
##' @param obj an object produced by a call to lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars orlgcpPredictMultitypeSpatialPlusPars
##' @param expon whether to exponentiate the results, so that the parameters beta haev the interpretation of "relative risk per unit increase in the covariate" default is TRUE
##' @param LaTeX whether to print paramter names using LaTeX symbols (if the table is later to be exported to a LaTeX document)
##' @param ... other arguments
##' @return a data frame containing the median, 0.025 and 0.975 quantiles.
##' @seealso \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export

parsummary <- function(obj,expon=TRUE,LaTeX=FALSE,...){

    mode <- "spatial"
    
    if(expon){
        if(LaTeX){
            nms <- c("$\\sigma$","$\\phi$","$\\theta$")
            if(inherits(obj,"lgcpPredictMultitypeSpatialPlusParameters")){
                nms <- c("$\\sigma","$\\phi","$\\theta")
            }
        }
        else{
            nms <- c("sigma","phi","theta")
        }
        latexnms <- c("$\\sigma$","$\\phi$","$\\theta$")
    }
    else{
        if(LaTeX){
            nms <- c("$\\log(\\sigma)$","$\\log(\\phi)$","$\\log(\\theta)$")
            if(inherits(obj,"lgcpPredictMultitypeSpatialPlusParameters")){
                nms <- c("$\\log(\\sigma","$\\log(\\phi","$\\log(\\theta")
            }
        }
        else{
            nms <- c("log(sigma)","log(phi)","log(theta)")
        } 
        latexnms <- nms <- c("$\\log(\\sigma)$","$\\log(\\phi)$","$\\log(\\theta)$")        
    }    
    cnams <- c("median" , "lower 95\\% CRI" , "upper 95\\% CRI")
    
    
    fun <- function(x){
        return(quantile(x,c(0.5,0.025,0.975)))
    }
    
    verbify <- function(x){
        return(paste("\\verb=",x,"=",sep=""))
    }
        
    info <- c()
    rnams <- c() 
    
    parnames <- c()
    modno <- c()
    
    if(!inherits(obj$priors,"list")){
        if(inherits(obj$priors$etaprior,"LogGaussianPrior")){
            for (i in 1:length(obj$priors$etaprior$mean)){
                ylab <- paste(nms[i],sep="")            
                rnams <- c(rnams,ylab)
                if(expon){
                    info <- rbind(info,fun(exp(obj$etarec[,i])))
                }
                else{
                    info <- rbind(info,fun(obj$etarec[,i]))
                }
                parnames <- c(parnames,latexnms[i])
                if(i==3){
                    mode <- "spatiotemporal"
                }
            }
        }
        else{
            stop("prior for eta currently unsupported by parsummary")
        }
        
        if(inherits(obj$priors$betaprior,"GaussianPrior")){
            for (i in 1:length(obj$priors$betaprior$mean)){                
                if(expon){
                    if(LaTeX){
                        ylab <- paste("$\\exp(\\beta_{",names(coefficients(obj$glmfit))[i],"})$",sep="")
                    }
                    else{
                        ylab <- paste("exp(beta_",names(coefficients(obj$glmfit))[i],")",sep="")
                    }
                    parnames <- c(parnames,verbify(names(coefficients(obj$glmfit))[i]))
                    rnams <- c(rnams,ylab)
                    info <- rbind(info,fun(exp(obj$betarec[,i])))
                }
                else{
                    if(LaTeX){
                        ylab <- paste("$\\beta_{",names(coefficients(obj$glmfit))[i],"}$",sep="")
                    }
                    else{
                        ylab <- paste("beta_",names(coefficients(obj$glmfit))[i],sep="")
                    }
                    parnames <- c(parnames,verbify(names(coefficients(obj$glmfit))[i]))
                    rnams <- c(rnams,ylab)
                    info <- rbind(info,fun(obj$betarec[,i]))
                }
            }
        }
        else{
            stop("prior for beta currently unsupported by parsummary")
        }
    }
    else{ # multivariate
        mode <- "multivariate"
        etaidx <- 1        
        for(k in 1:length(obj$priors)){
            if(inherits(obj$priors[[k]]$etaprior,"LogGaussianPrior")){
                for (i in 1:length(obj$priors[[k]]$etaprior$mean)){
                    ylab <- paste(nms[i],"_",k,ifelse(expon,"$",")$"),sep="")
                    rnams <- c(rnams,ylab)
                    if(expon){
                        info <- rbind(info,fun(exp(obj$etarec[,etaidx])))
                    }
                    else{
                        info <- rbind(info,fun(obj$etarec[,etaidx]))
                    }
                    parnames <- c(parnames,latexnms[i])
                    etaidx <- etaidx + 1
                    modno <- c(modno,k)
                }
            }
            else{
                stop("prior for eta currently unsupported by parsummary")
            }
        }
        
        betaidx <- 1     
        for(k in 1:(length(obj$priors)-1)){    
            if(inherits(obj$priors[[k]]$betaprior,"GaussianPrior")){
                for (i in 1:length(obj$priors[[k]]$betaprior$mean)){                    
                    if(expon){
                        if(LaTeX){
                            ylab <- paste("$\\exp[\\beta(",k,")(",names(coefficients(obj$glmfit[[k]]))[i],")]$",sep="")
                        }
                        else{
                            ylab <- paste("exp[beta(",k,")(",names(coefficients(obj$glmfit[[k]]))[i],")]",sep="")
                        }
                        parnames <- c(parnames,verbify(names(coefficients(obj$glmfit[[k]]))[i]))
                        rnams <- c(rnams,ylab)
                        info <- rbind(info,fun(exp(obj$betarec[,betaidx])))
                    }
                    else{
                        if(LaTeX){
                            ylab <- paste("$\\beta(",k,")(",names(coefficients(obj$glmfit[[k]]))[i],")$",sep="")
                        }
                        else{
                            ylab <- paste("beta(",k,")(",names(coefficients(obj$glmfit[[k]]))[i],")",sep="")
                        }
                        parnames <- c(parnames,verbify(names(coefficients(obj$glmfit[[k]]))[i]))
                        rnams <- c(rnams,ylab)
                        info <- rbind(info,fun(obj$betarec[,betaidx]))
                    }
                    betaidx <- betaidx + 1
                    modno <- c(modno,k)
                }
            }
            else{
                stop("prior for beta currently unsupported by parsummary")
            }
        }
    }
    
    rownames(info) <- rnams
    colnames(info) <- cnams
    
    attr(info,"parnames") <- parnames
    attr(info,"mode") <- mode
    attr(info,"modno") <- modno     

    return(info)       
}


##' textsummary function
##'
##' A function to print a text description of the inferred paramerers beta and eta from a call to the function lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars or lgcpPredictMultitypeSpatialPlusPars
##'
##' @param obj an object produced by a call to lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars orlgcpPredictMultitypeSpatialPlusPars 
##' @param digits see the option "digits" in ?format 
##' @param scientific see the option "scientific" in ?format
##' @param inclIntercept logical: whether to summarise the intercept term, default is FALSE.
##' @param ... other arguments passed to the function "format" 
##' @return A text summary, that can be pasted into a LaTeX document and later edited.
##' @seealso \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export

textsummary <- function(obj,digits=3,scientific=-3,inclIntercept=FALSE,...){
    psu <- parsummary(obj)
    parnames <- attr(psu,"parnames") 
    mode <- attr(psu,"mode")
    modno <- attr(psu,"modno")

    form <- function(x,...){
        xtxt <- format(x,digits=digits,scientific=scientific,...)
        if(length(grep("e",xtxt))>0){
            spl <- unlist(strsplit(xtxt,"e"))
            xtxt <- paste(spl[1],"$\\times10^{",as.character(as.numeric(spl[2])),"}$",sep="")
        }    
        return(xtxt)
    }
    
    parvals <- matrix(sapply(psu,form),nrow=nrow(psu),ncol=ncol(psu)) 
    
    sigf <- apply(psu[,2:3],1,function(x){ifelse((x[1]<1 & x[2]<1)|(x[1]>1 & x[2]>1),TRUE,FALSE)})
    
       
     
    nspar <- 2
    if(mode=="spatiotemporal"){
        nspar <- 3
    }
    np <- nrow(psu) - nspar
    
    if(nspar==2){
        lfsen <- paste("A summary of the parameters of the latent field is as follows. The parameter ",parnames[1]," had median ",
                        parvals[1,1]," (95\\% CRI ",parvals[1,2]," to ",parvals[1,3],") and the parameter ",parnames[2]," had median ",
                        parvals[2,1]," (95\\% CRI ",parvals[2,2]," to ",parvals[2,3],").",sep="")
    }
    else if(nspar==3){
        lfsen <- paste("A summary of the parameters of the latent field is as follows. The parameter ",parnames[1]," had median ",
                        parvals[1,1]," (95\\% CRI ",parvals[1,2]," to ",parvals[1,3],"); the parameter ",parnames[2]," had median ",
                        parvals[2,1]," (95\\% CRI ",parvals[2,2]," to ",parvals[2,3],"); and the parameter ",parnames[3]," had median ",
                        parvals[3,1]," (95\\% CRI ",parvals[3,2]," to ",parvals[3,3],").",sep="")
    }
    
    
    tab <- parvals[(nspar+1):nrow(parvals),]
    tsi <- sigf[(nspar+1):nrow(parvals)]
    efd <- sapply(psu[,1],function(x){ifelse(x>1,1,-1)})[(nspar+1):nrow(parvals)] 
    prn <- parnames[(nspar+1):nrow(parvals)]
    
    if(!inclIntercept){
        if(any(prn=="\\verb=(Intercept)=")){
            idx <- which(prn=="\\verb=(Intercept)=")
            tab <- tab[-idx,,drop=FALSE]
            tsi <- tsi[-idx] # significant?
            efd <- efd[-idx] # effect direction
            prn <- prn[-idx]
        }
    }
    
    sumrow <- function(infor){
        i <- infor[1]
        dir <- infor[2]
        direc <- "increase"
        if(dir==-1){
            direc <- "reduction"
        }
        return(paste("each unit increase in ",pr[i]," led to a ",direc," in relative risk with median ",tt[i,1]," (95\\% CRI ",tt[i,2]," to ",tt[i,3],")",sep=""))
    }
    
    pardesc <- c()
    pardesc2 <- c()
    
    if(all(tsi)){
        pardesc <- "All of the main effects were found to be significant: "
        tt <- tab[tsi,drop=FALSE,]
        ef <- efd[tsi]
        pr <- prn[tsi]
        ord <- order(ef)
        tt <- tt[ord,,drop=FALSE]
        ef <- ef[ord]
        pr <- pr[ord]
        apply(cbind(1:nrow(tt),ef),1,function(x){pardesc <<- paste(pardesc,sumrow(x),ifelse(x[1]==nrow(tt),". ","; "),sep="")})
    }
    else{
        if(any(tsi)){
            pardesc <- "The following effects were found to be significant: "
            tt <- tab[tsi,,drop=FALSE]
            ef <- efd[tsi]
            pr <- prn[tsi]
            ord <- order(ef)
            tt <- tt[ord,,drop=FALSE]
            ef <- ef[ord]
            pr <- pr[ord]
            apply(cbind(1:nrow(tt),ef),1,function(x){pardesc <<- paste(pardesc,sumrow(x),ifelse(x[1]==nrow(tt),". ","; "),sep="")})
        }
        
        
        
        if(all(!tsi)){
            pardesc2 <- "None of the main effects were found to be significant: "
        }
        else{
            pardesc2 <- "The remainder of the main effects were not found to be significant: "
        }
        tt <- tab[!tsi,,drop=FALSE]
        ef <- efd[!tsi]
        pr <- prn[!tsi]
        ord <- order(ef)
        tt <- tt[ord,,drop=FALSE]
        ef <- ef[ord]
        pr <- pr[ord]
        apply(cbind(1:nrow(tt),ef),1,function(x){pardesc2 <<- paste(pardesc2,sumrow(x),ifelse(x[1]==nrow(tt),". ","; "),sep="")})
    }      
    #browser()

    cat("\n")
    cat(lfsen)
    cat("\n")
    cat("\n")
    cat(pardesc)
    cat("\n")
    cat("\n")
    cat(pardesc2)
    cat("\n")

}





##' etavals function
##'
##' A function to return the sampled eta from a call to the function lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars or lgcpPredictMultitypeSpatialPlusPars
##'
##' @param lg an object produced by a call to lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars orlgcpPredictMultitypeSpatialPlusPars
##' @return the posterior sampled eta
##' @seealso \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}
##' @export
etavals <- function(lg){
    return(lg$etarec)
}


##' betavals function
##'
##' A function to return the sampled beta from a call to the function lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars or lgcpPredictMultitypeSpatialPlusPars 
##'
##' @param lg an object produced by a call to lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars orlgcpPredictMultitypeSpatialPlusPars 
##' @return the posterior sampled beta
##' @seealso \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{etavals}
##' @export
betavals <- function(lg){
    return(lg$betarec)
}

##' ltar function
##'
##' A function to return the sampled log-target from a call to the function lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, 
##' lgcpPredictSpatioTemporalPlusPars or lgcpPredictMultitypeSpatialPlusPars. This is used as a convergence diagnostic.
##'
##' @param lg an object produced by a call to lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars orlgcpPredictMultitypeSpatialPlusPars 
##' @return the log-target from each saved iteration of the MCMC chain.
##' @seealso \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export
ltar <- function(lg){
    return(lg$tarrec)
}



##' postcov function
##'
##' Generic function for producing plots of the posterior covariance function from a call to the function lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars,
##' lgcpPredictSpatioTemporalPlusPars or lgcpPredictMultitypeSpatialPlusPars.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method postcov
##' @seealso
##' \link{postcov.lgcpPredictSpatialOnlyPlusParameters},\link{postcov.lgcpPredictAggregateSpatialPlusParameters}, \link{postcov.lgcpPredictSpatioTemporalPlusParameters}, \link{postcov.lgcpPredictMultitypeSpatialPlusParameters},  
##'\link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{exceedProbs}, \link{betavals}, \link{etavals}  
##' @export

postcov <- function(obj,...){
    UseMethod("postcov")
}



##' postcov.lgcpPredictSpatioTemporalPlusParameters function
##'
##' A function for producing plots of the posterior spatiotemporal covariance function. 
##'
##' @method postcov lgcpPredictSpatioTemporalPlusParameters
##' @param obj an lgcpPredictSpatioTemporalPlusParameters object
##' @param qts vector of quantiles of length 3, default is 0.025, 0.5, 0.975
##' @param covmodel the assumed covariance model. NULL by default, this information is read in from the object obj, so generally does not need to be set.
##' @param ask parameter "ask", see ?par
##' @param ... additional arguments
##' @return a plot of the posterior spatial covariance function and temporal correlation function.
##' @usage "postcov(obj,qts=c(0.025,0.5,0.975),covmodel=NULL,ask=TRUE,...)"
##' @seealso \link{postcov.lgcpPredictSpatialOnlyPlusParameters}, \link{postcov.lgcpPredictAggregateSpatialPlusParameters}, \link{postcov.lgcpPredictSpatioTemporalPlusParameters}, \link{postcov.lgcpPredictMultitypeSpatialPlusParameters},
##' \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export

postcov.lgcpPredictSpatioTemporalPlusParameters <- function(obj,qts=c(0.025,0.5,0.975),covmodel=NULL,ask=TRUE,...){
    postcov.lgcpPredictSpatialOnlyPlusParameters(obj=obj,qts=qts,covmodel=covmodel,ask=ask,...)   
   
    xrg <- seq(0,diff(range(obj$aggtimes))+1,length.out=100)
    n <- nrow(obj$etarec) 
    if(length(qts)!=3){
        stop("Argument qts must have length 3")
    }
    qts <- sort(qts)
    postc <- t(sapply(1:n,function(i){exp(-xrg*exp(obj$etarec[i,3]))}))
    qtil <- apply(postc,2,quantile,qts)
    plot(xrg,qtil[1,],type="l",ylim=c(0,max(qtil[3,],na.rm=TRUE)),xlab="Time",ylab="Correlation",lty="dashed")
    lines(xrg,qtil[2,])
    lines(xrg,qtil[3,],lty="dashed")
    legend("topright",lty=c("solid","dashed"),legend=c(paste(qts[2],"quantile"),paste(qts[1],"-",qts[3]," CRI",sep="")))
}



##' postcov.lgcpPredictSpatialOnlyPlusParameters function
##'
##' A function for producing plots of the posterior spatial covariance function. 
##'
##' @method postcov lgcpPredictSpatialOnlyPlusParameters
##' @param obj an lgcpPredictSpatialOnlyPlusParameters object
##' @param qts vector of quantiles of length 3, default is 0.025, 0.5, 0.975
##' @param covmodel the assumed covariance model. NULL by default, this information is read in from the object obj, so generally does not need to be set.
##' @param ask parameter "ask", see ?par
##' @param ... additional arguments
##' @return a plot of the posterior covariance function.
##' @usage "postcov(obj,qts=c(0.025,0.5,0.975),covmodel=NULL,ask=TRUE,...)"
##' @seealso \link{postcov.lgcpPredictSpatialOnlyPlusParameters}, \link{postcov.lgcpPredictAggregateSpatialPlusParameters}, \link{postcov.lgcpPredictSpatioTemporalPlusParameters}, \link{postcov.lgcpPredictMultitypeSpatialPlusParameters},
##' \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export

postcov.lgcpPredictSpatialOnlyPlusParameters <- function(obj,qts=c(0.025,0.5,0.975),covmodel=NULL,ask=TRUE,...){
    par(ask=ask)
    if(is.null(covmodel)){
        covmodel <- obj$covFct
        if(is.null(covmodel)){
            covmodel <- exponentialCovFct
            warning("Exponential Covariance function assumed.",immediate.=TRUE)    
        }
    }
    xrg <- seq(0,3*max(exp(obj$etarec[,2])),length.out=100)
    n <- nrow(obj$etarec) 
    if(length(qts)!=3){
        stop("Argument qts must have length 3")
    }
    qts <- sort(qts)
    postc <- t(sapply(1:n,function(i){covmodel(xrg,CovParameters=list(sigma=exp(obj$etarec[i,1]),phi=exp(obj$etarec[i,2])))$eval}))
    qtil <- apply(postc,2,quantile,qts)
    plot(xrg,qtil[1,],type="l",ylim=c(0,max(qtil[3,],na.rm=TRUE)),xlab="Range",ylab="Covariance",lty="dashed")
    lines(xrg,qtil[2,])
    lines(xrg,qtil[3,],lty="dashed")
    legend("topright",lty=c("solid","dashed"),legend=c(paste(qts[2],"quantile"),paste(qts[1],"-",qts[3]," CRI",sep="")))

}



##' postcov.lgcpPredictMultitypeSpatialPlusParameters function
##'
##' A function for producing plots of the posterior covariance function. 
##'
##' @method postcov lgcpPredictMultitypeSpatialPlusParameters
##' @param obj an lgcpPredictMultitypeSpatialPlusParameters object
##' @param qts vector of quantiles of length 3, default is 0.025, 0.5, 0.975
##' @param covmodel the assumed covariance model. NULL by default, this information is read in from the object obj, so generally does not need to be set.
##' @param ask parameter "ask", see ?par
##' @param ... additional arguments
##' @return plots of the posterior covariance function for each type.
##' @usage "postcov(obj,qts=c(0.025,0.5,0.975),covmodel=NULL,ask=TRUE,...)"
##' @seealso \link{postcov.lgcpPredictSpatialOnlyPlusParameters}, \link{postcov.lgcpPredictAggregateSpatialPlusParameters}, \link{postcov.lgcpPredictSpatioTemporalPlusParameters}, \link{postcov.lgcpPredictMultitypeSpatialPlusParameters},
##' \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export

postcov.lgcpPredictMultitypeSpatialPlusParameters <- function(obj,qts=c(0.025,0.5,0.975),covmodel=NULL,ask=TRUE,...){
    par(ask=ask)
    plotfun <- function(i,obj,qts,covmodel,ask,...){
        if(is.null(covmodel)){
            covmodel <- obj$covFct[[i]]
            if(is.null(covmodel)){
                covmodel <- exponentialCovFct
                warning("Exponential Covariance function assumed.",immediate.=TRUE)    
            }
        }
        sidx <- 2*i-1
        pidx <- 2*i
        xrg <- seq(0,3*max(exp(obj$etarec[,pidx])),length.out=100)
        n <- nrow(obj$etarec) 
        if(length(qts)!=3){
            stop("Argument qts must have length 3")
        }
        qts <- sort(qts)
        postc <- t(sapply(1:n,function(i){covmodel(xrg,CovParameters=list(sigma=exp(obj$etarec[i,sidx]),phi=exp(obj$etarec[i,pidx])))$eval}))
        qtil <- apply(postc,2,quantile,qts)
        plot(xrg,qtil[1,],type="l",ylim=c(0,max(qtil[3,],na.rm=TRUE)),xlab="Range",ylab="Covariance",lty="dashed")
        lines(xrg,qtil[2,])
        lines(xrg,qtil[3,],lty="dashed")
        legend("topright",lty=c("solid","dashed"),legend=c(paste(qts[2],"quantile"),paste(qts[1],"-",qts[3]," CRI",sep="")))
    }
    
    sapply(1:length(obj$covFct),plotfun,obj=obj,qts=qts,covmodel=covmodel,ask=ask,...)
}



##' postcov.lgcpPredictAggregateSpatialPlusParameters function
##'
##' A function for producing plots of the posterior covariance function. 
##'
##' @method postcov lgcpPredictAggregateSpatialPlusParameters
##' @param obj an lgcpPredictAggregateSpatialPlusParameters object
##' @param qts vector of quantiles of length 3, default is 0.025, 0.5, 0.975
##' @param covmodel the assumed covariance model. NULL by default, this information is read in from the object obj, so generally does not need to be set.
##' @param ask parameter "ask", see ?par
##' @param ... additional arguments
##' @return ...
##' @usage "postcov(obj,qts=c(0.025,0.5,0.975),covmodel=NULL,ask=TRUE,...)"
##' @seealso \link{postcov.lgcpPredictSpatialOnlyPlusParameters}, \link{postcov.lgcpPredictAggregateSpatialPlusParameters}, \link{postcov.lgcpPredictSpatioTemporalPlusParameters}, \link{postcov.lgcpPredictMultitypeSpatialPlusParameters},
##' \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export

postcov.lgcpPredictAggregateSpatialPlusParameters <- function(obj,qts=c(0.025,0.5,0.975),covmodel=NULL,ask=TRUE,...){
    postcov.lgcpPredictSpatialOnlyPlusParameters(obj=obj,qts=qts,covmodel=covmodel,ask=ask,...)
}


##' condProbs function
##'
##' A function to compute the conditional type-probabilities from a multivariate LGCP. See the vignette "Bayesian_lgcp" for a full explanation of this.\cr
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
##' The term 'conditional probability of type k' means the probability that at a particular location there 
##' will be an event of type k, which denoted p_k.
##'
##' @param obj an lgcpPredictMultitypeSpatialPlusParameters object
##' @return an lgcpgrid object containing the consitional type-probabilities for each type
##' @seealso \link{segProbs}, \link{postcov.lgcpPredictSpatialOnlyPlusParameters}, \link{postcov.lgcpPredictAggregateSpatialPlusParameters}, \link{postcov.lgcpPredictSpatioTemporalPlusParameters}, \link{postcov.lgcpPredictMultitypeSpatialPlusParameters},
##' \link{ltar}, \link{autocorr}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export

condProbs <- function(obj){

    cins <- obj$cellInside
    Nfields <- length(obj$Zlist)
    betalenlist <- sapply(obj$Zlist,ncol)
    betaidx <- c(0,cumsum(betalenlist))
    nits <- nrow(obj$etarec)
    bottomleft <- matrix(FALSE,obj$ext*obj$M,obj$ext*obj$N)
    bottomleft[1:obj$M,1:obj$N] <- TRUE
    ZML <- lapply(obj$Zlist,function(x){x[bottomleft,,drop=FALSE]})

    getField <- function(obj,itno,fno){
        fn <- paste(obj$gridfunction$dirname,"simout.nc",sep="")
        ncdata <- nc_open(fn)
        ans <- ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(1,1,1,fno,itno), count=c(-1,-1,-1,1,1))
        nc_close(ncdata)
        return(ans)
    }
    
    getRisk <- function(obj,itno,fno){
        if(fno>Nfields){
            stop("Error in function condProbs, getRisk subroutine")
        }
        S <- getField(obj,itno,fno)
        T <- getField(obj,itno,Nfields+1)
        betacols <- (betaidx[fno]+1):(betaidx[fno+1])
        beta <- obj$betarec[itno,betacols,drop=FALSE]
        return(cins*exp(matrix(ZML[[fno]]%*%t(beta),obj$M,obj$N)+S+T))
    }
    
    pb <- txtProgressBar(min=0,max=nits,style=3)
    cprob <- array(0,dim=c(obj$M,obj$N,Nfields))
    for(i in 1:nits){
        lsum <- 0
        for(j in 1:Nfields){
           assign(paste("l",j,sep=""),getRisk(obj,i,j))
           lsum <- lsum + get(paste("l",j,sep=""))
        }
        for(j in 1:Nfields){
           cprob[,,j] <- cprob[,,j] + get(paste("l",j,sep=""))/lsum
        } 
        setTxtProgressBar(pb,i)
    }
    close(pb)
    for(j in 1:Nfields){
       cprob[,,j] <- cins*cprob[,,j]/nits
    }
    return(lgcpgrid(cprob,xvals=obj$mcens,yvals=obj$ncens))
}


##' segProbs function
##'
##' A function to compute segregation probabilities from a multivariate LGCP. See the vignette "Bayesian_lgcp" for a full explanation of this.\cr
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
##' The term 'conditional probability of type k' means the probability that at a particular location, x, there 
##' will be an event of type k, we denote this p_k(x).
##' 
##' It is also of interest to scientists to be able to illustrate spatial regions where a genotype
##' dominates a posteriori. We say that type k dominates at position x if p_k(x)>c, where c (the parameter domprob) is a 
##' threshold is a threshold set by the user. Let A_k(c,q) denote the set of locations x for which P[p_k(x)>c|X] > q.
##' 
##' As the quantities c and q tend to 1 each area A_k(c,p) shrinks towards the empty set; this
##' happens more slowly in a highly segregated pattern compared with a weakly segregated one.
##'
##' The function segProbs computes P[p_k(x)>c|X] for each type, from which plots of P[p_k(x)>c|X] > q can be produced.
##'
##' @param obj an lgcpPredictMultitypeSpatialPlusParameters object
##' @param domprob the threshold beyond which we declare a type as dominant e.g. a value of 0.8 would mean we would consider each type to be dominant if the
##' conditional probability of an event of a given type at that location exceeded 0.8. 
##' @return an lgcpgrid object contatining the segregation probabilities.
##' @export

segProbs <- function(obj,domprob){

    cins <- obj$cellInside
    Nfields <- length(obj$Zlist)
    betalenlist <- sapply(obj$Zlist,ncol)
    betaidx <- c(0,cumsum(betalenlist))
    nits <- nrow(obj$etarec)
    bottomleft <- matrix(FALSE,obj$ext*obj$M,obj$ext*obj$N)
    bottomleft[1:obj$M,1:obj$N] <- TRUE
    ZML <- lapply(obj$Zlist,function(x){x[bottomleft,,drop=FALSE]})

    getField <- function(obj,itno,fno){
        fn <- paste(obj$gridfunction$dirname,"simout.nc",sep="")
        ncdata <- nc_open(fn)
        ans <- ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(1,1,1,fno,itno), count=c(-1,-1,-1,1,1))
        nc_close(ncdata)
        return(ans)
    }
    
    getRisk <- function(obj,itno,fno){
        if(fno>Nfields){
            stop("Error in function condProbs, getRisk subroutine")
        }
        S <- getField(obj,itno,fno)
        T <- getField(obj,itno,Nfields+1)
        betacols <- (betaidx[fno]+1):(betaidx[fno+1])
        beta <- obj$betarec[itno,betacols,drop=FALSE]
        return(cins*exp(matrix(ZML[[fno]]%*%t(beta),obj$M,obj$N)+S+T))
    }
    
    pb <- txtProgressBar(min=0,max=nits,style=3)
    segrec <- array(0,dim=c(obj$M,obj$N,Nfields))
    for(i in 1:nits){
        lsum <- 0
        for(j in 1:Nfields){
           assign(paste("l",j,sep=""),getRisk(obj,i,j))
           lsum <- lsum + get(paste("l",j,sep=""))
        }        
        for(j in 1:Nfields){
           cprob <- get(paste("l",j,sep=""))/lsum
           segrec[,,j] <- segrec[,,j] + matrix(as.vector(cprob>domprob)/nits,obj$M,obj$N)   
        } 
        setTxtProgressBar(pb,i)
    }
    close(pb)

    segrec <- lgcpgrid(segrec,xvals=obj$mcens,yvals=obj$ncens)

    attr(segrec,"domprob") <- domprob
    
    return(segrec)
}

##' addTemporalCovariates function
##'
##' A function to 'bolt on' temporal data onto a spatial covariate design matrix. The function takes a spatial design matrix, Z(s) and 
##' converts it to a spatiotemporal design matrix Z(s,t) when the effects can be separably decomposed i.e.,\cr
##' Z(s,t)beta = Z_1(s)beta_1 + Z_2(t)beta_2\cr
##' \cr
##' An example of this function in action is given in the vignette "Bayesian_lgcp", in the section on spatiotemporal data.
##' 
##' The main idea of this function is: having created a spatial Z(s) using getZmat, to create a dummy dataset tdata and temporal
##' formula corresponding to the temporal component of the separable effects. The entries in the model matrix Z(s,t) corresponsing to
##' the time covariates are constant over the observation window in space, but in general vary from time-point to time-point. 
##'
##' Note that if there is an intercept in the spatial part of the model e.g., X ~ var1 + var2, then in the temporal model,
##' the intercept should be removed i.e., t ~ tvar1 + tvar2 - 1  
##'
##' @param temporal.formula a formula of the form t ~ tvar1 + tvar2 etc. Where the left hand side is a "t". Note there should 
##' not be an intercept term in both of the the spatial and temporal components.
##' @param T the time point of interest
##' @param laglength  the number of previous time points to include in the analysis
##' @param tdata a data frame with variable t minimally including times (T-laglength):T and var1, var2 etc.
##' @param Zmat the spatial covariates Z(s), obtained by using the getZmat function.
##' @return A list of design matrices, one for each time, Z(s,t) for t in (T-laglength):T
##' @seealso \link{minimum.contrast}, \link{minimum.contrast.spatiotemporal}, \link{chooseCellwidth}, \link{getpolyol}, \link{guessinterp}, \link{getZmat},
##' \link{lgcpPrior}, \link{lgcpInits}, \link{CovFunction}
##' \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}, \link{lgcpPredictSpatioTemporalPlusPars}, 
##' \link{lgcpPredictMultitypeSpatialPlusPars}
##' @export

addTemporalCovariates <- function(temporal.formula,T,laglength,tdata,Zmat){

    if(attr(terms(temporal.formula),"intercept")==1){
        stop("temporal.formula must not include an intercept term, please specify as t ~ ... - 1, where ... are the existiing independent vasriable names")
    }
    
    varn <- variablesinformula(temporal.formula)[-1] #independent variable names      
    cidx <- match(varn,names(tdata))
    cidx <- cidx[!is.na(cidx)]
    
    temporal.formula <- as.formula(paste("t ~", paste(varn,collapse=" + "))) # trick to deal with factor variables ... removed below ... if this is not done, a problem occurs later in the call to glm in MALAlgcpSpatioTemporal.PlusPars 

    if (!inherits(T,"integer")){
	    warning("Converting T into integer value, see ?as.integer",immediate.=TRUE)
	    T <- as.integer(T) 
	}
	if (!inherits(laglength,"integer")){
	    warning("Converting laglength into integer values, see ?as.integer",immediate.=TRUE)
	    laglength <- as.integer(laglength) 
	}
	aggtimes <- T - laglength:0
	
	idx <- sapply(aggtimes,match,tdata$t)
    
    tdata <- tdata[idx,]	

    dmat <- tdata
    modmat <- model.matrix(temporal.formula,data=tdata)[,-1,drop=FALSE] # now remove intercept term, see above	
    
    dmat <- dmat[,cidx,drop=FALSE]    
    
    cins <- as.vector(attr(Zmat,"cellInside"))
    
    ZmatList <- list()
    for(i in 1:length(aggtimes)){
        zmtemp <- Zmat
        zmtemp <- cbind(zmtemp,matrix(modmat[i,],nrow=nrow(Zmat),ncol=ncol(modmat),byrow=TRUE))
        colnames(zmtemp) <- c(colnames(Zmat),colnames(modmat))
        
        dmattemp <- cbind(attr(Zmat,"data.frame"),dmat[i,])
        colnames(dmattemp) <- c(colnames(attr(Zmat,"data.frame")),colnames(dmat))
        
        attr(zmtemp,"data.frame") <- dmattemp
        attr(zmtemp,"cellInside") <- attr(Zmat,"cellInside")
        attr(zmtemp,"anymiss") <- attr(Zmat,"anymiss")
        attr(zmtemp,"mcens") <- attr(Zmat,"mcens")
        attr(zmtemp,"ncens") <- attr(Zmat,"ncens")
        attr(zmtemp,"M") <- attr(Zmat,"M")
        attr(zmtemp,"N") <- attr(Zmat,"N")
        attr(zmtemp,"mcens") <- attr(Zmat,"mcens")
        attr(zmtemp,"ncens") <- attr(Zmat,"ncens")
        attr(zmtemp,"polygonOverlay") <- attr(Zmat,"polygonOverlay")
        attr(zmtemp,"pixelOverlay") <- attr(Zmat,"pixelOverlay")
        attr(zmtemp,"FORM") <- as.formula(paste("X ~ ",as.character(attr(Zmat,"FORM"))[3],paste("+",varn,collapse=" ")))
        attr(zmtemp,"gridobj") <- attr(Zmat,"gridobj")
	    attr(zmtemp,"inclusion") <- attr(Zmat,"inclusion")
	    attr(zmtemp,"ext") <- attr(Zmat,"ext")
	    attr(zmtemp,"cellwidth") <- attr(Zmat,"cellwidth") 
              
        ZmatList[[i]] <- zmtemp
    }	
	
	return(ZmatList)   
}