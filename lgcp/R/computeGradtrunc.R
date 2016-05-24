##' computeGradtruncSpatial function
##'
##' \bold{Advanced use only.} A function to compute a gradient truncation parameter for 'spatial only' MALA via simulation. The function
##' requires an FFT 'grid' to be pre-computed, see \link{fftgrid}.
##'
##' @param nsims The number of simulations to use in computation of gradient truncation.
##' @param scale multiplicative scaling constant, returned value is scale (times) max(gradient over simulations). Default scale is 1.
##' @param nis cell counts on the extended grid
##' @param mu parameter of latent field, mu 
##' @param rootQeigs root of eigenvalues of precision matrix of latent field 
##' @param invrootQeigs reciprocal root of eigenvalues of precision matrix of latent field
##' @param scaleconst expected number of cases, or ML estimate of this quantity
##' @param spatial spatial at risk interpolated onto grid of requisite size
##' @param cellarea cell area
##' @return gradient truncation parameter
##' @seealso \link{fftgrid}
##' @export 

computeGradtruncSpatial <- function(nsims=100,
                                    scale=1,
                                    nis,
                                    mu,
                                    rootQeigs,
                                    invrootQeigs,
                                    scaleconst,
                                    spatial,
                                    cellarea){
    
    cat("Computing gradient truncation ...\n")
    grds <- c()
    pb <- txtProgressBar(min=1,max=nsims,style=3)
    M <- nrow(nis)
    N <- ncol(nis)
    for (i in 1:nsims){
        Gamma <- matrix(rnorm(M*N),M,N)
        Y <- YfromGamma(Gamma=Gamma,invrootQeigs=invrootQeigs,mu=mu)
        expY <- exp(Y)                  		                    
        grds[i] <- max(expY) ####suppressWarnings(max((-1)*Gamma +(1/length(Y))*Re(fft(fft(nis-scaleconst*spatial*expY*cellarea,inverse=TRUE)*rootQeigs,inverse=TRUE))))
        setTxtProgressBar(pb,i)                                                        		                    
    }
    close(pb)
    gt <- floor(scale*max(grds,na.rm=TRUE))
    if (is.na(gt) | is.infinite(gt) | is.nan(gt)){
        stop("Could not compute gradient truncation. To set manually, see gradtrunc argument of ?lgcpPredictSpatial") 
    }
    cat(paste("Using gradient truncation of ",gt,"\n",sep="")) 
    return(gt)
}


##' computeGradtruncSpatioTemporal function
##'
##' \bold{Advanced use only.} A function to compute a gradient truncation parameter for 'spatial only' MALA via simulation. The function
##' requires an FFT 'grid' to be pre-computed, see \link{fftgrid}.
##'
##' @param nsims The number of simulations to use in computation of gradient truncation.
##' @param scale multiplicative scaling constant, returned value is scale (times) max(gradient over simulations). Default scale is 1.
##' @param nis cell counts on the extended grid
##' @param mu parameter of latent field, mu 
##' @param rootQeigs root of eigenvalues of precision matrix of latent field 
##' @param invrootQeigs reciprocal root of eigenvalues of precision matrix of latent field
##' @param spatial spatial at risk interpolated onto grid of requisite size
##' @param temporal fitted temporal values
##' @param bt vectoer of variances b(delta t) in Brix and Diggle 2001
##' @param cellarea cell area
##' @return gradient truncation parameter
##' @seealso \link{fftgrid}
##' @export 

computeGradtruncSpatioTemporal <- function( nsims=100,
                                            scale=1,
                                            nis,
                                            mu,
                                            rootQeigs,
                                            invrootQeigs,
                                            spatial,
                                            temporal,
                                            bt,
                                            cellarea){
    
    cat("Computing gradient truncation ...\n")
    grds <- c()
    pb <- txtProgressBar(min=1,max=nsims,style=3)
    M <- nrow(nis[[1]])
    N <- ncol(nis[[1]])
    for (i in 1:nsims){
        Gamma <- list()                           
        
        lapply(1:length(nis),function(i){Gamma[[i]]<<-matrix(rnorm(M*N),M,N)})
        Y <- lapply(Gamma,YfromGamma,invrootQeigs=invrootQeigs,mu=mu)
        expY <- lapply(Y,exp)                 		                    
        
        grds[i] <- suppressWarnings(max(sapply(expY,max)))        
        
        setTxtProgressBar(pb,i)                                                        		                    
    }
    close(pb)
    gt <- floor(scale*max(grds,na.rm=TRUE))
    if (is.na(gt) | is.infinite(gt) | is.nan(gt)){
        stop("Could not compute gradient truncation. To set manually, see gradtrunc argument of ?lgcpPredictSpatial") 
    }
    cat(paste("Using gradient truncation of ",gt,"\n",sep="")) 
    return(gt)
}

