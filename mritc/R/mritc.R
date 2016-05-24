mritc <- function(intarr, mask, method=c("EM", "ICM", "HMRFEM", "MCMC",
                                  "MCMCsub", "PVHMRFEM", "MCMCsubbias")){
    if (length(dim(intarr)) != 3)
        stop("The intensity values of an MR image has to be in a 3D array.")

    if(length(dim(mask)) != 3)
        stop("The 'mask' has to be of dimension 3.")
    if(! all(unique(as.vector(mask)) %in% 0:1))
        stop("The value of mask has to be either 1 or 0.")
    if(all(mask==0))
        stop("All voxels are outside the mask.")
    
    method <- match.arg(method)
    
    y <- intarr[mask == 1]
    if(method == "MCMCsub") sub <- TRUE
    else sub <- FALSE
    mrispatial <- makeMRIspatial(mask, nnei=6, sub)
    if(method=="MCMCsubbias"){
        mrispatial26 <- makeMRIspatial(mask, nnei=26, sub=FALSE, bias=TRUE)
    }
    init <- initOtsu(y, 2)
    prop <- init$prop
    mu <- init$mu
    sigma <- init$sigma
  
    result <- switch(method,
              EM = mritc.em(y, prop, mu, sigma, verbose=TRUE),
              ICM = mritc.icm(y, mrispatial$neighbors, mrispatial$blocks,
                          mu=mu, sigma=sigma, verbose=TRUE),
              HMRFEM = mritc.hmrfem(y, mrispatial$neighbors, mrispatial$blocks,
                             mu=mu, sigma=sigma, verbose=TRUE),
              MCMC = mritc.bayes.nobias(y, mrispatial$neighbors, mrispatial$blocks,
                         mrispatial$sub, mrispatial$subvox,
                         spatialMat=diag(1,3), beta=0.7, mu, sigma,
                         niter=100, verbose=TRUE),
              MCMCsub = mritc.bayes.nobias(y, mrispatial$neighbors, mrispatial$blocks,
                            mrispatial$sub, mrispatial$subvox,
                            spatialMat=matrix(c(2,0,-1,0,2,0,-1,0,2), nrow=3),
                            beta=0.3, mu, sigma, niter=100, verbose=TRUE),
              MCMCsubbias = mritc.bayes.bias(y, mrispatial$neighbors,
                                mrispatial$blocks, mrispatial$subvox,
                                mrispatial26$neighbors, mrispatial26$blocks,
                                mrispatial26$weineighbors, mrispatial26$weights,  
                                spatialMat=matrix(c(2,0,-1,0,2,0,-1,0,2), nrow=3),
                                beta=0.3, mu, sigma, niter=1000, verbose=TRUE),
              PVHMRFEM = mritc.pvhmrfem(y, mrispatial$neighbors, mrispatial$blocks,
                                  mu=mu, sigma=sigma, verbose=TRUE))
    
    
    class(result) <- "mritc"
    result$method <- method
    result$mask <- mask

    result
}

