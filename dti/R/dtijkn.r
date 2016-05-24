dtiTensJKN <- function( object, what=c("FA","GA","MD"), coverage=c(.5,.8,.9,.95),mc.cores = setCores( , reprt = FALSE)) {
  ##
  ##  Estimate bias and variance of tensor characteristics by Jackknife   
  ##  object need to be of S4-class dtiData 
  #
  ##  First tensor estimates
  #
  
  if(is.null(mc.cores)) mc.cores <- 1
  mc.cores <- min(mc.cores,detectCores())
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  grad <- object@gradient
  ddim <- object@ddim
  ntotal <- prod(ddim)
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  sdcoef <- object@sdcoef
  z <- sioutlier(object@si,s0ind,mc.cores=mc.cores)
  #
  #  this does not scale well with openMP
  #
  cat("sioutlier completed\n")
  si <- array(z$si,c(ngrad,ddim))
  index <- z$index
  ngrad0 <- ngrad - length(s0ind)
  s0 <- si[s0ind,,,]
  si <- si[-s0ind,,,]
  if(ns0>1) {
    dim(s0) <- c(ns0,prod(ddim))
    s0 <- rep(1/ns0,ns0)%*%s0
    dim(s0) <- ddim
  }
  mask <- s0 > object@level
  mask <- connect.mask(mask)
  dim(si) <- c(ngrad0,prod(ddim))
  ttt <- array(0,dim(si))
  ttt[,mask] <- -log1p(sweep(si[,mask],2,as.vector(s0[mask]),"/")-1)
  #  suggestion by B. Ripley
  #   idsi <- 1:length(dim(si))
  #   ttt <- -log(sweep(si,idsi[-1],s0,"/"))
  ttt[is.na(ttt)] <- 0
  ttt[(ttt == Inf)] <- 0
  ttt[(ttt == -Inf)] <- 0
  dim(ttt) <- c(ngrad0,prod(ddim))
  cat("Data transformation completed ",format(Sys.time()),"\n")
  
  btbsvd <- svd(object@btb[,-s0ind])
  solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
  D <- solvebtb%*% ttt
  cat("Diffusion tensors (linearized model) generated ",format(Sys.time()),"\n")
  D[c(1,4,6),!mask] <- 1e-6
  D[c(2,3,5),!mask] <- 0
  D <- dti3Dreg(D,mc.cores=mc.cores)
  dim(D) <- c(6,ntotal)
  th0 <- array(s0,ntotal)
  th0[!mask] <- 0
  nvox <- sum(mask)
  ##
  ##  nonlinear regression
  ## 
  si <- array(z$si,c(ngrad,ddim))
  rm(z)
  ngrad0 <- ngrad
  df <- sum(table(object@replind)-1)
  param <- matrix(0,7,nvox)
  ms0 <- mean(s0[mask])
  mbv <- mean(object@bvalue[-object@s0ind])
  ##  use ms0 and mbv to rescale parameters such that they are of comparable magnitude
  param[1,] <- s0[mask]/ms0
  param[-1,] <- D2Rall(D[,mask]*mbv)
  ## use reparametrization D = R^T R
  sdcoef[-2] <- sdcoef[-2]/ms0# effect of rescaling of signal
  cat("start nonlinear regression",format(Sys.time()),"\n")
  if(mc.cores==1){
    param <- matrix(.C("dtens",
                       as.integer(nvox),
                       param=as.double(param),
                       as.double(matrix(si,c(ngrad,ntotal))[,mask]/ms0),
                       as.integer(ngrad),
                       as.double(object@btb/mbv),
                       as.double(sdcoef),
                       as.double(rep(0,ngrad)),#si
                       as.double(rep(1,ngrad)),#var                         
                       as.integer(1000),#maxit
                       as.double(1e-7),#reltol
                       PACKAGE="dti")$param,7,nvox)
  } else {
    x <- matrix(0,ngrad+7,nvox)
    x[1:7,] <- param
    x[-(1:7),] <- matrix(si,c(ngrad,ntotal))[,mask]/ms0
    param <- plmatrix(x,ptensnl,ngrad=ngrad,btb=object@btb/mbv,
                      sdcoef=sdcoef,maxit=1000,reltol=1e-7)
  }
  th0[mask] <- param[1,]*ms0
  D[,mask] <- R2Dall(param[-1,])/mbv
  dim(th0) <- ddim
  dim(D) <- c(6,ddim)
  dim(mask) <-  ddim
  gc()
  ## Regularizations
  if(any(is.na(D))|any(abs(D)>1e10)){
    dim(D) <- c(6,ntotal)
    indD <- (1:ntotal)[apply(is.na(D),2,any)]
    cat("NA's in D in ", length(indD),"voxel:",indD,"\n")
    D[,indD] <- c(1,0,0,1,0,1)
    mask[indD] <- FALSE
    indD <- (1:ntotal)[apply(abs(D)>1e10,2,any)]
    cat("Inf's in D in", length(indD)," voxel:",indD,"\n")
    D[,indD] <- c(1,0,0,1,0,1)
    mask[indD] <- FALSE
  }
  cat("Obtained tensor estimates",format(Sys.time()),"\n")
  ##
  ##   now prepare for Jackknife estimation of Bias and Variance
  ##
  z <- dtiind3D(D,mask,mc.cores=mc.cores, verbose = FALSE)
  if("FA" %in% what){
    fa <- pmax(z$fa,1e-100)
    tfa <- log(-log(fa))
    tfa1 <- tfa2 <- 0
  }
  if("GA" %in% what){
    ga <- pmax(z$ga,1e-100)
    tga <- log(ga)
    tga1 <- tga2 <- 0
  }
  if("MD" %in% what){
    md <- pmax(z$md,1e-100)
    tmd <- log(md)
    tmd1 <- tmd2 <- 0
  }
  ##
  ##  Leave-1-out jackknife estimation
  ##
  ngrad0 <- ngrad-ns0
  siind <- (1:ngrad)[-s0ind]
  dim(D) <- c(6,ntotal)
  for(i in 1:(ngrad0)){
    Djkn <- array(0,c(6,ntotal))
    j <- siind[i] # only use sb images
    param[1,] <- th0[mask]/ms0
    param[-1,] <- D2Rall(D[,mask]*mbv)
    if(mc.cores==1){
      param <- matrix(.C("dtens",
                         as.integer(nvox),
                         param=as.double(param),
                         as.double(matrix(si[-j,,,],c(ngrad-1,ntotal))[,mask]/ms0),
                         as.integer(ngrad-1),
                         as.double(object@btb[,-j]/mbv),
                         as.double(sdcoef),
                         as.double(rep(0,ngrad-1)),#si
                         as.double(rep(1,ngrad-1)),#var                         
                         as.integer(1000),#maxit
                         as.double(1e-7),#reltol
                         PACKAGE="dti")$param,7,nvox)
    } else {
      x <- matrix(0,ngrad+6,nvox)
      x[1:7,] <- param
      x[-(1:7),] <- matrix(si[-j,,,],c(ngrad-1,ntotal))[,mask]/ms0
      param <- plmatrix(x,ptensnl,ngrad=ngrad-1,btb=object@btb[,-j]/mbv,
                        sdcoef=sdcoef,maxit=1000,reltol=1e-7)
    }
    Djkn[,mask] <- R2Dall(param[-1,])/mbv
    z <- dtiind3D(Djkn,mask,mc.cores=mc.cores, verbose = FALSE)
    if("FA" %in% what){
      zfa <- pmax(z$fa,1e-100)
      zfa[is.na(zfa)] <- 1e-100
      zfa <- log(-log(zfa))
      tfa1 <- tfa1 + zfa
      tfa2 <- tfa2 + zfa*zfa
    }
    if("GA" %in% what){
      zga <- pmax(z$ga,1e-100)
      zga[is.na(zga)] <- 1e-100
      zga <- log(zga)
      tga1 <- tga1 + zga
      tga2 <- tga2 + zga*zga
    }
    if("MD" %in% what){
      zmd <- pmax(z$md,1e-100)
      zmd[is.na(zmd)] <- 1e-100
      zmd <- log(zmd)
      tmd1 <- tmd1 + zmd
      tmd2 <- tmd2 + zmd*zmd
    }
    cat("Step ",i," out of ",ngrad0," completed ",format(Sys.time()),"\n")
  }
  cat("We have all estimates, preparing results ",format(Sys.time()),"\n")
  ergs <- list(coverage=coverage)
  cq <- qnorm(.5+coverage/2)
  lcq <- length(cq)
  if("FA" %in% what){
    ##
    ##  Bias and variance under log(-log(FA)) transform
    ##
    tfa1 <- tfa1/ngrad0
    ergs$bfa <- array((ngrad0-1)*(tfa1-tfa),ddim)
    ergs$vfa <- array(pmax(0,(ngrad0-1)*(tfa2/ngrad0-tfa1*tfa1)),ddim)
    sdfa <- sqrt(ergs$vfa)
    faquant <- array(0,c(lcq,2,ddim))
    for(i in 1:lcq){
      faquant[i,1,,,] <- exp(-exp(array(tfa,ddim)-ergs$bfa+cq[i]*sdfa))
      faquant[i,2,,,] <- exp(-exp(array(tfa,ddim)-ergs$bfa-cq[i]*sdfa))
    }
    ergs$faquant <- faquant
    ergs$fa <- exp(-exp(array(tfa,ddim)))
    ergs$bcfa <- exp(-exp(array(tfa,ddim)-ergs$bfa))
  }
  if("GA" %in% what){
    ##
    ##  Bias and variance under log(GA) transform
    ##
    tga1 <- tga1/ngrad0
    ergs$bga <- array((ngrad0-1)*(tga1-tga),ddim)
    ergs$vga <- array(pmax(0,(ngrad0-1)*(tga2/ngrad0-tga1*tga1)),ddim)
    sdga <- sqrt(ergs$vga)
    gaquant <- array(0,c(lcq,2,ddim))
    for(i in 1:lcq){
      gaquant[i,1,,,] <- exp(array(tga,ddim)-ergs$bga-cq[i]*sdga)
      gaquant[i,2,,,] <- exp(array(tga,ddim)-ergs$bga+cq[i]*sdga)
    }
    ergs$gaquant <- gaquant
    ergs$ga <- exp(array(tga,ddim))
    ergs$bcga <- exp(array(tga,ddim)-ergs$bga)
  }
  if("MD" %in% what){
    ##
    ##  Bias and variance under log(GA) transform
    ##
    tmd1 <- tmd1/ngrad0
    ergs$bmd <- array((ngrad0-1)*(tmd1-tmd),ddim)
    ergs$vmd <- array(pmax(0,(ngrad0-1)*(tmd2/ngrad0-tmd1*tmd1)),ddim)
    sdmd <- sqrt(ergs$vmd)
    mdquant <- array(0,c(lcq,2,ddim))
    for(i in 1:lcq){
      mdquant[i,1,,,] <- exp(array(tmd,ddim)-ergs$bmd-cq[i]*sdmd)
      mdquant[i,2,,,] <- exp(array(tmd,ddim)-ergs$bmd+cq[i]*sdmd)
    }
    ergs$mdquant <- mdquant
    ergs$md <- exp(array(tmd,ddim))
    ergs$bcmd <- exp(array(tmd,ddim)-ergs$bmd)
  }
  ##  xx - contains estimate of xx from tensor model
  ##  bcxx - contains bias-corrected of xx 
  ##  bxx - contains bias estimate under transformation
  ##  vxx - contains variance estimate under transformation
  ##  xxquant[i,1,,,] - contains lower confidence limit for coverage probability
  ##              coverage[i]
  ##  xxquant[i,2,,,] - contains upper confidence limit for coverage probability
  ##              coverage[i]
  invisible(ergs)
}
