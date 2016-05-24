################################################################
#                                                              #
# Section for DTI functions                                    #
#                                                              #
################################################################

dtiTensor <- function(object,  ...) cat("No DTI tensor calculation defined for this class:",class(object),"\n")

setGeneric("dtiTensor", function(object,  ...) standardGeneric("dtiTensor"))

setMethod( "dtiTensor", "dtiData",
           function( object, 
                     method = c( "nonlinear", "linear", "quasi-likelihood"),
                     sigma = NULL, L = 1, 
                     mc.cores = setCores( , reprt = FALSE)) {
             
             ## check method! available are: 
             ##   "linear" - use linearized model (log-transformed)
             ##   "nonlinear" - use nonlinear model with parametrization according to Koay et.al. (2006)
             ##   "quasi-likelihood" - use nonlinear model Gaussian Approximation to \chi 
             ##
             ## in case of method="quasi-likelihood" we need estimates of sigma
             ##    possible formats: 
             ##       - scalar value    (global sigma, not recomended)
             ##       - array dim=ddim  (same sigma for all b-vaues (including 0), recommended for
             ##                         single shell experiments)
             ##       - array dim=c(ddim,nbv) (separate sigma for unique bvalues as given by 
             ##         function unifybvalues )
             ##       - array dim=c(ddim,ngrad)  
             method <- match.arg(method)
             
             if(is.null(mc.cores)) mc.cores <- 1
             mc.cores <- min(mc.cores,detectCores())
             args <- sys.call(-1)
             args <- c(object@call,args)
             ngrad <- object@ngrad
             grad <- object@gradient
             btb <- object@btb
             ddim <- object@ddim
             ntotal <- prod(ddim)
             s0ind <- object@s0ind
             ns0 <- length(s0ind)
             sdcoef <- object@sdcoef
             ngrad0 <- ngrad-ns0
             z <- sioutlier1(object@si,s0ind,object@level,mc.cores=mc.cores)
             #
             #  this does not scale well with openMP
             #
             cat("sioutlier completed\n")
             mask <- z$mask
             nvox <- sum(mask)
             ttt <- array(0,c(ngrad0,nvox))
             ttt <- -log1p(sweep(z$si[-s0ind,],2,as.vector(z$s0),"/")-1)
             #  suggestion by B. Ripley
             #   idsi <- 1:length(dim(si))
             #   ttt <- -log(sweep(si,idsi[-1],s0,"/"))
             ttt[is.na(ttt)] <- 0
             ttt[(ttt == Inf)] <- 0
             ttt[(ttt == -Inf)] <- 0
             cat("Data transformation completed ",format(Sys.time()),"\n")
             
             D <- matrix(0,6,ntotal)
             btbsvd <- svd(btb[,-s0ind])
             solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
             D[,mask] <- solvebtb%*% ttt
             cat("Diffusion tensors (linearized model) generated ",format(Sys.time()),"\n")
             rm(ttt)
             D[c(1,4,6),!mask] <- 1e-6
             D[c(2,3,5),!mask] <- 0
             D <- dti3Dreg(D,mc.cores=mc.cores)
             dim(D) <- c(6,ntotal)
             th0 <- array(0,ntotal)
             th0[mask] <- z$s0
             index <- z$index
             if(method %in%  c("nonlinear","quasi-likelihood")){
               #  method == "nonlinear" uses linearized model for initial estimates
               ngrad0 <- ngrad
               df <- sum(table(object@replind)-1)
               param <- matrix(0,7,nvox)
               ms0 <- mean(z$s0)
               mbv <- mean(object@bvalue[-object@s0ind])
               ##  use ms0 and mbv to rescale parameters such that they are of comparable magnitude
               param[1,] <- z$s0/ms0
               param[-1,] <- D2Rall(D[,mask]*mbv)
               ## use reparametrization D = R^T R
               sdcoef[-2] <- sdcoef[-2]/ms0# effect of rescaling of signal
               cat("start nonlinear regression",format(Sys.time()),"\n")
               if(mc.cores==1){
                 param <- matrix(.C("dtens",
                                    as.integer(nvox),
                                    param=as.double(param),
                                    as.double(z$si/ms0),
                                    as.integer(ngrad),
                                    as.double(btb/mbv),
                                    as.double(sdcoef),
                                    as.double(rep(0,ngrad)),#si
                                    as.double(rep(1,ngrad)),#var                         
                                    as.integer(1000),#maxit
                                    as.double(1e-7),#reltol
                                    PACKAGE="dti")$param,7,nvox)
               } else {
                 x <- matrix(0,ngrad+7,nvox)
                 x[1:7,] <- param
                 x[-(1:7),] <- z$si/ms0
                 param <- plmatrix(x,ptensnl,ngrad=ngrad,btb=btb/mbv,
                                   sdcoef=sdcoef,maxit=1000,reltol=1e-7)
               }
               th0[mask] <- pmax(1,param[1,]*ms0)
               D[,mask] <- R2Dall(param[-1,])/mbv
               cat("successfully completed nonlinear regression ",format(Sys.time()),"\n")
             }
             res <- matrix(0,ngrad,ntotal)
             zz <- .Fortran("tensres",
                            as.double(th0[mask]),
                            as.double(D[,mask]),
                            as.double(z$si),
                            as.integer(nvox),
                            as.integer(ngrad),
                            as.double(btb),
                            res = double(ngrad*nvox),
                            rss = double(nvox),
                            PACKAGE="dti")[c("res","rss")]
             res <- matrix(0,ngrad,ntotal)
             res[,mask] <- zz$res
             rss <- numeric(ntotal)
             rss[mask] <- zz$rss/(ngrad0-6)
             rm(zz)
             dim(th0) <- ddim
             dim(D) <- c(6,ddim)
             dim(res) <- c(ngrad,ddim)
             dim(rss) <- ddim
             dim(mask) <-  ddim
             gc()
             #
             #   get spatial correlation
             #
             if(any(is.na(res))){
               dim(res) <- c(ngrad,ntotal)
               indr <- (1:ntotal)[apply(is.na(res),2,any)]
               cat("NA's in res in voxel",indr,"\n")
               res[,indr] <- 0
             }
             if(any(is.na(D))|any(abs(D)>1e10)){
               dim(D) <- c(6,ntotal)
               indD <- (1:ntotal)[apply(is.na(D),2,any)]
               cat("NA's in D in ", length(indD),"voxel:",indD,"\n")
               D[,indD] <- c(1,0,0,1,0,1)
               mask[indD] <- FALSE
               z$si <- z$si[,-indD]
               indD <- (1:ntotal)[apply(abs(D)>1e10,2,any)]
               cat("Inf's in D in", length(indD)," voxel:",indD,"\n")
               D[,indD] <- c(1,0,0,1,0,1)
               mask[indD] <- FALSE
               z$si <- z$si[,-indD]
               nvox <- sum(mask)
               cat("voxel in mask",nvox,"\n")
               ## need to readjust if we had NA's in nonlinear regression
             }
             scorr <- mcorr(res,mask,ddim,ngrad0,lags=c(5,5,3),mc.cores=mc.cores)
             if(method=="quasi-likelihood"){
               if(is.null(sigma)||any(sigma<=0)){
                 cat("Please specify a valid estimate of sigma for quasi-likelihood\n 
            returning result for method='nonlinear' \n")
                 method <- "nonlinear"
               }
             }
             if(method=="quasi-likelihood"){
               ubv <- unifybvals(object@bvalue)
               uniquebv <- unique(ubv)
               nbv <- length(uniquebv)
               if(length(sigma)==1){
                 sigma <- array(sigma,ddim)[mask]
                 sigma <- array(sigma,c(nvox,ngrad))
               }
               else if(length(dim(sigma))==length(ddim)&&all(dim(sigma)==ddim)){
                 sigma <- array(sigma,ddim)[mask]
                 sigma <- array(sigma,c(nvox,ngrad))
               }
               else if(all(dim(sigma)==c(ddim,ngrad))){
                 sigma <- array(sigma,c(prod(ddim),ngrad))[mask,]
               } 
               else if(all(dim(sigma)==c(ddim,nbv))){
                 sigmain <- array(sigma,c(prod(ddim),nbv))
                 sigma <- array(0,c(nvox,ngrad))
                 for(i in 1:nbv) sigma[,ubv==uniquebv[i]] <- sigmain[mask,i]
               } else {
                 cat("Please specify a compatible estimate of sigma for quasi-likelihood\n 
            returning result for method='nonlinear' \n")
                 method <- "nonlinear"
               }
             }
             if(method=="quasi-likelihood"){
               param <- matrix(0,7,nvox)
               ## get sigma into the correct shape:
               cat("starting quasi-likelihood ",format(Sys.time()),"\n")
               dim(D) <- c(6,ntotal)
               param[1,] <- pmin(log(th0[mask]),10)
    ##  th0 = exp(param[1,]) needs to be positive
               param[-1,] <- D2Rall(D[,mask]*mbv)
               CL <- sqrt(pi/2)*gamma(L+1/2)/gamma(L)/gamma(3/2)
##               btb[,-s0ind] <- 0
### don't like this but otherwise we may see th0 to diverge if the model is inadequate ...
               if(mc.cores==1){
                  for(i in 1:nvox){
         ## z$si was created in sioutlier1, mask is set by thresholding with object$level + call to connect.mask
         ## dim(z$si) is c(ngrad,nvox)
                    param[,i] <- optim(param[,i],tchi,si=z$si[,i],
                                       sigma=sigma[i,],btb=btb/mbv,L=L,CL=CL,method="BFGS",
                                       control=list(reltol=1e-8,maxit=100))$par
                   if(i%/%1000*1000==i) cat(i,"voxel processed. Time:",format(Sys.time()),"\n")
                 }
               } else {
                 setCores(mc.cores)
                 x <- matrix(0,2*ngrad+7,nvox)
                 x[1:7,] <- param
                 x[(1:ngrad)+7,] <- z$si
                 x[(1:ngrad)+7+ngrad,] <- t(sigma)
                 ## implicit replication for skalar or two dimensional sigma
                 param <- plmatrix(x,ptenschi,fn=tchi,btb=btb/mbv,L=L,CL=CL)
                 cat(nvox,"voxel processed. Time:",format(Sys.time()),"\n")
               }
               D[,mask] <- R2Dall(param[-1,])/mbv
               th0[mask] <- exp(param[1,])
               cat("successfully completed quasi-likelihood ",format(Sys.time()),"\n")
             }
             ev <- dti3Dev(D,mask,mc.cores=mc.cores)
             dim(ev) <- c(3,ddim)   
             dim(D) <- c(6,ddim)   
             scale <- quantile(ev[3,,,][mask],.95,na.rm=TRUE)
             cat("estimated scale information",format(Sys.time()),"\n")  
             invisible(new("dtiTensor",
                           call  = args,
                           D     = D,
                           th0   = th0,
                           sigma = rss,
                           scorr = scorr$scorr, 
                           bw = scorr$bw, 
                           mask = mask,
                           hmax = 1,
                           gradient = object@gradient,
                           bvalue = object@bvalue,
                           btb   = btb,
                           ngrad = ngrad, # = dim(btb)[2]
                           s0ind = object@s0ind,
                           replind = object@replind,
                           ddim  = object@ddim,
                           ddim0 = object@ddim0,
                           xind  = object@xind,
                           yind  = object@yind,
                           zind  = object@zind,
                           voxelext = object@voxelext,
                           level = object@level,
                           orientation = object@orientation,
                           rotation = object@rotation,
                           source = object@source,
                           outlier = index,
                           scale = scale,
                           method = method)
             )
           })

D2Rall <- function(D){
  #
  #  used for initial values
  # in case of negative eigenvalues this uses c(1,0,0,1,0,1) for initial rho
  #
  nvox <- dim(D)[2]
  matrix(.Fortran("D2Rall",
                  as.double(D),
                  rho=double(6*nvox),
                  as.integer(nvox),
                  PACKAGE="dti")$rho,6,nvox)
}
R2Dall <- function(R){
  #
  #  used for initial values
  # in case of negative eigenvalues this uses c(1,0,0,1,0,1) for initial rho
  #
  nvox <- dim(R)[2]
  matrix(.Fortran("R2Dall",
                  as.double(R),
                  D=double(6*nvox),
                  as.integer(nvox),
                  PACKAGE="dti")$D,6,nvox)
}
#############
#############

dtiIndices <- function(object, ...) cat("No DTI indices calculation defined for this class:",class(object),"\n")

setGeneric("dtiIndices", function(object, ...) standardGeneric("dtiIndices"))

setMethod("dtiIndices","dtiTensor",function(object, mc.cores=setCores(,reprt=FALSE)) {
  args <- sys.call(-1)
  args <- c(object@call,args)
  ddim <- object@ddim
  n <- prod(ddim)
  z <- dtiind3D(object@D,object@mask,mc.cores=mc.cores)
  invisible(new("dtiIndices",
                call = args,
                fa = array(z$fa,ddim),
                ga = array(z$ga,ddim),
                md = array(z$md,ddim),
                andir = array(z$andir,c(3,ddim)),
                bary = array(z$bary,c(3,ddim)),
                gradient = object@gradient,
                bvalue = object@bvalue,
                btb   = object@btb,
                ngrad = object@ngrad, # = dim(btb)[2]
                s0ind = object@s0ind,
                ddim  = ddim,
                ddim0 = object@ddim0,
                voxelext = object@voxelext,
                orientation = object@orientation,
                rotation = object@rotation,
                xind  = object@xind,
                yind  = object@yind,
                zind  = object@zind,
                method = object@method,
                level = object@level,
                source= object@source)
  )
})

##
## Parallel version for quasi-likelihood
##
ptenschi <- function(x,fn,btb,L,CL){
  nvox <- dim(x)[2]
  param <- matrix(0,7,nvox)
  ngrad <- (dim(x)[1]-7)/2
  indsi <- (1:ngrad)+7
  indsg <- indsi+ngrad
  for(i in 1:nvox){
    param[,i] <-
      optim(x[1:7,i],fn,si=x[indsi,i],sigma=x[indsg,i],
            btb=btb,L=L,CL=CL,method="BFGS",
            control=list(reltol=1e-8,maxit=100))$par
  }
  param
}


tchi <- function(param,si,sigma,btb,L,CL){
  ##
  ##  Risk function for Diffusion Tensor model with
  ##  Gauss-approximation for noncentral chi
  ##
  ##   si are the original observations
  ##   si/sigma is assumed to follow a noncentral chi_{2L} distribution
  ##   sigma should be of length ng here
  ng <- dim(btb)[2]
#  cat("param",signif(param,4),"value")
  D <- .Fortran("rho2D0",
                as.double(param[-1]),
                D=double(6),
                PACKAGE="dti")$D
  logth <- param[1]
  if(logth>12) logth <- 12+log(logth/12)
  # logth should be < 12 for the solution
  # gDg <- D%*%btb ## b_i*g_i^TD g_i (i=1,ngrad)
  gvalue <- exp(logth-D%*%btb)/sigma
  #    mgvh <- -gvalue*gvalue/2
  #    muL <- CL*.C("hyperg_1F1_e",
  #               as.double(rep(-.5,ng)),
  #               as.double(rep(L,ng)),
  #               as.double(mgvh),
  #               as.integer(ng),
  #               val=as.double(mgvh),
  #               err=as.double(mgvh),
  #               status=as.integer(0*mgvh),
  #               PACKAGE="gsl")$val
  muL <- CL*hg1f1(rep(-.5,ng),rep(L,ng),-gvalue*gvalue/2)
#  vL <- 2*L+gvalue^2-muL^2
  vL <- pmax(1e-8,2*L+gvalue^2-muL^2)
  ## avoid negative variances that may result due to approximations
  ## within the iteration process 

  ## factor sigma in muL and sigma^2 in vL cancels in the quotient
  sum((si/sigma-muL)^2/vL)
}

