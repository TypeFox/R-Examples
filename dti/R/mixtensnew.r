dwiMixtensnew <- ## function(object, ...) cat("No dwiMixtensor calculation defined for this class:",class(object),"\n")
##
## setGeneric("dwiMixtensnew", function(object,  ...) standardGeneric("dwiMixtensnew"))
##
## setMethod("dwiMixtensnew","dtiData",
function(object, maxcomp=3, 
          model=c("MTiso","MTisoFA","MTisoEV"),
          fa=NULL, lambda=NULL, reltol=1e-10, maxit=5000, ngc=1000, 
          nguess=100*maxcomp^2, msc=c("BIC","AIC","AICC","none"),
          mc.cores = setCores(,reprt=FALSE)){
  #
  #  uses  S(g) = w_0 exp(-b*l_1) +\sum_{i} w_i exp(-b*l_2-b*(l_1-l_2)(g^T d_i)^2)
  #  w_i corresponds to th0* volume fraction of compartment i
  #  FA or (FA and l_2)  may be fixed  
  #  Optimization method: L-BFGS-B for tensor mixture models with isotropic compartment

  ## check model
  model <- match.arg(model)
  ## check msc
  msc <- match.arg(msc)
  factr <- reltol/1e-14 ## this is 1e6 
  set.seed(1)
  bvalue <- object@bvalue
  maxbv <- max(bvalue)
  bvalue <- bvalue/maxbv
  maxc <- .866
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  ddim <- object@ddim
  nvox <- prod(ddim)
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  if(5*(2+3*maxcomp)>ngrad){
    #     maxcomp <- max(1,trunc((ngrad-5)/15))
    cat("Maximal number of components reduced to", maxcomp,"due to insufficient
           number of gradient directions\n")
  }
  #
  #   which model should be used
  #
  if(model=="MTisoEV"&&(is.null(fa)||is.null(lambda))){
    cat("No eigenvalues specified with model=='MTisoEV'\n")
    if(is.null(fa)){
      cat("setting model='MTiso'\n")
      model<-"MTiso"
    } else {
      cat("setting model='MTisoFA'\n")
      model<-"MTisoFA"
    }
  }
  imodel <- switch(model,MTisoEV=0,MTisoFA=1,MTiso=2)
  alpha <- switch(model,MTisoEV=(fa^2+sqrt(fa^2*(3-2*fa^2)))/(1-fa^2),
                  MTisoFA=(fa^2+sqrt(fa^2*(3-2*fa^2)))/(1-fa^2),
                  MTiso=0)# not used
  lambda <- switch(model,MTisoEV=lambda*maxbv/(1+alpha),
                   MTisoFA=0,# will be adjusted later
                   MTiso=0)#  will not be used
  #
  #  First tensor estimates to generate eigenvalues and -vectors
  #
  prta <- Sys.time()
  cat("Start tensor estimation at",format(prta),"\n")
  tensorobj <- dtiTensor(object, mc.cores = mc.cores)
  cat("Start evaluation of eigenstructure at",format(Sys.time()),"\n")
  z <- dtieigen(tensorobj@D, tensorobj@mask, mc.cores = mc.cores)
  rm(tensorobj)
  gc()
  fa <- array(z$fa,ddim)
  ev <- array(z$ev,c(3,ddim))*maxbv
  #
  #  rescale by bvalue to go to implemented scale
  #
  andir <- array(z$andir,c(3,2,ddim))
  rm(z)
  gc()
  #
  #  prepare parameters for searching initial estimates
  #
    lambdahat <- ev[3,,,] 
    # use third ev instead of (ev[2,,,]+ev[3,,,])/2 to avoid effects from mixtures
    alphahat <- if(imodel==2) (ev[1,,,]-lambdahat)/lambdahat else alpha
    lambdahat <- .5*median(lambdahat[!is.na(alphahat)&fa>.3])
    if(imodel==2) alphahat <- median(alphahat[!is.na(alphahat)&fa>.3])
    fahat <- alphahat/sqrt(3+2*alphahat+alphahat^2)
    cat("Using lambda_2=",lambdahat,"fa=",fahat," and alpha=",alphahat,"in initial estimates\n")
  
  cat("Start search outlier detection at",format(Sys.time()),"\n")
  #
  #  replace physically meaningless S_i by mena S_0 values
  #
  z <- sioutlier(object@si,s0ind,mc.cores=mc.cores)
  #
  #  this does not scale well with openMP
  #
  cat("End search outlier detection at",format(Sys.time()),"\n")
  si <- array(z$si,c(ngrad,ddim))
  index <- z$index
  rm(z)
  gc()
  cat("Start generating auxiliary objects",format(Sys.time()),"\n")
  #
  #  compute mean S_0, s_i/S_0 (siq), var(siq) and mask
  #
  nvox <- prod(ddim[1:3])
  cat("means0:")
  t1 <- Sys.time()
  z <- .Fortran("means0",# mixtensbv.f
                  as.double(si[s0ind,,,,drop=FALSE]),
                  as.integer(nvox),
                  as.integer(ns0),
                  as.integer(object@level),
                  s0=double(nvox),
                  mask=logical(nvox),
                  PACKAGE="dti")[c("s0","mask")]
  s0 <- array(z$s0,ddim[1:3])
  mask <- array(z$mask,ddim[1:3])
  means0 <- mean(s0[mask]) 
  #
  #  rescale s0 for numerical reasons
  #
  si <- si/means0
  npar <- switch(model,
                 MTisoEV=1+3*(0:maxcomp),
                 MTisoFA=2+3*(0:maxcomp),
                 MTiso=3+3*(0:maxcomp))
  #
  #   compute penalty for model selection, default BIC
  #
  penIC <- switch(msc,AIC=2*npar/ngrad,BIC=log(ngrad)*npar/ngrad,
                  AICC=(1+npar/ngrad)/(1-(npar+2)/ngrad),
                  None=log(ngrad)-log(ngrad-npar))
  cat("End generating auxiliary objects",format(Sys.time()),"\n")
  #
  #  avoid situations where si's are larger than s0
  #
  grad <- t(object@gradient)
  #
  #   determine initial estimates for orientations 
  #
  cat("Start search for initial directions at",format(Sys.time()),"\n")
  data("polyeders", envir = environment())
  polyeder <- icosa3 <- icosa3 ## WORKAROUND to make "polyeders" not global variable (for R CMD check)
  vert <- polyeder$vertices
  # remove redundant directions
  vind <- rep(TRUE,dim(vert)[2])
  vind[vert[1,]<0] <- FALSE
  vind[vert[1,]==0 & vert[2,] <0] <- FALSE
  vind[vert[1,]==0 & vert[2,] == 0 &vert[3,]<0] <- FALSE
  vert <- vert[,vind]
  #
  #  compute initial estimates (EV from grid and orientations from icosa3$vertices)
  #
  si <- aperm(si,c(4,1:3))
  dim(si) <- c(ngrad,nvox)
    siind <- wi <- matrix(0,maxcomp+1,nvox)
    siind[1,!mask] <- -1
    krit <- numeric(nvox)
    nvoxm <- sum(mask)
    if(mc.cores<=1){
      z <-   getsiindbv(si[,mask],grad,bvalue,t(vert),alphahat,lambdahat,
                       maxcomp,maxc=maxc,nguess=nguess)
      krit[mask] <- z$krit # sqrt(sum of squared residuals) for initial estimates
      siind[,mask] <- z$siind # components 1: model order 2:
      wi[,mask] <- z$wi # weigths
      # grid index for EV 2+(1:m) index of orientations
    } else {
      mc.cores.old <- setCores(,reprt=FALSE)
      setCores(mc.cores)
      x <- si[,mask]
      nvico <- dim(vert)[2]
      dgrad <- matrix(abs(grad%*%vert),ngrad,nvico)
      dgrad <- dgrad/max(dgrad)
      dgradi <- matrix(abs(t(vert)%*%vert),nvico,nvico)
      dgradi <- dgradi/max(dgradi)
      isample <- selisample(nvico,maxcomp,nguess,dgradi,maxc)
      nguess <- length(isample)/maxcomp
      cat("using ",nguess,"guesses for initial estimates\n")
      z <- plmatrix(x,pgetsiindbv,grad=grad,bv=bvalue,nvico=nvico,
                    dgrad=dgrad,dgradi=dgradi,isample=isample,alpha=alphahat,
                    lambda=lambdahat,maxcomp=maxcomp,maxc=maxc,nguess=nguess)
      setCores(mc.cores.old,reprt=FALSE)
      krit[mask] <- z[1,] # risk
      siind[,mask] <- z[2:(maxcomp+2),] # siind 
      wi[,mask] <- z[-(1:(maxcomp+2)),] # wi 
      
    }
  cat("Model orders for initial estimates")
  print(table(siind[1,]))
  cat("End search for initial values at",format(Sys.time()),"\n")
  #  logarithmic eigen values
  orient <- array(0,c(2,maxcomp,ddim))
  prt0 <- Sys.time()
    siind <- siind[-1,,drop=FALSE]
    if(mc.cores<=1){
      cat("Starting parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
      z <- switch(imodel+1,.C("mixtrl0b", 
                              as.integer(nvoxm),#nvoxm
                              as.integer(siind[,mask]),#siind 
                              as.double(wi[,mask]),# wi
                              as.integer(ngrad),#ngrad
                              as.integer(maxcomp),#maxcomp
                              as.integer(maxit),#maxit
                              as.double(t(grad)),#grad_in
                              as.double(bvalue),#bv_in
                              lambda  = as.double(lambda),#lambda_in
                              alpha   = as.double(alpha),#alpha_in
                              as.double(factr),#factr
                              as.double(penIC),#penIC
                              as.double(krit[mask]),#best rss
                              as.double(vert),#vert
                              as.double(si[,mask]),#si_in
                              sigma2  = double(nvoxm),#sigma2_ret error variance 
                              orient  = double(2*maxcomp*nvoxm),#orient_ret phi/theta for all mixture tensors
                              order   = integer(nvoxm),#order_ret selected order of mixture
                              mix     = double((maxcomp+1)*nvoxm),#mixture weights
                              PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")],
                  .C("mixtrl1b", 
                     as.integer(nvoxm),#n1
                     as.integer(siind[,mask]),#siind 
                     as.double(wi[,mask]),# wi
                     as.integer(ngrad),#ngrad
                     as.integer(maxcomp),#maxcomp
                     as.integer(maxit),#maxit
                     as.double(t(grad)),#grad_in
                     as.double(bvalue),#bv_in
                     as.double(lambdahat),#lambda_in
                     alpha = as.double(alpha),#alpha_in
                     as.double(factr),#factr
                     as.double(penIC),#penIC
                     as.double(krit[mask]),#best rss
                     as.double(vert),#vert
                     as.double(si[,mask]),#si_in
                     sigma2  = double(nvoxm),#sigma2_ret error variance 
                     orient  = double(2*maxcomp*nvoxm),#orient_ret phi/theta for all mixture tensors
                     order   = integer(nvoxm),#order_ret selected order of mixture
                     lambda  = double(nvoxm),#lambda_ret lambda_2 
                     mix     = double((maxcomp+1)*nvoxm),#mixture weights
                     PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")],
                  .C("mixtrl2b", 
                     as.integer(nvoxm),#n1
                     as.integer(siind[,mask]),#siind 
                     as.double(wi[,mask]),# wi
                     as.integer(ngrad),#ngrad
                     as.integer(maxcomp),#maxcomp
                     as.integer(maxit),#maxit
                     as.double(t(grad)),#grad_in
                     as.double(bvalue),#bv_in
                     as.double(lambdahat),#lambda_in
                     as.double(alphahat),#alpha_in
                     as.double(factr),#factr
                     as.double(penIC),#penIC
                     as.double(krit[mask]),#best rss
                     as.double(vert),#vert
                     as.double(si[,mask]),#si_in
                     sigma2  = double(nvoxm),#sigma2_ret error variance 
                     orient  = double(2*maxcomp*nvoxm),#orient_ret phi/theta for all mixture tensors
                     order   = integer(nvoxm),#order_ret selected order of mixture
                     alpha   = double(nvoxm),#alpha_ret alpha=(lambda_1-lambda_2)/lambda_2 
                     lambda  = double(nvoxm),#lambda_ret lambda_2 
                     mix     = double((maxcomp+1)*nvoxm),#mixture weights
                     PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")])
      cat("End parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
      sigma2 <-  array(0,ddim)
      sigma2[mask] <- z$sigma2
      orient <- matrix(0,2*maxcomp,nvox)
      orient[,mask] <- z$orient
      dim(orient) <- c(2, maxcomp, ddim)
      order <- array(0, ddim)
      order[mask] <- z$order
      lev <- matrix(0,2,nvox)
      lev[2,mask] <- z$lambda
      lev[1,mask] <- (z$alpha+1)*z$lambda
      dim(lev) <- c(2,ddim)
      mix <- matrix(0,maxcomp+1,nvox)
      mix[,mask] <- z$mix
      dim(mix) <- c(maxcomp+1, ddim)
    } else {
      cat("Starting parameter estimation and model selection (C-code) on",mc.cores," cores",format(Sys.time()),"\n")
      x <- matrix(0,ngrad+2+2*maxcomp,sum(mask))
      dim(si) <- c(ngrad,nvox)
      x[1:ngrad,] <- si[,mask]
      x[ngrad+1,] <- krit[mask]
      dim(siind) <- c(maxcomp,nvox)
      x[ngrad+2:(1+maxcomp),] <- siind[,mask] 
      x[-(1:(1+maxcomp)),] <- wi[,mask] 
      res <- matrix(0,5+3*maxcomp,nvox)
      res[,mask] <- switch(imodel+1,
                           plmatrix(x,pmixtn0b,ngrad=ngrad,maxcomp=maxcomp,maxit=maxit,
                                    grad=grad,bv=bvalue,lambda=lambda,alpha=alpha,factr=factr,
                                    penIC=penIC,vert=vert,mc.cores=mc.cores),# model=0
                           plmatrix(x,pmixtn1b,ngrad=ngrad,maxcomp=maxcomp,maxit=maxit,
                                    grad=grad,bv=bvalue,lambda=lambdahat,alpha=alpha,factr=factr,
                                    penIC=penIC,vert=vert,mc.cores=mc.cores),# model=1
                           plmatrix(x,pmixtn2b,ngrad=ngrad,maxcomp=maxcomp,maxit=maxit,
                                    grad=grad,bv=bvalue,lambda=lambdahat,alpha=alphahat,factr=factr,
                                    penIC=penIC,vert=vert,mc.cores=mc.cores))# model=2
      cat("End parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
      rm(x)
      gc()
      sigma2 <-  array(res[2,],ddim)
      orient <- array(res[maxcomp+4+1:(2*maxcomp),], c(2, maxcomp, ddim))
      order <- array(as.integer(res[1,]), ddim)
      lev <- array(res[3:4,], c(2,ddim))
      mix <- array(res[4+(1:maxcomp),], c(maxcomp, ddim))
    }
    th0 <- apply(mix,2:4,sum)
    mix <- sweep(mix[-1,,,],2:4,th0,"/")
    th0 <- th0*means0
    save(krit,wi,mask,z,mix,th0,file="tmp.rsc")
    method <- switch(imodel+1,"MTisoEV","MTisoFA","MTiso")
    model <- switch(imodel+1,"iso-prolate-fixedev","iso-prolate-fixedfa","iso-prolate")
##  invisible(new("dwiMixtensor",
  list(              model = model,
                call   = args,
                ev     = lev/maxbv,
                mix    = mix,
                orient = orient,
                order  = order,
                p      = 0,
                s0=s0,
                z=z,
                th0    = th0,
                sigma  = sigma2*means0^2,
                scorr  = array(1,c(1,1,1)), 
                bw     = c(0,0,0), 
                mask   = mask,
                hmax   = 1,
                gradient = object@gradient,
                bvalue = object@bvalue,
                btb    = object@btb,
                ngrad  = object@ngrad, # = dim(btb)[2]
                s0ind  = object@s0ind,
                replind = object@replind,
                ddim   = object@ddim,
                ddim0  = object@ddim0,
                xind   = object@xind,
                yind   = object@yind,
                zind   = object@zind,
                voxelext = object@voxelext,
                level = object@level,
                orientation = object@orientation,
                rotation = object@rotation,
                source = object@source,
                outlier = index,
                scale = 1,
                method = method)
##  )
}
##)

imtfbv <- function(par,si,grad,bv,rho){
##  RSS for
##  Tensor mixture model with isotrop compartment  
##  (and multiple b-values)
   npar <- length(par)
   nc <-npar/3-1
   n <- length(si)
   .Fortran("imtfunbv",
            as.double(par),
            as.double(si),
            as.double(grad),
            as.double(bv),
            as.integer(npar),
            as.integer(nc),
            as.integer(n),
            as.double(rho),
            double(3*nc),#evc
            double(nc),#w
            rss=double(1),#residual sum of squares
            PACKAGE="dti")$rss
   }
imtgrbv <- function(par,si,grad,bv,rho){
##  Gradient for
##  Tensor mixture model with isotrop compartment  
##  (and multiple b-values)
   npar <- length(par)
   nc <-npar/3-1
   n <- length(si)
   .Fortran("imtgrdbv",
            as.double(par),
            as.double(si),
            as.double(grad),
            as.double(bv),
            as.integer(npar),
            as.integer(nc),
            as.integer(n),
            as.double(rho),
            double(3*nc),#evc
            double(nc),#w
            double(n),#f
            double(n*nc),#evg
            double(n),#el1
            double(n*nc),#el2k
            grad=double(npar),#gradient
            PACKAGE="dti")$grad
   }

imtfb1 <- function(par,alpha,si,grad,bv,rho){
##  RSS for
##  Tensor mixture model with isotrop compartment  
##  (and multiple b-values)
##  fixed FA (alpha=fa^2+sqrt(fa^2*(3-2*fa^2))/(1-fa^2)
##  l12 = l2*alpha   l1 = l2*(alpha+1)
   npar <- length(par)
   nc <- (npar-2)/3
   n <- length(si)
   .Fortran("imtfunb1",
            as.double(par),
            as.double(alpha),
            as.double(si),
            as.double(grad),
            as.double(bv),
            as.integer(npar),
            as.integer(nc),
            as.integer(n),
            as.double(rho),
            double(3*nc),#evc
            double(nc),#w
            rss=double(1),#residual sum of squares
            PACKAGE="dti")$rss
   }
imtgrb1 <- function(par,alpha,si,grad,bv,rho){
##  Gradient for
##  Tensor mixture model with isotrop compartment  
##  (and multiple b-values)
##  fixed FA (alpha=fa^2+sqrt(fa^2*(3-2*fa^2))/(1-fa^2)
##  l12 = l2*alpha   l1 = l2*(alpha+1)
   npar <- length(par)
   nc <- (npar-2)/3
   n <- length(si)
   .Fortran("imtgrdb1",
            as.double(par),
            as.double(alpha),
            as.double(si),
            as.double(grad),
            as.double(bv),
            as.integer(npar),
            as.integer(nc),
            as.integer(n),
            as.double(rho),
            double(3*nc),#evc
            double(nc),#w
            double(n),#f
            double(n*nc),#evg
            double(n),#el1
            double(n*nc),#el2k
            grad=double(npar),#gradient
            PACKAGE="dti")$grad
   }
imtfb0 <- function(par,alpha,l2,si,grad,bv,rho){
##  RSS for
##  Tensor mixture model with isotrop compartment  
##  (and multiple b-values)
##  fixed FA and l2 (alpha=fa^2+sqrt(fa^2*(3-2*fa^2))/(1-fa^2)
##  l12 = l2*alpha   l1 = l2*(alpha+1)
   npar <- length(par)
   nc <- (npar-1)/3
   n <- length(si)
   .Fortran("imtfunb0",
            as.double(par),
            as.double(alpha),
            as.double(l2),
            as.double(si),
            as.double(grad),
            as.double(bv),
            as.integer(npar),
            as.integer(nc),
            as.integer(n),
            as.double(rho),
            double(3*nc),#evc
            double(nc),#w
            rss=double(1),#residual sum of squares
            PACKAGE="dti")$rss
   }
imtgrb0 <- function(par,alpha,l2,si,grad,bv,rho){
##  Gradient for
##  Tensor mixture model with isotrop compartment  
##  (and multiple b-values)
##  fixed FA and l2 (alpha=fa^2+sqrt(fa^2*(3-2*fa^2))/(1-fa^2)
##  l12 = l2*alpha   l1 = l2*(alpha+1)
   npar <- length(par)
   nc <- (npar-1)/3
   n <- length(si)
   .Fortran("imtgrdb0",
            as.double(par),
            as.double(alpha),
            as.double(l2),
            as.double(si),
            as.double(grad),
            as.double(bv),
            as.integer(npar),
            as.integer(nc),
            as.integer(n),
            as.double(rho),
            double(3*nc),#evc
            double(nc),#w
            double(n),#f
            double(n*nc),#evg
            double(n),#el1
            double(n*nc),#el2k
            grad=double(npar),#gradient
            PACKAGE="dti")$grad
   }
mtfbv <- function(par,si,grad,bv,rho){
##  RSS for
##  Tensor mixture model without isotrop compartment  
##  (and multiple b-values)
   npar <- length(par)
   nc <-(npar-2)/3
   n <- length(si)
   .Fortran("mtfunbv",
            as.double(par),
            as.double(si),
            as.double(grad),
            as.double(bv),
            as.integer(npar),
            as.integer(nc),
            as.integer(n),
            as.double(rho),
            double(3*nc),#evc
            double(nc),#w
            rss=double(1),#residual sum of squares
            PACKAGE="dti")$rss
   }
mtgrbv <- function(par,si,grad,bv,rho){
##  Gradient for
##  Tensor mixture model without isotrop compartment  
##  (and multiple b-values)
   npar <- length(par)
   nc <- (npar-2)/3
   n <- length(si)
   .Fortran("mtgrdbv",
            as.double(par),
            as.double(si),
            as.double(grad),
            as.double(bv),
            as.integer(npar),
            as.integer(nc),
            as.integer(n),
            as.double(rho),
            double(3*nc),#evc
            double(nc),#w
            double(n),#f
            double(n*nc),#evg
            double(n*nc),#el2k
            grad=double(npar),#gradient
            PACKAGE="dti")$grad
   }
mtfb1 <- function(par,alpha,si,grad,bv,rho){
##  RSS for
##  Tensor mixture model without isotrop compartment  
##  (and multiple b-values)
##  fixed FA (alpha=fa^2+sqrt(fa^2*(3-2*fa^2))/(1-fa^2)
##  l12 = l2*alpha   l1 = l2*(alpha+1)
   npar <- length(par)
   nc <- (npar-1)/3
   n <- length(si)
   .Fortran("mtfunb1",
            as.double(par),
            as.double(alpha),
            as.double(si),
            as.double(grad),
            as.double(bv),
            as.integer(npar),
            as.integer(nc),
            as.integer(n),
            as.double(rho),
            double(3*nc),#evc
            double(nc),#w
            rss=double(1),#residual sum of squares
            PACKAGE="dti")$rss
   }
mtgrb1 <- function(par,alpha,si,grad,bv,rho){
##  Gradient for
##  Tensor mixture model without isotrop compartment  
##  (and multiple b-values)
##  fixed FA (alpha=fa^2+sqrt(fa^2*(3-2*fa^2))/(1-fa^2)
##  l12 = l2*alpha   l1 = l2*(alpha+1)
   npar <- length(par)
   nc <- (npar-1)/3
   n <- length(si)
   .Fortran("mtgrdb1",
            as.double(par),
            as.double(alpha),
            as.double(si),
            as.double(grad),
            as.double(bv),
            as.integer(npar),
            as.integer(nc),
            as.integer(n),
            as.double(rho),
            double(3*nc),#evc
            double(nc),#w
            double(n),#f
            double(n*nc),#evg
            double(n*nc),#el2k
            grad=double(npar),#gradient
            PACKAGE="dti")$grad
   }
mtfb0 <- function(par,alpha,l2,si,grad,bv,rho){
##  RSS for
##  Tensor mixture model without isotrop compartment  
##  (and multiple b-values)
##  fixed FA (alpha=fa^2+sqrt(fa^2*(3-2*fa^2))/(1-fa^2)
##  l12 = l2*alpha   l1 = l2*(alpha+1)
   npar <- length(par)
   nc <- npar/3
   n <- length(si)
   .Fortran("mtfunb0",
            as.double(par),
            as.double(alpha),
            as.double(l2),
            as.double(si),
            as.double(grad),
            as.double(bv),
            as.integer(npar),
            as.integer(nc),
            as.integer(n),
            as.double(rho),
            double(3*nc),#evc
            double(nc),#w
            rss=double(1),#residual sum of squares
            PACKAGE="dti")$rss
   }
mtgrb0 <- function(par,alpha,l2,si,grad,bv,rho){
##  Gradient for
##  Tensor mixture model without isotrop compartment  
##  (and multiple b-values)
##  fixed FA and l2 (alpha=fa^2+sqrt(fa^2*(3-2*fa^2))/(1-fa^2)
##  l12 = l2*alpha   l1 = l2*(alpha+1)
   npar <- length(par)
   nc <- npar/3
   n <- length(si)
   .Fortran("mtgrdb0",
            as.double(par),
            as.double(alpha),
            as.double(l2),
            as.double(si),
            as.double(grad),
            as.double(bv),
            as.integer(npar),
            as.integer(nc),
            as.integer(n),
            as.double(rho),
            double(3*nc),#evc
            double(nc),#w
            double(n),#f
            double(n*nc),#evg
            double(n*nc),#el2k
            grad=double(npar),#gradient
            PACKAGE="dti")$grad
   }
getsiindbv <- function(si,grad,bv,vico,alpha,lambda,maxcomp=3,maxc=.866,nguess=100){
  # assumes dim(grad) == c(ngrad,3)
  # assumes dim(si) == c(ngrad,n1,n2,n3)
  # SO removed
  ngrad <- dim(grad)[1]
  nvico <- dim(vico)[1]
  nsi <- dim(si)[1]
  dgrad <- matrix(abs(grad%*%t(vico)),ngrad,nvico)
  dgrad <- dgrad/max(dgrad)
  dgradi <- matrix(abs(vico%*%t(vico)),nvico,nvico)
  dgradi <- dgradi/max(dgradi)
  isample <- selisample(nvico,maxcomp,nguess,dgradi,maxc)
  nguess <- length(isample)/maxcomp
  #
  #  eliminate configurations with close directions 
  #
  # this provides configurations of initial estimates with minimum angle between 
  # directions > acos(maxc)
  nvoxel <- prod(dim(si)[-1])
  cat("using ",nguess,"guesses for initial estimates\n")
  z <- .Fortran("getsiibv",
                as.double(si),
                as.integer(nsi),
                as.integer(nvoxel),
                as.integer(maxcomp),
                as.double(dgrad),
                as.double(bv),
                as.integer(nvico),
                as.double(alpha),
                as.double(lambda),
                double(nsi*nvico),
                as.integer(isample),
                as.integer(nguess),
                double(nsi),
                double(nsi),#z0
                double(nsi*(maxcomp+1)),
                siind=integer((maxcomp+1)*nvoxel),
                wi=double((maxcomp+1)*nvoxel),
                krit=double(nvoxel),
                as.integer(maxcomp+1),
                PACKAGE="dti")[c("siind","krit","wi")]
  dim(z$siind) <- dim(z$wi) <- c(maxcomp+1,nvoxel)
  siind <- z$siind
  wi <- z$wi
  krit <- z$krit
  # now voxel where first tensor direction seems important
  list(siind=array(siind,c(maxcomp+1,dim(si)[-1])),
       wi=array(wi,c(maxcomp+1,dim(si)[-1])),
       krit=array(krit,dim(si)[-1]))
}

rskmixb2 <- function(par,si,grad,bv){
   npar <- length(par)
   ng <- dim(grad)[2]
   .Fortran("rskmixb2",
            as.double(par),
            as.integer(npar),
            as.double(si),
            as.double(grad),
            as.double(bv),
            as.integer(ng),
            rsk=double(1),
            PACKAGE="dti")$rsk
}
drskmb2 <- function(par,si,grad,bv){
   npar <- length(par)
   ng <- dim(grad)[2]
   .Fortran("drskmb2",
            as.double(par),
            as.integer(npar),
            as.double(si),
            as.double(grad),
            as.double(bv),
            as.integer(ng),
            drsk=double(npar),
            PACKAGE="dti")$drsk
}
fmixureb <- function(par,grad,bv){
  npar <- length(par)
  ng <- dim(grad)[2]
  fval <- numeric(ng)
  for(i in 1:ng) {
     fval[i] <- .Fortran("fmixturb",
           as.double(par[1:(npar-3)]),
           as.integer(npar-3),
           as.double(par[npar-2]),
           as.double(par[npar-1]),
           as.double(par[npar]),
           as.double(grad[,i]),
           as.double(bv[i]),
           fval=double(1),##
           PACKAGE="dti")$fval
  }
  fval
}
