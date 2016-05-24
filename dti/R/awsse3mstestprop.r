# This file contains the implementation of dti.smooth() for 
# "dtiData" Adaptive smoothing in SE(3) considering b=0 as an individual shell

dwi.smooth.testprop <-  function(spatialdim,ngrad,bv,D0=1.3e-3,th0,df,kstar,lambda=15,kappa0=.9,ncoils=df/2,sigma=NULL,ws0=1,usemaxni=TRUE,seed=1,minlevel=1e-6,maxz=25,diffz=.5){
  #  spatialdim - spatial dimention of images
  #  ngrad      - number of gradients per shell if length(bv) != ngrad, otherwise total number of gradients
  #  
  #
  #
  data("optgradients", envir = environment())
  grad <- dti::optgrad[[ngrad-5]]
  ns0 <- 1
  lbv <- length(bv)
  if(lbv!=ngrad){
    bv <- unique(bv)
    lbv <- length(bv)
    grad0 <- c(0,0,0)
    bv0 <- 0
    for(i in 1:lbv){
      grad0 <- cbind(grad0,grad)
      bv0 <- c(bv0,rep(bv[i],ngrad))
    }
    ngrad <- lbv*ngrad
    grad <- grad0[,-1]
    bv <- bv0[-1]
  }
  #
  #  prepare for criterion related statistics
  #
  vext <- c(1,1)
  nvox <- prod(spatialdim)
  set.seed(seed)
  ze <- seq(0,maxz,diffz)
  nz <- length(ze)
  exceedence0  <- exceedenceb  <- matrix(0,nz,kstar) # this is used to store   
  if(minlevel < 5/nvox) {
    minlevel <- 5/nvox
    cat("minlevel reset to",minlevel,"due to insufficient size of test sample\n")
  }
  elevel <- trunc(log(1e-6,10))
  levels <- as.vector(outer(c(.5,.2,.1),10^(-0:elevel)))
  levels <- levels[levels>=minlevel]
  par(mfrow=c(1,2),mar=c(3,3,3,3),mgp=c(2,1,0))
  
  kldistnorm1 <- function(th1,y,df){
    #require(gsl)
    L <- df/2
    m1 <- sqrt(pi/2)*gamma(L+1/2)/gamma(1.5)/gamma(L)*hg1f1(-0.5,L, -th1^2/2)
    (m1-y)^2/2/(2*L+th1^2-m1^2)
  }
  #
  #   generate data
  #
  s0 <- array(sqrt(rchisq(nvox,df,th0^2)),spatialdim)
  mask <- array(TRUE,spatialdim)
  if(is.null(sigma)) sigma <- awssigmc(s0,12,mask,ncoils,vext,h0=1.25,verbose=FALSE)$sigma
  thb <- th0*exp(-bv*D0)
  sb <- aperm(array(sqrt(rchisq(nvox*ngrad,df,thb^2)),c(ngrad,spatialdim)),c(2:4,1))
  level <- 0
  varstats <- sofmchi(ncoils)
  ddim <- spatialdim
  msstructure <- getnext3g(grad,bv)
  nshell <- msstructure$nbv
  cat("generated nc-chi data: image size:",spatialdim,"number of gradients:",ngrad,
      "number of shells",nshell,"degrees of freedom:",df,"\n")
  cat("Noncentrality parameters for shells (bv=0 first)",th0,unique(thb),"\n")
  cat("specified effective number of coils",ncoils,"estimated sigma",sigma,"\n")    
  cat("parameters: lambda=",lambda,"kappa0=",kappa0,"ncoils=",ncoils,"kstar=",kstar,"\n")
  if(is.null(kappa0)){
    #  select kappa based on variance reduction on the sphere
    warning("You need to specify  kappa0  ")
    stop()
  }
  #
  #  rescale so that we have Chi-distributed values
  #
  sb <- sb/sigma
  s0 <- s0/sigma
  #  th0 <- array(s0,c(ddim,nshell+1))
  ni0 <- array(1,ddim)
  gradstats <- getkappasmsh3(grad, msstructure)
  #     save(gradstats,file="gradstats.rsc")
  hseq <- gethseqfullse3msh(kstar,gradstats,kappa0,vext=vext)
  nind <- as.integer(hseq$n*1.25)
  # make it nonrestrictive for the first step
  z <- list(th=array(1,dim(sb)), th0=array(1,dim(s0)), ni = array(1,dim(sb)), ni0 = array(1,dim(s0)))
  if(usemaxni){
    ni <- array(1,dim(sb))
    ni0 <- array(1,dim(s0))
  }
  minlevel <- gamma(ncoils+0.5)/gamma(ncoils)*sqrt(2)
  #  thats the mean of the central chi distribution with 2*ncoils df
  prt0 <- Sys.time()
  cat("adaptive smoothing in SE3, kstar=",kstar,"\n")
  kinit <- if(lambda<1e10) 0 else kstar
  mc.cores <- setCores(,reprt=FALSE)
  for(k in kinit:kstar){
    gc()
    hakt <- hseq$h[,k+1]
    hakt0 <- mean(hakt)
    t0 <- Sys.time()
    thnimsh <- interpolatesphere0(z$th,z$th0,ni,ni0,msstructure,mask)
    t1 <- Sys.time()
    param <- lkfullse3msh(hakt,kappa0/hakt,gradstats,vext,nind) 
    param0 <- lkfulls0(hakt0,vext,nind) 
    vs2 <- varstats$s2[findInterval(thnimsh$mstheta, varstats$mu, all.inside = TRUE)]/2
    vs02 <- varstats$s2[findInterval(thnimsh$msth0, varstats$mu, all.inside = TRUE)]/2
    t2 <- Sys.time()
    z <- .Fortran("adsmse3s",
                  as.double(sb),#y
                  as.double(s0),#y0
                  as.double(thnimsh$mstheta),#th
                  as.double(thnimsh$msni),#ni/si^2
                  as.double(thnimsh$msth0),#th0
                  as.double(thnimsh$msni0),#ni0/si^2
                  as.double(vs2),#si^2/2
                  as.double(vs02),#si^2/2 for s0
                  as.logical(mask),#mask
                  as.integer(nshell+1),#ns number of shells
                  as.integer(ddim[1]),#n1
                  as.integer(ddim[2]),#n2
                  as.integer(ddim[3]),#n3
                  as.integer(ngrad),#ngrad
                  as.double(lambda),#lambda
                  as.double(ws0),# wghts0 rel. weight for s0 image
                  as.integer(mc.cores),#ncores
                  as.integer(param$ind),#ind
                  as.double(param$w),#w
                  as.integer(param$n),#n
                  as.integer(param0$ind),#ind0
                  as.double(param0$w),#w0
                  as.integer(param0$n),#n0
                  th=double(prod(ddim)*ngrad),#thn
                  ni=double(prod(ddim)*ngrad),#nin
                  th0=double(prod(ddim)),#th0n
                  ni0=double(prod(ddim)),#ni0n
                  double(ngrad*mc.cores),#sw
                  double(ngrad*mc.cores),#swy
                  double((nshell+1)*mc.cores),#thi
                  double((nshell+1)*mc.cores),#nii
                  double((nshell+1)*mc.cores),#fsi2                
                  PACKAGE="dti")[c("ni","th","ni0","th0")]
    t3 <- Sys.time()
    if(usemaxni){
      ni <- z$ni <- if(usemaxni) pmax(ni,z$ni)
      ni0 <- z$ni0 <- if(usemaxni) pmax(ni0,z$ni0)
      nie <- max(ni)
      nie0 <- max(ni0)
    } 
    nie <- max(z$ni)
    nie0 <- max(z$ni0)
    dim(z$th) <- c(ddim,ngrad)
    dim(z$th0) <- c(ddim)
    dim(z$ni) <- c(ddim,ngrad)
    dim(z$ni0) <- c(ddim)
    gc()
    cat("Step",k,"completet at",format(t3),"\n")
    s0hat <- pmax(z$th0,minlevel)*sigma
    sbhat <- aperm(pmax(z$th,minlevel)*sigma,c(4,1:3))# bring spherical component to front
    nihat <- aperm(z$ni,c(4,1:3)) # bring spherical component to front
    kldist0 <- kldistnorm1(th0,s0hat,df)
    kldistb <- kldistnorm1(thb,sbhat,df)
    exceedence0[,k] <- .Fortran("exceed",
                                as.double(kldist0),
                                as.integer(length(kldist0)),
                                as.double(ze/nie0),
                                as.integer(nz),
                                exprob=double(nz),
                                PACKAGE="dti")$exprob
    exceedenceb[,k] <- .Fortran("exceed",
                                as.double(kldistb),
                                as.integer(length(kldistb)),
                                as.double(ze/nie),
                                as.integer(nz),
                                exprob=double(nz),
                                PACKAGE="dti")$exprob
    if(k>2){
      contour(ze,1:k,exceedence0[,1:k],levels=levels,ylab="step",xlab="z",
              main=paste("S0-image lambda=",lambda," Exceed. Prob."))
      contour(ze,1:k,exceedenceb[,1:k],levels=levels,ylab="step",xlab="z",
              main=paste("Sb-images lambda=",lambda," Exceed. Prob."))
    }
  }
  list(exceedence0=exceedence0,exceedenceb=exceedenceb,levels=levels,z=ze)
}


