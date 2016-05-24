likesigmaf <- function(sigma,wj,Sj,L){
  ni <- sum(wj)
  ksi <- sum(wj*Sj^2)/ni
  if(ksi-2*L*sigma^2 > 0){
    z <- .Fortran("lncchi",as.double(sigma),
                  as.double(ni),
                  as.double(ksi),
                  as.double(wj),
                  as.double(Sj),
                  as.double(L),
                  as.integer(length(Sj)),
                  double(floor(L+10)),
                  ergs=double(1),
                  PACKAGE="dti")$ergs-L
  } else {# central case
    #z <- -(2*L-1)*sum(wj*log(Sj))/ni+ksi/2/sigma^2+2*L*log(sigma)+(L-1)*log(2)+lgamma(L)
    z <- -(L-1)*sum(wj*log(Sj))/ni+ksi/2/sigma^2+2*L*log(sigma)+(L-1)*log(2)+lgamma(L)
  }
  z
}


mlikesigmaf <- function(sigma,wj,Sj,L){
  .Fortran("localmin",as.double(sigma/10),
           as.double(sigma*10),
           as.double(wj),
           as.double(Sj),
           as.double(L),
           as.integer(length(Sj)),
           as.double(1e-8),
           as.integer(100),
           double(floor(L+10)),
           xmin=double(1),
           fmin=double(1),
           PACKAGE="dti")[c("xmin","fmin")]
}
#
#
#      estimate variance parameter in a multicoil system
#
#
awslsigmc <- function(y,                 # data
                      steps,             # number of iteration steps for PS
                      mask = NULL,       # data mask, where to do estimation
                      ncoils = 1,        # number of coils for parallel MR image acquisition
                      vext = c( 1, 1),   # voxel extensions
                      lambda = 5,       # adaptation parameter for PS
                      minni = 2,         # minimum sum of weights for estimating local sigma
                      hsig = 5,          # bandwidth for median smoothing local sigma estimates
                      sigma = NULL,
                      family = c("NCchi"),
                      verbose = FALSE,
                      u=NULL#,
                      #bc=FALSE # bias correction ...
) {
  ## some functions for pilot estimates
  
  #   IQQ <- function (x, q = .25, na.rm = FALSE, type = 7) 
  #     diff(quantile(as.numeric(x), c(q, 1-q), na.rm = na.rm, names = FALSE, type = type))
  
  family <- match.arg(family)
  
  IQQ <- function (x, q = .25, na.rm = FALSE, type = 7){ 
    cqz <- qnorm(.05)/qnorm(q)
    x <- as.numeric(x)
    z <- diff(quantile(x, c(q, 1-q), na.rm = na.rm, names = FALSE, type = type))
    z0 <- 0
    while(abs(z-z0)>1e-5*z0){
      # outlier removal
      z0 <- z
      #       cat(sum(x>(z0*cqz))," ")
      x <- x[x<(z0*cqz)]
      z <- diff(quantile(x, c(q, 1-q), na.rm = na.rm, names = FALSE, type = type))  
      #       cat(z0,z,"\n")
    }
    z
  }
  
  IQQdiff <- function(y, mask, q = .25, verbose = FALSE) {
    cq <- qnorm(1-q)*sqrt(2)*2
    sx <- IQQ( diff(y[mask]), q)/cq
    sy <- IQQ( diff(aperm(y,c(2,1,3))[aperm(mask,c(2,1,3))]), q)/cq
    sz <- IQQ( diff(aperm(y,c(3,1,2))[aperm(mask,c(3,1,2))]), q)/cq
    if(verbose) cat( "Pilot estimates of sigma", sx, sy, sz, "\n")
    min( sx, sy, sz)
  }
  
  estsigma <- function(y, mask, q, L, sigma, verbose=FALSE){
    meany <- mean(y[mask]^2)
    eta <- sqrt(max( 0, meany/sigma^2-2*L))
    m <- sqrt(pi/2)*gamma(L+1/2)/gamma(L)/gamma(3/2)*hg1f1(-1/2,L,-eta^2/2)
    v <- max( .01, 2*L+eta^2-m^2)
    if(verbose) cat( eta, m, v, "\n")
    IQQdiff( y, mask, q)/sqrt(v)
  }
  varstats <- sofmchi(ncoils)
  if("NCchi" == family){
    th <- seq(0,30,.01)
    z <- fncchiv(th,varstats)
    minz <- min(z)
    th <- th[z>min(z)]
    minth <- min(th)
    z <- z[z>min(z)]
    nz <- nls(z~(a*th+b*th*th+c*th*th*th+d)/(a*th+b*th*th+c*th*th*th+1+d),data=list(th=th,z=z),start=list(a=1,b=1,c=1,d=1))
    vpar <- c(minth,minz,coef(nz))
    ##  this provides an excellent approximation for the variance reduction in case of low ncp 
  }
  
  if(length(vext)==3) vext <- vext[2:3]/vext[1]
  ## dimension and size of cubus
  ddim <- dim(y)
  n <- prod(ddim)
  
  ## test dimension
  if (length(ddim) != 3) stop("first argument should be a 3-dimentional array")
  
  ## check mask
  if (is.null(mask)) mask <- array(TRUE, ddim)
  if(length(mask) != n) stop("dimensions of data array and mask should coincide")
  
  ## initial value for sigma_0 
  if(is.null(sigma)){
    # sigma <- sqrt( mean( y[mask]^2) / 2 / ncoils)
    sigma <- IQQdiff( y, mask, .25, verbose=verbose)
    if("NCchi" == family){
      sigma <- estsigma( y, mask, .25, ncoils, sigma)
      sigma <- estsigma( y, mask, .25, ncoils, sigma)
      sigma <- estsigma( y, mask, .25, ncoils, sigma)
    }
  }
  ##
  ##   Prepare for diagnostics plots
  ##
  if(verbose){
    mslice <-  (ddim[3]+1)/2
    ymslice <- y[,,mslice]
    ymslice[!mask[,,mslice]] <- 0
    if(!is.null(u)&&"NCchi" == family){
      par(mfrow=c(2,4),mar=c(3,3,3,1),mgp=c(2,1,0))
    } else {
      par(mfrow=c(2,3),mar=c(3,3,3,1),mgp=c(2,1,0))
    }
  } else {
    cat("step")
  }
  ## define initial arrays for parameter estimates and sum of weights (see PS)
  th <- y
  ksi <- array( y^2, ddim)
  ni <- array( 1, ddim)
  sigma <- array(sigma, ddim)
  # initialize array for local sigma by global estimate
  mc.cores <- setCores(,reprt=FALSE)
  ## preparations for median smoothing
  nwmd <- (2*as.integer(hsig)+1)^3
  parammd <- .Fortran("paramw3",
                      as.double(hsig),
                      as.double(c(1,1)),
                      ind=integer(3*nwmd),
                      w=double(nwmd),
                      n=as.integer(nwmd),
                      PACKAGE = "dti")[c("ind","w","n")]
  nwmd <- parammd$n
  parammd$ind <- parammd$ind[1:(3*nwmd)]
  dim(parammd$ind) <- c(3,nwmd)
  ## iterate PS starting with bandwidth h0
  for (i in 1:steps) {
    
    h <- 1.25^((i-1)/3)
    nw <- prod(2*as.integer(h/c(1,vext))+1)
    #    cat("nw=",nw,"h/vext=",h/c(1,1/vext),"ih=",as.integer(h/c(1,1/vext)),"\n")
    param <- .Fortran("paramw3",
                      as.double(h),
                      as.double(vext),
                      ind=integer(3*nw),
                      w=double(nw),
                      n=as.integer(nw),
                      PACKAGE = "dti")[c("ind","w","n")]
    nw <- param$n
    param$ind <- param$ind[1:(3*nw)]
    dim(param$ind) <- c(3,nw)
    param$w   <- param$w[1:nw]    
      ## perform one step PS with bandwidth h
      z <- .Fortran("awslchi2",
                    as.double(y),        # data
                    as.double(ksi),# \sum w_ij S_j^2
                    ni = as.double(ni),
                    as.double(sigma),
                    as.double(vpar),# parameters for var. reduction
                    as.double(ncoils),
                    as.logical(mask),
                    as.integer(ddim[1]),
                    as.integer(ddim[2]),
                    as.integer(ddim[3]),
                    as.integer(param$ind),
                    as.double(param$w),
                    as.integer(nw),
                    as.double(minni),
                    double(nw*mc.cores), # wad(nw,nthreds)
                    double(nw*mc.cores), # sad(nw,nthreds)
                    as.double(lambda),
                    as.integer(mc.cores),
                    as.integer(floor(ncoils)),
                    double(floor(ncoils)*mc.cores), # work(L,nthreds)
                    th = double(n),
                    sigman = double(n),
                    ksi = double(n),
                    PACKAGE = "dti")[c("ni","ksi","th","sigman")]
      thchi <- z$th
      ksi <- z$ksi
      thchi[!mask] <- 0
      ## extract sum of weigths (see PS) and consider only voxels with ni larger then mean
    th <- array(z$th,ddim)
    ni <- array(z$ni,ddim)
    mask[z$sigman==0] <- FALSE
    nmask <- sum(mask)
    if(verbose) cat("local estimation in step ",i," completed",format(Sys.time()),"\n") 
    ##
    ##  nonadaptive smoothing of estimated standard deviations
    ##
    if(any(ni[mask]>minni)){
       z$sigman[mask] <- z$sigman[mask]*(1+runif(nmask,-.0001,.0001))
    ##
    ##  avoid ties in local neighborhood cubes of z$sigman
    ##
       sigma <- .Fortran("mediansm",
                      as.double(z$sigman),
                      as.logical(mask),
                      as.integer(ddim[1]),
                      as.integer(ddim[2]),
                      as.integer(ddim[3]),
                      as.integer(parammd$ind),
                      as.integer(nwmd),
                      double(nwmd*mc.cores), # work(nw,nthreds)
                      as.integer(mc.cores),
                      sigman = double(n),
                      PACKAGE = "dti")$sigman
    } 
    dim(sigma) <- ddim
    mask[sigma==0] <- FALSE
    if(verbose) cat("local median smoother in step ",i," completed",format(Sys.time()),"\n") 
    ##
    ##  diagnostics
    ##
    if(verbose){
      meds <- median(sigma[mask])
      means <- mean(sigma[mask])
      image(ymslice,col=grey(0:255/255))
      title(paste("S  max=",signif(max(y[mask]),3)," median=",signif(median(y[mask]),3)))
      image(th[,,mslice],col=grey(0:255/255))
      title(paste("E(S)  max=",signif(max(th[mask]),3)," median=",signif(median(th[mask]),3)))
      image(sigma[,,mslice],col=grey(0:255/255),zlim=c(0,max(sigma[mask])))
      title(paste("sigma max=",signif(max(sigma[mask]),3)," median=",signif(meds,3)))
      image(ni[,,mslice],col=grey(0:255/255))
      title(paste("Ni    max=",signif(max(ni[mask]),3)," median=",signif(median(ni[mask]),3)))
      plot(density(sigma[mask]),main="density of sigma")
      plot(density(ni[mask]),main="density of Ni")
      cat("mean sigma",means,"median sigma",meds,"sd sigma",sd(sigma[mask]),"\n")
      if(!is.null(u)&&"NCchi" == family){
        thchims <- fncchir(th/sigma,varstats)*sigma
        thchims[!mask] <- u[!mask]
        image(abs(thchims-u)[,,mslice],col=grey(0:255/255))
        title("abs Error in thchims")
        plot(density((thchims-u)[mask]),main="density of thchims-u")
        cat("MAE(th)",mean(abs(thchi-u)[mask]),"RMSE(th)",sqrt(mean((thchi-u)[mask]^2)),"MAE(thms)",mean(abs(thchims-u)[mask]),"RMSE(thms)",sqrt(mean((thchims-u)[mask]^2)),"\n")
      }
    } else {
      cat(" ",i)
    }
  }
  ## END PS iteration
  if(!verbose) cat("\n")
  thchi <- fncchir(th/sigma,varstats)*sigma
  thchi[!mask] <- 0
  ## this is the result (th is expectation, not the non-centrality parameter !!!)
  invisible(list(sigma = sigma,
                 sigmal = array(z$sigman,ddim),
                 theta = th, 
                 thchi = thchi,
                 ni  = ni,
                 mask = mask))
}
#
#
#      estimate variance parameter in a multicoil system
#
#
awssigmc <- function(y,                 # data
                     steps,             # number of iteration steps for PS
                     mask = NULL,       # data mask, where to do estimation
                     ncoils = 1,        # number of coils for parallel MR image acquisition
                     vext = c( 1, 1),   # voxel extensions
                     lambda = 20,       # adaptation parameter for PS
                     h0 = 2,            # initial bandwidth for first step in PS
                     verbose = FALSE, 
                     sequence = FALSE,  # return estimated sigma for intermediate steps of PS?
                     hadj = 1,          # adjust parameter for density() call for mode estimation
                     q = .25,  # for IQR
                     qni = .8,
                     method=c("VAR","MAD")  # for variance, alternative "MAD" for mean absolute deviation
) {
  method <- match.arg(method)
  ## some functions for pilot estimates
  IQQ <- function (x, q = .25, na.rm = FALSE, type = 7) 
    diff(quantile(as.numeric(x), c(q, 1-q), na.rm = na.rm, names = FALSE, type = type))
  
  IQQdiff <- function(y, mask, q = .25, verbose = FALSE) {
    cq <- qnorm(1-q)*sqrt(2)*2
    sx <- IQQ( diff(y[mask]), q)/cq
    sy <- IQQ( diff(aperm(y,c(2,1,3))[aperm(mask,c(2,1,3))]), q)/cq
    sz <- IQQ( diff(aperm(y,c(3,1,2))[aperm(mask,c(3,1,2))]), q)/cq
    if(verbose) cat( "Pilot estimates of sigma", sx, sy, sz, "\n")
    min( sx, sy, sz)
  }
  
  estsigma <- function(y, mask, q, L, sigma, verbose=FALSE){
    meany <- mean(y[mask]^2)
    eta <- sqrt(max( 0, meany/sigma^2-2*L))
    m <- sqrt(pi/2)*gamma(L+1/2)/gamma(L)/gamma(3/2)*hg1f1(-1/2,L,-eta^2/2)
    v <- max( .01, 2*L+eta^2-m^2)
    if(verbose) cat( eta, m, v, "\n")
    IQQdiff( y, mask, q)/sqrt(v)
  }
  
  varstats <- sofmchi(ncoils)
  if(length(vext)==3) vext <- vext[2:3]/vext[1]
  ## dimension and size of cubus
  ddim <- dim(y)
  n <- prod(ddim)
  
  ## test dimension
  if (length(ddim) != 3) stop("first argument should be a 3-dimentional array")
  
  ## check mask
  if (is.null(mask)) mask <- array(TRUE, ddim)
  if(length(mask) != n) stop("dimensions of data array and mask should coincide")
  
  ## initial value for sigma_0 
  # sigma <- sqrt( mean( y[mask]^2) / 2 / ncoils)
  sigma <- IQQdiff( y, mask, q, verbose=verbose)
  #  cat( "sigmahat1", sigma, "\n")
  sigma <- estsigma( y, mask, q, ncoils, sigma)
  #  cat( "sigmahat2", sigma, "\n")
  sigma <- estsigma( y, mask, q, ncoils, sigma)
  #  cat( "sigmahat3", sigma, "\n")
  sigma <- estsigma( y, mask, q, ncoils, sigma, verbose=verbose)
  #  cat( "sigmahat4", sigma,"\n")
  
  ## define initial arrays for parameter estimates and sum of weights (see PS)
  th <- array( 1, ddim)
  ni <- array( 1, ddim)
  #  y <- y/sigma # rescale to avoid passing sigma to awsvchi
  minlev <- sqrt(2)*gamma(ncoils+.5)/gamma(ncoils)
  if (sequence) sigmas <- lind <- minni <- numeric(steps)
  mc.cores <- setCores(,reprt=FALSE)
  ## iterate PS starting with bandwidth h0
  for (i in 1:steps) {
    
    h <- h0 * 1.25^((i-1)/3)
    nw <- prod(2*as.integer(h/c(1,vext))+1)
    #    cat("nw=",nw,"h/vext=",h/c(1,1/vext),"ih=",as.integer(h/c(1,1/vext)),"\n")
    param <- .Fortran("paramw3",
                      as.double(h),
                      as.double(vext),
                      ind=integer(3*nw),
                      w=double(nw),
                      n=as.integer(nw),
                      PACKAGE = "dti")[c("ind","w","n")]
    nw <- param$n
    param$ind <- param$ind[1:(3*nw)]
    dim(param$ind) <- c(3,nw)
    param$w   <- param$w[1:nw]    
    fncchi <- fncchiv(th/sigma,varstats)
    ## correction factor for variance of NC Chi distribution
    ## perform one step PS with bandwidth h
    if(method=="VAR"){
      z <- .Fortran("awsvchi",
                    as.double(y),        # data
                    as.double(th),       # previous estimates
                    ni = as.double(ni),
                    as.double(fncchi/2),
                    as.logical(mask),
                    as.integer(ddim[1]),
                    as.integer(ddim[2]),
                    as.integer(ddim[3]),
                    as.integer(param$ind),
                    as.double(param$w),
                    as.integer(nw),
                    as.double(lambda),
                    as.double(sigma),
                    th = double(n),
                    sy = double(n),
                    PACKAGE = "dti")[c("ni","th","sy")]
    } else {
      #    cat("n",n,prod(ddim[1:3]),"ly",length(y),"lth",length(th),"lni",length(ni),"lfns",length(fncchi),"lmask",length(mask),"nw",nw,"lind",length(param$ind),"lw",length(param$w),"mc.cores",mc.cores,"\n")
      z <- .Fortran("awsadchi",
                    as.double(y),        # y(n1,n2,n3)
                    as.double(th),       # th(n1,n2,n3)
                    ni = as.double(ni),  # ni(n1,n2,n3)
                    as.double(fncchi/2), # fns(n1,n2,n3)
                    as.logical(mask),    # mask(n1,n2,n3)
                    as.integer(ddim[1]), # n1
                    as.integer(ddim[2]), # n2
                    as.integer(ddim[3]), # n3
                    as.integer(param$ind), # ind(3,nw)
                    as.double(param$w), # w(nw)
                    as.integer(nw), # nw
                    as.double(lambda), # lambda
                    as.double(sigma), # sigma
                    double(nw*mc.cores), # wad(nw,nthreds)
                    as.integer(mc.cores), # nthreds
                    th = double(n), # thn(n1*n2*n3)
                    sy = double(n), # sy(n1*n2*n3)
                    PACKAGE = "dti")[c("ni","th","sy")]       
    }
    ## extract sum of weigths (see PS) and consider only voxels with ni larger then mean
    th <- z$th
    ni <- z$ni
    ni[!mask]<-1
    ind <- (ni > .9999*quantile(ni[ni>1],qni))#&(z$th>sigma*minlev)
    ## use correction factor for sd of NC Chi distribution 
    sy1 <- z$sy[ind]
    th1 <- th[ind]
    sy1 <- sy1/fncchis(th1/sigma,varstats)
    ## use the maximal mode of estimated local sd parameters, exclude largest values for better precision
    dsigma <- density( sy1[sy1>0], n = 4092, adjust = hadj, to = min( max(sy1[sy1>0]), median(sy1[sy1>0])*5) )
    sigma <- dsigma$x[dsigma$y == max(dsigma$y)][1]
    
    if (sequence) {
      sigmas[i] <- sigma
      lind[i] <- sum(ind)
      minni[i] <- min(ni[ind])
    }
    
    if (verbose) {
      plot(dsigma, main = paste( "estimated sigmas step", i, "h=", signif(h,3)))
      cat( "step", i, "h=", signif( h, 3), "quantiles of ni", signif( quantile(ni[ind]), 3), "mean", signif( mean(ni[ind]), 3), "\n")
      cat( "quantiles of sigma", signif( quantile(sy1[sy1>0]), 3), "mode", signif( sigma, 3), "\n")
    }
  }
  ## END PS iteration
  
  
  ## this is the result (th is expectation, not the non-centrality parameter !!!)
  invisible(list(sigma = if(sequence) sigmas else sigma,
                 theta = th, 
                 lind  = if(sequence) lind else sum(ind), 
                 minni  = if(sequence) minni else min(ni[ind])))
}

aflsigmc <- function(y,ncoils,level=NULL,mask=NULL,h=2,hadj=1,vext = c( 1, 1)){
  ##
  ##   estimate effective sigma and effective ncoils (L) according to Aja-Fernandez 2013
  ##
  #    sh2B stands for \hat{\sigma}^2_{n,B}
  #    sh2L stands for \hat{\sigma}^2_{nL}
  #    sh2eB stands for \hat{\sigma}^2_{eff,B}
  #    sh2eS stands for \hat{\sigma}^2_{eff,S}
  #    LeB stands for \hat{L}_{eff,B}
  #    vmlb  -  Var(M_L(x)|x \in B)
  #    vmlbx -  Var(M_L(x)|x \in B)_x
  #    phix  -  \hat{\Phi}_n(x)
  #    
  ddim <- dim(y)
  n <- prod(ddim)
  if(is.null(level)&is.null(mask)){
    warning("no background definition, need either level or mask")
    return(invisible(NULL))
  }
  if(is.null(level)) indB <- !mask else indB <- y<level
  sh2Lsimple <- mean(y[indB]^2)/2
  mask1 <- array(TRUE,ddim)
  ##
  ##  local variance estimates in vx
  ## 
  vx <- .Fortran("afmodevn",
                 as.double(y),
                 as.integer(ddim[1]),
                 as.integer(ddim[2]),
                 as.integer(ddim[3]),
                 as.logical(mask1),
                 as.double(h),
                 as.double(vext),
                 sigma = double(n),
                 PACKAGE = "dti")$sigma
  dim(vx) <- ddim
  vxb <- vx[indB]
  vxs <- vx[!indB]
  dv2b <- density( vxb[vxb>0], n = 4092, adjust = hadj, to = min( max(vxb[vxb>0]), median(vxb[vxb>0])*5) )
  sh2B <- dv2b$x[dv2b$y == max(dv2b$y)][1]
  dv2s <- density( vxs[vxs>0], n = 4092, adjust = hadj, to = min( max(vxs[vxs>0]), median(vxs[vxs>0])*5) )
  sh2eS <- dv2s$x[dv2s$y == max(dv2s$y)][1]
  m2 <- .Fortran("afmodem2",
                 as.double(y),
                 as.integer(ddim[1]),
                 as.integer(ddim[2]),
                 as.integer(ddim[3]),
                 as.logical(mask1),
                 as.double(h),
                 as.double(vext),
                 sm = double(n),
                 PACKAGE = "dti")$sm
  dim(m2) <- ddim
  m2b <- m2[indB]
  dm2b <- density( m2b[m2b>0], n = 4092, adjust = hadj, to = min( max(m2b[m2b>0]), median(m2b[m2b>0])*5) )
  sh2L <- dm2b$x[dm2b$y == max(dm2b$y)][1]/2
  cat("sh2B",sh2B,"sh2Lsimple",sh2Lsimple,"sh2L",sh2L,"\n")
  ##
  ##  now get sh2eB and LeB
  ##
  LeBn <- ncoils
  sh2eBn <- sh2B/2/(LeBn-gamma(LeBn+.5)^2/gamma(LeBn)^2)
  LeB <- 1
  sh2eB <- 0
  while(abs(sh2eB-sh2eBn)>1e-5||abs(LeB-LeBn)>1e-5){
    LeB <- LeBn
    sh2eB <- sh2eBn
    LeBn <- sh2L/sh2eB
    sh2eBn <- sh2B/2/(LeBn-gamma(LeBn+.5)^2/gamma(LeBn)^2)
    cat("LeB",LeB,"LeBn",LeBn,"sh2eB",sh2eB,"sh2eBn",sh2eBn,"\n")
  }
  sh2eB <- sh2eBn
  LeB <- sh2L/sh2eB
  phix <- pmax(0,pmin(1,sh2L/(m2-sh2L)))
  dim(phix) <- dim(m2)
  seff <- sqrt((1-phix)*sh2eS+phix*sh2eB)
  Leff <- sh2L/seff^2
  list(seff=seff,Leff=Leff,sh2L=sh2L,sh2B=sh2B,sh2eS=sh2eS,sh2eB=sh2eB,LeB=LeB,m2=m2,vx=vx,indB=indB)
}

afsigmc <- function(y,                 # data
                    level = NULL,             # threshold for background separation
                    mask = NULL,       # data mask, where to do estimation, needs to refer to background if level == NULL
                    ncoils = 1,        # number of coils for parallel MR image acquisition
                    vext = c( 1, 1),   # voxel extensions
                    h = 2,             # bandwidth for local averaging
                    verbose = FALSE,
                    hadj = 1,           # adjust parameter for density() call for mode estimation
                    method = c("modevn","modem1chi","bkm2chi","bkm1chi")# methods according to table 2 in Aja-Ferbnandez (2009)
) {
  #  for method=="modevn"  mask should refer to voxel within the head
  #  for all other methods to the background
  #
  method <- tolower(method)
  method <- match.arg(method)
  if(method!="modevn"&is.null(level)&is.null(mask)) {
    stop("need information on background using either level or mask")
  }
  ## dimension and size of cubus
  ddim <- dim(y)
  n <- prod(ddim)
  
  ## test dimension
  if (length(ddim) != 3) stop("first argument should be a 3-dimentional array")
  
  ## check mask
  if (is.null(mask)) mask <- array(TRUE, ddim)
  if(!is.null(level)){
    if (method=="modevn"){
      mask[y<level] <- FALSE
    } else {
      mask[y>level] <- FALSE
    }
  }
  if(length(mask) != n) stop("dimensions of data array and mask should coincide")
  
  ## let FORTRAN do the calculation
  if(method%in%c("modevn","modem1chi")){
    if(method=="modevn"){
      sigma <- .Fortran("afmodevn",
                        as.double(y),
                        as.integer(ddim[1]),
                        as.integer(ddim[2]),
                        as.integer(ddim[3]),
                        as.logical(mask),
                        as.double(h),
                        as.double(vext),
                        sigma = double(n),
                        PACKAGE = "dti")$sigma
      sigma <- sigma/2/(ncoils-gamma(ncoils+.5)^2/gamma(ncoils)^2)
      sigma <- array( sqrt(sigma), ddim)
    } else {
      afactor <- sqrt(1/2)*gamma(ncoils)/gamma(ncoils+.5)
      sigma <- .Fortran("afmodem1",
                        as.double(y),
                        as.integer(ddim[1]),
                        as.integer(ddim[2]),
                        as.integer(ddim[3]),
                        as.logical(mask),
                        as.double(h),
                        as.double(vext),
                        sigma = double(n),
                        PACKAGE = "dti")$sigma
      sigma <- array( afactor*sigma, ddim)
    }
    ##  use the maximal mode of estimated local variance parameters, exclude largest values for better precision
    dsigma <- density( sigma[sigma>0], n = 4092, adjust = hadj, to = min( max(sigma[sigma>0]), median(sigma[sigma>0])*5) )
    sigmag <- dsigma$x[dsigma$y == max(dsigma$y)][1]
    if(verbose){
      plot(dsigma, main = paste( "estimated sigmas h=", signif( h, 3)))
      cat("quantiles of sigma", signif( quantile(sigma[sigma>0]), 3), "mode", signif( sigmag, 3), "\n")
    }
  } else {
    if(method=="bkm2chi"){
      sigmag <- sqrt(mean(y[mask]^2)/2/ncoils)
    } else {
      sigmag <- mean(y[mask])*sqrt(ncoils/2)*gamma(ncoils)/gamma(ncoils+.5)/sqrt(ncoils)
    }
  }
  ## this is the estimate
  sigmag
}




#
#    R - function  aws  for likelihood  based  Adaptive Weights Smoothing (AWS)
#    for local constant Gaussian, Bernoulli, Exponential, Poisson, Weibull and  
#    Volatility models                                                         
#
#    emaphazises on the propagation-separation approach 
#
#    Copyright (C) 2006 Weierstrass-Institut fuer
#                       Angewandte Analysis und Stochastik (WIAS)
#
#    Author:  Joerg Polzehl
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
#  USA.
#
#     default parameters:  see function setawsdefaults
#       
awslinsd <- function(y,hmax=NULL,hpre=NULL,h0=NULL,mask=NULL,
                     ladjust=1,wghts=NULL,varprop=.1,A0,A1)
{
  #
  #    first check arguments and initialize
  #
  homogen <- TRUE
  wghts <- NULL
  
  args <- match.call()
  dy<-dim(y)
  if(length(dy)!=3) stop("Image should be 3D")
  #
  #   set appropriate defaults
  #
  if(is.null(wghts)) wghts <- c(1,1,1)
  wghts <- wghts[1]/wghts[2:3]
  cpar<-setawsdefaults(mean(y),ladjust,hmax,wghts)
  if(is.null(mask)) {
    if(length(dy)==0) mask <- rep(TRUE,length(y)) else mask <- array(TRUE,dy)
  }
  lambda <- cpar$lambda
  maxvol <- cpar$maxvol
  k <- cpar$k
  kstar <- cpar$kstar
  hmax <- cpar$hmax
  n<-length(y)
  # 
  #   family dependent transformations 
  #
  zfamily <- awsgfamily(y,h0,3)
  sigma2 <- zfamily$sigma2
  h0 <- zfamily$h0
  rm(zfamily)
  # now check which procedure is appropriate
  ##  this is the version on a grid
  n <- length(y)
  n1 <- dy[1]
  n2 <- dy[2]
  n3 <- dy[3]
  #
  #    Initialize  for the iteration
  #  
  zobj<-list(ai=y, bi= rep(1,n), theta= y, fix=!mask)
  mae<-NULL
  lambda0<-1e50 # that removes the stochstic term for the first step, initialization by kernel estimates
  #
  #   produce a presmoothed estimate to stabilze variance estimates
  #
  if(is.null(hpre)) hpre<-20^(1/3)
  dlw<-(2*trunc(hpre/c(1,wghts))+1)[1:3]
  hobj <- .Fortran("caws03d",as.double(y),
                   as.logical(mask),
                   as.integer(n1),
                   as.integer(n2),
                   as.integer(n3),
                   as.double(hpre),
                   theta=as.double(zobj$theta),
                   bi=as.double(zobj$bi),
                   double(prod(dlw)),
                   as.double(wghts),
                   PACKAGE="dti")[c("bi","theta")]
  dim(hobj$theta) <- dim(hobj$bi) <- dy
  #
  #   iteratate until maximal bandwidth is reached
  #
  #  cat("Progress:")
  #  total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
  pb <- txtProgressBar(0, kstar, style = 3)
  while (k<=kstar) {
    hakt0 <- gethani(1,10,1.25^(k-1),c(1,0,0,1,0,1),c(1,wghts),1e-4)
    hakt <- gethani(1,10,1.25^k,c(1,0,0,1,0,1),c(1,wghts),1e-4)
    dlw<-(2*trunc(hakt/c(1,wghts))+1)[1:3]
    if(any(h0>0)) lambda0<-lambda0*Spatialvar.gauss(hakt0/0.42445/4,h0,3)/Spatialvar.gauss(hakt0/0.42445/4,1e-5,3)
    # Correction for spatial correlation depends on h^{(k)} 
    hakt0<-hakt
    # heteroskedastic Gaussian case
    zobj <- .Fortran("cgaws",as.double(y),
                     as.logical(mask),
                     as.double(sigma2),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(n3),
                     hakt=as.double(hakt),
                     hhom=as.double(rep(1,n)),
                     as.double(lambda0),
                     as.double(zobj$theta),
                     bi=as.double(zobj$bi),
                     gi=double(n),
                     gi2=double(n),
                     theta=double(n),
                     double(prod(dlw)),
                     as.double(wghts),
                     PACKAGE="dti")[c("bi","hhom","theta","gi","gi2","hakt")]
    dim(zobj$theta)<-dim(zobj$gi)<-dim(zobj$gi2)<-dim(zobj$bi)<-dy
    hhom <- zobj$hhom
    #
    #    Calculate MAE and MSE if true parameters are given in u 
    #    this is for demonstration and testing for propagation (parameter adjustments) 
    #    only.
    #
    #   Prepare for next iteration
    #
    #
    #   Create new variance estimate
    #
    vobj <- awsgsigma2(y,mask,hobj,zobj,varprop,h0,A0,A1)
    sigma2 <- vobj$sigma2inv
    coef <- vobj$coef
    rm(vobj)
    lambda0<-lambda
    setTxtProgressBar(pb, k)
    #    if (max(total) >0) {
    #      cat(signif(total[k],2)*100,"% . ",sep="")
    #    }
    k <- k+1
    gc()
  }
  close(pb)
  # cat("\n")
  ###                                                                       
  ###            end iterations now prepare results                                                  
  ###                                 
  list(theta = zobj$theta, vcoef=coef, mask=mask)
}
###########################################################################
#
#   Auxialiary functions
#
############################################################################
#
#   transformations for Gaussian case with variance modelling
#
############################################################################
awsgfamily <- function(y,h0,d){
  sigma2 <- max(1,IQRdiff(as.vector(y))^2)
  if(any(h0)>0) sigma2<-sigma2*Varcor.gauss(h0)
  #  cat("Estimated variance: ", signif(sigma2,4),"\n")
  sigma2 <- rep(sigma2, length(y))
  dim(sigma2) <- dim(y)
  sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes 
  list(sigma2=sigma2,h0=h0)
}
############################################################################
#
#  estimate inverse of variances
#
############################################################################
awsgsigma2 <- function(y,mask,hobj,tobj,varprop,h0,thmin,thmax){
  ## specify sd to be linear to mean
  thrange <- range(y[mask])
  #    thmax <- thrange[2]-.2*diff(thrange)
  #    thmin <- thrange[1]+.1*diff(thrange)
  ind <- tobj$gi>1.5&mask&tobj$theta>thmin&tobj$theta<thmax
  absresid <- abs((y-tobj$theta)[ind]*tobj$gi[ind]/sqrt(tobj$gi[ind]^2-tobj$gi2[ind]))/.8
  #     absresid <- abs((y-tobj$theta)[ind])/.8
  theta <- tobj$theta[ind]
  wght <- (tobj$gi[ind]^2-tobj$gi2[ind])/tobj$gi[ind]^2
  coef <- coefficients(lm(absresid~theta,weights=wght^2))
  # force positive variance for positive mean by increasing variance estimate
  if(coef[2] < 0){
    coef[1] <- coefficients(lm(absresid~1,weights=wght^2))
    coef[2] <- 0
  }
  if(coef[1] <0.5){
    coef[2] <- coefficients(lm(absresid~theta-1,weights=wght^2))
    coef[1] <- 0
  }
  gamma <- pmin(tobj$gi/hobj$bi,1)
  theta <- gamma*tobj$theta+(1-gamma)*hobj$theta
  sigma2 <- (coef[1]+coef[2]*theta)^2
  varquantile <- varprop*mean(sigma2[mask])
  sigma2 <- pmax(sigma2,varquantile)
  #  cat("Estimated mean variance",signif(mean(sigma2),3)," Variance parameters:",signif(coef,3),"\n")
  list(sigma2inv=1/sigma2,coef=coef)
}
#######################################################################################
#
#        Auxilary functions
#
#######################################################################################
#
#        Set default values
#
# default values for lambda and tau are chosen by propagation condition
# (strong version) with alpha=0.05 (Gaussian) and alpha=0.01 (other family models)
# see script aws_propagation.r
#
#######################################################################################
setawsdefaults <- function(meany,ladjust,hmax,wghts){
  qlambda <- .98
  if(is.null(hmax)) hmax <- 5
  lambda <- ladjust*qchisq(qlambda,1)*2
  maxvol <- getvofh(hmax,c(1,0,0,1,0,1),c(1,wghts))
  kstar <- as.integer(log(maxvol)/log(1.25))
  k <- 6
  cat("Estimating variance model using PS with with lambda=",signif(lambda,3)," hmax=",hmax,"number of iterations:",kstar-k+1,"\n")
  list(lambda=lambda,hmax=hmax,kstar=kstar,maxvol=maxvol,k=k,wghts=wghts)
}
#######################################################################################
#
#    IQRdiff (for robust variance estimates
#
#######################################################################################
IQRdiff <- function(y) IQR(diff(y))/1.908
