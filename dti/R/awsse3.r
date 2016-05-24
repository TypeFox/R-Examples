# This file contains the implementation of dti.smooth() for 
# "dtiData" Adaptive smoothing in SE(3)

dwi.smooth <- function(object, ...) cat("No DTI smoothing defined for this class:",class(object),"\n")

setGeneric("dwi.smooth", function(object, ...) standardGeneric("dwi.smooth"))

setMethod("dwi.smooth", "dtiData", function(object,kstar,lambda=20,kappa0=NULL,ncoils=1,sigma=NULL,level=NULL,
                                            vred=4,xind=NULL,yind=NULL,zind=NULL,verbose=FALSE,dist=1,model=
                                              c("Gapprox","Gapprox2","Chi","Chi2")){
  ## check model
  model <- match.arg(model)  
  args <- sys.call(-1)
  args <- c(object@call,args)
  sdcoef <- object@sdcoef
  level <- object@level
  vext <- object@voxelext[2:3]/object@voxelext[1]
  if(length(sigma)==1) {
    cat("using supplied sigma",sigma,"\n")
  } else {
    mask <- getmask(object,level)$mask
    sigma <- numeric(object@ngrad)
    for(i in 1:object@ngrad){
      sigma[i] <- awssigmc(object@si[,,,i],12,mask,ncoils,vext,h0=1.25,verbose=verbose)$sigma
      cat("image ",i," estimated sigma",sigma[i],"\n")
    }
    cat("quantiles of estimated sigma values",quantile(sigma),"\n")
    sigma <- median(sigma)
    cat("using median estimated sigma",sigma,"\n")
  }
  model <- switch(model,"Chi2"=1,"Chi"=0,"Gapprox2"=2,"Gapprox"=3,1)
  if(model>=2) varstats <- sofmchi(ncoils)
  #
  #  Chi2 uses approx of noncentral Chi^2 by rescaled central Chi^2 with adjusted DF 
  #        and smoothes Y^2
  #  Chi uses approx of noncentral Chi^2 by rescaled central Chi^2 with adjusted DF  
  #        and smoothes Y
  #   uses approximation of noncentral Chi by a Gaussian (with correct mean and variance)
  #        and smoothes Y
  
  #
  if(!(is.null(xind)&is.null(yind)&is.null(zind))){
    if(is.null(xind)) xind <- 1:object@ddim[1]
    if(is.null(yind)) yind <- 1:object@ddim[2]
    if(is.null(zind)) zind <- 1:object@ddim[3]
    object <- object[xind,yind,zind]
  }
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  ngrad <- ngrad - ns0
  grad <- object@gradient[,-s0ind]
  bvalues <- object@bvalue[-s0ind]
  multishell <- sd(bvalues) > mean(bvalues)/50
  if(multishell) {
    msstructure <- getnext3g(grad,bvalues)
    #     save(grad,bvalues,msstructure,file="msstructure.rsc")
    model <- 3
    nshell <- msstructure$nbv
  }
  sb <- object@si[,,,-s0ind]
  s0 <- object@si[,,,s0ind]
  if(is.null(kappa0)){
    #  select kappa based on variance reduction on the sphere
    if(is.null(vred)) {
      warning("You need to specify either kappa0 or vred\n returning unsmoothed object")
      return(object)
    }
    if(!is.numeric(vred)|vred<1){
      warning("vred needs to be >= 1\n returning unsmoothed object")
      return(object)
    }
    kappa0 <- suggestkappa(grad,vred,dist)$kappa
  }
  #
  #  rescale so that we have Chi-distributed values
  #
  if(ns0>1){
    dim(s0) <- c(prod(ddim),ns0)
    if(model%in%c(1,2)){
      s0 <- sqrt(s0^2%*%rep(1,ns0))        
      mask <- s0 > sqrt(ns0)*level
    } else {
      s0 <- s0%*%rep(1,ns0)
      mask <- s0 > ns0*level
    }
    # S0 will be noncentral Chi with 2*ns0*ncoils DF for model 1 and 2
    dim(s0) <- ddim
  } else {
    mask <- s0 > level
  }
  sb <- sb/sigma
  s0 <- s0/sigma
  if(model==1){
    #
    #   use squared values for Chi^2
    #
    sb <- sb^2
    s0 <- s0^2
    sigma <- sigma^2
  }
  th0 <- s0
  ni0 <- array(1,ddim)
  if(multishell){
    gradstats <- getkappasmsh(grad, msstructure,dist=dist)
    #     save(gradstats,file="gradstats.rsc")
    hseq <- gethseqfullse3msh(kstar,gradstats,kappa0,vext=vext)
  } else {
    gradstats <- getkappas(grad,dist=dist)
    hseq <- gethseqfullse3(kstar,gradstats,kappa0,vext=vext)
    hseq$h <- cbind(rep(1,ngrad),hseq$h)
  }
  nind <- as.integer(hseq$n*1.25)#just to avoid n being to small due to rounding
  hseq <- hseq$h
  # make it nonrestrictive for the first step
  ni <- array(1,dim(sb))
  minlevel <- if(model==1) 2*ncoils else sqrt(2)*gamma(ncoils+.5)/gamma(ncoils)  
  minlevel0 <- if(model==1) 2*ns0*ncoils else sqrt(2)*gamma(ns0*ncoils+.5)/gamma(ns0*ncoils)
  #  z <- list(th=pmax(sb,minlevel), ni = ni)
  z <- list(th=array(minlevel,dim(sb)), ni = ni)
  #  th0 <- pmax(s0,minlevel0)
  th0 <- array(minlevel0,dim(s0))
  prt0 <- Sys.time()
  cat("adaptive smoothing in SE3, kstar=",kstar,if(verbose)"\n" else " ")
  kinit <- if(lambda<1e10) 0 else kstar
  mc.cores <- setCores(,reprt=FALSE)
  for(k in kinit:kstar){
    gc()
    hakt <- hseq[,k+1]
    if(multishell){
      thmsh <- interpolatesphere(z$th,msstructure)
      param <- lkfullse3msh(hakt,kappa0/hakt,gradstats,vext,nind) 
      #       save(z,thmsh,param,msstructure,file=paste("thmsh",k,".rsc",sep=""))
      if(length(sigma)==1) {
        z <- .Fortran("adsmse3m",
                      as.double(sb),#y
                      as.double(thmsh),#th
                      ni=as.double(z$ni),#ni
                      as.double(fncchiv(thmsh,varstats)),#sthi
                      as.logical(mask),#mask
                      as.integer(nshell),# number of shells
                      as.integer(ddim[1]),#n1
                      as.integer(ddim[2]),#n2
                      as.integer(ddim[3]),#n3
                      as.integer(ngrad),#ngrad
                      as.double(lambda),#lambda
                      as.integer(ncoils),#ncoils
                      as.integer(mc.cores),#ncores
                      as.integer(param$ind),#ind
                      as.double(param$w),#w
                      as.integer(param$n),#n
                      th=double(prod(ddim)*ngrad),#thn
                      double(ngrad*mc.cores),#sw
                      double(ngrad*mc.cores),#swy
                      double(nshell*mc.cores),#si
                      double(nshell*mc.cores),#thi
                      PACKAGE="dti")[c("ni","th")]
        dim(z$th) <- dim(z$ni) <- c(ddim,ngrad)
        gc()
      } else {
        warning("not yet implemented for heterogenious variances\n
             returning original object")
        return(object)
      }
    } else {
      param <- lkfullse3(hakt,kappa0/hakt,gradstats,vext,nind) 
      if(length(sigma)==1) {
        z <- .Fortran("adsmse3p",
                      as.double(sb),
                      as.double(z$th),
                      ni=as.double(z$ni/if(model>=2) fncchiv(z$th,varstats) else 1),
                      as.logical(mask),
                      as.integer(ddim[1]),
                      as.integer(ddim[2]),
                      as.integer(ddim[3]),
                      as.integer(ngrad),
                      as.double(lambda),
                      as.integer(ncoils),
                      as.integer(mc.cores),
                      as.integer(param$ind),
                      as.double(param$w),
                      as.integer(param$n),
                      th=double(prod(ddim)*ngrad),
                      double(prod(ddim)*ngrad),#ldf (to precompute lgamma)
                      double(ngrad*mc.cores),
                      double(ngrad*mc.cores),
                      as.integer(model),
                      as.double(minlevel),
                      PACKAGE="dti")[c("ni","th")]
        #     save(z,sb,thk,ni,param,model,ddim,lambda,ngrad,model,minlevel,ncoils,mask,file="tmpi.rsc")
        gc()
      } else {
        warning("not yet implemented for heterogenious variances\n
             returning original object")
        return(object)
      }
    }
    if(verbose){
      dim(z$ni)  <- c(prod(ddim),ngrad)
      cat("k:",k,"h_k:",signif(max(hakt),3)," quartiles of ni",signif(quantile(z$ni[mask,]),3),
          "mean of ni",signif(mean(z$ni[mask,]),3),
          " time elapsed:",format(difftime(Sys.time(),prt0),digits=3),"\n")
    } else {
      cat(".")
    }
    param <- reduceparam(param)
    z0 <- .Fortran("asmse30p",
                   as.double(s0),
                   as.double(th0),
                   ni=as.double(ni0/if(model==2) fncchiv(th0,varstats) else 1),
                   as.logical(mask),
                   as.integer(ddim[1]),
                   as.integer(ddim[2]),
                   as.integer(ddim[3]),
                   as.double(lambda),
                   as.integer(ncoils*ns0),
                   as.integer(param$ind),
                   as.double(param$w),
                   as.integer(param$n),
                   as.integer(param$starts),
                   as.integer(param$nstarts),
                   th0=double(prod(ddim)),
                   double(prod(ddim)),#ldf (to precompute lgamma)
                   double(param$nstarts),#swi
                   as.integer(model),
                   as.double(minlevel0),
                   PACKAGE="dti")[c("ni","th0")]
    #     save(z0,s0,th0,ni0,param,model,ddim,lambda,ns0,model,minlevel,ncoils,mask,file="tmp.rsc")
    th0 <- z0$th0
    ni0 <- z0$ni
    rm(z0)
    gc()
    if(verbose){
      cat("End smoothing s0: quartiles of ni",signif(quantile(ni0[mask]),3),
          "mean of ni",signif(mean(ni0[mask]),3),
          " time elapsed:",format(difftime(Sys.time(),prt0),digits=3),"\n")
    }
  }
  ngrad <- ngrad+1
  si <- array(0,c(ddim,ngrad))
  #
  #  back to original scale
  #
  s0factor <- switch(model+1,ns0,ns0,sqrt(ns0),ns0)
  #cat("sigma",sigma,"s0factor",s0factor,"minlevel0",minlevel,"maxth0",max(th0),"lth0",length(th0),"\n")
  #  si[,,,1] <-  pmax(th0,minlevel0)*sigma/s0factor
  #  si[,,,-1] <- pmax(z$th,minlevel)*sigma
  bvalue <- c(0,object@bvalue[-object@s0ind])
  si[,,,1] <-  th0*sigma/s0factor
  si[,,,-1] <- z$th*sigma
  object@si <- if(model==1) sqrt(si) else si
  object@gradient <- grad <- cbind(c(0,0,0),grad)
  object@bvalue <- bvalue
  object@btb <- sweep(create.designmatrix.dti(grad), 2, bvalue, "*")
  object@s0ind <- as.integer(1)
  object@replind <- as.integer(1:ngrad)
  object@ngrad <- as.integer(ngrad)
  object@call <- args
  object
}
)


lkfullse3 <- function(h,kappa,gradstats,vext,n){
  ngrad <- dim(gradstats$bghat)[2]
  if(length(h)<ngrad) h <- rep(h[1],ngrad)
  z <- .Fortran("lkfulse3",
                as.double(h),
                as.double(kappa),
                as.double(gradstats$k456),
                as.integer(ngrad),
                as.double(vext),
                ind=integer(5*n),
                w=double(n),
                n=as.integer(n),
                as.integer(gradstats$dist),
                PACKAGE="dti")[c("ind","w","n")]
  dim(z$ind) <- c(5,n)
  list(h=h,kappa=kappa,ind=z$ind[,1:z$n],w=z$w[1:z$n],nind=z$n)
}


