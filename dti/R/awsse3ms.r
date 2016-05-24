# This file contains the implementation of dti.smooth() for 
# "dtiData" Adaptive smoothing in SE(3) considering b=0 as an individual shell

dwi.smooth.ms <- function(object, ...) cat("No DTI smoothing defined for this class:",class(object),"\n")

setGeneric("dwi.smooth.ms", function(object, ...) standardGeneric("dwi.smooth.ms"))

setMethod("dwi.smooth.ms", 
          "dtiData", 
          function(object,
                   kstar,
                   lambda = 12,
                   kappa0 = .5,
                   ncoils = 1,
                   sigma = NULL,
                   ws0 = 1,
                   level = NULL,
                   mask = NULL,
                   xind = NULL,
                   yind = NULL,
                   zind = NULL,
                   verbose = FALSE,
                   usemaxni = TRUE,
                   memrelease = TRUE){
            
            # make the call part of the object
            args <- sys.call(-1)
            args <- c(object@call,args)
            
            # we need a lot of object properties
            ddim <- object@ddim
            ngrad <- object@ngrad
            s0ind <- object@s0ind
            ns0 <- length(s0ind)
            ngrad <- ngrad - ns0
            grad <- object@gradient[,-s0ind]
            bvalues <- object@bvalue[-s0ind]
            sdcoef <- object@sdcoef
            level <- object@level
            vext <- object@voxelext[2:3]/object@voxelext[1]
            
            # determine the mesh of gradient triangles for interpolation
            # determine the number of shells, hmmm: difficult task!
            # this requires unique values for the bvalues
            # cf. the corresponding code in multishell.r:getnext3g
            msstructure <- getnext3g(grad, bvalues)
            bv <- msstructure$bv
            nshell <- as.integer(msstructure$nbv)
            ubv <- msstructure$ubv
            
            # consider the noise
            varstats <- sofmchi(ncoils)
            if (length(sigma) == 1) {
              sigmacase <- 1
              if (verbose) cat("using supplied sigma ", sigma, "\n")
            } else if (identical(dim(sigma), ddim)) {
              sigmacase <- 2
              if (verbose) cat("using supplied array of sigma\n")
            } else if (identical(dim(sigma), c(ddim, nshell+1L))) {
              sigmacase <- 3
              if (is.null(mask)) mask <- getmask(object, level)$mask
              if (verbose) cat("using supplied array of sigma per shell\n")
            } else {
              if (is.null(mask)) mask <- getmask(object, level)$mask
              sigma <- numeric(object@ngrad)
              for(i in 1:object@ngrad){
                ## FOR FUTURE USE:
                # awslsigmc(object@si[ , , , i], 12, ncoils = ncoils, lambda = 5, verbose = verbose, hsig = 5, mask = mask)
                ## END
                sigma[i] <- awssigmc(object@si[,,,i], 12, mask, ncoils, vext, h0=1.25, verbose=verbose)$sigma
                if (verbose) cat("image ", i," estimated sigma", sigma[i], "\n")
              }
              if (verbose) cat("quantiles of estimated sigma values", quantile(sigma), "\n")
              sigma <- median(sigma)
              if (verbose) cat("using median estimated sigma", sigma,"\n")
              sigmacase <- 1
            }
            
            
            # 
            sb <- object@si[,,,-s0ind]
            s0 <- object@si[,,,s0ind]
            if(is.null(kappa0)){
              #  select kappa based on variance reduction on the sphere
              warning("You need to specify  kappa0  returning unsmoothed object")
              return(object)
            }
            #
            #  rescale so that we have Chi-distributed values
            #
            if (sigmacase == 1) {
              sb <- sb/sigma
              s0 <- s0/sigma
            } else if (sigmacase == 2) {
              sb <- sweep(sb, 1:3, sigma, "/")
              s0 <- sweep(s0, 1:3, sigma, "/")
            } else { # now is per shell
              s0 <- sweep(s0, 1:3, sigma[, , , 1], "/")
              for (shnr in 1:nshell) {
                indbv <- (1:length(bv))[bv == ubv[shnr]]
                sb[, , , indbv] <- sweep(sb[, , , indbv], 1:3, sigma[, , , shnr+1], "/")
              }
            }
            if(ns0>1){
              dim(s0) <- c(prod(ddim),ns0)
              s0 <- s0%*%rep(1/sqrt(ns0),ns0)
              #  make sure linear combination of s0 has same variance as original 
              dim(s0) <- ddim
            }
            if (is.null(mask)) mask <- getmask(object, level)$mask
            # determine minimal subcube that contains the mask ore use supplied information
            if(is.null(xind)) xind <- (1:ddim[1])[apply(mask,1,any)]
            if(is.null(yind)) yind <- (1:ddim[2])[apply(mask,2,any)]
            if(is.null(zind)) zind <- (1:ddim[3])[apply(mask,3,any)]
            # reduce everything to this subcube 
            ##
            ##  check memory size needed for largest vector
            ##
            vsize <- (nshell+1)*ngrad*length(xind)*length(yind)*length(zind)
            if(vsize>2^31-1){
              cat("region specified is to large, a vector of size",vsize,"needs to be passed through
      the .Fortran, this is limited to 2^31-1 in R\n reducing the zind")
              maxzind <- (2^31-1)/(nshell+1)/ngrad/length(xind)/length(yind)
              nzind <- apply(mask,3,sum)
              
            }
            sb <- sb[xind, yind, zind,]
            s0 <- s0[xind, yind, zind]
            mask <- mask[xind, yind, zind]
            if (length(sigma) > 1) {
              if(length(dim(sigma)) == 3) {
                sigma <- sigma[xind, yind, zind]
              } else {
                sigma <- sigma[xind, yind, zind, ]
              }
            }
            gc()
            ddim0 <- dim(s0)
            ddimb <- dim(sb)
            gradstats <- getkappasmsh3(grad, msstructure)
            hseq <- gethseqfullse3msh(kstar,gradstats,kappa0,vext=vext)
            nind <- as.integer(hseq$n*1.25)
            # make it nonrestrictive for the first step
            z <- list(th=array(1,ddimb), th0=array(1,ddim0), ni = array(1,ddimb), ni0 = array(1,ddim0))
            if(usemaxni){
              ni <- array(1,ddimb)
              ni0 <- array(1,ddim0)
            }
            prt0 <- Sys.time()
            cat("adaptive smoothing in SE3, kstar=",kstar,if(verbose)"\n" else " ")
            kinit <- if(lambda<1e10) 0 else kstar
            mc.cores <- setCores(,reprt=FALSE)
            gc()
            for(k in kinit:kstar){
              hakt <- hseq$h[,k+1]
              t0 <- Sys.time()
              thnimsh <- interpolatesphere1(z$th,z$th0,z$ni,z$ni0,msstructure,mask)
              if(memrelease) rm(z)
              gc()
              t1 <- Sys.time()
              param <- lkfullse3msh(hakt,kappa0/hakt,gradstats,vext,nind) 
              hakt0 <- mean(hakt)
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
                            as.double(vs2),#var/2
                            as.double(vs02),#var/2 for s0
                            as.logical(mask),#mask
                            as.integer(nshell+1),#ns number of shells
                            as.integer(ddim0[1]),#n1
                            as.integer(ddim0[2]),#n2
                            as.integer(ddim0[3]),#n3
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
                            th=double(prod(ddimb)),#thn
                            ni=double(prod(ddimb)),#nin
                            th0=double(prod(ddim0)),#th0n
                            ni0=double(prod(ddim0)),#ni0n
                            double(ngrad*mc.cores),#sw
                            double(ngrad*mc.cores),#swy
                            double((nshell+1)*mc.cores),#thi
                            double((nshell+1)*mc.cores),#nii
                            double((nshell+1)*mc.cores),#fsi2
                            PACKAGE="dti")[c("ni","th","ni0","th0")]
              t3 <- Sys.time()
              gc()
              if(memrelease) rm(thnimsh,vs2,vs02)
              gc()
              if(usemaxni){
                ni <- z$ni <- if(usemaxni) pmax(ni,z$ni)
                ni0 <- z$ni0 <- if(usemaxni) pmax(ni0,z$ni0)
              }
              dim(z$th) <- ddimb
              dim(z$th0) <- ddim0
              dim(z$ni) <- ddimb
              dim(z$ni0) <- ddim0
              if(verbose){
                dim(z$ni) <- c(prod(ddim0),ngrad)
                cat("k:",k,"h_k:",signif(max(hakt),3)," quartiles of ni",signif(quantile(z$ni[mask,]),3),
                    "mean of ni",signif(mean(z$ni[mask,]),3),"\n              quartiles of ni0",signif(quantile(z$ni0[mask]),3),
                    "mean of ni0",signif(mean(z$ni0[mask]),3),
                    " time elapsed:",format(difftime(Sys.time(),prt0),digits=3),"\n")
                cat("interpolation:",format(difftime(t1,t0),digits=3),
                    "param:",format(difftime(t2,t1),digits=3),
                    "smoothing:",format(difftime(t3,t2),digits=3),"\n")
              } else {
                cat(".")
              }
            }
            ngrad <- ngrad+1
            #
            #  one s0 image only
            #
            si <- array(object@si,c(ddim,ngrad))
            #
            #  back to original scale
            #
            if (length(sigma) > 1) {
              if(length(dim(sigma)) == 3) {
                si[xind, yind, zind, 1] <-  z$th0/sqrt(ns0)*sigma
              } else {
                si[xind, yind, zind, 1] <-  z$th0/sqrt(ns0)*sigma[, , , 1]
              }
            } else {
              si[xind, yind, zind, 1] <-  z$th0/sqrt(ns0)*sigma
            }
            #  go back to original s0 scale
            #  for the DWI we need to scale back differently
            if (sigmacase == 1) {
              si[xind, yind, zind, -1] <- z$th*sigma
            } else if (sigmacase == 2) {
              si[xind, yind, zind, -1] <- sweep(z$th, 1:3, sigma, "*")
            } else if (sigmacase == 3) { # now sigma is an array with (identical(dim(sigma), ddim) == TRUE)
              for (shnr in 1:nshell) {
                indth <- (1:length(bv))[bv == ubv[shnr]]
                indsi <- indth+1
                si[xind, yind, zind, indsi] <- sweep(z$th[, , , indth], 1:3, sigma[, , , shnr+1], "*")
                }
            } else {
                si[xind, yind, zind, -1] <- z$th*sigma[,,,-1]
            }
            object@si <-  si
            object@gradient <- grad <- cbind(c(0,0,0),grad)
            object@bvalue <- bvalue <- c(0,object@bvalue[-object@s0ind])
            object@btb <- sweep(create.designmatrix.dti(grad), 2, bvalue, "*")
            object@s0ind <- as.integer(1)
            object@replind <- as.integer(1:ngrad)
            object@ngrad <- as.integer(ngrad)
            object@call <- args
            object
          }
)


getkappas3 <- function(grad, trace = 0){
  #
  #  dist: acos(g_i%*%g_j)
  #
  ngrad <- dim(grad)[2]
  kappa456 <- array(0, c(3, ngrad, ngrad))
  bghat <- betagamma(grad)
  for(i in 1:ngrad) kappa456[1,i,] <- acos(pmin(1,abs(grad[,i]%*%grad)))
  list(k456 = kappa456, bghat = bghat, dist=3)
}
getkappasmsh3 <- function(grad, msstructure, trace = 0){
  ngrad <- dim(grad)[2]
  nbv <- msstructure$nbv
  bv <- msstructure$bv
  ubv <- msstructure$ubv
  bvind <- k456 <- bghat <- list(NULL)
  for(i in 1:nbv){
    #
    #   collect information for spherical distances on each schell separately
    #
    ind <- (1:ngrad)[bv==ubv[i]]
    z <- getkappas3(grad[,ind], trace = trace)
    bvind[[i]] <- ind
    k456[[i]] <- z$k456
    bghat[[i]] <- z$bghat
  }
  list(k456 = k456, bghat = bghat, bvind = bvind, dist=3, nbv = nbv, ngrad = ngrad)
}
interpolatesphere0 <- function(theta,th0,ni,ni0,n3g,mask){
  ##  interpolate estimated thetas to get values on all spheres
  ##  n3g  generated by function  getnext3g
  ##  dim(theta) = c(n1,n2,n3,ngrad)
  nbv <- n3g$nbv
  bv <- n3g$bv
  ubv <- n3g$ubv
  dtheta <- dim(theta)
  dth0 <- dim(th0)
  ng <- dtheta[4]
  n <- prod(dtheta[1:3])
  nmask <- sum(mask)
  dim(theta)  <- dim(ni)  <- c(n,ng)
  mstheta <- msni <- array(0,c(nbv+1,n,ng))
  t1 <- Sys.time()
  z <- .Fortran("ipolsp",
                as.double(theta[mask,]),
                as.double(th0[mask]),
                as.double(ni[mask,]),
                as.double(ni0[mask]),
                as.integer(nmask),
                as.integer(ng),
                as.integer(n3g$ind),
                as.double(n3g$w),
                as.integer(nbv),
                as.integer(nbv+1),
                msth=double((nbv+1)*nmask*ng),
                msni=double((nbv+1)*nmask*ng),
                PACKAGE="dti")[c("msth","msni")]
  cat("time for sb-interpolation", format(difftime(Sys.time(),t1),digits=3),"\n")
  mstheta[,mask,] <- z$msth
  msni[,mask,] <- z$msni
  #  now fill vector for s0
  msth0 <- msni0 <- array(0,c(nbv+1,n))
  msth0[1,] <- th0
  msni0[1,] <- ni0
  for(i in 1:nbv){
    indi <- (1:ng)[bv==ubv[i]]
    lindi <- length(indi)
    msth0[i+1,mask] <- theta[mask,indi]%*%rep(1/lindi,lindi)
    #  correct value would be 
    #  msni0[i+1,] <- 1/(1/(ni[,indi])%*%rep(1/lindi,lindi)^2)
    #  try to be less conservative by ignorin squares in w
    msni0[i+1,mask] <- 1/((1/ni[mask,indi])%*%rep(1/lindi,lindi))   
  }
  dim(msth0) <- dim(msni0) <- c(nbv+1,dth0)
  dim(mstheta) <- dim(msni) <- c(nbv+1,dtheta)
  list(mstheta=mstheta,msni=msni,msth0=msth0,msni0=msni0)
}

interpolatesphere1 <- function(theta,th0,ni,ni0,n3g,mask){
  ##  interpolate estimated thetas to get values on all spheres
  ##  n3g  generated by function  getnext3g
  ##  dim(theta) = c(n1,n2,n3,ngrad)
  nbv <- n3g$nbv
  bv <- n3g$bv
  ubv <- n3g$ubv
  dtheta <- dim(theta)
  dth0 <- dim(th0)
  ng <- dtheta[4]
  n <- prod(dtheta[1:3])
  #t1 <- Sys.time()
  z <- .Fortran("ipolsp1",
                as.double(theta),
                as.double(th0),
                as.double(ni),
                as.double(ni0),
                as.logical(mask),
                as.integer(n),
                as.integer(ng),
                as.integer(n3g$ind),
                as.double(n3g$w),
                as.integer(nbv),
                as.integer(nbv+1),
                msth=double((nbv+1)*n*ng),
                msni=double((nbv+1)*n*ng),
                PACKAGE="dti")[c("msth","msni")]
  #cat("time for sb-interpolation", format(difftime(Sys.time(),t1),digits=3),"\n")
  #  now fill vector for s0
  dim(theta)  <- dim(ni)  <- c(n,ng)
  nmask <- sum(mask)
  msth0 <- msni0 <- array(0,c(nbv+1,n))
  msth0[1,] <- th0
  msni0[1,] <- ni0
  for(i in 1:nbv){
    indi <- (1:ng)[bv==ubv[i]]
    lindi <- length(indi)
    #   msth0[i+1,mask] <- theta[mask,indi]%*%rep(1/lindi,lindi)
    msth0[i+1,mask] <- .Fortran("getmsth0",
                                as.double(theta[mask,indi]),
                                as.integer(nmask),
                                as.integer(lindi),
                                msth0=double(nmask),
                                PACKAGE="dti")$msth0
    #  correct value would be 
    #  msni0[i+1,] <- 1/(1/(ni[,indi])%*%rep(1/lindi,lindi)^2)
    #  try to be less conservative by ignorin squares in w
    #   msni0[i+1,mask] <- 1/((1/ni[mask,indi])%*%rep(1/lindi,lindi))   
    msni0[i+1,mask] <- .Fortran("getmsni0",
                                as.double(ni[mask,indi]),
                                as.integer(nmask),
                                as.integer(lindi),
                                msni0=double(nmask),
                                PACKAGE="dti")$msni0
  }
  dim(msth0) <- dim(msni0) <- c(nbv+1,dth0)
  dim(z$msth) <- dim(z$msni) <- c(nbv+1,dtheta)
  list(mstheta=z$msth,msni=z$msni,msth0=msth0,msni0=msni0)
}

lkfulls0 <- function(h,vext,n){
  z <- .Fortran("lkfuls0",
                as.double(h),
                as.double(vext),
                ind=integer(3*n),
                w=double(n),
                n=as.integer(n),
                PACKAGE="dti")[c("ind","w","n")]
  dim(z$ind) <- c(3,n)
  list(h=h,ind=z$ind[,1:z$n],w=z$w[1:z$n],nind=z$n)
}


