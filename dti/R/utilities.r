################################################################
#                                                              #
# Section for Utility functions                                #
#                                                              #
################################################################ 

sdpar <- function(object,  ...) cat("No method defined for class:",class(object),"\n")

setGeneric("sdpar", function(object,  ...) standardGeneric("sdpar"))

setMethod("sdpar", "dtiData", function(object,
                                       level = NULL,
                                       sdmethod = "none",
                                       interactive = TRUE,
                                       threshfactor = 1) {
  # determine interval of linearity
  if(!(sdmethod%in%c("none","sd","mad"))){
    warning("sdmethod needs to be either 'none','sd' or 'mad'")
    return(object)
  }
  if(prod(object@ddim)==1){
    warning("you need more than one voxel to model variances")
    return(object)
  }
  sdcoef <- object@sdcoef
  level0 <- if(is.null(level)) object@level else max(0,level)
  s0ind<-object@s0ind
  s0 <- object@si[,,,s0ind,drop=FALSE]
  ls0ind <- length(s0ind)
  A0 <- level0
  if(ls0ind>1) {
    dim(s0) <- c(prod(object@ddim),ls0ind)
    s0mean <- s0%*%rep(1/ls0ind,ls0ind)
    A1 <- quantile(s0mean[s0mean>0],.98)
    dim(s0mean) <- object@ddim
  } else {
    dim(s0) <- object@ddim
    A1 <- quantile(s0[s0>0],.98)
  }
  if(interactive) {
    oldpar <- par( mfrow = c( 1, 3), mar = c( 3, 3, 3, 1), mgp = c( 2, 1, 0))
    img <- if(ls0ind>1) s0mean[,,(object@ddim[3]-1)%/%2+1] else s0[,,(object@ddim[3]-1)%/%2+1]
    maximg <- max(img)
    accept <- FALSE
    ddim <- object@ddim
    ddm1 <- ddim-1
    bw <- min(bw.nrd(if(ls0ind>1) s0mean[s0mean>min(s0mean)] else s0[s0>0]),diff(range(if(ls0ind>1) s0mean else s0))/256)
    z <- density(if(ls0ind>1) s0mean[s0mean>0&s0mean<A1] else s0[s0>0&s0<A1],bw = max(bw,.01),,n=1024)
    indx1 <- trunc(0.05*ddm1[1]):trunc(0.95*ddm1[1])+1
    indx2 <- trunc(0.1*ddm1[1]):trunc(0.9*ddm1[1])+1
    indx3 <- trunc(0.15*ddm1[1]):trunc(0.85*ddm1[1])+1
    indy1 <- trunc(0.05*ddm1[2]):trunc(0.95*ddm1[2])+1
    indy2 <- trunc(0.1*ddm1[2]):trunc(0.9*ddm1[2])+1
    indy3 <- trunc(0.15*ddm1[2]):trunc(0.85*ddm1[2])+1
    indz1 <- trunc(0.05*ddm1[3]):trunc(0.95*ddm1[3])+1
    indz2 <- trunc(0.1*ddm1[3]):trunc(0.9*ddm1[3])+1
    indz3 <- trunc(0.15*ddm1[3]):trunc(0.85*ddm1[3])+1
    z1 <- density(if(ls0ind>1) s0mean[indx1,indy1,indz1][s0mean[indx1,indy1,indz1]>0] else s0[indx1,indy1,indz1][s0[indx1,indy1,indz1]>0],bw=bw,n=1024)
    z2 <- density(if(ls0ind>1) s0mean[indx2,indy2,indz2][s0mean[indx2,indy2,indz2]>0] else s0[indx2,indy2,indz2][s0[indx2,indy2,indz2]>0],bw=bw,n=1024)
    z3 <- density(if(ls0ind>1) s0mean[indx3,indy3,indz3][s0mean[indx3,indy3,indz3]>0] else s0[indx3,indy3,indz3][s0[indx3,indy3,indz3]>0],bw=bw,n=1024)
    n <- prod(ddim)
    n1 <- length(indx1)*length(indy1)*length(indz1)
    n2 <- length(indx2)*length(indy2)*length(indz2)
    n3 <- length(indx3)*length(indy3)*length(indz3)
    ylim <- range(z$y,z1$y*n1/n,z2$y*n2/n,z3$y*n3/n)
    while(!accept){
      plot(z,type="l",main="Density of S0 values and cut off point",ylim=ylim)
      lines(z1$x,z1$y*n1/n,col=2)
      lines(z2$x,z2$y*n2/n,col=3)
      lines(z3$x,z3$y*n3/n,col=4)
      lines(c(A0,A0),c(0,max(z$y)/2),col=2,lwd=2)
      legend(min(A0,0.25*max(z$x)),ylim[2],c("Full cube",paste("Central",(n1*100)%/%n,"%"),
                                             paste("Central",(n2*100)%/%n,"%"),paste("Central",(n3*100)%/%n,"%")),col=1:4,lwd=rep(1,4))
      cat("A good cut off point should be left of support of the density of grayvalues within the head\n")
      show.image(make.image(img/maximg))
      title("Central slice: Intensity values")
      show.image(make.image((img<A0)))
      title("Central slice: voxel not in mask")
      a <- readline(paste("Accept current cut off point",A0," (Y/N):"))
      if (toupper(a) == "N") {
        cutpoint <-  readline("Provide value for cut off point:")
        cutpoint <- if(!is.null(cutpoint)) as.numeric(cutpoint) else A0
        if(!is.na(cutpoint)) {
          level0 <- A0 <- cutpoint
        }
      } else {
        accept <- TRUE
      }
    }
    par( oldpar)
  } else {
    if(is.null(level)){
      ddim <- object@ddim
      indx1 <- trunc(0.4*ddim[1]):trunc(0.6*ddim[1])
      indy1 <- trunc(0.4*ddim[2]):trunc(0.6*ddim[2])
      indz1 <- trunc(0.7*ddim[3]):trunc(0.7*ddim[3])
      A0a <- quantile(if(ls0ind>1) s0mean[indx1,indy1,indz1][s0mean[indx1,indy1,indz1]>1] else s0[indx1,indy1,indz1][s0[indx1,indy1,indz1]>1],.01)/(1+1/length(object@s0ind))
      #  A0a provides a guess for a threshold based on lower quantiles of intensities
      #  in a central cube (probably contained within the head)
      #  the last factor adjusts for increased accuracy with replicated s0-values
      indx1 <- c(1:trunc(0.15*ddim[1]),trunc(0.85*ddim[1]):ddim[1])
      indy1 <- c(1:trunc(0.15*ddim[2]),trunc(0.85*ddim[2]):ddim[2])
      indz1 <- c(1:trunc(0.15*ddim[3]),trunc(0.85*ddim[3]):ddim[3])
      A0b <- quantile(if(ls0ind>1) s0mean[indx1,indy1,indz1] else s0[indx1,indy1,indz1],.99)
      #  A0a provides a guess for a threshold based on upper quantiles of intensities
      #  in cubes located at the edges (probably only containing noise
      level0 <- A0 <- min(A0a,A0b)*threshfactor
    } 
  }
  A0 <- max(level0,A1/200,1)
  # avoid A0=0 since this may lead to Inf weights in dtiTensor
  # determine parameters for linear relation between standard deviation and mean
  if(sdmethod=="none"){
     sdcoef0 <- c(1,0) ## use OLSE
  } else {
    if(ls0ind>1) {
       s0sd <- apply(s0,1,sdmethod)
       ind <- s0mean>A0&s0mean<A1
       if(length(ind)<2){
         warning("you need more than one voxel to model variances choice of A0/A1 to restrictive")
         return(object)
       }
       sdcoef0 <- coefficients(lm(s0sd[ind]~s0mean[ind]))
       if(sdcoef0[1]<0){
          sdcoef0 <- numeric(2)
          sdcoef0[1] <- .25  # this is an arbitrary (small) value to avaoid zero variances
          sdcoef0[2] <- coefficients(lm(s0sd[ind]~s0mean[ind]-1))
       }
       if(sdcoef0[2]<0){
          sdcoef0 <- numeric(2)
          sdcoef0[1] <- max(0.25,mean(s0sd[ind]))
          sdcoef0[2] <- 0
       }
     } else {
       sdcoef0 <- awslinsd(s0,hmax=5,mask=NULL,A0=A0,A1=A1)$vcoef
     }
     cat("Estimated parameters:",signif(sdcoef0[1:2],3),"Interval of linearity",signif(A0,3),"-",signif(A1,3),"\n")
  }
  object@level <- level0
  object@sdcoef[1:4] <- c(sdcoef0,A0,A1)
  object
})

getsdofsb <- function(object,  ...) cat("No method defined for class:",class(object),"\n")

setGeneric("getsdofsb", function(object,  ...) standardGeneric("getsdofsb"))

setMethod("getsdofsb","dtiData", function(object,qA0=.1,qA1=.98,nsb=NULL,level=NULL){
  # determine interval of linearity
  if(prod(object@ddim)==1){
    warning("you need more than one voxel to model variances")
    return(object)
  }
  if(length(object@sdcoef)==4) object@sdcoef<- c(object@sdcoef,rep(0,4))
  ngrad <- object@ngrad
  if(is.null(level)) level <- object@level
  s0ind<-object@s0ind
  s0 <- object@si[,,,s0ind]
  ls0ind <- length(s0ind)
  A0 <- level
  if(ls0ind>1) {
    dim(s0) <- c(prod(object@ddim),ls0ind)
    s0 <- s0%*%rep(1/ls0ind,ls0ind)
  }
  dim(s0) <- object@ddim
  mask <- s0 > level
  vext <- object@voxelext
  ns0 <- length(s0ind)
  ngrad0 <- ngrad - ns0
  if(is.null(nsb)) nsb <- ngrad0
  set.seed(1)
  snsb <- sample(ngrad0,nsb)
  sb <- extract(object,what="sb")$sb[,,,snsb,drop=FALSE]
  A0 <- quantile(sb,qA0)
  A1 <- quantile(sb,qA1)
  sdcoef1 <- coef1 <- coef2 <- numeric(nsb)
  for(i in 1:nsb) {     
    z <- awslinsd(sb[,,,i],hmax=5,mask=mask,A0=A0,A1=A1)$vcoef
    cat("standard deviation parameters trial i",z,"\n")
    coef1[i] <- z[1]
    coef2[i] <- z[2]
    sdcoef1[i] <- z[1]+z[2]*mean(sb[,,,i][mask])
  }
  # determine parameters for linear relation between standard deviation and mean
  #  object@sdcoef[5:8] <- c(median(sdcoef1),A0,A1,nsb)
  object@sdcoef[5:8] <- c(median(coef1),median(coef2),A0,A1)
  cat("Estimated parameters:",object@sdcoef[5:6],signif(mean(sdcoef1),3),"Interval used",signif(A0,3),"-",signif(A1,3),"\n")
  object
})
############### [

setMethod("[","dtiData",function(x, i, j, k, drop=FALSE){
  args <- sys.call(-1)
  args <- c(x@call,args)
  if (missing(i)) i <- TRUE
  if (missing(j)) j <- TRUE
  if (missing(k)) k <- TRUE
  if (is.logical(i)) ddimi <- x@ddim[1] else ddimi <- length(i)
  if (is.logical(j)) ddimj <- x@ddim[2] else ddimj <- length(j)
  if (is.logical(k)) ddimk <- x@ddim[3] else ddimk <- length(k)
  swap <- rep(FALSE,3)
  if (!is.logical(i)) swap[1] <- i[1] > i[length(i)]
  if (!is.logical(j)) swap[2] <- j[1] > j[length(j)]
  if (!is.logical(k)) swap[3] <- k[1] > k[length(k)]
  orientation <- x@orientation
  gradient <- x@gradient
  if(swap[1]) {
    orientation[1] <- (orientation[1]+1)%%2
    gradient[1,] <- -gradient[1,]
  }
  if(swap[2]) {
    orientation[2] <- (orientation[2]+1)%%2+2
    gradient[2,] <- -gradient[2,]
  }
  if(swap[3]) {
    orientation[3] <- (orientation[3]+1)%%2+4
    gradient[3,] <- -gradient[3,]
  }
  invisible(new("dtiData",
                call   = args,
                si     = x@si[i,j,k,,drop=FALSE],
                gradient = gradient,
                bvalue = x@bvalue,
                btb    = x@btb,
                ngrad  = x@ngrad,
                s0ind  = x@s0ind,
                replind = x@replind,
                ddim   = c(ddimi,ddimj,ddimk),
                ddim0  = x@ddim0,
                xind   = x@xind[i],
                yind   = x@yind[j],
                zind   = x@zind[k],
                sdcoef = x@sdcoef,
                level  = x@level,
                voxelext = x@voxelext,
                orientation = as.integer(orientation),
                rotation = x@rotation,
                source = x@source)
  )
})


subsetg <- function(x, ind){
  args <- sys.call(-1)
  args <- c(x@call,args)
  if (missing(ind)) ind <- 1:x@ngrad
  s0ind <- (1:length(ind))[ind%in%x@s0ind]  
  invisible(new("dtiData",
                call   = args,
                si     = x@si[,,,ind,drop=FALSE],
                gradient = x@gradient[,ind],
                bvalue = x@bvalue[ind],
                btb    = x@btb[,ind],
                ngrad  = as.integer(length(ind)),
                s0ind  = s0ind,
                replind = replind(x@gradient[,ind]),
                ddim   = x@ddim,
                ddim0  = x@ddim0,
                xind   = x@xind,
                yind   = x@yind,
                zind   = x@zind,
                sdcoef = x@sdcoef,
                level  = x@level,
                voxelext = x@voxelext,
                orientation = as.integer(x@orientation),
                rotation = x@rotation,
                source = x@source)
  )
}

combineDWIdata <- function(x1, x2, s0strategy="first"){
  # s0 strategies: "first", "second", "both", "rfirst", "rsecond", "rboth"
  args <- sys.call(-1)
  args <- c(x1@call,x2@call,args)
  if(any(x1@ddim!=x2@ddim)) return(warning("Incompatible dimensions"))
  s01 <- x1@si[,,,x1@s0ind]
  s02 <- x2@si[,,,x2@s0ind]
  ls01 <- length(x1@s0ind)
  ls02 <- length(x2@s0ind)
  if(s0strategy%in%c("both","rfirst","rboth")) dim(s01) <- c(prod(x1@ddim),ls01)
  if(s0strategy%in%c("both","rsecond","rboth")) dim(s02) <- c(prod(x2@ddim),ls02)
  grad1 <- x1@gradient[,-x1@s0ind]
  grad2 <- x2@gradient[,-x2@s0ind]
  s0 <- switch(s0strategy,"first" = s01, 
               "second" = s02, 
               "both" = array(c(s01,s02),c(x1@ddim,ls01+ls02)),
               "rfirst" = s01%*%rep(1/ls01,ls01), 
               "rsecond" = s02%*%rep(1/ls02,ls02), 
               "rboth" = cbind(s01,s02)%*%rep(1/(ls01+ls02),ls01+ls02))
  ls0 <- switch(s0strategy,"first" = ls01, 
                "second" = ls02, 
                "both" = ls01+ls02,
                "rfirst" = 1, 
                "rsecond" = 1, 
                "rboth" = 1)
  ns1 <- x1@ngrad-ls01
  ns2 <- x2@ngrad-ls02
  si <- array(c(s0,x1@si[,,,-x1@s0ind],x2@si[,,,-x2@s0ind]),c(x1@ddim,ls0+ns1+ns2))
  gradient <- matrix(c(rep(0,3*ls0),x1@gradient[,-x1@s0ind],
                       x2@gradient[,-x2@s0ind]),c(3,ls0+ns1+ns2))
  bvalue <- c(rep(0,ls0),x1@bvalue[-x1@s0ind],x2@bvalue[-x2@s0ind])
  btb <- matrix(c(rep(0,6*ls0),x1@btb[,-x1@s0ind],
                  x2@btb[,-x2@s0ind]),c(6,ls0+ns1+ns2))
  invisible(new("dtiData",
                call   = args,
                si     = si,
                gradient = gradient,
                bvalue = bvalue,
                btb    = btb,
                ngrad  = as.integer(ls0+ns1+ns2),
                s0ind  = 1:ls0,
                replind = replind(gradient),
                ddim   = x1@ddim,
                ddim0  = as.integer(c(x1@ddim0[1:3],ls0+ns1+ns2)),
                xind   = x1@xind,
                yind   = x1@yind,
                zind   = x1@zind,
                sdcoef = x1@sdcoef,
                level  = x1@level,
                voxelext = x1@voxelext,
                orientation = as.integer(x1@orientation),
                rotation = x1@rotation,
                source = x1@source)
  )
  
}

##############

setMethod("[","dtiTensor",function(x, i, j, k, drop=FALSE){
  args <- sys.call(-1)
  args <- c(x@call,args)
  if (missing(i)) i <- TRUE
  if (missing(j)) j <- TRUE
  if (missing(k)) k <- TRUE
  if (is.logical(i)) ddimi <- x@ddim[1] else ddimi <- length(i)
  if (is.logical(j)) ddimj <- x@ddim[2] else ddimj <- length(j)
  if (is.logical(k)) ddimk <- x@ddim[3] else ddimk <- length(k)
  swap <- rep(FALSE,3)
  if (!is.logical(i)) swap[1] <- i[1] > i[length(i)]
  if (!is.logical(j)) swap[2] <- j[1] > j[length(j)]
  if (!is.logical(k)) swap[3] <- k[1] > k[length(k)]
  orientation <- x@orientation
  gradient <- x@gradient
  btb <- x@btb
  D <- x@D
  if(swap[1]) {
    orientation[1] <- (orientation[1]+1)%%2
    gradient[1,] <- -gradient[1,]
    btb[2:3,] <- - btb[2:3,]
    D[2:3,,,] <- - D[2:3,,,]
  }
  if(swap[2]) {
    orientation[2] <- (orientation[2]+1)%%2+2
    gradient[2,] <- -gradient[2,]
    btb[c(2,5),] <- - btb[c(2,5),]
    D[c(2,5),,,] <- - D[c(2,5),,,]
  }
  if(swap[3]) {
    orientation[3] <- (orientation[3]+1)%%2+4
    gradient[3,] <- -gradient[3,]
    btb[c(3,5),] <- - btb[c(3,5),]
    D[c(3,5),,,] <- - D[c(3,5),,,]
  }
  ind <- 1:prod(x@ddim)
  if(length(x@outlier)>0){
    ind <- rep(FALSE,prod(x@ddim))
    ind[x@outlier] <- TRUE
    dim(ind) <- x@ddim
    ind <- ind[i,j,k]
    outlier <- (1:length(ind))[ind]
  } else {
    outlier <- numeric(0)
  }
  
  invisible(new("dtiTensor",
                call  = args, 
                D     = D[,i,j,k,drop=FALSE],
                th0   = x@th0[i,j,k,drop=FALSE],
                sigma = if(x@method=="linear") x@sigma[i,j,k,drop=FALSE] else array(1,c(1,1,1)),
                scorr = x@scorr, 
                bw = x@bw,
                mask = x@mask[i,j,k,drop=FALSE],
                hmax = x@hmax,
                gradient = gradient,
                bvalue = x@bvalue,
                btb   = btb,
                ngrad = x@ngrad,
                s0ind = x@s0ind,
                replind = x@replind,
                ddim  = c(ddimi,ddimj,ddimk),
                ddim0 = x@ddim0,
                xind  = x@xind[i],
                yind  = x@yind[j],
                zind  = x@zind[k],
                voxelext = x@voxelext,
                level = x@level,
                orientation = as.integer(orientation),
                rotation = x@rotation,
                outlier = outlier,
                scale = x@scale,
                source = x@source,
                method = x@method)
  )
})

setMethod("[", "dkiTensor",
          function(x, i, j, k, drop = FALSE) {
            
            args <- sys.call(-1)
            args <- c(x@call, args)
            
            if (missing(i)) i <- TRUE
            if (missing(j)) j <- TRUE
            if (missing(k)) k <- TRUE
            
            if (is.logical(i)) ddimi <- x@ddim[1] else ddimi <- length(i)
            if (is.logical(j)) ddimj <- x@ddim[2] else ddimj <- length(j)
            if (is.logical(k)) ddimk <- x@ddim[3] else ddimk <- length(k)
            
            swap <- rep(FALSE, 3)
            if (!is.logical(i)) swap[1] <- (i[ 1] > i[length(i)])
            if (!is.logical(j)) swap[2] <- (j[ 1] > j[length(j)])
            if (!is.logical(k)) swap[3] <- (k[ 1] > k[length(k)])
            orientation <- x@orientation
            gradient <- x@gradient
            btb <- x@btb
            D <- x@D[, i, j, k, drop = FALSE]
            W <- x@W[, i, j, k, drop = FALSE]
            if( swap[1]) {
              orientation[1] <- (orientation[1] + 1) %% 2
              gradient[1, ] <- - gradient[1, ]
              btb[2:3, ] <- - btb[2:3, ]
              D[2:3, , , ] <- - D[2:3, , , ]
              warning("[: kurtosis tensor not correctly transformed for reverse order!")
              ## TODO: determine elements to change!
              ## W[ c(?), , , ] <- - W[ c(?), , , ] 
            }
            if(swap[2]) {
              orientation[2] <- (orientation[2] + 1) %% 2 + 2
              gradient[2, ] <- - gradient[2, ]
              btb[c(2, 5), ] <- - btb[c(2, 5), ]
              D[c(2, 5), , , ] <- - D[c(2, 5), , , ]
              warning("[: kurtosis tensor not correctly transformed for reverse order!")
              ## TODO: determine elements to change!
              ## W[ c(?), , , ] <- - W[ c(?), , , ] 
            }
            if(swap[3]) {
              orientation[3] <- (orientation[3] + 1) %% 2 + 4
              gradient[3, ] <- -gradient[3, ]
              btb[c(3, 5), ] <- - btb[c(3, 5), ]
              D[c(3, 5), , , ] <- - D[c(3, 5), , , ]
              warning( "[: kurtosis tensor not correctly transformed for reverse order!")
              ## TODO: determine elements to change!
              ## W[ c(?), , , ] <- - W[ c(?), , , ] 
            }
            
            ind <- 1:prod(x@ddim)
            if (length(x@outlier) > 0) {
              ind <- rep(FALSE, prod(x@ddim))
              ind[x@outlier] <- TRUE
              dim(ind) <- x@ddim
              ind <- ind[i, j, k]
              outlier <- (1:length(ind))[ind]
            } else {
              outlier <- numeric(0)
            }
            
            invisible(new("dkiTensor",
                          call        = args, 
                          D           = D,
                          W           = W,
                          th0         = x@th0[i, j, k, drop=FALSE],
                          ## TODO: determination of sigma!
                          sigma       = if(x@method == "linear") x@sigma[i, j, k, drop=FALSE] else array(1, c(1, 1, 1)),
                          scorr       = x@scorr, 
                          bw          = x@bw,
                          mask        = x@mask[ i, j, k, drop=FALSE],
                          hmax        = x@hmax,
                          gradient    = gradient,
                          bvalue      = x@bvalue,
                          btb         = btb,
                          ngrad       = x@ngrad,
                          s0ind       = x@s0ind,
                          replind     = x@replind,
                          ddim        = c(ddimi, ddimj, ddimk),
                          ddim0       = x@ddim0,
                          xind        = x@xind[i],
                          yind        = x@yind[j],
                          zind        = x@zind[k],
                          voxelext    = x@voxelext,
                          level       = x@level,
                          orientation = as.integer(orientation),
                          rotation    = x@rotation,
                          outlier     = outlier,
                          scale       = x@scale,
                          source      = x@source,
                          method      = x@method))
          })

#############

setMethod("[", "dkiIndices",
          function(x, i, j, k, drop=FALSE){
            
            args <- sys.call(-1)
            args <- c(x@call,args)
            
            if (missing(i)) i <- TRUE
            if (missing(j)) j <- TRUE
            if (missing(k)) k <- TRUE
            
            if (is.logical(i)) ddimi <- x@ddim[1] else ddimi <- length(i)
            if (is.logical(j)) ddimj <- x@ddim[2] else ddimj <- length(j)
            if (is.logical(k)) ddimk <- x@ddim[3] else ddimk <- length(k)
            
            swap <- rep(FALSE,3)
            if (!is.logical(i)) swap[1] <- i[1] > i[length(i)]
            if (!is.logical(j)) swap[2] <- j[1] > j[length(j)]
            if (!is.logical(k)) swap[3] <- k[1] > k[length(k)]
            
            orientation <- x@orientation
            gradient <- x@gradient
            btb <- x@btb
            andir <- x@andir
            
            if(swap[1]) {
              orientation[1] <- (orientation[1]+1)%%2
              gradient[1, ] <- -gradient[1, ]
              btb[2:3, ] <- - btb[2:3, ]
              andir[1, , , ] <- - andir[1, , , ]
            }
            if(swap[2]) {
              orientation[2] <- (orientation[2]+1)%%2+2
              gradient[2, ] <- -gradient[2, ]
              btb[c(2, 5), ] <- - btb[c(2, 5), ]
              andir[2, , , ] <- - andir[2, , , ]
            }
            if(swap[3]) {
              orientation[3] <- (orientation[3]+1)%%2+4
              gradient[3, ] <- -gradient[3, ]
              btb[c(3, 5), ] <- - btb[c(3, 5), ]
              andir[3, , , ] <- - andir[3, , , ]
            }
            
            invisible(new("dkiIndices",
                          call = args,
                          fa = x@fa[i, j, k, drop=FALSE],
                          ga = x@ga[i, j, k, drop=FALSE],
                          md = x@md[i, j, k, drop=FALSE],
                          andir = andir[, i, j, k, drop=FALSE],
                          bary = x@bary[, i, j, k, drop=FALSE],
                          k1 = x@k1[i, j, k, drop=FALSE],
                          k2 = x@k2[i, j, k, drop=FALSE],
                          k3 = x@k3[i, j, k, drop=FALSE],
                          mk = x@mk[i, j, k, drop=FALSE],
                          mk2 = x@mk2[i, j, k, drop=FALSE],
                          kaxial = x@kaxial[i, j, k, drop=FALSE],
                          kradial = x@kradial[i, j, k, drop=FALSE],
                          fak = x@fak[i, j, k, drop=FALSE],
                          gradient = gradient,
                          bvalue = x@bvalue,
                          btb   = btb,
                          ngrad = x@ngrad,
                          s0ind = x@s0ind,
                          ddim  = c(ddimi,ddimj,ddimk),
                          ddim0 = x@ddim0,
                          voxelext = x@voxelext,
                          orientation = as.integer(orientation),
                          rotation = x@rotation,
                          xind  = x@xind[i],
                          yind  = x@yind[j],
                          zind  = x@zind[k],
                          method = x@method,
                          level = x@level,
                          source= x@source)
            )
          })

#############
setMethod("[","dwiMixtensor",function(x, i, j, k, drop=FALSE){
  args <- sys.call(-1)
  args <- c(x@call,args)
  if (missing(i)) i <- TRUE
  if (missing(j)) j <- TRUE
  if (missing(k)) k <- TRUE
  if (is.logical(i)) ddimi <- x@ddim[1] else ddimi <- length(i)
  if (is.logical(j)) ddimj <- x@ddim[2] else ddimj <- length(j)
  if (is.logical(k)) ddimk <- x@ddim[3] else ddimk <- length(k)
  
  ind <- 1:prod(x@ddim)
  if(length(x@outlier)>0){
    ind <- rep(FALSE,prod(x@ddim))
    ind[x@outlier] <- TRUE
    dim(ind) <- x@ddim
    ind <- ind[i,j,k]
    outlier <- (1:length(ind))[ind]
  } else {
    outlier <- numeric(0)
  }
  #  cat("indix i",i,"\n")
  #  cat("indix j",j,"\n")
  #  cat("indix k",k,"\n")
  swap <- rep(FALSE,3)
  if (!is.logical(i)) swap[1] <- i[1] > i[length(i)]
  if (!is.logical(j)) swap[2] <- j[1] > j[length(j)]
  if (!is.logical(k)) swap[3] <- k[1] > k[length(k)]
  #  cat("swap",swap,"\n")
  orientation <- x@orientation
  gradient <- x@gradient
  btb <- x@btb
  orient <- x@orient
  if(swap[1]) {
    orientation[1] <- (orientation[1]+1)%%2
    gradient[1,] <- -gradient[1,]
    btb[2:3,] <- - btb[2:3,]
    orient <- - orient
  }
  if(swap[2]) {
    orientation[2] <- (orientation[2]+1)%%2+2
    gradient[2,] <- -gradient[2,]
    btb[c(2,5),] <- - btb[c(2,5),]
    orient[2,,,,] <- - orient[2,,,,]
  }
  if(swap[3]) {
    orientation[3] <- (orientation[3]+1)%%2+4
    gradient[3,] <- -gradient[3,]
    btb[c(3,5),] <- - btb[c(3,5),]
    orient[1,,,,] <- pi - orient[1,,,,]
  }
  #  cat("new orientation",orientation,"\n")
  #  cat("indix i",i,"\n")
  #  cat("indix j",j,"\n")
  #  cat("indix k",k,"\n")
  invisible(new("dwiMixtensor",
                call  = args, 
                ev     = x@ev[,i,j,k,drop=FALSE],
                mix    = x@mix[,i,j,k,drop=FALSE],
                orient = orient[,,i,j,k,drop=FALSE],
                order  = x@order[i,j,k,drop=FALSE],
                p      = x@p,
                th0   = x@th0[i,j,k,drop=FALSE],
                sigma = x@sigma[i,j,k,drop=FALSE],
                scorr = x@scorr, 
                bw = x@bw,
                mask = x@mask[i,j,k,drop=FALSE],
                hmax = x@hmax,
                gradient = gradient,
                bvalue = x@bvalue,
                btb   = btb,
                ngrad = x@ngrad,
                s0ind = x@s0ind,
                replind = x@replind,
                ddim  = c(ddimi,ddimj,ddimk),
                ddim0 = x@ddim0,
                xind  = x@xind[i],
                yind  = x@yind[j],
                zind  = x@zind[k],
                voxelext = x@voxelext,
                level = x@level,
                orientation = as.integer(orientation),
                rotation = x@rotation,
                outlier = outlier,
                scale = x@scale,
                source = x@source,
                method = x@method)
  )
})

#############

setMethod("[","dtiIndices",function(x, i, j, k, drop=FALSE){
  args <- sys.call(-1)
  args <- c(x@call,args)
  if (missing(i)) i <- TRUE
  if (missing(j)) j <- TRUE
  if (missing(k)) k <- TRUE
  if (is.logical(i)) ddimi <- x@ddim[1] else ddimi <- length(i)
  if (is.logical(j)) ddimj <- x@ddim[2] else ddimj <- length(j)
  if (is.logical(k)) ddimk <- x@ddim[3] else ddimk <- length(k)
  swap <- rep(FALSE,3)
  if (!is.logical(i)) swap[1] <- i[1] > i[length(i)]
  if (!is.logical(j)) swap[2] <- j[1] > j[length(j)]
  if (!is.logical(k)) swap[3] <- k[1] > k[length(k)]
  orientation <- x@orientation
  gradient <- x@gradient
  btb <- x@btb
  andir <- x@andir
  if(swap[1]) {
    orientation[1] <- (orientation[1]+1)%%2
    gradient[1,] <- -gradient[1,]
    btb[2:3,] <- - btb[2:3,]
    andir[1,,,] <- - andir[1,,,]
  }
  if(swap[2]) {
    orientation[2] <- (orientation[2]+1)%%2+2
    gradient[2,] <- -gradient[2,]
    btb[c(2,5),] <- - btb[c(2,5),]
    andir[2,,,] <- - andir[2,,,]
  }
  if(swap[3]) {
    orientation[3] <- (orientation[3]+1)%%2+4
    gradient[3,] <- -gradient[3,]
    btb[c(3,5),] <- - btb[c(3,5),]
    andir[3,,,] <- - andir[3,,,]
  }
  
  invisible(new("dtiIndices",
                call = args,
                fa = x@fa[i,j,k,drop=FALSE],
                ga = x@ga[i,j,k,drop=FALSE],
                md = x@md[i,j,k,drop=FALSE],
                andir = andir[,i,j,k,drop=FALSE],
                bary = x@bary[,i,j,k,drop=FALSE],
                gradient = gradient,
                bvalue = x@bvalue,
                btb   = btb,
                ngrad = x@ngrad,
                s0ind = x@s0ind,
                ddim  = c(ddimi,ddimj,ddimk),
                ddim0 = x@ddim0,
                voxelext = x@voxelext,
                orientation = as.integer(orientation),
                rotation = x@rotation,
                xind  = x@xind[i],
                yind  = x@yind[j],
                zind  = x@zind[k],
                method = x@method,
                level = x@level,
                source= x@source)
  )
})

###########

setMethod("[","dwiQball",function(x, i, j, k, drop=FALSE){
  args <- sys.call(-1)
  args <- c(x@call,args)
  if (missing(i)) i <- TRUE
  if (missing(j)) j <- TRUE
  if (missing(k)) k <- TRUE
  if (is.logical(i)) ddimi <- x@ddim[1] else ddimi <- length(i)
  if (is.logical(j)) ddimj <- x@ddim[2] else ddimj <- length(j)
  if (is.logical(k)) ddimk <- x@ddim[3] else ddimk <- length(k)
  swap <- rep(FALSE,3)
  if (!is.logical(i)) swap[1] <- i[1] > i[length(i)]
  if (!is.logical(j)) swap[2] <- j[1] > j[length(j)]
  if (!is.logical(k)) swap[3] <- k[1] > k[length(k)]
  if(any(swap)) {
    warning("can't reverse order of indices")
    return(invisible(x))
  }
  ind <- 1:prod(x@ddim)
  if(length(x@outlier)>0){
    ind <- rep(FALSE,prod(x@ddim))
    ind[x@outlier] <- TRUE
    dim(ind) <- x@ddim
    ind <- ind[i,j,k]
    outlier <- (1:length(ind))[ind]
  } else {
    outlier <- numeric(0)
  }
  
  invisible(new("dwiQball",
                call  = args, 
                order = x@order,
                forder = x@forder,
                lambda = x@lambda,
                sphcoef = if(x@what%in%c("sqrtODF")) x@sphcoef[,,i,j,k,drop=FALSE] else  x@sphcoef[,i,j,k,drop=FALSE],
                varsphcoef = x@varsphcoef,
                th0   = x@th0[i,j,k,drop=FALSE],
                sigma = x@sigma[i,j,k,drop=FALSE],
                scorr = x@scorr, 
                bw = x@bw,
                mask = x@mask[i,j,k,drop=FALSE],
                hmax = x@hmax,
                gradient = x@gradient,
                btb   = x@btb,
                bvalue = x@bvalue,
                ngrad = x@ngrad,
                s0ind = x@s0ind,
                replind = x@replind,
                ddim  = c(ddimi,ddimj,ddimk),
                ddim0 = x@ddim0,
                xind  = x@xind[i],
                yind  = x@yind[j],
                zind  = x@zind[k],
                voxelext = x@voxelext,
                level = x@level,
                orientation = x@orientation,
                rotation = x@rotation,
                outlier = outlier,
                scale = x@scale,
                source = x@source,
                what = x@what)
  )
})


########## extract()

#extract <- function(x, ...) cat("Data extraction not defined for this class:",class(x),"\n")

#setGeneric("extract", function(x, ...) standardGeneric("extract"))
#setGeneric("extract", function(x, ...) cat("Data extraction not defined for this class:",class(x),"\n"), package="dti")

setMethod("extract","dtiData",function(x,
                                       what=c("data","gradient","btb","s0","sb","siq"), 
                                       xind=TRUE, yind=TRUE, zind=TRUE){
  what <- tolower(what) 
  ## check what
  what <- match.arg(what, several.ok = TRUE)
  swap <- rep(FALSE,3)
  if(is.numeric(xind)) swap[1] <- xind[1]>xind[length(xind)]
  if(is.numeric(yind)) swap[2] <- yind[1]>yind[length(yind)]
  if(is.numeric(zind)) swap[3] <- zind[1]>zind[length(zind)]
  if(any(swap)) {
    warning("can't reverse order of indices ")
    return(NULL)
  }
  x <- x[xind,yind,zind]
  z <- list(NULL)
  if("gradient" %in% what) z$gradient <- x@gradient
  if("btb" %in% what) z$btb <- x@btb
  if("s0" %in% what) z$s0 <- x@si[,,,x@s0ind,drop=FALSE]
  if("sb" %in% what) z$sb <- x@si[,,,-x@s0ind,drop=FALSE]
  if("siq" %in% what) {
    S0 <- x@si[,,,x@s0ind,drop=FALSE]
    Si <- x@si[,,,-x@s0ind,drop=FALSE]
    z$siq <- sweep(Si,1:3,apply(S0,1:3,mean),"/")
    z$siq[is.na(z$siq)] <- 0
  }
  if("data" %in% what) z$data <- x@si
  invisible(z)
})

#############

setMethod("extract","dwiMixtensor",function(x, 
                                            what=c("andir","order","ev","mix","s0","mask","fa","eorder","bic","aic"), 
                                            xind=TRUE, yind=TRUE, zind=TRUE){
  what <- tolower(what) 
  ## check what
  what <- match.arg(what, several.ok = TRUE)
  swap <- rep(FALSE,3)
  if(is.numeric(xind)) swap[1] <- xind[1]>xind[length(xind)]
  if(is.numeric(yind)) swap[2] <- yind[1]>yind[length(yind)]
  if(is.numeric(zind)) swap[3] <- zind[1]>zind[length(zind)]
  if(any(swap)){
    warning("can't reverse order of indices ")
    return(NULL)
  }
  x <- x[xind,yind,zind]
  n1 <- x@ddim[1]
  n2 <- x@ddim[2]
  n3 <- x@ddim[3]
  z <- list(NULL)
  if("order" %in% what) z$order <- x@order
  if("ev" %in% what) { 
    ev <- array(0,c(3,dim(x@ev)[-1]))
    ev[1:2,,,] <- x@ev
    ev[3,,,] <- x@ev[2,,,]
    z$ev <- ev
  }
  if("mix" %in% what) z$mix <- x@mix
  if("andir" %in% what) {
    orient <- x@orient
    andir <- array(0,c(3,prod(dim(orient))/2))
    dim(orient) <- c(2,prod(dim(orient))/2)
    sth <- sin(orient[1,])
    andir[1,] <- sth*cos(orient[2,])
    andir[2,] <- sth*sin(orient[2,])
    andir[3,] <- cos(orient[1,])
    z$andir <- array(andir,c(3,dim(x@orient)[-1]))
  }
  if("s0" %in% what) z$s0 <- x@th0
  if("mask" %in% what) z$mask <- x@mask
  if("fa" %in% what){
    alpha <- (x@ev[1,,,]-x@ev[2,,,])/x@ev[2,,,]
    fa <- alpha/sqrt(3+2*alpha+alpha^2)
    fa[x@order==0] <- 0
    dim(fa) <- x@ddim
    z$fa <- fa
  }
  if("eorder" %in% what) {
    maxorder <- dim(x@mix)[1]
    mix <- x@mix
    dim(mix) <- c(maxorder,n1*n2*n3)     
    smix <- rep(1,maxorder)%*%mix
    mix <- sweep(mix,2,smix,"/")
    # the last two lines are needed for models with isotropic compartment
    z$eorder <- array((2*(1:maxorder)-1)%*%mix,x@ddim)
  }
  if("bic" %in% what) {
    ngrad <- x@ngrad      
    ns0 <- length(x@s0ind)
    iso <- apply(x@mix,-1,sum)
    iso <- iso>0&&iso<1e0-1e-8
    penBIC <- log(ngrad-ns0)/(ngrad-ns0)*(iso+1+2*x@order)
    z$bic <- array(log(pmax(1e-10,x@sigma))+penBIC,dim(x@sigma))
  }
  if("aic" %in% what) {
    ngrad <- x@ngrad      
    ns0 <- length(x@s0ind)
    iso <- apply(x@mix,-1,sum)
    iso <- iso>0&&iso<1e0-1e-8
    penAIC <- 2/(ngrad-ns0)*(iso+1+2*x@order)
    z$aic <- array(log(pmax(1e-10,x@sigma))+penAIC,dim(x@sigma))
  }
  invisible(z)
})



setMethod("extract","dtiTensor",function(x, 
                                         what=c("tensor","fa","ga","md","evalues","andir","s0","mask","bic","aic","outlier"), 
                                         xind=TRUE, yind=TRUE, zind=TRUE,mc.cores=setCores(,reprt=FALSE)){
  what <- tolower(what) 
  ## check what
  what <- match.arg(what, several.ok = TRUE)
  swap <- rep(FALSE,3)
  if(is.numeric(xind)) swap[1] <- xind[1]>xind[length(xind)]
  if(is.numeric(yind)) swap[2] <- yind[1]>yind[length(yind)]
  if(is.numeric(zind)) swap[3] <- zind[1]>zind[length(zind)]
  if(any(swap)) {
    warning("can't reverse order of indices ")
    return(NULL)
  }
  
  x <- x[xind,yind,zind]
  n1 <- x@ddim[1]
  n2 <- x@ddim[2]
  n3 <- x@ddim[3]
  ddim <- x@ddim
  nvox <- prod(ddim)
  needev <- ("fa" %in% what) || ("ga" %in% what) || ("md" %in% what) || ("evalues" %in% what)
  needall <- needev && ("andir" %in% what)
  
  z <- list(NULL)
  if(needall){
    erg <- dti3Dall(x@D,x@mask,mc.cores=mc.cores)
    if("fa" %in% what) z$fa <- array(erg$fa,x@ddim)
    if("ga" %in% what) z$ga <- array(erg$ga,x@ddim)
    if("md" %in% what) z$md <- array(erg$md,x@ddim)
    if("evalues" %in% what) z$evalues <- array(erg$ev,c(3,n1,n2,n3))
    if("andir" %in% what) z$andir <- array(erg$andir,c(3,n1,n2,n3))
  } else {
    if(needev){
      ev <- array(dti3Dev(x@D,x@mask,mc.cores=mc.cores),c(3,n1,n2,n3))
      ev[ev<1e-12] <- 1e-12 
      if("fa" %in% what) {
        dd <- apply(ev^2,2:4,sum)
        md <- (ev[1,,,]+ev[2,,,]+ev[3,,,])/3
        sev <- sweep(ev,2:4,md)
        z$fa <- array(sqrt(1.5*apply(sev^2,2:4,sum)/dd),x@ddim)
      }
      if("ga" %in% what) {
        sev <- log(ev)
        md <- (sev[1,,,]+sev[2,,,]+sev[3,,,])/3
        sev <- sweep(sev,2:4,md)
        ga <- sqrt(apply(sev^2,2:4,sum))
        ga[is.na(ga)] <- 0
        z$ga <- array(ga,x@ddim)
      }
      if("md" %in% what) z$md <- array((ev[1,,,]+ev[2,,,]+ev[3,,,])/3,x@ddim)
      if("evalues" %in% what) z$evalues <- array(ev,c(3,x@ddim))
    }
    if("andir" %in% what){
      ev <- array(dti3Dand(x@D,x@mask,mc.cores=mc.cores),c(3,nvox))
    }
  }
  if("tensor" %in% what) z$tensor <- array(x@D,c(6,x@ddim))
  if("s0" %in% what) z$s0 <- array(x@th0,x@ddim)
  if("mask" %in% what) z$mask <- x@mask
  if("bic" %in% what) {
    ngrad <- x@ngrad      
    ns0 <- length(x@s0ind)
    penBIC <- log(ngrad-ns0)/(ngrad-ns0)*6
    z$bic <- array(log(pmax(1e-10,x@sigma))+penBIC,dim(x@sigma))
  }
  if("aic" %in% what) {
    ngrad <- x@ngrad      
    ns0 <- length(x@s0ind)
    penAIC <- 12/(ngrad-ns0)
    z$aic <- array(log(pmax(1e-10,x@sigma))+penAIC,dim(x@sigma))
  }
  if("outlier" %in% what) {
    ind <- 1:prod(x@ddim)
    ind <- rep(FALSE,prod(x@ddim))
    if(length(x@outlier)>0) ind[x@outlier] <- TRUE
    dim(ind) <- x@ddim
  }
  invisible(z)
})

##############

setMethod("extract","dtiIndices",function(x, 
                                          what=c("fa","andir","ga","md","bary"), xind=TRUE, yind=TRUE, zind=TRUE){
  what <- tolower(what)
  ## check what
  what <- match.arg(what, several.ok = TRUE)
  swap <- rep(FALSE,3)
  if(is.numeric(xind)) swap[1] <- xind[1]>xind[length(xind)]
  if(is.numeric(yind)) swap[2] <- yind[1]>yind[length(yind)]
  if(is.numeric(zind)) swap[3] <- zind[1]>zind[length(zind)]
  if(any(swap)) {
    warning("can't reverse order of indices ")
    return(NULL)
  }
  
  x <- x[xind,yind,zind]
  n1 <- x@ddim[1]
  n2 <- x@ddim[2]
  n3 <- x@ddim[3]
  
  z <- list(NULL)
  if("fa" %in% what) z$fa <- x@fa
  if("ga" %in% what) z$ga <- x@ga
  if("md" %in% what) z$md <- x@md
  if("andir" %in% what) z$andir <- x@andir
  if("bary" %in% what) z$bary <- x@bary
  invisible(z)
})

##############

setMethod("extract","dwiQball",function(x,     
                                        what=c("sphcoef","s0","mask","bic","aic","outlier"),
                                        xind=TRUE, yind=TRUE, zind=TRUE){
  what <- tolower(what) 
  ## check what
  what <- match.arg(what, several.ok = TRUE)
  swap <- rep(FALSE,3)
  if(is.numeric(xind)) swap[1] <- xind[1]>xind[length(xind)]
  if(is.numeric(yind)) swap[2] <- yind[1]>yind[length(yind)]
  if(is.numeric(zind)) swap[3] <- zind[1]>zind[length(zind)]
  if(any(swap)) {
    warning("can't reverse order of indices ")
    return(NULL)
  }
  
  x <- x[xind,yind,zind]
  n1 <- x@ddim[1]
  n2 <- x@ddim[2]
  n3 <- x@ddim[3]
  
  z <- list(NULL)
  if("sphcoef" %in% what) z$sphcoef <- x@sphcoef
  if("s0" %in% what) z$s0 <- x@th0
  if("mask" %in% what) z$mask <- x@mask
  if("bic" %in% what) {
    ngrad <- x@ngrad      
    ns0 <- length(x@s0ind)
    ord <- x@order
    penBIC <- log(ngrad-ns0)/(ngrad-ns0)*(ord+1)*(ord+2)/2
    z$bic <- array(log(pmax(1e-10,x@sigma))+penBIC,dim(x@sigma))
  }
  if("aic" %in% what) {
    ngrad <- x@ngrad      
    ns0 <- length(x@s0ind)
    ord <- x@order
    penAIC <- (ord+1)*(ord+2)/(ngrad-ns0)
    z$aic <- array(log(pmax(1e-10,x@sigma))+penAIC,dim(x@sigma))
  }
  if("outlier" %in% what) {
    ind <- 1:prod(x@ddim)
    ind <- rep(FALSE,prod(x@ddim))
    if(length(x@outlier)>0) ind[x@outlier] <- TRUE
    dim(ind) <- x@ddim
  }
  invisible(z)
})

getmask <- function(object,  ...) cat("No method defined for class:",class(object),"\n")

setGeneric("getmask", function(object,  ...) standardGeneric("getmask"))

setMethod("getmask","dtiData",function(object, level=NULL, prop=.4, size=3){
  if(is.null(level)) level <- object@level
  if(!is.null(level)){ 
    z <- .Fortran("getmask",
                  as.double(object@si[,,,object@s0ind]),
                  as.integer(object@ddim[1]),
                  as.integer(object@ddim[2]),
                  as.integer(object@ddim[3]),
                  as.integer(length(object@s0ind)),
                  as.double(level),
                  as.integer(size),
                  as.double(prop),
                  s0=double(prod(object@ddim)),
                  mask=logical(prod(object@ddim)),
                  PACKAGE="dti")[c("s0","mask")]
  } else {
    z <- list(s0=object@si[,,,object@s0ind],mask=array(TRUE,object@ddim))
  }
  dim(z$s0) <- dim(z$mask) <- object@ddim
  z
}
)
setMethod("getmask","array",function(object, level=NULL, prop=.4, size=3){
  if(length(dim(object))!=3){
    level <- NULL
    warning("array dimension need to be of length 3, returning trivial mask")
  } else {
    ddim <- dim(object)
  }
  if(!is.null(level)){ 
    z <- .Fortran("getmask",
                  as.double(object),
                  as.integer(ddim[1]),
                  as.integer(ddim[2]),
                  as.integer(ddim[3]),
                  as.integer(1),
                  as.double(level),
                  as.integer(size),
                  as.double(prop),
                  s0=double(prod(ddim)),
                  mask=logical(prod(ddim)),
                  PACKAGE="dti")[c("s0","mask")]
  } else {
    z <- list(s0=object,mask=array(TRUE,dim(object)))
  }
  dim(z$s0) <- dim(z$mask) <- dim(object)
  z
}
)
selectCube <- function(xind,yind,zind,ddim,maxobj){ 
  if(is.null(xind)) xind <- 1:ddim[1]
  if(is.null(yind)) yind <- 1:ddim[2]
  if(is.null(zind)) zind <- 1:ddim[3]
  n1 <- length(xind)
  n2 <- length(yind)
  n3 <- length(zind)
  if(any(c(xind[1],yind[1],zind[1]) < 1 || any(c(xind[n1]-ddim[1],yind[n2]-ddim[2],zind[n3]-ddim[3])>0))){
    stop("Error in index specification, specified cube exceeds dimensions of object")
  }
  l1 <- (xind[n1]-xind[1]+1-n1) == 0
  l2 <- (yind[n2]-yind[1]+1-n2) == 0
  l3 <- (zind[n3]-zind[1]+1-n3) == 0
  if(!(l1&&l1&&l3)){
    stop("Error in index specification, xind, yind and zind need to be increasing
         and consecutive")
  }
  n <- n1*n2*n3
  if(n>maxobj) {
    cat("size of data cube",n," exceeds maximum of",maxobj,"\n")
    mod2 <- function(x) x==x%/%2*2
    nz <- max(1,maxobj%/%(n1*n2))
    nx <- min(n1,floor(sqrt(maxobj/nz)))
    ny <- min(n2,maxobj%/%(nx*nz))
    xm <- (xind[1]+xind[n1])%/%2
    xd <- nx%/%2
    xind <- if(mod2(nx)) (xm-xd+1):(xm+xd) else (xm-xd):(xm+xd)
    ym <- (yind[1]+yind[n2])%/%2
    yd <- ny%/%2
    yind <- if(mod2(ny)) (ym-yd+1):(ym+yd) else (ym-yd):(ym+yd) 
    zm <- (zind[1]+zind[n3])%/%2
    zd <- nz%/%2
    zind <- if(mod2(nz)) (zm-zd+1):(zm+zd) else (zm-zd):(zm+zd)
    n1 <- length(xind)
    n2 <- length(yind)
    n3 <- length(zind)
    n <- n1*n2*n3
  }
  if(n==0) stop("Empty cube specified")
  list(xind=xind,yind=yind,zind=zind, n=n, n1=n1, n2=n2, n3=n3)
}

