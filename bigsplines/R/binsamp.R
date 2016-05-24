binsamp <-
  function(x,xrng=NULL,nmbin=11,nsamp=1,alg=c("new","old")){
    # bin-sampled knots
    
    x <- as.matrix(x)
    xdim <- dim(x)
    if(is.null(xrng)){
      if(xdim[2]>1){xrng <- apply(x,2,range)} else{xrng <- matrix(range(x),2,1)}
    }
    mysamp <- function(z){ if(length(z)==1L){z} else {sample(z,size=min(nsamp,length(z)))} }
    nmbin <- as.integer(nmbin)
    if(length(nmbin)!=xdim[2]){nmbin <- rep(nmbin[1],xdim[2])}
    if(any(nmbin<2L)){stop("Must set input 'nmbin' >= 2 for all predictors.")}
    gvec <- matrix(1,xdim[1],1)
    kconst <- 1
    if(alg[1]=="old"){
      for(kk in 1:xdim[2]){
        gvec  <-  gvec + kconst*round((nmbin[kk]-1L)*((x[,kk]-xrng[1,kk])/(xrng[2,kk]-xrng[1,kk])))
        kconst  <-  kconst*nmbin[kk]
      }
    } else {
      for(kk in 1:xdim[2]){
        gvec  <-  gvec + kconst*pmin(round((nmbin[kk]-1)*((x[,kk]-xrng[1,kk])/(xrng[2,kk]-xrng[1,kk]))),nmbin[kk]-1L)
        kconst  <-  kconst*nmbin[kk]
      }
    }
    gvec <- as.factor(gvec)
    return(unlist(tapply(1:xdim[1],gvec,mysamp)))
    
  }
