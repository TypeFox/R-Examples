KMC <- function(cl, con, rates, dx, dy, dz, type = "BKL", nSimu = 10, nHops = 1e+07, seed = NULL){
  if(missing(cl))
    stop("Please provide a 'cluster' object (cl). See 'makeCluster'")
  if(missing(con))
    stop("Please specify the connectivity of the percolation network (con)")
  if(missing(rates))
    stop("Please specify the charge transfer rates (rates)")
  if(missing(dx))
    stop("Please specify the inter-site distances (dx)")
  if(missing(dy))
    stop("Please specify the inter-site distances (dy)")
  if(missing(dz))
    stop("Please specify the inter-site distances (dz)")
  
  if(is.list(rates)){
    warning("'rates' is a list. Only the first element has been used")
    rates <- rates[[1]]
  }
  if(is.null(dim(rates)) ) rates <- t(as.matrix(rates))
  if(is.data.frame(rates)) rates <-   as.matrix(rates)
  if(!is.matrix(rates)) stop("'rates' must be a matrix or a data.frame")
  
  if(is.list(dx)){
    warning("'dx' is a list. Only the first element has been used")
    dx <- dx[[1]]
  }
  if(is.null(dim(dx)) ) dx <- t(as.matrix(dx))
  if(is.data.frame(dx)) dx <-   as.matrix(dx)
  if(!is.matrix(dx)) stop("'dx' must be a matrix or a data.frame")
  
  if(is.list(dy)){
    warning("'dy' is a list. Only the first element has been used")
    dy <- dy[[1]]
  }
  if(is.null(dim(dy)) ) dy <- t(as.matrix(dy))
  if(is.data.frame(dy)) dy <-   as.matrix(dy)
  if(!is.matrix(dy)) stop("'dy' must be a matrix or a data.frame")
  
  if(is.list(dz)){
    warning("'dz' is a list. Only the first element has been used")
    dz <- dz[[1]]
  }
  if(is.null(dim(dz)) ) dz <- t(as.matrix(dz))
  if(is.data.frame(dz)) dz <-   as.matrix(dz)
  if(!is.matrix(dz)) stop("'dz' must be a matrix or a data.frame")
  
  dims <- list(dim(rates),dim(dx),dim(dy),dim(dz))
  if(any(!sapply(dims,identical,dims[[1]])))
    stop("'rates', 'dx', 'dy' and 'dz' must have the same dimensions")
  
  if(ncol(con) != 2)
    stop("'con' must be a two-column matrix or data.frame")
  if(nrow(con) != ncol(rates))
    stop("'con' and rates mismatch: nrow(con) != ncol(rates)")
  if(is.data.frame(con))
    con <- as.matrix(con)
  
  if(!(type%in%c("BKL","FRM")))
    stop("Unrecognized type of simulation. 'type' must be equal to 'BKL' or 'FRM'")
  
  if(is.null(seed)){
    useless <- runif(1)
    seed <- sample(.Random.seed[-(1:2)],1)
  }
  
  set.seed(seed)
  seed <- matrix(as.integer(floor(runif(12*nSimu,-1E7,1E7))), ncol=12, nrow=nSimu)
  simu.names <- paste("Simu",gsub(" ","0",format(1:nSimu)),sep="")
  seed <- split(seed, simu.names)
  
  nFrames <- nrow(rates)
  nDimers <- ncol(rates)
  
  frame.names <- rownames(rates)
  dimer.names <- colnames(rates)
  
  site.levels <- as.factor(c(con[,1],con[,2]))
  levels(site.levels) <- 1:nlevels(site.levels)
  con <- t(matrix(as.integer(site.levels), ncol=2))
  
  inputs <- list(
    nHops=as.numeric(nHops),
    nFrames=as.integer(nFrames),
    nDimers=as.integer(nDimers),
    rates=rates,
    dx=dx,
    dy=dy,
    dz=dz,
    con=con
  )

  KMCSingle <- function(seed, inputs, type){
    outputs <- list(
      distx = array(as.numeric(0.0E0),dim=inputs$nFrames),
      disty = array(as.numeric(0.0E0),dim=inputs$nFrames),
      distz = array(as.numeric(0.0E0),dim=inputs$nFrames),
      time  = array(as.numeric(0.0E0),dim=inputs$nFrames),
      nhop  = array(as.numeric(0.0E0),dim=c(inputs$nFrames,inputs$nDimers))
      )
    fortran.args <- c(list(seed), inputs, outputs)
    if(type=="FRM"){
      OUT <- do.call(.Fortran,c("frm",fortran.args))
      OUT <- OUT[c("distx","disty","distz","time","nhop")]
    }
    else if(type=="BKL"){
      OUT <- do.call(.Fortran,c("bkl",fortran.args))
      OUT <- OUT[c("distx","disty","distz","time","nhop")]
    }
    return(OUT)
  }

  clusterCall(cl, library.dynam, chname="ChargeTransport", package="ChargeTransport", lib.loc=.libPaths())
  
#   loc <- system.file(package="ChargeTransport")
#   lib <- file.path(loc, "libs",paste("ChargeTransport",.Platform$dynlib.ext, sep=""))
#   clusterCall(cl, dyn.load, lib)
#   clusterCall(cl, dyn.load, "/home/jide/Dev/R/Packages/ChargeTransport/1.0/src/BKL.so")
#   clusterCall(cl, dyn.load, "/home/jide/Dev/R/Packages/ChargeTransport/1.0/src/FRM.so")

#   OUT <- parLapply(cl, seed, KMCSingle, inputs, type)
  tm <- system.time(
    OUT <- parLapply(cl, 1:length(seed),
    function(n, seed, inputs, type){
      to.return <- KMCSingle(seed[[n]], inputs, type)
      cat("Simulation", format(n, width=nchar(length(seed))),
          "/", length(seed), "done...\n")
      return(to.return)
    }, seed, inputs, type))
  print(tm)
  distx <- sapply(OUT, function(x) return(x$distx))
  disty <- sapply(OUT, function(x) return(x$disty))
  distz <- sapply(OUT, function(x) return(x$distz))
  time  <- sapply(OUT, function(x) return(x$time ))
  nhop  <- sapply(OUT, function(x) return(x$nhop ))
  distx <- do.call(array, list(distx, dim=c(nFrames, nSimu)))
  disty <- do.call(array, list(disty, dim=c(nFrames, nSimu)))
  distz <- do.call(array, list(distz, dim=c(nFrames, nSimu)))
  time  <- do.call(array, list(time , dim=c(nFrames, nSimu)))
  nhop  <- do.call(array, list(nhop , dim=c(nFrames, nDimers, nSimu)))
  
  dimnames(distx) <- list(frame.names, simu.names)
  dimnames(disty) <- list(frame.names, simu.names)
  dimnames(distz) <- list(frame.names, simu.names)
  dimnames(time ) <- list(frame.names, simu.names)
  dimnames(nhop ) <- list(frame.names, dimer.names, simu.names)
  
  to.return <- list(distx = distx, disty = disty, distz = distz, time = time, nhop = nhop)
  attr(to.return, "class") <- "KMC"

  return(to.return)
}
