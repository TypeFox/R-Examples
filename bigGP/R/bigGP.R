bigGP.init <- function(P = NULL, parallelRNGpkg = "rlecuyer", seed = 0){
  if(.bigGP$P == 0) { # .bigGP is P=0 by default 
    cat("Initializing processes:")

    if(mpi.comm.size() == 0){
      if(is.null(P)){
        stop("bigGP.init: You must specify the number of worker processes to initiate (potentially each having multiple cores to allow for threading) and which should be equal to D(D+1)/2 for some integer D.")
      } else{
        mpi.spawn.Rslaves(nslaves = P)
      }
    } else{
      if(!is.null(P) && P != mpi.comm.size() - 1)
        warning(paste("bigGP.init: Number of worker processes requested, ", P, " does not match number of active MPI slave processes, ", mpi.comm.size()-1, "; ignoring 'P' argument.", sep = ""))      
    }
    
    .bigGP.fill()
    
    if(sum(unlist(mpi.remote.exec(require, "bigGP", ret = TRUE))) != .bigGP$P) # i.e., this package
      stop("bigGP.init: error in loading bigGP on slaves.")
    mpi.remote.exec(.Call("init_comms", as.integer(mpi.comm.c2f()), PACKAGE="bigGP"), ret = FALSE)
    out <- .Call("init_comms", as.integer(mpi.comm.c2f()), PACKAGE="bigGP")
    
    if(out != .bigGP$D)
      stop("bigGP.init: number of processes may not be consistent with the partition number, D.")
    
    mpi.bcast.cmd(.bigGP.fill())
        
    cat("... Done.\n")
    
    cat("Using ", .bigGP$P, " processes with a partition size (D) of ", .bigGP$D, ".\n", sep = "")
    
    if(is.null(parallelRNGpkg) || parallelRNGpkg == ""){
      warning("bigGP.init: Initializing process seeds sequentially; not guaranteed to give independent streams; please use rlecuyer or rsprng to be certain.")
      mpi.bcast.Robj2slave(seed)
      mpi.bcast.cmd(set.seed(mpi.comm.rank() + seed)) 
    } else{
      if(parallelRNGpkg == "rlecuyer"){
        if(sum(unlist(mpi.remote.exec(requireNamespace, "rlecuyer", ret = TRUE))) != .bigGP$P)
          stop("bigGP.init: error in using rlecuyer on slaves.")
        RNGkind("L'Ecuyer-CMRG")
        mpi.setup.rngstream(iseed = seed) # first value indicates the RNGkind
      } else{
        if(parallelRNGpkg == "rsprng"){
            stop("bigGP.init: the rsprng package is no longer available on CRAN. Advanced users who install rsprng from the CRAN archived packages can uncomment the appropriate lines in bigGP.init() below this stop() call and rebuild the bigGP package.")
          #if(sum(unlist(mpi.remote.exec(requireNamespace, "rsprng", ret = TRUE))) != .bigGP$P)
          # stop("bigGP.init: error in loading rsprng on slaves.")
          #mpi.bcast.Robj2slave(seed)
          #mpi.bcast.cmd(rsprng::init.sprng(nstream = .bigGP$P, streamno = mpi.comm.rank()-1, seed = seed))
        } else {
          warning('bigGP.init: parallelRNGpkg ', parallelRNGpkg, ' not recognized. Initializing process seeds sequentially; not guaranteed to give independent streams; please use rlecuyer or rsprng to be certain.')
          mpi.bcast.Robj2slave(seed)
          mpi.bcast.cmd(set.seed(mpi.comm.rank() + seed)) 
        }
      }
    }
  }
  invisible(NULL)
}



calcIJ <- function(D) {
  'Finds row and column indices of subblock assigned to the process within the h=1th part of the overall matrix'
  rank <- mpi.comm.rank()
  if( rank < 1 || rank > D * (D + 1) / 2 )
    warning(paste("calcIJ: Invalid rank", rank, ".", sep = " "))
  I <- 0
  J <- 0
  Dd <- D
  rd <- rank - 1 # -1 accounts for shift to master not being a worker process
  while( rd >= Dd ) {
    J <- J + 1
    rd <- rd - Dd
    Dd <- Dd - 1
  }
  I <- J+rd
  return(list(I = as.integer(I), J = as.integer(J)))
}

calcD <- function(P){
  ' Calculates the partition number D given number of processes P; i.e., D such that D*(D+1)/2 <= P '
  D <- floor((sqrt(1+8*P) - 1) / 2)
  if(P != D*(D+1)/2)
    stop(paste("calcD: Number of slave processes, ", P, " is not equal to D(D+1)/2 for integer D.", sep = ""))
  return(as.integer(D))
}
                
getDistributedVectorLength <- function(n, h = 1){
  if(.bigGP$I == .bigGP$J) return(h * ceiling(n / (.bigGP$D * h))) else return(0)
}

getDistributedTriangularMatrixLength <- function(n, h = 1){
  if(.bigGP$I == .bigGP$J){
    return(h * (h+1) / 2 * (ceiling(n / (.bigGP$D * h)))^2)
  } else{
    return(h^2 * (ceiling(n / (.bigGP$D * h)))^2)
  }
}

getDistributedRectangularMatrixLength <- function(n1, n2, h1 = 1, h2 = 1){
  len <- h1 * h2 * ceiling(n1 / (.bigGP$D * h1)) * ceiling(n2 / (.bigGP$D * h2))
  if(.bigGP$I == .bigGP$J) return(len) else return(2*len)
}

remoteGetIndices <- function(type = "vector", objName, objPos = ".GlobalEnv", n1, n2 = NULL, h1 = 1, h2 = 1) {
  if(!is.element(type, c("vector", "triangular", "symmetric", "rectangular")))
    stop("remoteGetIndices: type must be one of 'vector', 'triangular', 'symmetric', 'rectangular'.")
  .n1 <- n1; .n2 <- n2; .h1 <- h1; .h2 <- h2
  mpi.bcast.Robj2slave(.n1)
  mpi.bcast.Robj2slave(.n2)
  mpi.bcast.Robj2slave(.h1)
  mpi.bcast.Robj2slave(.h2)
  if(type == "vector")
    mpi.bcast.cmd(.tmp <- localGetVectorIndices(.n1, .h1))
  if(type == "symmetric" || type == "triangular")
    mpi.bcast.cmd(.tmp <- localGetTriangularMatrixIndices(.n1, .h1))
  if(type == "rectangular")
    mpi.bcast.cmd(.tmp <- localGetRectangularMatrixIndices(.n1, .n2, .h1, .h2))
  mpi.remote.exec(localAssign, objName, ".tmp", objPos)
  remoteRm(.tmp)
  return(NULL)
}
   

localGetVectorIndices <- function(n, h = 1){
  ' Finds the indices of the entries this processor owns of a vector'
  if(.bigGP$I == .bigGP$J){
    bs <- (n + .bigGP$D * h - 1)%/%(.bigGP$D *h)
    ind <- matrix(nrow = bs*h, ncol = 1)
    for( JJ in 0:(h-1) ) {
      ind[JJ*bs+(0:(bs-1))+1, 1] <- JJ*bs*.bigGP$D + .bigGP$J*bs + (0:(bs-1))
    }
    ind[ind >= n] <- 0
  } else{
    ind <- numeric(0)
                                        # ind <- matrix(nrow = 0, ncol = 1) # this was what Ben had, but I don't think he uses the index vector anyway; it seems only to be used as input to the user-defined mean/cov fxns and for my post-processing
  }
  return(ind + 1) # +1 for R 1-based indexing
}

localGetTriangularMatrixIndices <- function(n, h = 1){
  ' Finds the indices of the entries this processor owns of a lower-triangular matrix, as a two-column matrix'
  bs <- (n+.bigGP$D*h-1) %/% (.bigGP$D*h)
  if( .bigGP$I == .bigGP$J ) {
    ind <- matrix(nrow = bs*bs*h*(h+1)/2, ncol = 2)
  } else {
    ind <- matrix(nrow = bs*bs*h*h, ncol = 2)
  }
  start <- 1
  for( JJ in 0:(h-1) ) {
    for( II in JJ:(h-1) ) {
      if( (II == JJ) || (.bigGP$I == .bigGP$J) ) { #one block
        ind[start:(start+(bs*bs-1)), 1] <- II*bs*.bigGP$D+.bigGP$I*bs + rep(0:(bs-1), bs)
        ind[start:(start+(bs*bs-1)), 2] <- JJ*bs*.bigGP$D+.bigGP$J*bs + rep(0:(bs-1), each = bs)
        start <- start + bs*bs
      } else { # two blocks
        ind[start:(start+(bs*bs-1)), 1] <- II*bs*.bigGP$D+.bigGP$I*bs + rep(0:(bs-1), bs)
        ind[start:(start+(bs*bs-1)), 2] <- JJ*bs*.bigGP$D+.bigGP$J*bs + rep(0:(bs-1), each = bs)
        start <- start + bs*bs
        ind[start:(start+(bs*bs-1)), 1] <- II*bs*.bigGP$D+.bigGP$J*bs + rep(0:(bs-1), bs)
        ind[start:(start+(bs*bs-1)), 2] <- JJ*bs*.bigGP$D+.bigGP$I*bs + rep(0:(bs-1), each = bs)
        start <- start + bs*bs
      }
    }
  }                                      
  ind[ind >= n] <- 0
  return(ind + 1)
}
                                    
localGetRectangularMatrixIndices <- function(n1, n2, h1 = 1, h2 = 1) {
  bsr <- (n2 + .bigGP$D*h2 - 1) %/% (.bigGP$D*h2)
  bsc <- (n1 + .bigGP$D*h1 - 1) %/% (.bigGP$D*h1)
  
  if( .bigGP$I == .bigGP$J ) {
    ind <- matrix(nrow=bsr*bsc*h2*h1,ncol=2)
  } else {
    ind <- matrix(nrow=2*bsr*bsc*h2*h1,ncol=2)
  }
  start <- 1
  for( JJ in 0:(h1-1) ) {
    for( II in 0:(h2-1) ) {
      if( (.bigGP$I == .bigGP$J) ) { #one block
        ind[start:(start+bsc*bsr-1),2] = II*bsr*.bigGP$D+.bigGP$I*bsr + rep(0:(bsr-1) , bsc)
        ind[start:(start+bsc*bsr-1),1] = JJ*bsc*.bigGP$D+.bigGP$J*bsc + rep(0:(bsc-1) , each = bsr)
        start = start + bsc*bsr 
      } else { # two blocks
        ind[start:(start+bsc*bsr-1),2] = II*bsr*.bigGP$D+.bigGP$I*bsr + rep(0:(bsr-1) , bsc)
        ind[start:(start+bsc*bsr-1),1] = JJ*bsc*.bigGP$D+.bigGP$J*bsc + rep(0:(bsc-1) , each = bsr)
        start = start + bsc*bsr 
        
        ind[start:(start+bsc*bsr-1),2] = II*bsr*.bigGP$D+.bigGP$J*bsr + rep(0:(bsr-1) , bsc)
        ind[start:(start+bsc*bsr-1),1] = JJ*bsc*.bigGP$D+.bigGP$I*bsc + rep(0:(bsc-1) , each = bsr)
        start = start + bsc*bsr 
      }
    }
  }
  ind[ind[ , 2] >= n2, ] <- 0
  ind[ind[ , 1] >= n1, ] <- 0
  return(ind + 1)
}



alloc <- function(input, inputPos = '.GlobalEnv'){
  # forces allocation of memory to avoid issues with promises being passed into C code and overwriting of input
  if(is.numeric(input)){
    tmp <- 0; length(tmp) <- input
    return(tmp)
  } else{
    if(is.character(input)){
      return(1*get(input, pos = eval(as.name(inputPos))))
    } else stop("alloc: 'input' must be a numeric value or character string.")
  } 
}

.bigGP.fill <- function(init = FALSE) {
  'initializes .bigGP object holding the information on the distributed setup'
  if(init) {
    .bigGP$P <- 0
    .bigGP$D <- 0
    .bigGP$I <- -1
    .bigGP$J <- -1
  } else {
    .bigGP$P <- mpi.comm.size() - 1
    .bigGP$D <- calcD(.bigGP$P)
    if(mpi.comm.rank()) {
      tmp <- calcIJ(.bigGP$D)
      .bigGP$I <- tmp$I
      .bigGP$J <- tmp$J
    }
  }
  invisible(NULL)
}


".bigGP" <- new.env()  # initialize so that it is in the bigGP namespace and can be modified via .bigGP.fill()
.bigGP.fill(init = TRUE)



".onAttach" <- function (lib, pkg) {
  # the following allows bigGP.init() to assign D and P into .bigGP; I'm mimicing spam's use of .Spam
  # not needed now that .bigGP is an environment
#   unlockBinding(".bigGP", asNamespace("bigGP"))
  packageStartupMessage("
=========================================================================================
Loading bigGP.\n
Warning: before using bigGP, you must initialize the slave processes using bigGP.init() \n
(which can also be done indirectly via initializing a krigeProblem object).\n
If R was started through mpirun/orterun/mpiexec, please quit by using bigGP.quit().
=========================================================================================\n
")
}

bigGP.quit <- function(save = "no"){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0) {
      mpi.close.Rslaves()
    }
    mpi.quit(save)
  }
}


bigGP.exit <- function(){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0) {
      mpi.close.Rslaves()
    }
    mpi.exit()
    detach(package:bigGP, unload = TRUE)
  }
}

# removed this for 0.1-3 as it was causing an error in R CMD check:
#* checking whether the namespace can be unloaded cleanly ... WARNING
#---- unloading
#[1] "Detaching Rmpi. Rmpi cannot be used unless relaunching R."
#Warning message:
#.onUnload failed in unloadNamespace() for 'bigGP', details:
#  call: detach(package:Rmpi)
#  error: invalid 'name' argument 

#.onUnload <- function(libpath){
#  if (is.loaded("mpi_initialize")){
#    if (mpi.comm.size(1) > 0) {
#      mpi.close.Rslaves()
#    }
#    mpi.exit()
#  }
#}

if(FALSE) {
  .Last.lib <- function(libpath){
   if (is.loaded("mpi_initialize")){
      if (mpi.comm.size(1) > 0) {
         print("Please use mpi.close.Rslaves() to close slaves.")
         mpi.close.Rslaves()
      }
      print("Please use mpi.quit() to quit R")
      .Call("mpi_finalize")
   }
 }
}
