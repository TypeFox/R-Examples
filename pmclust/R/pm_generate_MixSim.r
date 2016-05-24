### This file contains a data generation by MixSim.

generate.MixSim.spmd <- function(N, p, K, MixSim.obj = NULL,
    MaxOmega = NULL, BarOmega = NULL, PiLow = 1.0, sph = FALSE, hom = FALSE){
  ### Obtain a set of parameters from MixSim.
  if(spmd.comm.rank() == 0){
    if(is.null(MixSim.obj)){
      MixSim.obj <- MixSim::MixSim(BarOmega, MaxOmega, K = K, p = p,
                                   PiLow = PiLow, sph = sph, hom = hom)
    }
    if(class(MixSim.obj) != "MixSim"){
      stop("MixSim.obj is not a MixSim class.")
    }
  } else{
    MixSim.obj <- NULL
  }
  MixSim.obj <- spmd.bcast.object(MixSim.obj)

  ### Arrange samples for each processors.
  N.allspmds <- unlist(lapply(get.jid(N, all = TRUE), length),
                       use.names = FALSE)
  N.spmd <- N.allspmds[spmd.comm.rank() + 1]
  ret.dataset <- MixSim::simdataset(N.spmd, MixSim.obj$Pi, MixSim.obj$Mu,
                                    MixSim.obj$S)

  data.simu <- ret.dataset$X
  data.class <- ret.dataset$id
  data.n.class <- tabulate(data.class, nbins = K)
  data.n.class <- spmd.allreduce.double(as.double(data.n.class), double(K),
                                        op = "sum")

  ret <- list(K = K, p = p, N = N, N.allspmds = N.allspmds,
              N.spmd = N.spmd,
              X.spmd = data.simu, CLASS.spmd = data.class,
              N.CLASS.spmd = data.n.class,
              MixSim.obj = MixSim.obj)
  ret
} # End of generate.MixSim.spmd().

generate.MixSim <- generate.MixSim.spmd
