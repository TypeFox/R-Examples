###############################################################
#
# check out how much bias protection is obtainable for
#     95% efficiency in ideal model in particular models
#
###############################################################

require(ROptEst)

# some check function
checkOut <- function(L2M, nbd, extraICs = NULL, biastype = symmetricBias(),
                     normtype= NormType()){
  ## L2M : L2 model
  ## nbd : neighborhood
  ## extraICs : list of further ICs)
  ## biastype and normtype (...)
  
  force(normtype)
  force(biastype)
  RobM <- InfRobModel(center = L2M, neighbor = nbd)
  ## rmx - IC
  IC.rmx <- radiusMinimaxIC(L2Fam=L2M, neighbor=nbd,
                            risk=asMSE(biastype=biastype, normtype=normtype),
                            upRad=15,loRad=0.005)
  ## Anscombe - IC
  IC.ans <- optIC(model = RobM, risk = asAnscombe(biastype=biastype,
                                                  normtype=normtype))
  ## Lower case / Most Bias robust Estimator
  IC.mbe <- optIC(model = RobM, risk = asBias(biastype=biastype,
                                              normtype=normtype), tol = 1e-10)

  todo <- list(rmx=IC.rmx,ans=IC.ans,mbe=IC.mbe)
  if(!is.null(extraICs)){
     namICs <- names(extraICs)
     extraICs <- lapply(extraICs, function(x) makeIC(x,L2M))
     names(extraICs) <- namICs
     todo <- c(todo,extraICs)
  }
  ie <- 1/c(unlist(lapply(todo,getMaxIneff, neighbor=nbd)))
}

contnb <- ContNeighborhood(radius = 0.5)
totvnb <- TotalVarNeighborhood(radius = 0.5)
medianmad <- list(function(x)sign(x),function(x)sign(abs(x)-qnorm(.75)))
### Normal location and scale --- takes ~2min:
system.time({print(round(mineff.ls <- checkOut(L2M = NormLocationScaleFamily(), nbd = contnb,
         extraICs = list(medmad=medianmad)),3))})
### Normal location
system.time(print(round(mineff.l <- checkOut(L2M = NormLocationFamily(), nbd = contnb),3)))
### Normal scale convex contamination:
system.time(print(round(mineff.sc <- checkOut(L2M = NormScaleFamily(), nbd = contnb),3)))
### Normal scale total variation
system.time(print(round(mineff.sv <- checkOut(L2M = NormScaleFamily(), nbd = totvnb),3)))
### Poisson(lambda=1) convex contamination:
system.time(print(round(mineff.pc <- checkOut(L2M = PoisFamily(lambda = 1), nbd = contnb),3)))
### Poisson(lambda=1) scale convex contamination:
system.time(print(round(mineff.pv <- checkOut(L2M = PoisFamily(lambda = 1), nbd = totvnb),3)))



