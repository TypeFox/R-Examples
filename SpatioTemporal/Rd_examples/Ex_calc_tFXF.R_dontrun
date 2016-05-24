##create a trend
trend <- cbind(1:5,sin(1:5))
##an index of locations
idx <- c(rep(1:3,3),1:2,2:3)
##a list of time points for each location/observation
T <- c(rep(1:3,each=3),4,4,5,5)

##expand the F matrix to match the locations/times in idx/T.
F <- trend[T,]

##create a covariance matrix for locations and each of 
C <- makeSigmaNu(c(1,1), as.matrix(dist(1:max(idx))),
            blocks1=c(3,3,3,2,2), ind1=idx)

##compute F' %*% C %*% F
tFmatF <- calc.tFXF(F, C, idx, block.sizes = c(3,3,3,2,2),
                    n.blocks = 5)
##which is equivalent of
tFmatF.alt <- calc.tFX(F, t(calc.tFX(F, C, idx)), idx)

range(tFmatF-tFmatF.alt)
\dontshow{
  if( max(abs(tFmatF-tFmatF.alt)) > 1e-13 ){
    stop("calc.tFXF 2: Results not equal")
  }
}
