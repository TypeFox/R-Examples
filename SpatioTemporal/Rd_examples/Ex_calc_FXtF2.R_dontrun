##create a trend
trend <- cbind(1:5,sin(1:5))
##an index of locations
idx <- c(rep(1:3,3),1:2,2:3)
idx2 <- c(rep(1:2,3),2,2)
##a list of time points for each location/observation
T <- c(rep(1:3,each=3),4,4,5,5)
T2 <- c(rep(1:3,each=2),4,5)

##expand the F matrix to match the locations/times in idx/T.
F <- trend[T,]
F2 <- trend[T2,]

##first column gives time and second location for each observation
cbind(T, idx)
##...and for the second set
cbind(T2, idx2)

##create a cross covariance matrix
C <- makeSigmaB(list(c(1,1),c(1,.5)), crossDist(1:max(idx),1:max(idx2)))

##compute F %*% X %*% F2'
FXtF2 <- calc.FXtF2(F, C, loc.ind=idx, F2=F2, loc.ind2=idx2)

##which is equivalent to
FXtF2.alt <- expandF(F, idx) %*% C %*% t( expandF(F2, idx2) )

range(FXtF2 - FXtF2.alt)
\dontshow{
  if( max(abs(FXtF2 - FXtF2.alt)) > 1e-13 ){
    stop("calc.FXtF2: Results not equal")
  }
}
