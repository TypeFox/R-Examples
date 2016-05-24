MMRtime <- function(x, d, age) {
#### estimate Mean/Median Residual lifetime over age. 
temp <- WKM( x=x , d=d )
tivec <- temp$times
pivec <- temp$jump

if( age >= tivec[length(tivec)] ) stop("age too large")
if( age < tivec[1] ) warning("age smaller than first event time")

pivec[ tivec < age ] <- 0
Sage <- sum( pivec )

fenzi <- sum( (tivec - age)*pivec )
MRtime <- fenzi/Sage

Ptheta <- Sage/2
Cprob <- cumsum(pivec)
posi <- sum(Cprob < Ptheta)
theta <- tivec[posi+1]

list(MeanResidual = MRtime, MedianResidual = theta - age)
}
