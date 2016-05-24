sparseSignal <- function(N,s,b=1.0,delta=1e-7,nlev=0.05,slev=0.9) {
 # Part of R1Magic by mehmet.suzen@physics.org
 noise <- sample(seq(0,b*nlev,delta),N,replace=TRUE)
 spikeValues <- sample(seq(b*slev,b,delta),s,replace=TRUE)
 sparseS <- matrix(noise,N,1)
 spikeLocations <- sample(1:N,s)
 sparseS[spikeLocations] <- spikeValues
 return(sparseS)
}

