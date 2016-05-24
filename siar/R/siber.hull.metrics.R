siber.hull.metrics <- function(X,Y,G,R=10^4){

reps <- 10^4

## now loop through the data and calculate the ellipses
M <- length(unique(G))

# split the isotope data based on group
spx <- split(X,G)
spy <- split(Y,G)


# some matrices and vectors in which to store the results
sim.mu.X <- matrix(data=0,nrow=reps,ncol=M)
sim.mu.Y <- matrix(data=0,nrow=reps,ncol=M)

for (i in 1:M){

  #-----------------------------------------------------------------------------

  mymodel <- bayesMVN(spx[[i]],spy[[i]])
  #mymodel <- bayestwoNorm(x,y)

  #-----------------------------------------------------------------------------


  if (i == 1){
    simC <- mymodel$b[,1]
    simN <- mymodel$b[,2]
  }
  else{
    simC <- cbind(simC,mymodel$b[,1])
    simN <- cbind(simN,mymodel$b[,2])
  }
  rm(mymodel)
}

# now loop through all the simulated group means
# and calculate the layman metrics
nr <- nrow(simC) # the number of reps.. same as reps defined above

dNr <- numeric(nr)
dCr <- numeric(nr)
TA <- numeric(nr)
CD <- numeric(nr)
MNND <- numeric(nr)
SDNND <- numeric(nr)

for (i in 1:nr) {

  layman <- laymanmetrics( simC[i,], simN[i,] )

  # extract the metrics into their named vectors
  dNr[i] <- layman$dN_range
  dCr[i] <- layman$dC_range
  TA[i] <- layman$hull$TA
  CD[i] <- layman$CD
  MNND[i] <- layman$MNND
  SDNND[i] <- layman$SDNND

}

metrics <- cbind( dNr, dCr, TA, CD, MNND, SDNND)


return(metrics)

}