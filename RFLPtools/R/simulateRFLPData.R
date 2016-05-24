###############################################################################
## R Code to simulate reference data
###############################################################################

## N: number of samples for each of the number of bands
## nrBands: number of bands for which data is simulated
##          total number of samples: N*length(nrBands)
## BandCenters: molecular weights are randomly generated around "band centers"
## delta: use uniform distribution with +/- delta around the band centers
## refData: if TRUE also Taxonname and Accession are added
simulateRFLPdata <- function(N = 10, nrBands = 3:12, 
                             bandCenters = seq(100, 800, by = 100),
                             delta = 50, refData = FALSE){
  if(length(N) > 1){
    N <- N[1]
    warning("Only the first element of 'N' is used.")
  }
  if(N <= 0) stop("'N' has to be a positive integer!")
  N <- trunc(N)
  
  if(any(nrBands <= 0)) stop("'nrBands' has to be a vector of positive integer!")
  if(any(bandCenters <= 0)) stop("'bandCenters' has to be a vector of positive reals!")
  if(length(delta) > 1){
    delta <- delta[1]
    warning("Only the first element of 'delta' is used.")
  }
  if(delta <= 0) stop("'delta' has to be a positive real!")
  
  ## data matrix
  simData <- matrix(NA, nrow = sum(N*nrBands), ncol = 3)
  colnames(simData) <- c("Sample", "Band", "MW")
  
  row.count <- 0
  sample.count <- 1
  for(i in nrBands){
    for(j in 1:N){
      ## randomly select "band centers" (with replacement!)
      Bcent <- sample(bandCenters, i, replace = TRUE)
      ## simulate molecular weights
      simData[row.count+(1:i),] <- c(rep(sample.count, i),
                                     1:i,
                                     sort(runif(i, min = Bcent-delta, max = Bcent+delta)))
      row.count <- row.count + i
      sample.count <- sample.count + 1
    }
  }
  
  ## Generate data.frame and add column Enzyme 
  simData <- data.frame(simData, Enzyme = "Enzyme 1")
  simData$Sample <- paste("Sample", simData$Sample)
  if(refData){
    simData$Taxonname <- simData$Sample
    simData$Accession <- simData$Sample
  }
  simData
}
