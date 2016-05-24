buildsurfs <-
function(reps=1000, verbose=TRUE) {

  #--------------------------------------------------------------
  # 
  # TITLE:     buildsurfs()
  # AUTHOR:    TARMO REMMEL
  # DATE:      23 JULY 2013
  # CALLS:     CARsimu(), ClassStat()
  # CALLED BY: NA
  # NEEDS:     SDMTools LIBRARY
  # NOTES:     USED TO BUILD AN ARRAY OF EXPECTED CLASS METRIC
  #            RESULTS ALONG WITH THEIR VARIABILITY BASED ON
  #            1000 REALIZATIONS ACROSS 9 LEVELS OF THEMATIC
  #            PROPORTION, AND 11 LEVELS OF RHO.  RESULTS ARE
  #            STORED FOR 38 CLASS METRICS.
  #--------------------------------------------------------------

  library(SDMTools)

  # DEFINE BINARY PROPORTION INTERVALS (9 OF THEM)
  propvals <- seq(10,90,by=10)

  # DEFINE RHO INTERVALS (11 OF THEM)
  rhovals <- seq(0,0.2499999, by=0.2499999/10)

  # BUILD AN ARRAY TO STORE RESULTS DIMENSIONS:[METRIC,PROPORTION,RHO,REPLICATE]
  storage <- array(data=NA, dim=c(38,9,11,reps))
  
  # LOOP THROUGH COMBINATIONS AND SIMULATE REPLICATES COMPUTING CLASS METRICS FOR EACH
  for(prop in 1:9) {
    for(rho in 1:11) {
      for(replicate in 1:reps) {
        
        if(verbose) {
          # INDICATE WHICH COMBINATION AND REPLICATE IS CURRENTLY BEING PROCESSED
          cat(propvals[prop], rhovals[rho], replicate, "\n", sep=" ")
        }
        
        # PRODUCE SIMULATED REALIZATION WITH GIVEN RHO AND PROPORTION PARAMETERS
        realizationtemp <- CARsimu(rho = rhovals[rho], rajz = FALSE)
        realization <- quantile(realizationtemp, propvals[prop]/100)
        GARB <- realizationtemp > realization[1]
        GARB <- factor(GARB)
        GARB <- as.numeric(GARB)
        realization <- GARB
        dim(realization) <- c(64,64)

        # COMPUTE AND STORE CLASS METRICS   
        results <- ClassStat(realization)
                  
        # WRITE METRICS TO APPROPRIATE ARRAY AND LOCATION
        storage[,prop,rho,replicate] <- results[1,]
        dim(storage) <- c(38,9,11,reps)
              
        
      } # END FOR: REPLICATE
    } # ENF FOR: RHO 
  } # END FOR: PROP
  
  if(verbose) {
    cat("\nDone.\n")
  }

  storage <- as.numeric(storage)
  storage <- array(storage, c(38,9,11,reps))
  return(storage)
  
}
