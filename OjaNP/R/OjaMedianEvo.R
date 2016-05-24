`ojaMedianEvo` <-
function(X, control=ojaMedianControl(...), ...){
   x <- y <- 1
   
   rows <- dim(X)[1]
   cols <- dim(X)[2]

    solution <-  .Call("ojaEvo",data=X, initialSigma=as.numeric(control$sigmaInit), sigmaAdaptation=as.numeric(control$sigmaAda), adaptationFactor=control$adaFactor, iterations=as.numeric(control$iter), 
                        useAllSubsets=control$useAllSubsets, numberOfSubsetsUsed=as.numeric(control$nSubsetsUsed), sigmaLog10Decrease=as.numeric(control$sigmaLog10Dec),
                       storeSubdeterminants=control$storeSubDet);
     output <- solution$best
     output
  }
