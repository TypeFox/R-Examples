## "_helpers" files are not exported in namespace.

#########################
####################### not sure what below code is good for







###### this function is a waste of time. absolute wrong way.
###### You want to condition on each starting point and then paste together.
## estimtae E(fnc(N+) | data) etc for N-, Rt.
acceptcounts <- function(sims, data, fnc=function(x){x;}){
  N <- length(sims);
  indices <- condIndices(sims, data);
  sims <- sims[indices];
  condSimsCount <- length(sims);

  print(sum(indices))
  
  #Now calculate what you want; here we just calculate summary stat stuff.
  summStats <- sapply(sims, function(aList){BDsummaryStats(aList)}, simplify=TRUE); # need pass T
  results <- apply(summStats, 1, function(x){sum(fnc(x));});
#  for (a in 1:N){
#    if (sims[[a]]$states[ length(sims[[a]]$states)] == endstate){
#      endstateCount <- endstateCount+1;
#      waits <- waitTimes(sims[[a]]$states, sims[[a]]$times, T)
#      jumps <- NijBD(sims[[a]]$states);
#      maxState <- length(jumps[1,])-1;
#      Holdtime <- Holdtime + seq(0,maxState,1) %*% waits;
#      Nplus <- Nplus + sum(jumps[2,]);
   #   Nminus <- Nminus + sum(jumps[1,]);
  #  }
  #}
  
  names(results) <- c("Nplus", "Nminus", "Holdtime");
  fullresults <- list(results,  results / N, results / condSimsCount);
  names(fullresults) <- c("count", "joint", "cond");
  fullresults;

}


#accept-reject algorithm
##function for simulating birth death and counting wait times nad jump up
## sims is list of BDMCs
acceptcounts.1 <- function(sims,  endstate=1){
  waits <- NULL;  jumps <- NULL;
  endstateCount <- 0;
  Holdtime <- Nplus <- Nminus <- 0;
  N <- length(sims);
  endStates <- unlist(lapply(sims, function(aList){ getStates(aList)[length(getStates(aList))]; }));  
  sims <- sims[endStates == endstate];
  endstateCount <- length(sims);
  summStats <- sapply(sims, function(aList){BDsummaryStats(aList)}, simplify=TRUE); # need pass T
  results <- apply(summStats, 1, sum)
  names(results) <- c("Nplus", "Nminus", "Holdtime");
  fullresults <- list(results,  results / N, results / endstateCount);
  names(fullresults) <- c("count", "joint", "cond");
  fullresults;
}

