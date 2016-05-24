
#sims should be a list in which each element is of type CTMC, i.e. has
# states, times and T. (times includes 0 as first time).
#The CTMCs should fit together!  Because this stores only the times of jumps,
#i.e. in a "fully observed" format, regardless of how many observations are
#passed in
#Returns a CTMC from alliding these CTMCs together.
combineCTMC <- function(sims){
  currentTime <- 0;
  numSims <- length(sims);
  if(is(sims[[1]],"CTMC")){
    combinedCTMC <- list(states= getStates(sims[[1]])[1], times=c(0), T=0);
  }
  else {
    combinedCTMC <- list(states= sims[[1]]$states[1], times=c(0), T=0);
  }
  #combinedCTMC <- list(states= NULL, times=NULL, T=0);
  for (i in 1:numSims){
    currCTMC <- sims[[i]];

    if (is(currCTMC, "CTMC")){# S4 class, not just a list
      currCTMClength <- length(getStates(currCTMC));
      currStates <- getStates(currCTMC)[2:currCTMClength]
      currTimes <- getTimes(currCTMC)[2:currCTMClength];
      currT <- getT(currCTMC);
    }
    else{
      currCTMClength <- length(currCTMC$states);
      currStates <- currCTMC$states[2:currCTMClength]
      currTimes <- currCTMC$times[2:currCTMClength];
      currT <- currCTMC$T;
    }
    if (currCTMClength > 1){
      combinedCTMC$states <- c(combinedCTMC$states, currStates);
      combinedCTMC$times <- c(combinedCTMC$times, currentTime+currTimes);
    }
    currentTime <- currentTime + currT;
  }
  combinedCTMC$T <- currentTime;
  combinedCTMC;
}
