###Charles Doss

### These are functions that apply to Continuous Time Markov Chains (CTMCs).
###  They should eventually be moved into an OOP paradigm, but haven't yet.
### A "CTMC" is a list with components "times," "states," and "T."




#function for CTMC; count number jumps from i to j, for all i,j.
#CTMC should be a fully or partially observed ctstime real markovchain vector.
#For now it is presumed the states are integers, and are bounded below by zero.
#careful with indices, off-by-one: index1 is state0, 2 is 1 ,etc...
#returns square matrix
Nij <- function(CTMC){
  n <- length(CTMC);
  currIdx <- CTMC[1]+1;
  nextIdx <- NULL;
  result <- matrix(0, max(CTMC)+1, max(CTMC)+1); #Don't need diagonals, but oh well.
  
  for (count in 1:(n-1)){
    nextIdx <- CTMC[count+1]+1 ;
    result[currIdx,nextIdx] <- result[currIdx,nextIdx]+1;
    currIdx <- nextIdx;
  }
  result;
}

#####begin utility function: getPartialData ###########################

#This function takes a full CTMC history (ie simulated one) and reads off
#the value of the CTMC at fixed points. i.e. it turns full data into partial data.
#observeTimes is vector of times at which process is to be observed. INCREASING.
#CTMC is the actual markov chain; (components times,states, and T).

getPartialData <- function (observeTimes, CTMC){
  if (inherits(CTMC, "CTMC")){ 
    times <- getTimes(CTMC);
    states <- getStates(CTMC);
    T <- getT(CTMC);
  }
  else if (class(CTMC)[1] == "list"){
    times <- CTMC$times;
    states <- CTMC$states;
    T <- CTMC$T;
  }
  len <- length(observeTimes);
  observeStates <- vector("numeric", length=len);
  numJumps <- length(states);
  if (observeTimes[1] < times[1] || observeTimes[len] > T){
    print("Can't observe the CTMC at the given times. error.")
    stop();
  }
  j <- 2; #j is one larger than the state you'll be in
  for (i in 1:len){
    ##DONT need to reset j, because the observeTimes should be increasing.
    while (observeTimes[i] > times[j] && j <= numJumps){
      j <- j+1;
    }
    if (j<=numJumps){
      observeStates[i] <- states[j-1];
    }
    else{
      observeStates[i] <- states[numJumps];
    }
  }
  new("CTMC_PO_1", states=observeStates, times=observeTimes)
}


#####end utility function: getPartialData ###########################


setMethod("plot", signature(x="CTMC", y="missing"),
          function(x,y, xlab="Time", ylab="State", type="s", ... ){
            times <- getTimes(x)
            states <- getStates(x)
            n <- length(states)
            times <- c(times, getT(x))
            states <- c(states, states[n])
            plot(times, states, type=type, ylab=ylab, xlab=xlab, ...);
          });


setMethod("plot", signature(x="CTMC_PO_1", y="missing"),
          function(x, xlab="Time", ylab="State", type="l", ...){
            plot(getTimes(x), getStates(x), type=type,
                 ylab=ylab, xlab=xlab,
                 ...);
          })


#Get state of a markov chain stored in my format at increasing sorted times Ts
#For use with fully observed markov chains
getMCstate <- function(CTMC, Ts){
  if (is(CTMC, "CTMC")){
    theTimes <- getTimes(CTMC);
    theT <- getT(CTMC);
    theStates <- getStates(CTMC);
  }
  else {
    theTimes <- CTMC$times;
    theT <- CTMC$T
    theStates <- CTMC$states;
  }
  N <- length(theTimes);
  M <- length(Ts);
  if ( (Ts[1]< theTimes[1]) || (Ts[M]> theT) ) {
    print("getMCstate: bad t passed");
    stop();
  }
  theStates[findInterval(Ts, theTimes)];
}

#Get the ("sub") markov chains of argument CTMC up until time T
#For use with fully observed markov chains
getSubMC <- function(CTMC, T){
  theTimes <- CTMC$times;
  N <- length(theTimes);
#  M <- length(Ts);
#  if ( (Ts[1]< theTimes[1]) || (Ts[M]> CTMC$T) ) {
  if ( (T < theTimes[1]) || (T > CTMC$T) ) {
    print("getMCstate: bad t passed");
  }
  newCTMC <- list(states=-1,times=-1,T=-1);
  lastIdx <- findInterval(T, theTimes);
  newCTMC$states <- CTMC$states[1:lastIdx];
  newCTMC$times <- CTMC$times[1:lastIdx];
  newCTMC$T <- T;
  newCTMC
}


#above function returns the first list, NULLs included, if you want..
##NOTE: returns elt in list; so i is off by 1. For time of first jumpl, pass i=2!
getIthJumpTimes <- function(timesList, i){
  result<- sapply(timesList, function(timeVec){if (length(timeVec)>=i){timeVec[i]} }, simplify=FALSE);
  result <- result[sapply(result, function(x){!is.null(x)}, simplify=TRUE)];
  result <- simplify(result); #list to vector
  ##should remove uses of 'simplify' use 'unlist' function  
  ##  result <- unlist(result);
}

#Takes one ctmc rather than the above getIthJumpTimes which takes a whole list of times
# either ith jumptime or null if no ith jumptime.
##NOTE: returns elt in list; so i is off by 1. For time of first jumpl, pass i=2!
getIthJumpTime <- function(CTMC, i){
  if (is(CTMC, "CTMC")){
    myTimes <- getTimes(CTMC);
  }
  else {
    myTimes <- CTMC$times
  }
  if (length(myTimes) >= i ){ return(myTimes[i]); }
  else return(NULL);
}

getIthState <- function(CTMC, i){
  if (is(CTMC, "CTMC")){
    myStates <- getStates(CTMC);
  }
  else {
    myStates <- CTMC$states
  }
  if (length(myStates) >= i ){ return(myStates[i]); }
  else return(NULL);
}


