#Classes for use with birth-death-immigration, fully and partially observed.

#note: CTMC is integer valued
#times should be increasing.
setClass(Class="CTMC", representation(states="numeric", times="numeric", T="numeric"),
         validity=function(object){
           #Partially observed classes are free; just check it has the right
           #slots.  Should be built in S4 but documentation is nonexistent.
           n <- length(object@times);
           errorRet <- "";
           if (!(mode(object@states) == "numeric"))
             errorRet <- paste(errorRet,
                               "Need to pass in a 'states' numeric vector.");
           if (!(mode(object@times) == "numeric"))
             errorRet <- paste(errorRet, "Need to pass in a 'times' numeric vector.");
           if (!(mode(object@T)== "numeric" && (length(object@T)==1)) )
             errorRet <- paste(errorRet, "Need to pass in a 'T' length 1 numeric.");
           if (!(length(object@states) == length(object@times)))
             errorRet <- paste(errorRet, "Length of 'states' and length of 'times' should be equal.");
           if (n>1){
             diffs <- object@times[2:n]-object@times[1:(n-1)];
             if (!all(diffs>=0)){
               errorRet <- paste(errorRet, "Times should be an increasing sequence of (real) numbers.");
               print("You passed in for times")
               print(object@times); ### would prefer this to be part of errorRet ; need to paste a sequence into a string, not sure how

             }
           }
           if (!(object@T >= object@times[n])){
             errorRet <- paste(errorRet, "Need total observation time 'T' larger than final data observation time 'times[n]'.");
             print(paste("T is",object@T), digits=20)
             print("and times are")
             print(object@times, digits=20)
           }
           if (errorRet == "")
             {return(TRUE);}
           else {
             return(errorRet);
           }
         })


setClass(Class="BDMC", contains=c("CTMC"),
         validity=function(object){
           ret <- (getValidity(getClassDef("CTMC")))(object);
           n <- length(object@states); #CTMC checks this is equal lengthtimes
           if (n>=2){ #else it's automaticallyok
             diffs <- object@states[2:n] - object@states[1:(n-1)];
             ret <- ret&& ((sum(diffs==-1) + sum(diffs==1)) == (n-1));
           }
           ret;
         }) #same data types

setClass(Class="CTMC_PO_1", representation(states="numeric", times="numeric"),
         validity=function(object){
           #Partially observed classes are free; just check it has the right
           #slots.  Should be built in S4 but documentation is nonexistent.
           return(
                  (mode(object@states) == "numeric" &&
                  mode(object@times) == "numeric" &&
                  (length(object@states) == length(object@times)))
                  ||
                  (mode(object@states) == "NULL" &&
                   mode(object@times) == "NULL")
                  );
           })


setClass(Class="CTMC_PO_many", representation(BDMCsPO="list"),
         validity=function(object){
           BD1checker <- getValidity(getClassDef("CTMC_PO_1"));
           bools <- sapply(object@BDMCsPO, BD1checker);
           return(all(bools));
         })

###
setClass(Class="CTMC_many", representation(CTMCs="list"),
         validity=function(object){
           ctmc1Checker <- getValidity(getClassDef("CTMC"));
           bools <- sapply(object@CTMCs, ctmc1Checker);
           return(all(bools))
         })

setClass(Class="BDMC_many", contains=c("CTMC_many"),
         validity=function(object){
           ## DONT need to do CTMCmany check, b/c CTMC is
           ## checked by Validity of "BDMC" already.
           ##ret <- (getValidity(getClassDef("CTMC_many")))(object);
           bd1Checker <- getValidity(getClassDef("BDMC"));
           bools <- sapply(object@CTMCs, bd1Checker)
           return(all(bools))
         })



if (!isGeneric("getBDMCsPOlist")) {
  if (is.function("getBDMCsPOlist"))
    fun <- getBDMCsPOlist
  else fun <- function(object) standardGeneric("getBDMCsPOlist")
  setGeneric("getBDMCsPOlist", fun)
}



if (!isGeneric("getStates")) {
  if (is.function("getStates"))
    fun <- getStates
  else fun <- function(object) standardGeneric("getStates")
  setGeneric("getStates", fun)
}


if (!isGeneric("getTimes")) {
  if (is.function("getTimes"))
    fun <- getTimes
  else fun <- function(object) standardGeneric("getTimes")
  setGeneric("getTimes", fun)
}


if (!isGeneric("getT")) {
  if (is.function("getT"))
    fun <- getT
  else fun <- function(object) standardGeneric("getT")
  setGeneric("getT", fun)
}

if (!isGeneric("getTs")) {
  if (is.function("getTs"))
    fun <- getT
  else fun <- function(object) standardGeneric("getTs")
  setGeneric("getTs", fun)
}


setMethod("getStates", "CTMC", function(object) {object@states});
setMethod("getStates", "BDMC", function(object) {object@states});
setMethod("getStates", "CTMC_PO_1", function(object) {object@states});

setMethod("getTimes", "CTMC", function(object) {object@times});
setMethod("getTimes", "BDMC", function(object) {object@times});
setMethod("getTimes", "CTMC_PO_1", function(object) {object@times});

setMethod("getT", "CTMC", function(object) {object@T});
setMethod("getT", "BDMC", function(object) {object@T});
setMethod("getT", "CTMC_PO_1", function(object) {n <- length(object@times); object@times[n]-object@times[1];}); ## use tail(,1)
setMethod("getT", "CTMC_PO_many", function(object) {Ts <- sapply(object@BDMCsPO, getT);  sum(Ts); });
setMethod("getT", "CTMC_many", function(object) {Ts <- sapply(object@CTMCs, getT);  sum(Ts); });
setMethod("getT", "BDMC_many", function(object) {Ts <- sapply(object@CTMCs, getT);  sum(Ts); });

setMethod("getTs", "CTMC_PO_many", function(object) {sapply(object@BDMCsPO, getT);});
setMethod("getTs", "CTMC_many", function(object) {sapply(object@CTMCs, getT);});
setMethod("getTs", "BDMC_many", function(object) {sapply(object@CTMCs, getT);});


setMethod("[", "CTMC_PO_many", function(x,i,j="missing",...,drop=TRUE){
  new("CTMC_PO_many", BDMCsPO=x@BDMCsPO[i,...,drop=drop])
})
setMethod("[[", "CTMC_PO_many", function(x,i,j="missing",...,drop=TRUE){
  x@BDMCsPO[[i,...,drop=drop]];
})
setMethod("[", "CTMC_many", function(x,i,j="missing",...,drop=TRUE){
  new("CTMC_many", CTMCs=x@CTMCs[i,...,drop=drop]);
})
setMethod("[[", "CTMC_many", function(x,i,j="missing",...,drop=TRUE){
  x@CTMCs[[i,...,drop=drop]];
})

##accessor
setMethod("getBDMCsPOlist", "CTMC_PO_many", function(object){object@BDMCsPO})


###Methods for converting from back and forth between lists and classes
CTMC2list <- function(aCTMC){
  states <- getStates(aCTMC);
  times <- getTimes(aCTMC);
  T <- getT(aCTMC);
  list(states=states,times=times,T=T);
}

list2CTMC <- function(aCTMC){
  states <- aCTMC$states;
  times <- aCTMC$times;
  T <- aCTMC$T;
  new("CTMC", states=states,times=times,T=T);
}









#################################begin utility func "BDsummaryStats"


##Gets summary statistics for a single birth death markov chain passed in
##in format $states $times $T

if (!isGeneric("BDsummaryStats")) {
  fun <- function(sim) standardGeneric("BDsummaryStats")
  setGeneric("BDsummaryStats", fun)
}
setMethod("BDsummaryStats", "BDMC_many",
          definition=
          function(sim){
            res <- sapply(sim@CTMCs, BDsummaryStats) ##can i specify ".BDMC" version?
            apply(res,1, sum)            
          })

##### this gives a warning, "no definition for BDMC" upon package installation.
#####   so it's not defining the classes first, for some reason.  why would that be.
setMethod("BDsummaryStats", "BDMC",
          definition=
          function(sim){
            waits <- waitTimes(sim@states, sim@times, sim@T)
            jumps <- NijBD(sim@states);
            maxState <- length(jumps[1,])-1;
            Holdtime <- seq(0,maxState,1) %*% waits;
            Nplus <- sum(jumps[2,]);
            Nminus <-sum(jumps[1,]);
            results <- c(Nplus, Nminus, Holdtime);
            names(results) <- c("Nplus", "Nminus", "Holdtime");
            results;
          })

setMethod("BDsummaryStats", "list",
          definition=
          function(sim){
                                        #      endstateCount <- endstateCount+1;
            waits <- waitTimes(sim$states, sim$times, sim$T)
            jumps <- NijBD(sim$states);
            maxState <- length(jumps[1,])-1;
                                        #  maxState <- max(sim$states);
            Holdtime <- seq(0,maxState,1) %*% waits;
            Nplus <- sum(jumps[2,]);
            Nminus <-sum(jumps[1,]);
            results <- c(Nplus, Nminus, Holdtime);
            names(results) <- c("Nplus", "Nminus", "Holdtime");
            results;
          })



#################################end utility func "BDsummaryStats"




###################begin utility BDsummaryStats.PO 

if (!isGeneric("BDsummaryStats.PO")) {
  fun <- function(dat) standardGeneric("BDsummaryStats.PO")
  setGeneric("BDsummaryStats.PO", fun)
}

setMethod("BDsummaryStats.PO", "list",
          definition=
          function(dat){
            n <- length(dat$states);
            T <- dat$times[n];
            diffs <- dat$states[2:n] - dat$states[1:n-1]
            maxState <- max(dat$states);
            Holdtime <- seq(0,maxState,1) %*% waitTimes(dat$states, dat$times, T);
            results <- c(  sum(diffs[diffs>0]),
                         abs(sum(diffs[diffs<0])),
                         Holdtime);
            names(results) <- c("Nplus", "Nminus", "Holdtime");
            results;
          })

setMethod("BDsummaryStats.PO", "CTMC_PO_1",
          definition=
          function(dat){
            n <- length(dat@states);
            T <- dat@times[n];
            diffs <- dat@states[2:n] - dat@states[1:n-1]
            maxState <- max(dat@states);
            Holdtime <- seq(0,maxState,1) %*% waitTimes(dat@states, dat@times, T);
            results <- c(  sum(diffs[diffs>0]),
                         abs(sum(diffs[diffs<0])),
                         Holdtime);
            names(results) <- c("Nplus", "Nminus", "Holdtime");
            results;
          })

setMethod("BDsummaryStats.PO", "CTMC_PO_many",
          definition=
          function(dat){
            res <- sapply(dat@BDMCsPO, BDsummaryStats.PO)
            apply(res,1,sum)
          })

###################end utility BDsummaryStats.PO 
