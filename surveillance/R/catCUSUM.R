#########################################################################
# Categorical CUSUM for y_t \sim M_k(n_t, \pi_t) for t=1,...,tmax
# Workhorse function doing the actual computations - no semantic checks
# are performed here, we expect "proper" input.
#
# Params:
#  y - (k) \times tmax observation matrix for all categories
#  pi0 - (k) \times tmax in-control prob vector for all categories
#  pi1 - (k) \times tmax out-of-control prob vector for all categories
#  dfun - PMF function of the categorical response, i.e. multinomial, binomial,
#         beta-binom, etc.
#  n   - vector of dim tmax containing the varying sizes
#  h   - decision threshold of the Categorical CUSUM
#########################################################################


catcusum.LLRcompute <- function(y, pi0, pi1, h, dfun, n, calc.at=TRUE,...) {
  #Initialize variables
  t <- 0
  stopped <- FALSE
  S <- numeric(ncol(y)+1)
  U <- numeric(ncol(y)+1)

  ##Check if dfun is the binomial
  isBinomialPMF <- isTRUE(attr(dfun,which="isBinomialPMF"))

  #Run the Categorical LR CUSUM
  while (!stopped) {
    #Increase time
    t <- t+1

    #Compute log likelihood ratio
    llr <-  dfun(y=y[,t,drop=FALSE], size=n[t], mu=pi1[,t,drop=FALSE], log=TRUE,...) - dfun(y=y[,t,drop=FALSE], size=n[t], mu=pi0[,t,drop=FALSE], log=TRUE, ...)
    #Add to CUSUM
    S[t+1] <- max(0,S[t] + llr)

    #For binomial data it is also possible to compute how many cases it would take
    #to sound an alarm given the past.
    if (nrow(y) == 2 & calc.at) {
      ##For the binomial PMF it is possible to compute the number needed for an
      ##alarm exactly
      if (isBinomialPMF) {
        ##Calculations in ../maple/numberneededbeforealarm.mw.
        at <- (h - S[t] - n[t] * ( log(1 - pi1[1,t]) - log(1-pi0[1,t]))) / (log(pi1[1,t]) - log(pi0[1,t]) - log(1-pi1[1,t]) + log(1-pi0[1,t]))
        U[t+1] = ceiling(max(0,at))
        ##Note: U[t+1] Can be higher than corresponding n_t.
        if (U[t+1]>n[t]) U[t+1] <- NA
      } else {
        #Compute the value at by trying all values betweeen 0 and n_t. If
        #no alarm, then we know the value for an alarm must be larger than y_t
        if (S[t+1]>h) {
          ay <- rbind(seq(0,y[1,t],by=1),n[t]-seq(0,y[1,t],by=1))
        } else {
          ay <- rbind(seq(y[1,t],n[t],by=1),n[t]-seq(y[1,t],n[t],by=1))
        }

        llr <-  dfun(ay, size=n[t], mu=pi1[,t,drop=FALSE], log=TRUE,...) - dfun(ay, size=n[t], mu=pi0[,t,drop=FALSE], log=TRUE, ...)
        alarm <- llr > h-S[t]

        ##Is any a_t==TRUE?, i.e. does a y_t exist or is the set over which to
        ##take the minimum empty?
        if (any(alarm)) {
          U[t+1] <- ay[1,which.max(alarm)]
        } else {
          U[t+1] <- NA
        }
      }
    }

    ##Only run to the first alarm. Then reset.
    if ((S[t+1] > h) | (t==ncol(y))) { stopped <- TRUE}
  }

  ##If no alarm at the end put rl to end (its censored! hoehle: Actually it should be length+1!
  ##but the chopping is written such that copying occurs until the final index (hence we can't
  ##just do ncol(pi0)+1
  ##Hence, N is more like the last index investigated.
  if (any(S[-1]>h)) {
    t <- which.max(S[-1] > h)
  } else {
    t <- ncol(pi0) ##Last one
  }
  ##Missing: cases needs to be returned!
  return(list(N=t,val=S[-1],cases=U[-1]))
}

######################################################################
## Wrap function to process sts object by categoricalCUSUM (new S4
## style). Time varying number of counts is found in slot populationFrac.
##
## Params:
##  control - list with the following components
##    * range - vector of indices in disProgObj to monitor
##    * h     - threshold, once CUSUM > h we have an alarm
##    * pi0 - (k-1) \times tmax in-control prob vector for all but ref cat
##    * pi1 - (k-1) \times tmax out-of-control prob vector for all but ref cat
##    * dfun - PMF to use for the computations, dmultinom, dbinom, dBB, etc.
## ... - further parameters to be sent to dfun
######################################################################

categoricalCUSUM <- function(stsObj,
                               control = list(range=NULL,h=5,
                                 pi0=NULL, pi1=NULL, dfun=NULL, ret=c("cases","value")),...) {

  ##Set the default values if not yet set
  if(is.null(control[["pi0",exact=TRUE]])) {
    stop("Error: No specification of in-control proportion vector pi0!")
  }
  if(is.null(control[["pi1",exact=TRUE]])) {
    stop("Error: No specification of out-of-control proportion vector pi1!")
  }
  if(is.null(control[["dfun",exact=TRUE]])) {
    stop("Error: No specification of the distribution to use, e.g. dbinom, dmultinom or similar!")
  }

  if(is.null(control[["h",exact=TRUE]]))
    control$h <- 5
  if(is.null(control[["ret",exact=TRUE]]))
  	control$ret <- "value"

  ##Extract the important parts from the arguments
  range <- control$range
  y <- t(stsObj@observed[range,,drop=FALSE])
  pi0 <- control[["pi0",exact=TRUE]]
  pi1 <- control[["pi1",exact=TRUE]]
  dfun <- control[["dfun",exact=TRUE]]
  control$ret <- match.arg(control$ret, c("value","cases"))
  ##Total number of objects that are investigated. Note this
  ##can't be deduced from the observed y, because only (c-1) columns
  ##are reported so using: n <- apply(y, 2, sum) is wrong!
  ##Assumption: all populationFrac's contain n_t and we can take just one
  n <- stsObj@populationFrac[range,1,drop=TRUE]

  ##Semantic checks
  if ( ((ncol(y) != ncol(pi0)) | (ncol(pi0) != ncol(pi1))) |
      ((nrow(y) != nrow(pi0)) | (nrow(pi0) != nrow(pi1)))) {
    stop("Error: dimensions of y, pi0 and pi1 have to match")
  }
  if ((control$ret == "cases") & nrow(pi0) != 2) {
    stop("Cases can only be returned in case k=2.")
  }
  if (length(n) != ncol(y)) {
    stop("Error: Length of n has to be equal to number of columns in y.")
  }
  ##Check if all n entries are the same
  if (!all(apply(stsObj@populationFrac[range,],1,function(x) all.equal(as.numeric(x),rev(as.numeric(x)))))) {
    stop("Error: All entries for n have to be the same in populationFrac")
  }

  ##Reserve space for the results
  ##start with cusum[timePoint -1] = 0, i.e. set cusum[1] = 0
  alarm <- matrix(data = FALSE, nrow = length(range), ncol = nrow(y))
  upperbound <- matrix(data = 0, nrow = length(range), ncol = nrow(y))

  ##Small helper function to be used along the way --> move to other file!
  either <- function(cond, whenTrue, whenFalse) {
    if (cond) return(whenTrue) else return(whenFalse)
  }

  ##Setup counters for the progress
  doneidx <- 0
  N <- 1
  noofalarms <- 0
  noOfTimePoints <- length(range)

  #######################################################
  ##Loop as long as we are not through the entire sequence
  #######################################################
  while (doneidx < noOfTimePoints) {
    ##Run Categorical CUSUM until the next alarm
    res <- catcusum.LLRcompute(y=y, pi0=pi0, pi1=pi1, n=n, h=control$h, dfun=dfun,calc.at=(control$ret=="cases"),...)

    ##In case an alarm found log this and reset the chart at res$N+1
    if (res$N < ncol(y)) {
      ##Put appropriate value in upperbound
      upperbound[1:res$N + doneidx,]  <- matrix(rep(either(control$ret == "value", res$val[1:res$N] ,res$cases[1:res$N]),each=ncol(upperbound)),ncol=ncol(upperbound),byrow=TRUE)
      alarm[res$N + doneidx,] <- TRUE

      ##Chop & get ready for next round
      y <- y[,-(1:res$N),drop=FALSE]
      pi0 <- pi0[,-(1:res$N),drop=FALSE]
      pi1 <- pi1[,-(1:res$N),drop=FALSE]
      n <- n[-(1:res$N)]

      ##Add to the number of alarms
      noofalarms <- noofalarms + 1
    }

    doneidx <- doneidx + res$N
  }

  ##Add upperbound-statistic of last segment (note: an alarm might or might be reached here)
  upperbound[(doneidx-res$N+1):nrow(upperbound),]  <- matrix( rep(either(control$ret == "value", res$val, res$cases),each=ncol(upperbound)),ncol=ncol(upperbound),byrow=TRUE)
  ##Inherit alarms as well (last time point might contain an alarm!)
  alarm[(doneidx-res$N+1):nrow(upperbound),] <- matrix( rep(res$val > control$h,each=ncol(alarm)), ncol=ncol(alarm),byrow=TRUE)

  # Add name and data name to control object
  control$name <- "categoricalCUSUM"
  control$data <- NULL #not supported anymore

  #New direct calculations on the sts object
  stsObj@observed <- stsObj@observed[control$range,,drop=FALSE]
  stsObj@epoch <- stsObj@epoch[control$range,drop=FALSE]
  stsObj@state <- stsObj@state[control$range,,drop=FALSE]
  stsObj@populationFrac <- stsObj@populationFrac[control$range,,drop=FALSE]
  stsObj@alarm <- alarm
  stsObj@upperbound <- upperbound
  stsObj@control <- control

  #Fix the corresponding start entry
  if (stsObj@epochAsDate==FALSE){
    start <- stsObj@start
    new.sampleNo <- start[2] + min(control$range) - 1
    start.year <- start[1] + (new.sampleNo - 1) %/% stsObj@freq
    start.sampleNo <- (new.sampleNo - 1) %% stsObj@freq + 1
    stsObj@start <- c(start.year,start.sampleNo)
  } else {
    stsObj@start <- c(isoWeekYear(epoch(stsObj)[1])$ISOYear,isoWeekYear(epoch(stsObj)[1])$ISOWeek)
  }

  #Ensure dimnames in the new object ## THIS NEEDS TO BE FIXED!
  #stsObj <- fix.dimnames(stsObj)

  #Done
  return(stsObj)
}
