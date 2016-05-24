globalVariables(c("trialObj","futime","event","exit","trt","totInfec","entry","N",
                  "trialListCensor","out","n","Nvacc","approxOk","V","P","pwList",
                  "N.vax.arms"))

## Create a auxillary function "is.TRUE" to use for logical comparisons
## is.TRUE(x) returns TRUE for each element of 'x' that is TRUE, and
##   returns FALSE for all else (i.e. for FALSE and NA)
is.TRUE <- function(x) 
{
  if ( !is.logical(x) ) stop("Argument to 'is.TRUE' must be of type Logical")
  x & !is.na(x)
}


visitScheduleTruncated <-  function(visitWeeks, thruWeek=NULL, nVisits=NULL)
{
  if ( is.null( c(thruWeek, nVisits) ) )
    stop("Must specify one of the arguments: thruWeek or nVisits\n")
  
  ## determine how many more visit weeks need to be added to 'visitWeeks'
  if ( !is.null(nVisits) )
  {
    visitWeeks[ 1:nVisits ]
  } else
  {
    ## find values of 'visitWeeks' that are strictly greater than 'thruWeek'
    gt_thruWeek <- which( visitWeeks > thruWeek )
    
    if ( length(gt_thruWeek) > 0 ) {
      return( visitWeeks[ 1:gt_thruWeek[1] ] )
    } else {
      ## if we get here it means that thruWeek is the last week of our visit schedule
      ## and for some reason we still need to add another week - don't know why...
      ## I'm going to see if my code will work okay if I just return the entire visit
      ## vector and don't add another 'imaginary' visit after the last one. 
      return( visitWeeks )
    }
  }
}

## This function takes as input a vector of times (for 505 we are
## using time in weeks since first vaccination), and returns the
## time of the next or previous scheduled visit.  The returned times 
## are according to the visit schedule, not based off ppt data.
##
## The required argument 'protocolVstFunc', must be a function whose
## first argument takes a time (or vector of times) and that returns
## a vector of all scheduled visit times thru the largest time given
## (e.g. see function 'visitSchedule505')

getVisitWeek <- function( week, visitWeeks, whichVisit=c("next","previous"))
{
  
  whichVisit <- match.arg( whichVisit )
  
  ## check for NAs
  noNAs <- !is.na( week )
  
  wk <- week[ noNAs ]
  
  ## get scheduled visit weeks through max week in 'week'
  schedVstWeeks <- visitScheduleTruncated(visitWeeks, thruWeek=max(wk))
  
  ## if schedVstWeeks has only one value, then need to append
  ## on a 2nd, so avoid errors in 'cut'
  if (length(schedVstWeeks) == 1)
    schedVstWeeks <- c(schedVstWeeks, Inf)
  
  ## find the intervals that the values of 'wk' lie in
  interval <- cut( wk, breaks=schedVstWeeks, right=TRUE, labels=FALSE)
  
  ## If we want the "previous" visit, we return the lower bound of the
  ## interval (schedVstWeeks[interval]), if we want the "next" we return
  ## the upper bound of the interval (schedVstWeeks[interval+1]  )
  if (whichVisit == "previous")
  {
    week[ noNAs ] <- schedVstWeeks[ interval ]
  } else if (whichVisit == "next")
  {
    week[ noNAs ] <- schedVstWeeks[ interval + 1 ]
  }
  
  week 
}

## Function returning a list with components: enrollRate, dropRate, infecRate
## sampled from their joint prior distribution. 
##
## Right now we're assuming independence so we sample from each of three
## prior separately

## This function has a generic interface that must be maintained.
## However the guts of it can be altered as desired as can be the form
## that is assumed for 'paramList'
##
## n = number of sets of parameter values to sample
GammaDist <- function(n, paramList)
{
  do.call( rgamma, as.list( c(n=n, paramList) ) )
}


## function for inputing ...
gammaInput <- function(alpha, beta)
{ 
  list( shape = alpha, rate = beta )
}

## This prior distribution returns the same values passed to it (i.e., it is a non-random point-mass distribution)
Constant <- function(n, paramList)
{
  rep( paramList, times=n )
}

## 'sampleRates' samples 'n' values from each of the probability distributions specified in the 'from' list
## 'from' contains two sublists: 'Function' (distributions) and 'params' (parameters of the distributions)
sampleRates <- function(n, from)
{
  Function <- from$Function # a list
  Params <- from$params     # a list
  
  out <- lapply( 1:length(Params), function(i, n) Function[[i]](n, Params[[i]]), n=n ) 
  names( out ) <- names( Params )
  
  ## rearrange (if needed)
  out <- out[ c("enrollment", "dropout", "infection") ]
  
  ## now rename
  names( out ) <- c("enrollRate", "dropRate", "infecRate")
  
  out
}

## Function to pull out from the interim data the std. summary stats 
## we need for estimating the rates
extractSummaryMeasures <- function( obsEDIobj, enrollPeriod )
{
  
  enrollTime <- ceiling( obsEDIobj$enrollTime )
  nEnroll <- length( enrollTime )
  
  ## Time during which enrollment happened in interim data 
  enrollWeeks <- min( enrollTime) : max( enrollTime )
  wksPerPeriod <- table( findInterval( enrollWeeks, enrollPeriod$start) )
  
  ## which periods are represented in 'cntsByPeriod'
  w.periods <- as.integer( names( wksPerPeriod ) )
  
  relativeRates <- enrollPeriod$relativeRates
  weighted_nWeeksEnroll <- sum( relativeRates[w.periods] * wksPerPeriod )
  
  
  ## dropout info: number of dropouts, person weeks at risk
  nDropouts <- sum( !is.na( obsEDIobj$dropTime ) )
  personWeeksAtRisk_Dropout <- sum( obsEDIobj$obsTime )
  
  
  ## infection info: number infected, person weeks at risk
  ## For infected ppts, measure time at risk through midpoint of
  ## (lastNegTestTime, infecDxTime).  For others through lastNegTestTime.
  ## 
  nInfected <- sum( !is.na( obsEDIobj$infecDxTime ) )
  
  ## average the values: (lastNegTestTime, infecDxTime)
  aveLastNeg_DxTime <-
    rowMeans( obsEDIobj[, c("lastNegTestTime","infecDxTime")], na.rm=TRUE)
  personWeeksAtRisk_Infection <- sum( aveLastNeg_DxTime )
  
  list( enrollment = list( nEnroll, weighted_nWeeksEnroll ),
        dropout = list( nDropouts, personWeeksAtRisk_Dropout ),
        infection = list( nInfected, personWeeksAtRisk_Infection )
  )
}


enrollmentEstimationFunction <- function( obsEDIobj, enrollPeriod )
{
  enrollTime <- ceiling( obsEDIobj$enrollTime )
  nEnroll <- length( enrollTime )
  
  ## Time during which enrollment happened in interim data 
  enrollWeeks <- min( enrollTime) : max( enrollTime )
  wksPerPeriod <- table( findInterval( enrollWeeks, enrollPeriod$start) )
  
  ## which periods are represented in 'cntsByPeriod'
  w.periods <- as.integer( names( wksPerPeriod ) )
  
  relativeRates <- enrollPeriod$relativeRates
  weighted_nWeeksEnroll <- sum( relativeRates[w.periods] * wksPerPeriod )
  
  ## return estimate
  (nEnroll / weighted_nWeeksEnroll)
}

dropoutEstimationFunction <- function( obsEDIobj, enrollPeriod )
{
  ## dropout info: number of dropouts, person weeks at risk
  nDropouts <- sum( !is.na( obsEDIobj$dropTime ) )
  personWeeksAtRisk_Dropout <- sum( obsEDIobj$obsTime )
  
  ## return estimate
  ( nDropouts / personWeeksAtRisk_Dropout )
}


infectionUpdateFunction <- function( obsEDIobj, enrollPeriod, hyperParams )
{
  ## infection info: number infected, person weeks at risk
  ## For infected ppts, measure time at risk through midpoint of
  ## (lastNegTestTime, infecDxTime).  For others through lastNegTestTime.
  ## 
  nInfected <- sum( !is.na( obsEDIobj$infecDxTime ) )
  
  ## average the values: (lastNegTestTime, infecDxTime)
  aveLastNeg_DxTime <-
    rowMeans( obsEDIobj[, c("lastNegTestTime","infecDxTime")], na.rm=TRUE)
  personWeeksAtRisk_Infection <- sum( aveLastNeg_DxTime )
  
  mapply( sum, hyperParams, list(nInfected, personWeeksAtRisk_Infection),
          SIMPLIFY = FALSE )
}



## Function that produces a posterior object that is constant for some
## parameters and not for others.  
##
## If a function is provided in 'estimationFunction' for a parameter,
## then that parameter will be estimated from the data using that function
## and that value will be taken as fixed for PET prediction

## If a function is provided in 'updateFunction' for a parameter,
## then that parameter will be updated from the data using that function
##
## paramEstimationFunction should be a named list, with names matching parameters
constructPartialPosteriorObject <-
  function( obsEDIobj, enrollPeriod, priorObj, 
            estimationFunction=NULL, updateFunction=NULL )
  {
    
    if ( any( names(estimationFunction) %in% names(updateFunction)) )
      stop("Parameters cannot be listed in both estimationFunction and updateFunction\n")
    
    posteriorObj <- priorObj
    
    
    ## define constant function 'Const'
    Const <- function(n, p) {
      if (length(p)==1) return( rep(p, times=n) )
      stop("Constant value must have length==1\n")
    }
    
    ## if there are parameters to be estimated, then do:
    for ( param in names( estimationFunction ) )
    {
      ## estimate value for parameter
      posteriorObj$params[[ param ]] <- 
        estimationFunction[[ param ]]( obsEDIobj, enrollPeriod )
      
      ## assign new function to be constant
      posteriorObj$Function[[ param ]] <- Const 
    }
    
    ## if there are parameters to be updated, then do:
    for ( param in names( updateFunction ) )
    {
      ## update parameter based on observed data
      posteriorObj$params[[ param ]] <- 
        updateFunction[[ param ]]( obsEDIobj, enrollPeriod, 
                                   priorObj$param[[ param ]] )
    }
    
    ## return posterior object
    posteriorObj
  }            



constructPosteriorObject <- function( obsEDIobj, enrollPeriod, priorObj )
{
  
  ## Since we're using conjugate priors, the form of the posterior
  ## is same as the prior (just different param values)
  posteriorFunc <- priorObj$Function
  
  ## extract summary data from interim observed data object 
  updateList <- extractSummaryMeasures( obsEDIobj, enrollPeriod )
  
  posteriorParams <- priorObj$params
  ## Update the priors with observed information 
  for ( comp in names(posteriorParams) ) 
  {
    posteriorParams[[comp]] <- 
      mapply( sum, posteriorParams[[comp]], updateList[[comp]],
              SIMPLIFY = FALSE )
  }
  
  list( Function = posteriorFunc, params = posteriorParams )
}

## Generate enrollment times
## --------------------------
## Sub-function used in 'getEnrollment'.  This specifies the 
## distribution from which the enrollment data come.  This is
## where the actual sampling occurs, and the function to swap out
## if you want to make different distributional assumptions.

generateEnrollment <- function(start, end, rate)
{
  ## number of enrollees during period is poisson distributed with 
  ## rate = 'rate' * (end - start + 1), and the times are uniformly
  ## distributed on the interval (start,end)
  N <- rpois(1, lambda = rate * (end - (start-1)))
  list( N = N,  times = sort( runif(N, min=start-1, max=end) ) )
}


## Generate dropout times
## ----------------------
## Specify the number of participants to generate times for (N)
## and the dropout rate (rate) 
##
## Returns continuously valued times

generateDropout <- function(N, rate)
{
  ## changing the 'rate' parameter to a form that gives the
  ## prob(event in interval (0,1)) = rate, i.e., P(dropout time <=1) = 1 - exp(-lambda) = rate
  rexp(N, rate= -log( 1 - rate) )
}


## Generate infection times  
## ------------------------
## These should be taken as the time at which the infection is
## diagnosable (by Elisa), not when infection occurs.
## 
## Specify the number of participants to generate times for (N)
## and the infection rate (rate) 
##
## Returns continuously valued times

generateInfection <- function(N, rate)
{
  ## the argument 'rate' must either be a numeric vector of length 1 
  ## or a data.frame containing columns: start, end and rate
  if ( length(rate) == 1 && !is.data.frame(rate)) {
    
    ## changing the 'rate' parameter to a form that gives the
    ## prob(event in interval (0,1)) = rate
    rexp(N, rate= -log( 1 - rate) )
  } else
  {
    if ( is.data.frame(rate) )    
    {
      ## ensure the rows of data frame are ordered properly w.r.t time
      rate <- rate[ order( rate$end ), ]
      
      ## for each row of 'rate' generate 'N' infection times.  The final
      ## value will be the minimum of these times (or NA) - see below
      ## for details
      nRows <- nrow(rate)
      offsetTimes <- c(rate$start[1]-1, rate$end[ -nRows ]) 
      rowTimes <- 
        lapply(1:nRows, function(i, N, r, offset, max) {
          
          ## generate exponential with approp. rate and add offset
          times <- offset[i] + rexp(N, -log( 1 - r[i]) )
          
          ## values greater than 'max' get set to missing
          times[ times > max[i] ] <- NA 
          
          times 
        }, N=N, r=rate$rate, offset=offsetTimes, max=rate$end)
      
      do.call(pmin, c(rowTimes, na.rm=TRUE))
      
    } else
      stop("Argument 'rate' must be a data.frame containing both rate ",
           "and time period information in order to utilize multiple ",
           "rates\n")
  }
}

## we typically supply 'rate', 'enrollPeriod' data.frame and 'nEnroll', and keep 'nWeeks' and 'maxEnroll' =NULL

getEnrollment <-
  function( rate, enrollPeriod=NULL, nEnroll=NULL, nWeeks=NULL, 
            startWeek=1, maxEnroll=NULL)
  {
    if ( is.null( enrollPeriod ) )
      enrollPeriod <- data.frame( start=1, end=NA, relativeRate=1 )
    
    nPeriods <- nrow( enrollPeriod )
    
    ## all periods have defined 'start's, but not necessarily 'end's
    startPeriod <- max( which( startWeek >=  enrollPeriod$start ) )
    
    ## define 'endWeek'
    endWeek <- ifelse( !is.null(nWeeks), startWeek + nWeeks - 1, NA)
    
    ## use 'endWeek' if it's not NA, else use max period
    endPeriod <- ifelse( !is.na(endWeek), 
                         max( which(endWeek >= enrollPeriod$start) ),
                         nPeriods)
    
    ## if a single rate was specified, repeat it so length = nPeriods
    if ( length(rate) == 1 )
      rate <- rep(rate, nPeriods)
    
    ## initialize 'N' (number enrolled) and 'times' (enrollment times)
    N <- 0
    times <- NULL
    
    ## if 'nWeeks' is specified, generate all nWeeks of data,
    ## then truncate to 'maxEnroll' if necessary
    if ( !is.null(nWeeks) )
    {
      for ( i in startPeriod:endPeriod )
      {
        ## get weekly enrollment rate for period 'i'
        rate.i <- rate[i] * enrollPeriod$relativeRate[i] 
        
        start.i <- ifelse(i > startPeriod, enrollPeriod$start[i], startWeek)
        end.i <- ifelse(i < endPeriod, enrollPeriod$end[i], endWeek)
        
        ## generateEnrollment returns a list with components 'N' and 'times'
        out <- generateEnrollment(start.i, end.i, rate.i)
        
        N <- sum(N, out$N)
        times <- c( times, out$times )
      }
      if ( !is.null(maxEnroll) && N > maxEnroll )
      {
        N <- maxEnroll
        times <- times[ 1:N ]
      }
    } else {
      ## If 'nEnroll' is specified, then generate data until goal is reached
      for ( i in startPeriod:endPeriod )
      {
        ## get weekly enrollment rate for period 'i'
        rate.i <- rate[i] * enrollPeriod$relativeRate[i] 
        
        start.i <- ifelse(i > startPeriod, enrollPeriod$start[i], startWeek)
        end.i <- ifelse(i < endPeriod, enrollPeriod$end[i], endWeek)
        
        ## If end.i is 'NA' (i.e., not specified), then we set it to a value
        ## much larger than should be needed to generate remaining enrollees
        if ( is.na(end.i) )
          end.i <- start.i + ceiling( 10 * (nEnroll - N)/rate.i )
        
        ## generateEnrollment returns a list with components 'N' and 'times'
        out <- generateEnrollment(start.i, end.i, rate.i)
        
        N <- sum(N, out$N)
        times <- c( times, out$times )
        
        if ( N >= nEnroll ) break 
        
        if ( i == endPeriod )
          stop("Probable error in generateEnrollment(): ",
               " value of 'nEnroll' not reached\n\n")
      } 
      
      ## restrict to the first nEnroll enrollees
      N <- nEnroll
      times <- times[1:N]
    }
    
    list( N = N, times = times)
  }

getDropout <- function(N, rate)
{
  generateDropout(N, rate)
}

getInfection <- function(N, baseRate, relRates=NULL, trtAssgn=NULL)
{
  #infecRateTbl <- data.frame( trt     = c("C1", "T1", "T1"),
  #                            start   = c( 1,     1    27),
  #                            end     = c(NA,    26,   NA),
  #                            relRate = c( 1,     1,   VE)
  #                          )
  
  ## if no object is provided for relRates, then 'trtAssgn' isn't
  ## needed and all data is generated using the 'baseRate' value
  if ( is.null(relRates) )
  {
    ## need to fix this
    return( generateInfection(N, baseRate) )
    
  } else {
    
    nTrts  <- length( unique(relRates$trt) )
    nRates <- length( unique(relRates$relRate) )
    
    ## compute absolute rates from base-rate and relative rates
    relRates$rate <- relRates$relRate * baseRate
    
    ## If only one Trt or only one rate, then don't need treatment info
    if ( nTrts == 1 || nRates == 1 ) {
      ## get infection times for all ppt.s at once
      n <- N
      if (nRates == 1) {
        ## just one rate
        return( generateInfection(N, unique(relRates$rate) ) )
      } else {
        ## just one treatment group (and multiple rates)
        return( generateInfection(N, relRates[,c("rate","start","end")]) )
      }
    } else {
      
      ## More than one treatment and more than one rate - need treatment
      ## assignment info for this case
      if ( is.null(trtAssgn) )
        stop("Treatment assignment information is needed for simulation",
             "of infection times.\n", "Please provide this information ",
             "via argument 'trtAssgn'\n")
      
      ## create numeric vector to store infection times in, since we aren't
      ## able to generate them all at once and simply return them 
      infecTimes <- numeric(N) 
      
      ## get counts for all the treatments
      trtCnt <- table(trtAssgn)
      
      for ( trt in names(trtCnt) ) {
        ## number of ppt.s with treatment 'trt'
        n.trt <- trtCnt[ trt ]
        
        ## identify which ppt.s have treatment 'trt'
        w.trt <- which( trtAssgn == trt )
        
        ## pull out rates just for 'trt'
        rates.trt <- relRates[ relRates$trt == trt, c("rate","start","end")]
        
        ## generate infection times for treatment 'trt' and store
        infecTimes[ w.trt ] <- generateInfection(n.trt, rates.trt)
      }                
      return( infecTimes )
    }    
  }
}


## Function to compute the Greatest Common Divisor of a set of integers
## (used by function getBlockSize() )
##
##  Arguments:  
##      nvec:  a vector of integers
##             The argument value need not be of type integer - it will be coerced
##             to integer internally
##             
##  Return Value:
##      the greatest common divisor (integer)
##
gcd <- function(nvec) {

    nComp <- length(nvec)
    ## Ensure that nvec has at least 2 components
    if ( nComp < 2 )
        stop("Error(gcd): Argument 'nvec' must have length >= 2\n")

    ## ensure that the input vector contains only integers
    if ( any( abs(round(nvec) - nvec) > sqrt(.Machine$double.eps)) )
        stop("Error(gcd): Non-integer values found in input vector\n")

    ## coerce to integer (to avoid dealing with numerical fuzz in comparisons
    nvec <- as.integer( round(nvec) )

    ## Defines a recursive function that computes the gcd of a pair of integers
    ## (taken off an R-mailing list post).
    gcd2 <- function(a,b) ifelse(b==0, a, gcd2(b, a %% b) )

    ## get pairwise gcd between current gcd and the next component
    for (i in 1:(nComp-1)) {
        nvec[i+1] <- gcd2(nvec[i], nvec[i+1])
    }
    nvec[ nComp ]
}


## getBlockSize:  
## A function that computes possible block sizes to use in a blocked randomization.
## It either reports the minimum possible block size (any multiple of it also is
## a usable value), or the smallest block size that falls within a given range
## (e.g. [10, 30] ) if the user provides a value via the 'range' argument.
## If 'range' is specified and there is no block size within that range, the function
## returns NA.
##
##  Arguments:  
##      nvec:  a vector of counts of ppt.s allocated to various treatment arms
##     range:  a length 2 integer vector of form [min, max] that gives the 
##             upper and lower bounds of a region of values to consider for the
##             blockSize of a randomization (the bounds are considered part of the
##             region).  If 'range' is specified then the function returns the
##             smallest multiple of the minimum block size that lies in the region
##             (if any) - or NA if none do.
##             
getBlockSize <- function(nvec, range=c(0,Inf)) {

    minBlk <- sum( round(nvec) / gcd(nvec) )

    ## if 'range' was specified by user do some checking on the value
    if ( !missing(range) ) {
      if (length(range) != 2)
          stop("Error(getBlockSize): The 'range' parameter must be a vector",
               " of length 2\n")

      if (is.infinite(range[1]))
          stop("Error(getBlockSize): The first value of the 'range' vector is",
               " infinite - this is not allowed\n")

      if ( range[2] < range[1] )
          stop("Error(getBlockSize): The first value of the 'range' vector must",
               " be less-than-or-equal-to the second value\n")

      ## If the range minimum was set to less than 1, reset to 1.  The blockSize
      ## must be an integer and a value below 1 is nonsensical
      if (range[1] < 1 )  range[1] <- 1

      ## Minimum multiple of minBlk needed to meet/exceed the minimum range value
      minMult <- ceiling(range[1]/minBlk)
      ## Maximum multiple of minBlk usable to meet/fall-below the maximum range value
      maxMult <-   floor(range[2]/minBlk)

      ## if the minimum is larger than the maximum, then no value lies in the range
      if ( minMult > maxMult )
        return (NA)
      else
        return ( minMult * minBlk )
    } else {
        ## if range was not specified
        return (minBlk)
    }
}


## Treatment assignment ties in with infection, since infections rates
## will vary across treatments - unless all are equally (in)effective.
getTreatmentAssignment <- 
  function(n, prob=NULL, nPerTrt=NULL, blockSize=NULL, seed=NULL, approxOK=FALSE)
  {
    ## Arguments:  
    ##-----------
    ##   n     (integer) total number of assignments to generate
    ##
    ##   prob  (numeric vector) A (possibly named) vector of assignment
    ##         probabilities for all treatments.  If the elements are named,
    ##         the names are used in the returned object (see details in
    ##         "Return Object" below)
    ##
    ##   nPerTrt  Controls how many ppts are assigned to each trt group.
    ##            Valid values are:
    ##              (character) "fixed" - 
    ##                  Assigns round(n * prob[i]) ppts to the i-th treatment
    ##                  (adjusted if sum( round(n * prob) ) != n ). 
    ##              (character) "random" -
    ##                  Number assigned to each trt is random - distributed as
    ##                  multinomial with probabilities given by 'prob'
    ##              (integer vector)
    ##                  Direct specification of the number ppt.s for each trt
    ##
    ##   blockSize  (integer) Specifies that a blocked randomization approach 
    ##              to treatment assignment should be used, with block size
    ##              size as given by this argument. 
    ##
    ##   seed  (integer) <Optional> Random number seed that can be specified 
    ##         in order to produce reproducible results 
    ##
    ##   approxOK (logical) when nPerTrt=="fixed" is it okay for (prob * n) to
    ##                      not be integers (i.e. is it okay for the number of
    ##                      of ppts in each trt to approximately match 'prob'
    ##                      (rather than exactly)?
    ##
    ## Return Object:
    ## --------------
    ##   A factor vector of length 'n' containing trt assignments. 
    ##   If the argument 'prob' is a named vector then those names are used as
    ##   the levels of the factor that is returned.  Otherwise some default names
    ##   are assigned by attaching a prefix to the treatment's integer code.
    ##   Treatment codes are 1, 2, ..., length(prob), with integer i representing
    ##   treatment with assignment probability of prob[i]
    
    ##
    if ( !is.null(blockSize) ) 
      useBlock <- TRUE
    
    if ( !is.null(prob) ) {
      useProb <- TRUE
      if ( !is.numeric(prob) || any(prob<0 | prob>1) || 
             abs( sum(prob) - 1 ) > .Machine$double.eps ^ 0.5 )
        stop("Argument 'prob' must be a vector of probabilites summing to 1.\n")
    }
    
    if ( !is.null(nPerTrt) ) {
      useNPT <- TRUE
      if ( is.character(nPerTrt) && is.null( prob ) ) 
        stop("Argument 'prob' must be specified too, when 'nPerTrt' is ",
             "set to one of: ('fixed', 'random')\n") 
    }
    
    ## set seed if user provided one
    if ( !is.null(seed) && is.numeric(seed) ) 
      set.seed(seed)
    
    ## if block randomization chosen, do this:
    if ( useBlock ) {
      
      ## the argument 'prob' takes precedence if both it and nPerTrt were given
      if ( useProb ) {  

        ## number of each trt group in each block
        nPerBlock <- round( blockSize * prob )

        ## Verify that the 'blockSize' specified is compatible with the treatment
        ## prob.s given by 'prob'
        if ( sum( nPerBlock ) != blockSize ) {
          stop("\n  The given value of blockSize, ", blockSize,
               ", is not compatible with the specified treatment\n",
               "distribution.  Please use function 'getBlockSize' ",
               "to determine a suitable\n block size (see ",
               "help(getBlockSize) ), and then re-run\n\n")
        }

        ## We allow user to use a poor block size if they really want to (may be
        ## necessary in odd cases?), but issue a warning about it
        if ( max( abs( nPerBlock/sum(nPerBlock) - prob ) > sqrt(.Machine$double.eps) )) {
           warning(
             "\n  The given value of blockSize, ", blockSize,
             ", is not compatible with the value\n",
             "of argument 'prob'.  Blocks of the specified size cannot create a\n",
             "treatment distribution that matches 'probs'.  We suggest using\n",
             "function 'getBlockSize' to determine a suitable block size and then\n",
             "re-running (see help(getBlockSize))\n\n", immediate=TRUE)
        } 
      } else {
        if ( useNPT ) {  
          if ( sum(useNPT) != blockSize )
            stop("The values of 'nPerTrt' must sum to the given blockSize \n")
          nPerBlock <- nPerTrt
        } else { 
          stop("Must specify either 'prob' or 'nPerTrt'\n")
        }
      }
      
      ## do Block randomization and return
      nBlocks <- ceiling( n / blockSize)
      
      nTrt <- length(nPerBlock)
      
      ## vector of length 'blockSize' containing the correct number of each trt in it (Doug's version: 'times=nPerBlock')
      trtBlock <- rep(1:nTrt, times=nPerBlock)
      
      ## generates a vector of randomized blocks
      trtVec <- unlist( lapply(1:nBlocks, function(i, trtBlk, bS)
        trtBlk[ order( runif(bS) ) ],
                               trtBlk=trtBlock, bS=blockSize) )
      
      ## keep the first 'n' of them (in case there are more than 'n')
      if ( n < length(trtVec) )
        trtVec <- trtVec[1:n]
      
    } else {   
      
      ## block randomization NOT chosen:
      if ( useProb ) {
        ## 'prob' given and 'nPerTrt' one of "fixed" or "random"
        if ( !useNPT || !nPerTrt %in% c("fixed","random") ) 
          stop("Either 'fixed' or 'random' must be specified for 'nPerTrt' ",
               "when 'prob' is given and blockSize is not.\n") 
        
        nTrt <- length(prob)
        
        if ( nPerTrt == "random" ) {
          trtVec <- sample.int(n = nTrt, size = n, prob = prob, replace=TRUE)
        } else { ## nPerTrt == "fixed"
          np <- round( prob * n ) 
          if ( !approxOk && !all.equal( np, prob * n) )
            stop("The given values of 'n' and 'prob' don't produce an integer",
                 "number of ppt.s for each treatment\n")
          if ( sum(np) != n )
            stop("The sum of round('n'*'prob') must equal'n' \n")
          
          ## create a vector with the correct number of each trt in it, and
          ## then randomize it
          trtVec <- rep.int(1:nTrt, times=np)[ order(runif(n)) ]
        }
      } else {
        ## 'prob' not given, so 'nPerTrt' must be an integer vector
        if ( !all.equal( nPerTrt, as.integer(nPerTrt)) || sum(nPerTrt) != n )
          stop("Value of 'nPerTrt' is not integer and/or does not sum to 'n' \n") 
        nTrt <- length(prob)
        trtVec <- rep.int(1:nTrt, times=nPerTrt)[ order(runif(n)) ]
      }
    }
    
    ## create names for 'prob' if they weren't supplied by user
    if ( !useProb || is.null( names(prob) ) ) {
      prefix <- "Trt"
      trtNames <- paste(prefix, 1:nTrt, sep="")
    } else {
      trtNames <- names(prob)
    }
    
    ## Last step, convert to factor
    trtVec <- factor(trtVec, labels=trtNames)
    
    trtVec
  }

## Usage:  
## The function is set up so that the user must define 
##   *either* 'nWeeks' (number of weeks of data to simulate)
##   *or* 'nEnroll' (number of partipants to enroll).
##
## The argument 'maxEnrollment' can be used along with 'nWeeks' if
## desired, to put a cap on enrollment.
##
## The function uses three "get" functions (in 'getEDIfuncs.R')
## as well as "visit schedule" functions (in 'visitSchedule.R'),




## For simulating a certain amount of data from the beginning of a trial
simulateObservedEDIdata <-
  function(rateParams, trtAssignProb, 
           infecRates, protocolVisitFunc,
           enrollPeriod, startWeek=1, 
           nEnroll=NULL, nWeeks=NULL,
           maxEnrollment=NULL)
  {
    
    ## very basic sanity check on args 'nEnroll' and 'nWeeks' ##
    if ( is.null(nEnroll) && is.null(nWeeks) ) 
      stop("You must specify either nEnroll or nWeeks\n") 
    
    if ( !is.null(nEnroll) && !is.null(nWeeks) ) 
      stop("You must specify only *one* of nEnroll or nWeeks\n")
    
    ## ---------------- LET THE GAMES BEGIN ----------------
    
    ## get enrollment count and times
    enr <- getEnrollment(
      rate = rateParams$enrollRate, 
      enrollPeriod = enrollPeriod,
      nEnroll = nEnroll,
      nWeeks = nWeeks, 
      startWeek = startWeek,
      maxEnroll= maxEnrollment )
    
    ## extract enrollment count
    N <- enr$N
    
    ## generate dropout and data for the 'N' participants
    dropout <- getDropout(N, rateParams$dropRate)
    
    ## generate treatment assignments
    nTrt <- length(trtAssignProb)
    trtAssignments <- getTreatmentAssignment(N, prob = trtAssignProb, 
                                             blockSize = 10*nTrt )
    
    ## generate infection data for the 'N' participants
    infect <- getInfection(N, baseRate=rateParams$infecRate, 
                           relRates = infecRates, 
                           trtAssgn = trtAssignments)
    
    
    ## If 'nWeeks' was specified, then we've established a bound on our
    ## follow-up time for participants. So we need to censor the dropout
    ## and infection times using that info.
    if ( !is.null( nWeeks ) )
    {
      ## amount of trial time simulated within this invocation of simulateEDI()
      simTime <- (startWeek - 1 + nWeeks) - enr$times
      
      ## censor dropout and infection times by 'simTime'
      dropout[ is.TRUE(dropout > simTime) ] <- NA
      infect[ is.TRUE(infect > simTime) ] <- NA
      
    } else
    {
      simTime <- NA
    }
    
    
    ## Censor infections that occur after dropout, as they can't be observed
    infect[ is.TRUE( infect > dropout ) ] <- NA
    
    
    ## *** NOTE  ****************************************************
    ##
    ##  We still have records with both dropout and infection times
    ##  recorded, and those all have dropout > infect.
    ##
    ##  We need to keep both times until we calculate the "DX time"
    ##  for the infections, then we can compare dropout to the DX 
    ##  time and keep the earlier of the two.
    ## **************************************************************
    
    
    ## ------------------------------------------------------------ ##
    ## Now, we start to construct the "observed" data, by applying  ##
    ## the study visit Map to the simulated data.  This will allow  ##
    ## us to figure out when: (a) the infections are diagnosed,     ##
    ## (b) the last visit occured (during the simulated period).    ##
    ## ------------------------------------------------------------ ##
    
    ## create observed data object to fill in
    obsEDI <- data.frame( trt = trtAssignments,
                          enrollTime = enr$times,
                          dropTime  = dropout,
                          infecDxTime = NA,
                          lastNegTestTime = NA,
                          obsTime = NA )
    
    ## get infection diagnosis dates for each non-NA infection time
    nonNA.inf <- !is.na( infect )
    
    ## Do the following only if there's at least one infection
    ## (NOTE - to reduce complexity, this section has redundant computation)
    if ( !all( is.na(infect) ) )
    {
      infecDX <- getVisitWeek( infect, protocolVisitFunc, whichVisit = "next")
      
      ## Compute the minimum of the infecDX time, dropout time and simTime
      ## Any events occuring *strictly* after this time are censored
      minTime <- pmin( infecDX, dropout, simTime, na.rm=TRUE)
      
      ## censor infecDX and dropout times
      infecDX[ is.TRUE(infecDX > minTime) ] <- NA
      dropout[ is.TRUE(dropout > minTime) ] <- NA
      
      ## store info in obsEDI
      obsEDI$infecDxTime <- infecDX
      obsEDI$dropTime <- dropout
    }
    
    
    ## Now that we're done making changes to 'dropout' and 'infecDxTime',
    ## we compute the amount of "observation time" ('obsTime') for each
    ## ppt.  This is equal to the 'simTime' for ppt.s without events,
    ## and equal to the event time for participants with events.
    obsEDI$obsTime <- pmin( obsEDI$infecDxTime, 
                            obsEDI$dropTime, simTime, na.rm=TRUE)
    
    
    ## Fill in the time of the last Negative HIV test.  This will be the
    ## last visit prior to 'obsTime'.  Why?  Because 'obsTime' equals
    ## infecDxTime for infecteds, dropTime for dropouts, and simTime for
    ## everyone else, and this is what we want.
    obsEDI$lastNegTestTime <- getVisitWeek( obsEDI$obsTime, 
                                            protocolVisitFunc, "previous")
    
    ## return obsEDI
    obsEDI
  }


## Usage:  
##
## For enrollment, the user must define either:
##   'nWeeksEnroll' (number of weeks of enrollment data to simulate)
##   *or* 'nEnroll' (number of partipants for whom to simulate 
##                   enrollment data).
##
## The argument 'maxEnrollment' can be used along with 'nWeeksEnroll'
## (if desired), to put a cap on enrollment.
##
## The two arguments 'nWeeksFU' and 'nWeeksTrialTime' each specify the 
## amount of follow-up to be performed: 
##   'nWeeksFU' specifies the amount of follow-up for each ppt  
##   'nWeeksTrialTime' specifies total trial duration - each ppt's
##                     follow-up time will vary based on enrollment time
## Only one of these two arguments can be specified.
##
## The function uses three "get" functions (in 'getEDIfuncs.R')
## as well as "visit schedule" functions (in 'visitSchedule.R'),



## For simulating a certain amount of data from the beginning of a trial
simFullEDIdata <-
  function(rateParams, trtAssignProb, blockSize, infecRates, protocolVisits,
           enrollPeriod, startWeek=1, 
           nEnroll=NULL, nWeeksEnroll=NULL, maxEnrollment=NULL,
           nWeeksFU=NULL, nWeeksTrialTime=NULL)
  {
    ## very basic sanity check on arguments 'nEnroll' and 'nWeeksEnroll'
    if ( is.null(nEnroll) && is.null(nWeeksEnroll) ) 
      stop("You must specify either nEnroll or nWeeksEnroll\n") 
    
    if ( !is.null(nEnroll) && !is.null(nWeeksEnroll) ) 
      stop("You must specify only *one* of nEnroll or nWeeksEnroll\n")
    
    if ( !is.null(nWeeksFU) && !is.null(nWeeksTrialTime) ) 
      stop("You must specify only *one* of nWeeksFU or nWeeksTrialTime\n")
    
    if (!is.null(nWeeksEnroll) && !is.null(nWeeksTrialTime) && nWeeksEnroll > nWeeksTrialTime) 
      stop("You have requested ", nWeeksEnroll," weeks of enrollment (via argument",
           " nWeeksEnroll), but\n", "allowed only ", nWeeksTrialTime,
           " weeks of trial time (via argument nWeeksTrialTime).\n")
    
    ## get enrollment count and times
    enr <- getEnrollment(
      rate = rateParams$enrollRate,
      enrollPeriod = enrollPeriod,
      nEnroll = nEnroll,
      nWeeks = nWeeksEnroll,
      startWeek = startWeek,
      maxEnroll= maxEnrollment )
    
    ## do a check immediately on enrollment times *if* the user
    ## specified 'nWeeksTrialTime'
    if ( !is.null( nWeeksTrialTime ) && !is.null(nEnroll) &&
           ( lastWeek <- ceiling(max(enr$times)) ) > nWeeksTrialTime ) {
      
      on.exit( cat("NOTE:  Enrollment of the", nEnroll, "participants took",
                   lastWeek, "weeks time, but only", nWeeksTrialTime, "weeks",
                   "of trial time were allowed by the user (via",
                   "argument 'nWeeksTrialTime').", 
                   "The number of enrollees is being reduced to enforce",
                   "the trial time limit specified.",
                   fill = 70, labels= "simFullEDIdata(): ") )
      
      
      ## Modify 'enr' to enforce 'nWeeksTrialTime'
      enr$times <- enr$times[ ceiling( enr$times ) <= nWeeksTrialTime ]
      enr$N <- length( enr$times )
    }
    
    ## extract enrollment count
    N <- enr$N
    
    ## generate dropout times for the 'N' participants from an exponential distribution
    dropout <- getDropout(N, rateParams$dropRate)
    
    ## generate treatment assignments
    nTrt <- length(trtAssignProb)
    
    trtAssignments <- getTreatmentAssignment(N, prob = trtAssignProb, 
                                             blockSize = blockSize )
    
    ## generate infection data for the 'N' participants
    infect <- getInfection(N, baseRate=rateParams$infecRate, 
                           relRates = infecRates, 
                           trtAssgn = trtAssignments)
    
    ## If 'nWeeksFU' or 'nWeeksTrialTime was specified, then we've got an 
    ## established bound on follow-up time.  We need to censor the dropout
    ## and infection times using that info.
    if ( !is.null( nWeeksFU ) || !is.null( nWeeksTrialTime ) )
    {
      ## amount of trial time simulated within this invocation of simulateEDI()
      if ( !is.null(nWeeksFU) ) {
        fuTime <- nWeeksFU 
      } else {
        ## else nWeeksTrialTime must have been specified
        fuTime <- (startWeek - 1 + nWeeksTrialTime) - enr$times
      }
      
      ## censor dropout and infection times by 'fuTime'
      dropout[ is.TRUE(dropout > fuTime) ] <- NA
      infect[ is.TRUE(infect > fuTime) ] <- NA
      
    } else
    {
      fuTime <- NA
    }
    
    
    ## Censor infections that occur after dropout, as they can't be observed
    infect[ is.TRUE( infect > dropout ) ] <- NA
    
    
    ## *** NOTE  ****************************************************
    ##
    ##  We still have records with both dropout and infection times
    ##  recorded, and those all have dropout > infect.
    ##
    ##  We need to keep both times until we calculate the "DX time"
    ##  for the infections, then we can compare dropout to the DX 
    ##  time and keep the earlier of the two.
    ## **************************************************************
    
    
    ## ------------------------------------------------------------ ##
    ## Now, we start to construct the "observed" data, by applying  ##
    ## the study visit Map to the simulated data.  This will allow  ##
    ## us to figure out when: (a) the infections are diagnosed,     ##
    ## (b) the last visit occurred (during the simulated period).    ##
    ## ------------------------------------------------------------ ##
    
    ## create observed data object to fill in
    obsEDI <- data.frame( trt = trtAssignments,
                          enrollTime = enr$times,
                          dropTime = dropout,
                          infecDxTime = NA,
                          lastNegTestTime = NA,
                          futime = NA )
    
    ## get infection diagnosis dates for each non-NA infection time
    nonNA.inf <- !is.na( infect )
    
    ## Do the following only if there's at least one infection
    ## (NOTE - to reduce complexity, this section has redundant computation)
    if ( !all( is.na(infect) ) )
    {
      ## protocolVisitFunc = RSA_vstSch_3mo, which returns a vector of scheduled visit weeks assuming 3-monthly
      ## testing after the last vaccination visit
      infecDX <- getVisitWeek( infect, protocolVisits, whichVisit = "next")
      
      ## Compute the minimum of the infecDX time, dropout time and fuTime
      ## Any events occurring *strictly* after this time are censored
      minTime <- pmin( infecDX, dropout, fuTime, na.rm=TRUE)
      
      ## censor infecDX and dropout times
      infecDX[ is.TRUE(infecDX > minTime) ] <- NA
      dropout[ is.TRUE(dropout > minTime) ] <- NA
      
      ## store info in obsEDI
      obsEDI$infecDxTime <- infecDX
      obsEDI$dropTime <- dropout
    }
    
    
    ## Now that we're done making changes to 'dropout' and 'infecDxTime',
    ## we compute the amount of "follow-up time" ('futime') for each
    ## ppt.  This is equal to the 'fuTime' for ppt.s without events,
    ## and equal to the event time for participants with events.
    obsEDI$futime <- pmin( obsEDI$infecDxTime, 
                           obsEDI$dropTime, fuTime, na.rm=TRUE)
    
    
    ## Fill in the time of the last Negative HIV test.  This will be the
    ## last visit prior to 'obsTime'.  Why?  Because 'obsTime' equals
    ## infecDxTime for infecteds, dropTime for dropouts, and fuTime for
    ## everyone else, and this is what we want.
    obsEDI$lastNegTestTime <- getVisitWeek( obsEDI$futime, 
                                            protocolVisits, "previous")
    
    ## return obsEDI
    obsEDI
  }

summInterimData <- function( obsEDI )
{
  ## enrollment info: *weighted* number of enrollees, weeks enrollment
  enrollTime <- ceiling( obsEDI$enrollTime )
  nWeeksEnroll <- max( enrollTime )
  
  ## dropout info: number of dropouts, person weeks at risk
  nDropouts <- sum( !is.na( obsEDI$dropTime ) )
  personWeeksAtRisk_Dropout <- sum( obsEDI$obsTime )
  
  ## infection info: number infected, person weeks at risk
  ## For infected ppts, measure time at risk through midpoint of
  ## (lastNegTestTime, infecDxTime).  For others through lastNegTestTime.
  ## 
  nInfected <- sum( !is.na( obsEDI$infecDxTime ) )
  
  ## average the values: (lastNegTestTime, infecDxTime)
  aveLastNeg_DxTime <- 
    rowMeans( obsEDI[, c("lastNegTestTime","infecDxTime")], na.rm=TRUE)
  personWeeksAtRisk_Infection <- sum( aveLastNeg_DxTime )
  
  cat("nWeeksEnroll =", nWeeksEnroll,
      "nDropouts =", nDropouts,
      "nInfected =", nInfected, "\n")
}

applyStopRules <- function(d, infectionTotals, boundLabel="highEff", lowerHRnoneff=NULL, upperHRnoneff=NULL, stage1HR=NULL, upperHRuncPower=NULL, highHR=NULL, 
                          alphaNoneff=0.05, alphaStage1=0.05, alphaUncPower=0.05, alphaHigh=0.05, post6moCut=26, post6moMonitor, estimand) {
  
  ## This function apply the stopping rules for nonefficacy and high efficacy monitoring.
  ## A vector of events where nonefficacy monitoring
  ## will occur is inputed by argument 'infectionTotals', for example (50, 70, 90, 110, 130, 150).
  ## If estimand="combined", the nonefficacy bound is reached if both 95% CIs for HR and cumulative incidence ratio (FR) 
  ## (1) lie above alternative HR ('HaHR'); (2) not below 'NullHR' (usually 1).
  ##
  ## Argument 'd' should be a data.frame containing (at least) columns:
  ##   'entry' - the entry time (in trial time)
  ##   'exit'  - time of event, trial completion or dropout (trial time)
  ##   'event' - 0/1 indicator of event of interest
  ##   'trt'   - 0/1 indicator of vaccine receipt
  ## 
  ## Other arguments:
  ## ---------------
  ## 'infectionTotals' - a vector specifying the total number of infections
  ##                     at which analyses should take place.  The contents 
  ##                     of this vector must be ordered (smallest to largest)
  ## 
  ## 'boundLabel'  - can only take on values "HighEff" or "NonEff"
  ##
  ## 'post6moCut'  - cut off time (in weeks)
  ## 'estimand'    - a character string specifying the estimand of interest (can be one of "combined", "cox", and "cuminc")
  ##
  ## check contents of 'd'
  if ( !all( c("entry","exit","event","trt") %in% names(d) ) )
    stop("DataFrame 'd' must contain columns: entry, exit, event, and trt\n")
  
  ## check specification of bounds
  Lower <- ifelse( boundLabel=="NonEff", FALSE, TRUE )
  
  if (Lower){ # high efficacy monitoring
    if (any(is.null(highHR),is.null(alphaHigh))){ stop("The arguments 'highHR' and 'alphaHigh' must be specified for high efficacy monitoring.\n") }
  } else {    # non-efficacy monitoring
    if (any(is.null(lowerHRnoneff),is.null(upperHRnoneff),is.null(alphaNoneff))){ stop("The arguments 'lowerHRnoneff', 'upperHRnoneff', and 'alphaNoneff' must be specified for non-efficacy monitoring.\n") }
  }
  
  if (estimand=="combined" & post6moMonitor){ stop("Post-6 months monitorting for non-efficacy is not implemented for the 'combined' approach.") }
  
  ## store the length of the infections totals
  L <- length(infectionTotals)
  
  ## sanity check: make sure that 'infectionTotals', 'lowerBounds' are all 
  ## ordered correctly, and no dups in 'infectionTotals'
  ## (infectionTotals increasing)
  if (L > 1) {
    if ( any( diff(infectionTotals) <= 0 ) ){
      stop("Argument 'infectionTotals' must be ordered by increasing",
           "magnitude.\n", "Please fix and retry.\n")
    }     
  }
  
  ## Make sure we have only two trt groups and that they're coded as 0 and 1
  uniq.trt <- sort( unique(d$trt) )
  if ( length(uniq.trt)>2 ) {
    warning("The data set given to 'applyStopRules' contains", length(uniq.trt),
            "treatments - it should only contain 2.\n") 
  } else if ( any( uniq.trt != c(0,1) ) ) { 
    warning("The data set given to 'applyStopRules' contains values of 'trt'",
            "other than 0 and 1, this is probably an error.\n",
            "The values are: ", uniq.trt[1], " and ", uniq.trt[2], "\n\n")
  }
  
  options( stringsAsFactors = FALSE)
  
  firstEnrollTime <- min(d$entry, na.rm=TRUE)
  
  ## select out events
  eventDF <- subset(d, event == 1, select=c(exit,trt))
  
  ## order events by time of occurrence 
  eventDF <- eventDF[ order(eventDF$exit), ,drop=FALSE]
  
  ## add column indicating the total number of infections accrued at each event time
  n <- nrow(eventDF)
  eventDF$totInfec <- 1:n
  
  ## add a final test at n for non efficacy test
  #if (!Lower) {
  #   if (max(infectionTotals)<n)
  #     infectionTotals = c(infectionTotals, n)
  #}
  ## get the number of tests specified
  nTests <- length(infectionTotals) # duplicate information: 'L' defined above the same
    
  ## extract 'exit' times corresponding to the infection totals given
  ## in 'infectionTotals' 
  testTimes <- subset(eventDF, totInfec %in% infectionTotals)$exit
  
  ## count how many times we have in 'testTimes' (we may not have accrued enough
  ## infections to perform all the tests specified via the argument 'infectionTotals'
  K <- length( testTimes )
  
  ## if we didn't reach any of the infection totals at which we planned to test
  ## then exit out, returning basic info
  if (K == 0) {
    return( list( finished = FALSE,
                  boundHit = NA, 
                  altDetected = NA,
                  stopTime = NA, 
                  stopInfectCnt = NA,
                  totInfecCnt = n,
                  totInfecSplit= table(eventDF$trt),
                  lastExitTime = max(d$exit) - firstEnrollTime) )   ## stage 1 exit time
  }
  
  ## order the full data by entry times 
  d <- d[ order(d$entry), ]
  
  
  ## create an object to be filled in during our testing loop below
  ## and which will be output by the function.  Contains summary info
  ## about the status at each test.  
  summObj <- data.frame(test = (1:K), 
                        testTime = testTimes,
                        ppts = integer(K),
                        FU= numeric(K),
                        infectTotal = infectionTotals[1:K],
                        infectSplit = character(K), 
                        
                        pptsPost6mo = integer(K),
                        FUPost6mo= numeric(K),
                        infectPost6mo = integer(K),
                        infectSplitPost6mo = character(K),
                        infectPctPost6mo = numeric(K),
                        
                        EstCumulatIncid = numeric(K),
                        EstHazardRatio = numeric(K))
  
  ## create the model ('coxFormula') to be fit by coxph()
  coxFormula <- Surv(futime.i, event.i) ~ trt
  
  ## indicator of whether trial is done (stopping bound hit)
  done <- FALSE
  
  ## initialize var 'bound', which will indicate which bound was hit
  bound <- NA
  
  nVaxInfStage1 <- NA
  altDetected <- NA
  
  for ( i in 1:K ) {
    
    t.i <- testTimes[i]
    
    ## restrict data to what would be observed at time 't.i'
    D <- subset(d, entry < t.i )  
    D$event.i <- D$event & (D$exit <= t.i) # a logical vector
    D$futime.i <- pmin(D$exit, t.i) - D$entry 
    
    # indicator of a post-6 months infection
    D$event.i6 <- D$event.i & (D$futime.i>post6moCut)
    
    ## populate summary object with basic info
    #summObj$EstHazardRatio[i] <- hr.i
    
    summObj$ppts[i] <- nrow(D)
    summObj$pptsPost6mo[i] <- sum( D$futime.i > post6moCut, na.rm=TRUE )
    
    summObj$FU[i] <- sum( D$futime.i )
    summObj$FUPost6mo[i] <- sum( pmax(D$futime.i - post6moCut, 0), na.rm=TRUE )
    
    infectTbl <- table( D$trt[ D$event.i ] )
    
    ## record infection "splits" (number of placebo infections: number of vacc infections)
    if ( length(infectTbl)== 2) { 
      infectionsInOnlyOneGroup <- FALSE
      summObj$infectSplit[i] <- paste("Pl:Vx =", paste(infectTbl,collapse=":") )
    } else {
      infectionsInOnlyOneGroup <- TRUE
      tbl <- c(0, 0)
      tbl[ 1 + as.numeric( names(infectTbl ) ) ] <- as.numeric(infectTbl)
      summObj$infectSplit[i] <- paste("Pl:Vx =", paste(tbl,collapse=":") )
    }
    
    summObj$infectPost6mo[i] <- sum( D$event.i & D$futime.i>post6moCut)
    summObj$infectPctPost6mo[i] <- 100 * ( summObj$infectPost6mo[i] / summObj$infectTotal[i] )
    
    ## record infection "splits" for post6month infections
    infectPost6moTbl <- table( D$trt[ D$event.i & D$futime.i>post6moCut] )
    if (length(infectPost6moTbl) == 2) {
      summObj$infectSplitPost6mo[i] <- paste("Pl:Vx =", paste(infectPost6moTbl,collapse=":") )
    } else {
      tbl <- c(0, 0)
      if (length(infectPost6moTbl) == 1)
        tbl[ 1 + as.numeric( names(infectPost6moTbl) ) ] <- as.numeric(infectPost6moTbl)
      summObj$infectSplitPost6mo[i] <- paste("Pl:Vx =", paste(tbl,collapse=":") )
    }    
    
    ## We need infections in both trt groups to create (Wald) confidence intervals for
    ## the hazard ratio and for the other ratio being used. So if we encounter a case
    ## where there are zero infections in one group, we will not try to get estimates
    ## but will skip ahead to the next analysis. 
    ##
    ## This approach is fine because the group with zero infections has to be the 
    ## vaccine group due to the harm monitoring we have in place (if there were zero
    ## placeo infections we'd have hit harm bound before getting to the noneff 
    ## monitoring).  And if there are 0 infections in the vaccine group we will never
    ## stop for non-efficacy, as that's the best-case scenario for a vaccine.
    
    if ( !Lower && infectionsInOnlyOneGroup ) {
      
      ## check on the code (make sure infections are in placebo group)
      if ( unique( D$trt[D$event.i] ) == 1 )
        stop("All infections are in the vaccine group.\n",
             "Either harm monitoring is not being done or else it is",
             "broken\n")
      
      ## fill in point estimates in 'summObj' and go to next time point
      summObj$EstHazardRatio[i]  <- 0
      summObj$EstCumulatIncid[i] <- 0
      
      next
    }
    
    
    ## run cox model
    coxPH.i <- coxph( coxFormula, data=D )
    if (post6moMonitor){ coxPH.i6 <- coxph(Surv(futime.i,event.i6) ~ trt, data=D) }
    
    ## extract the hazard ratio for 'trt' and store         
    hr.i <- exp( coxPH.i$coef )
    summObj$EstHazardRatio[i] <- hr.i
    
    if (estimand %in% c("combined", "cuminc")){
      ## calculate cumulative incidence (Fv/Fp) and its CI 
      KM <- survfit(coxFormula, data=D, error="greenwood")
      KM.sum <- summary(KM)
            
      # Nelson-Aalen estimates
      na.0 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=0"]/KM.sum$n.risk[KM.sum$strata=="trt=0"])
      varna.0 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=0"]/KM.sum$n.risk[KM.sum$strata=="trt=0"]^2)
      na.0 <- na.0[length(na.0)]
      varna.0 <- varna.0[length(varna.0)]      
      
      na.1 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=1"]/KM.sum$n.risk[KM.sum$strata=="trt=1"])
      varna.1 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=1"]/KM.sum$n.risk[KM.sum$strata=="trt=1"]^2)
      na.1 <- na.1[length(na.1)]
      varna.1 <- varna.1[length(varna.1)]
            
      # survival estimates
      S.0 <- exp(-na.0)
      varS.0 <- ifelse(is.na(varna.0), NA, exp(-2*na.0) * varna.0)
      S.1 <- exp(-na.1)
      varS.1 <- ifelse(is.na(varna.1), NA, exp(-2*na.1) * varna.1)
      
      # cumulative incidence ratio
      F.0 <- 1 - S.0
      F.1 <- 1 - S.1
      FR.i <- F.1/F.0
      varlogFR <- ifelse(is.na(varS.0) | is.na(varS.1), NA, varS.1/(F.1^2) + varS.0/(F.0^2))
      
      summObj$EstCumulatIncid[i] <- FR.i
      
      if (post6moMonitor){
        ## calculate cumulative incidence (Fv/Fp) and its CI 
        KM <- survfit(Surv(futime.i,event.i6) ~ trt, data=D, error="greenwood")
        KM.sum <- summary(KM)
                
        # Nelson-Aalen estimates
        na.0 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=0"]/KM.sum$n.risk[KM.sum$strata=="trt=0"])
        varna.0 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=0"]/KM.sum$n.risk[KM.sum$strata=="trt=0"]^2)
        na.0 <- na.0[length(na.0)]
        varna.0 <- varna.0[length(varna.0)]      
        
        na.1 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=1"]/KM.sum$n.risk[KM.sum$strata=="trt=1"])
        varna.1 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=1"]/KM.sum$n.risk[KM.sum$strata=="trt=1"]^2)
        na.1 <- na.1[length(na.1)]
        varna.1 <- varna.1[length(varna.1)]
        
        # survival estimates
        S.0 <- exp(-na.0)
        varS.0 <- ifelse(is.na(varna.0), NA, exp(-2*na.0) * varna.0)
        S.1 <- exp(-na.1)
        varS.1 <- ifelse(is.na(varna.1), NA, exp(-2*na.1) * varna.1)
        
        # cumulative incidence ratio
        F.0 <- 1 - S.0
        F.1 <- 1 - S.1
        FR.i6 <- F.1/F.0
        varlogFR6 <- ifelse(is.na(varS.0) | is.na(varS.1), NA, varS.1/(F.1^2) + varS.0/(F.0^2))
      }
    }
    
    ## Compare HR (and FR) to the boundaries
    if ( Lower ){  # high efficacy monitoring
      if (estimand=="cox"){
        HRci.up <- exp(coxPH.i$coef+qnorm(1-alphaHigh/2)*sqrt(coxPH.i$var))        
        if ( HRci.up < highHR  ) {
          done <- TRUE
          bound <- "HighEff"
        } else {
          if (i==K) {   ## the end of high eff monitoring
            done = TRUE
            bound = "notHighEff"
          }
        }
      }
      
      if (estimand!="cox"){
        FRci.up <- ifelse(is.na(varlogFR), NA, exp(log(FR.i)+qnorm(1-alphaHigh/2)*sqrt(varlogFR)))        
        if ( FRci.up < highHR ) {
          done <- TRUE
          bound <- "HighEff"
        } else {
          if (i==K) {   ## the end of high eff monitoring
            done = TRUE
            bound = "notHighEff"
          }
        }
      }
    } else {  # either non-efficacy monitoring or end-of-Stage 1 test for advancement into Stage 2
      if (i == K && infectionTotals[i]==n ) { ##  test for positive efficacy at the end of Stage 1 for advancement into Stage 2      
        done <- TRUE
        nVaxInfStage1 <- sum(D$event.i[D$trt==1])
        
        if (infectionsInOnlyOneGroup==TRUE){  # if there are zero infections in the vaccine group
          bound <- "Eff"
        } else {
          if (estimand=="cox"){ 
            HRci.up <- exp(coxPH.i$coef+qnorm(1-alphaStage1/2)*sqrt(coxPH.i$var))            
            bound <- ifelse(HRci.up < stage1HR, "Eff", "NonEffFinal")
          }
        
          if (estimand!="cox"){ 
            FRci.up <- ifelse(is.na(varlogFR), NA, exp(log(FR.i)+qnorm(1-alphaStage1/2)*sqrt(varlogFR)))            
            bound <- ifelse(FRci.up < stage1HR, "Eff", "NonEffFinal")
          }          
        }       
      } else {  # non-efficacy interim analysis before the end of Stage 1
        if (post6moMonitor){  # estimand must be different from 'combined'          
          if (estimand=="cox"){
            HRci.lw <- exp(coxPH.i$coef-qnorm(1-alphaNoneff/2)*sqrt(coxPH.i$var))  
            HRci.up <- exp(coxPH.i$coef+qnorm(1-alphaNoneff/2)*sqrt(coxPH.i$var))
            HRci.lw6 <- exp(coxPH.i6$coef-qnorm(1-alphaNoneff/2)*sqrt(coxPH.i6$var))  
            HRci.up6 <- exp(coxPH.i6$coef+qnorm(1-alphaNoneff/2)*sqrt(coxPH.i6$var))
            if ( all(lowerHRnoneff < c(HRci.lw, HRci.lw6)) && all(upperHRnoneff < c(HRci.up, HRci.up6)) ) {
              done <- TRUE
              bound <- "NonEffInterim"
            }
          }
          
          if (estimand=="cuminc"){
            FRci.lw <- ifelse(is.na(varlogFR), NA, exp(log(FR.i)-qnorm(1-alphaNoneff/2)*sqrt(varlogFR))) 
            FRci.up <- ifelse(is.na(varlogFR), NA, exp(log(FR.i)+qnorm(1-alphaNoneff/2)*sqrt(varlogFR)))
            FRci.lw6 <- ifelse(is.na(varlogFR6), NA, exp(log(FR.i6)-qnorm(1-alphaNoneff/2)*sqrt(varlogFR6))) 
            FRci.up6 <- ifelse(is.na(varlogFR6), NA, exp(log(FR.i6)+qnorm(1-alphaNoneff/2)*sqrt(varlogFR6)))
            if ( all(lowerHRnoneff < c(FRci.lw, FRci.lw6)) && all(upperHRnoneff < c(FRci.up, FRci.up6)) ) {
              done <- TRUE
              bound <- "NonEffInterim"
            }
          }
        } else {
          if (estimand=="combined"){
            HRci.lw <- exp(coxPH.i$coef-qnorm(1-alphaNoneff/2)*sqrt(coxPH.i$var))
            HRci.up <- exp(coxPH.i$coef+qnorm(1-alphaNoneff/2)*sqrt(coxPH.i$var))
            FRci.lw <- ifelse(is.na(varlogFR), NA, exp(log(FR.i)-qnorm(1-alphaNoneff/2)*sqrt(varlogFR)))
            FRci.up <- ifelse(is.na(varlogFR), NA, exp(log(FR.i)+qnorm(1-alphaNoneff/2)*sqrt(varlogFR)))          
            if ( all( lowerHRnoneff < c(HRci.lw, FRci.lw) ) && all( upperHRnoneff < c(HRci.up, FRci.up) ) ) {
              done <- TRUE
              bound <- "NonEffInterim"
            }
          }
          
          if (estimand=="cox"){
            HRci.lw <- exp(coxPH.i$coef-qnorm(1-alphaNoneff/2)*sqrt(coxPH.i$var))  
            HRci.up <- exp(coxPH.i$coef+qnorm(1-alphaNoneff/2)*sqrt(coxPH.i$var))
            if ( lowerHRnoneff < HRci.lw && upperHRnoneff < HRci.up ) {
              done <- TRUE
              bound <- "NonEffInterim"
            }
          }
          
          if (estimand=="cuminc"){
            FRci.lw <- ifelse(is.na(varlogFR), NA, exp(log(FR.i)-qnorm(1-alphaNoneff/2)*sqrt(varlogFR))) 
            FRci.up <- ifelse(is.na(varlogFR), NA, exp(log(FR.i)+qnorm(1-alphaNoneff/2)*sqrt(varlogFR)))
            if ( lowerHRnoneff < FRci.lw && upperHRnoneff < FRci.up ) {
              done <- TRUE
              bound <- "NonEffInterim"
            }
          }
        }                
      }
    }
    
    # test for design alternative
    if (done) {  # i.e., when the trial stops
      if (!is.null(upperHRuncPower)){
        if (infectionsInOnlyOneGroup==TRUE){  # if there are zero infections in the vaccine group
          altDetected <- !logical(length(upperHRuncPower))
        } else {
          if (estimand=="cox"){
            HRci.up <- exp(coxPH.i$coef+qnorm(1-alphaUncPower/2)*sqrt(coxPH.i$var))            
            altDetected <- HRci.up < upperHRuncPower
          }
          
          if (estimand!="cox"){
            FRci.up <- ifelse(is.na(varlogFR), NA, exp(log(FR.i)+qnorm(1-alphaUncPower/2)*sqrt(varlogFR)))            
            altDetected <- FRci.up < upperHRuncPower
          }
        }
      }
    }
    
    if (done) break
  }  
  
  if (done) {
    ## Time since the first subject was enrolled and number of infections at which we stopped
    stopTime <- t.i - firstEnrollTime
    Ninfec <- summObj$infectTotal[i]
  } else {
    stopTime <- NA 
    Ninfec   <- NA 
  }
  
  ## get the time of the last ppt's exit from the study
  ## Really only needed when the trial doesn't stop via testing (i.e.
  ## not enough infections accrue), but return for all cases anyway.
  lastExitTime <- max(d$exit) - firstEnrollTime
  
  ## subset summary object to just rows for tests done
  summObjsub <- summObj[1:i, ]
  
  if(estimand!="cuminc" & !infectionsInOnlyOneGroup & !is.na(bound)){ 
    if (bound=="NonEffInterim"){ alpha <- alphaNoneff }
    if (bound %in% c("Eff","NonEffFinal")){ alpha <- alphaStage1 }
    if (bound %in% c("HighEff","notHighEff")){ alpha <- alphaHigh }
    HRci.lw <- exp(coxPH.i$coef-qnorm(1-alpha/2)*sqrt(coxPH.i$var))  
    HRci.up <- exp(coxPH.i$coef+qnorm(1-alpha/2)*sqrt(coxPH.i$var))
    CI.out <- c(HRci.lw,HRci.up) 
  } else { 
    CI.out <- NA 
  }
  
  return( list( finished = done,
                boundHit = bound, 
                altDetected = altDetected,
                stopTime = stopTime, 
                stopInfectCnt = Ninfec,
                totInfecCnt = n,
                totInfecSplit= table(eventDF$trt),
                lastExitTime = lastExitTime,
                finalHRci = CI.out,
                summObj = summObjsub,
                nVaxInfStage1 = nVaxInfStage1) )
}

## Function to determine the number (i.e. count) of the first event that
## meets both of the following criteria:
##   (a) the cumulate percentage of events having some property meets or
##       exceeds the value given in 'minPercent',  and
##   (b) the number/count of the event is at least as large as 'minCount'
##
## The arguments to the function are: 
## ----------------------------------
##   'x': an indicator vector (0/1 or FALSE/TRUE) reporting whether each
##        event has some property.  The vector should only contain info
##        on *EVENTS*. 
##
##   'minPercent': the threshold that must be met for the cumulative 
##                 percentage of "positive" entries (1s or TRUEs) in 'x'
##
##   'minCount':   the minimum event total that must be attained along
##                 with the minimum percentage threshold
##
## Example:  
##   We have defined the infection total at which we wish to perform the first
##   non-efficacy analysis for a vaccine as the smallest infection count for
##   which the percentage of infections occurring at least 6-months
##   post-enrollment reaches/exceeds 30% of infections, with the additional
##   criteria that the total must also be at least 50 (to be of sufficient 
##   size that the "large sample" normal-approximations used in the monitoring
##   criteria are reasonable.
##
get_firstTest_infecCnt <- function(x, minPercent = 0, minCount = 0) {
  
  len.x <- length(x)
  pct <- 100 * cumsum( x ) / (1:len.x )
  
  pctCriteria <- (pct >= minPercent)
  cntCriteria <- ( ( 1:len.x ) >= minCount )
  
  if ( any(both <- pctCriteria & cntCriteria) ) {
    ## if any values meet both criteria, return the first one
    return( which(both)[1] )
  } else {
    ## if no values meet both criteria, then either:
    ##   (a) we never exceeded minCount, or
    ##   (b) we never exceeded minPercent 
    ## (possibly both, but not likely)
    if ( len.x < minCount ) {
      outname <- "Too Few Events"
    } else if ( pctCriteria[len.x] < minPercent ) {
      outname <- "Pct Too Low" 
    } else {
      outname <- "Other"
    }
    
    ## return an 'NA' with a name attached giving info on the reason
    return( structure( NA, names = outname) )
  }
}


## Function to implement "harm" monitoring on data.frame 'd' using boundaries
## specified in data.frame 'bounds'.
##
## The 'bounds' data frame must contain the columns specified by the arguments
## 'totInfecVar' and 'vaccInfecVar'.  These columns should contain a
## running total of infections in both trt groups (column 'totInfecVar'), and a 
## running total of infections in the vaccinee group ('vaccInfecVar').
## Other columns are allowed, but will be ignored.
##
## Data.frame 'd' will contain one row per infection, a 'trt' variable and a
## variable 'exit' containing the study time that the infection is observed
## (i.e. diagnosed) - the infected ppt.s "exit time" from un-infected Follow-up
do_harm_monitoring <- function(d, bounds, totInfecVar="N", vaccInfecVar="V") {
  
  ## sanity checks on input
  if ( length(d)==0 )
    stop("Argument 'd' to function 'do_harm_monitoring' has length 0.\n")
  
  if ( length(bounds)==0 )
    stop("Argument 'bounds' to function 'do_harm_monitoring' has length 0.\n")
  
  ## check that 'd' is ordered - else order it
  if ( nrow(d) > 1 & !all( diff(d$exit) >= 0 ) ) {
    d <- d[ order(d$exit), ]
  }
  
  ## get cumulative count of vaccinee infections
  vaccInfecCnt <- cumsum( d$trt > 0 ) # cumulatively counts vaccinee infections
  
  ## get count of number of infections
  totInfec <- length( vaccInfecCnt )
  
  ## restrict 'bounds' to apply to only infection totals we have
  bounds <- bounds[ bounds[, totInfecVar] <= totInfec , ]
  
  ## subset out the components of 'vaccInfecCnt' that correspond to the 
  ## infections totals at which we can stop for harm (the values in 
  ## bounds[[ totInfecVar ]]).   Compare the subsetted values to the 
  ## values in: bounds[[ vaccInfecVar ]].  If any are equal we stop.
  vaccInfecSub <- vaccInfecCnt[ bounds[[totInfecVar ]] ] # matching 'vaccInfecCnt' to 'bounds'
  harmBoundsHit <- ( vaccInfecSub == bounds[[ vaccInfecVar ]])
  
  ## remove NA (added by Yu on 12/28/2011)
  if ( any( harmBoundsHit, na.rm=TRUE ) ) {
    
    w.hit <- which( harmBoundsHit )[1]
    
    ## 'N' is the infection count at the harm bound hit, 
    ## 'V' is the vaccine total at the harm bound hit
    N <- bounds[[totInfecVar]][ w.hit ]
    V <- bounds[[vaccInfecVar]][ w.hit ]
    
    list( isHarm = TRUE,
          stopTime = d$exit[ N ], # 'd' is ordered by 'exit' time
          stopInfectCnt = N,
          stopInfecSplit = c(Vacc = V, Plac = N - V))
    
  } else {
    list( isHarm = FALSE )
  }
}


## Function to implement "harm" monitoring on data.frame 'd' using boundaries
## specified in data.frame 'bounds'.
##
## The 'bounds' data frame must contain the columns specified by the arguments
## 'totInfecVar' and 'vaccInfecVar'.  These should columns should contain a
## running total of infections in both trt groups (column 'totInfecVar'), and a 
## running total of infections in the vaccinee group ('vaccInfecVar').
## Other columns are allowed, but will be ignored.
##
## Data.frame 'd' will contain one row per infection, a 'trt' variable and a
## variable 'exit' containing the study time that the infection is observed
## (i.e. diagnosed) - the infected ppt.s "exit time" from un-infected Follow-up
do_harm_monitoring2 <- function(d, bounds, d2, stage1 =78,  totInfecVar="N", vaccInfecVar="V") {
  
  ## sanity checks on input
  if ( length(d)==0 )
    stop("Argument 'd' to function 'do_harm_monitoring' has length 0.\n")
  
  if ( length(bounds)==0 )
    stop("Argument 'bounds' to function 'do_harm_monitoring' has length 0.\n")
  
  
  ## check that 'd' is ordered - else order it
  if ( nrow(d) > 1 & !all( diff(d$exit) >= 0 ) ) {
    d <- d[ order(d$exit), ]
  }
  
  ## get cumulative count of vaccinee infections
  vaccInfecCnt <- cumsum( d$trt > 0 )
  
  ## get count of number of infections
  totInfec <- length( vaccInfecCnt )
  
  ## restrict 'bounds' to apply to only infection totals we have
  bounds <- bounds[ bounds[, totInfecVar] <= totInfec , ]
  
  ## subset out the components of 'vaccInfecCnt' that correspond to the 
  ## infections totals at which we can stop for harm (the values in 
  ## bounds[[ totInfecVar ]]).   Compare the subsetted values to the 
  ## values in: bounds[[ vaccInfecVar ]].  If any are equal we stop.
  vaccInfecSub <- vaccInfecCnt[ bounds[[totInfecVar ]] ]
  harmBoundsHit <- ( vaccInfecSub == bounds[[ vaccInfecVar ]])
  
  ## remove NA (added by Yu on 12/28/2011)
  if ( any( harmBoundsHit, na.rm=TRUE ) ) {
    
    w.hit <- which( harmBoundsHit )[1]
    
    ## 'N' is the infection count at the harm bound hit, 
    ## 'V' is the vaccine total at the harm bound hit
    N <- bounds[[totInfecVar]][ w.hit ]
    V <- bounds[[vaccInfecVar]][ w.hit ]
    
    ## calculate the stop time at stage 1, i.e., 
    ## follow all the subjects enrolled before 'stopTime' until 18 month
    t.i = d$exit[ N ]
    D <- subset(d2, entry < t.i )  
    
    stopTimeStg1 = max(D$exit)
    
    list( isHarm = TRUE,
          stopTime = d$exit[ N ],
          stopTimeStg1 = stopTimeStg1,
          stopInfectCnt = N,
          stopInfecSplit = c(Vacc = V, Plac = N - V))
    
  } else {
    list( isHarm = FALSE )
  }
}

## Functions for harm boundary
## This function creates an object containing the values of the
## function P(n, s) which Breslow (1970, JASA) defined to be
## for 0 <= s <= n <= N), the probability of the binomial random
## walk S_n continuing to S_n = s without "absorption" in the
## rejection region (i.e. without hitting the stopping bounds).
##
## N = max number of samples
##
## S_n = sum{ x_i, i=1,..,n }, and the {x_i} are iid bernoulli(p)
##
## To calculate this, we need to specify a value of 'p', 'N', and
## a vector 'B' of length N that specifies the stopping boundares
## associated with the values of 'n' from 1:N (respectively).  So,
## B[1] will be the stopping value for n=1, B[2] the stopping value
## for n=2, ..., and B[N] the stopping value for n=N.
##
##
## The object is a list with components:
##
##   'p' - value of 'p' passed to function
##   'N' - value of 'N' passed to function
##   'Bounds' - data.frame with columns 'n' = 1:N and
##              'StoppingBound'= argument 'Bound'
##
##   'Pns' =  list of length N, each sublist containing a numeric vector.
##            Pns[[ i ]] is a numeric vector containing values of P(n,s)
##            for n=i (i.e. values of P(i, s)) for values 0<= s <= i
##
##            Note that values of s >= Bound[i] will be set to zero (see
##            defn of P(n,s) at top to see why).
##
##   'Stop' - numeric vector of length N, with Stop[i] giving the prob.
##            of hitting the stopping bound (for first time) at n=i


pNS <- function(Bound, p=.5, N=45)
{
  if( length(Bound) != N )
    stop("Length of vector 'Bound' must equal value of argument 'N'")
  
  
  ## create 'Pns' and 'Stop'
  Stop <- numeric(N)
  Pns <- vector("list", length=N)
  
  if ( is.na(Bound[1]) || Bound[1]>1 )
  {
    Pns[[ 1 ]] <- c("0" = (1-p), "1" = p)
    Stop[ 1 ] <- 0
  } else
    stop("Why are you allowing stopping when n=1!!!")
  
  for (i in 2:N)
  {
    pv <- numeric(i+1)
    names(pv) <- as.character(0:i)
    Pns[[ i ]] <- pv
    
    max.S <- min( i, Bound[i]-1, na.rm=TRUE)
    
    Pns[[i]]["0"] <- (1-p)*Pns[[i-1]]["0"]
    
    for (s in as.character(1:max.S) )
    {
      s.minus.1 <- as.character( as.numeric(s)-1)
      
      if (as.numeric(s) < i )
        Pns[[i]][ s ] <- p*Pns[[i-1]][s.minus.1] + (1-p)*Pns[[i-1]][s]
      else
        Pns[[i]][ s ] <- p*Pns[[i-1]][s.minus.1]
    }
    
    ## stopping prob. is the prob. we were one below the current bound at the
    ## last 'n', times prob. that we had another 'success'
    if ( !is.na(Bound[i]) )
      Stop[ i ] <- p*Pns[[ i-1 ]][ as.character(Bound[i]-1) ]
  }
  
  Bounds <- data.frame( n=1:N, StoppingBound=Bound )
  
  totalStopProb <- sum( Stop )
  ExpStopTime <- sum( (1:N)*(Stop) )
  
  
  list(p = p, N = N, Bounds = Bounds, Pns = Pns, Stop = Stop,
       totalStopProb = totalStopProb,  ExpStopTime = ExpStopTime)
}

####################### end of function 'pNS'#############################

### THIS USES 'nonConstBounds' framework to do constant-bounds.
### Actually, it's constant from 5 to 45, it's zero before that, so
### technically it is non-constant...


## 'x' is the total number of infections (vacc + placebo)
## 'alphaVals' is a vector of nominal (un-adjusted) p-value thresholds
##    to use in establishing cutoffs for harm-monitoring.  This vector
##    must have the same length as 'startHarmMonitor'.  The i-th value
##    of 'alphaVals' applies to the i-th interval defined by 'startHarmMonitor'
## 'startHarmMonitor' gives the endpoints of all intervals.
##    The starting point of the first interval is 1, and of the the i-th interval
##    (for i>1) is 1 + startHarmMonitor[i-1]
semiConstSpending <- function(x, alphaVals, startHarmMonitor )
{
  which.interval <- findInterval(x, c(1, startHarmMonitor),
                                 rightmost.closed=TRUE)
  alphaVals[ which.interval ]
}

## function to create CSV filenames
fileNameFunc <- function( p, N, null.p)
{
  paste0("harmBounds_N=", N, "_alphaPerTest=", p, "_pVacc=", null.p, ".csv")
}

getAlphaPerTest <- function(harmMonitorRange, null.p){
  getCumAlpha <- function(alphaPerTest, harmMonitorRange, null.p){ 
    harmBounds <- getHarmBound(N=harmMonitorRange[2], per.test=alphaPerTest, harmMonitorRange=harmMonitorRange, null.p=null.p) 
    return(harmBounds$cumStopProb[NROW(harmBounds)] - 0.05)
  }
  return(uniroot(getCumAlpha, c(0,0.05), harmMonitorRange=harmMonitorRange, null.p=null.p)$root)
}

getHarmBound <- function(N,  ##Total number of infections desired for harm monitoring
                         per.test, ## value for per-test alpha level
                         harmMonitorRange,
                         null.p,
                         dataDir = NULL,
                         verbose = TRUE){
  ## Note: 
  ##   'null.p' = the probability that an infection occurs in a vaccinee, under the null 
  ##              hypothesis that infection is equally likely in vaccinees and placebo
  ##              recipients.  Hence 'null.p' equals the fraction of the populations that
  ##              has received vaccine.  This would be 0.5 under a 1:1 randomization, or
  ##              11/29 under a 11:18 randomization (V:P).
  
  ## Storage for output of per.test.alpha and
  summ <- data.frame(alpha = per.test, totAlpha = NA)
  
  ## We consider the total number of infected participants to be 'N' and
  ## assume apriori and equal likelihood of infection for vaccinees as for
  ## placebos.  So the number of infected vaccinees should follow a binomial
  ## distribution with size N and probability p=null.p (where 'null.p' is 
  ## the proportion of vaccinees in the trial).  We wish to have a total
  ## probability of a "type I error" (stopping the trial for 'harm' when there
  ## is no real difference between vaccine and placebo) of .05.  Our approach
  ## will be to test for harm after each new infection, and we wish to find a
  ## fixed 'alpha' value to use for each test, so that the overall prob. of a
  ## false positive over the course of the trial is .05.
  ##
  ## To do this, we will iteratively choose a value alpha, generate a set of
  ## 'stopping bounds' corresponding to that alpha, and then estimate the
  ## overall type I error rate for those bounds.  We then go back and adjust
  ## our alpha value (up or down) depending on whether our estimated type I
  ## error is too high or too low.  Repeat until desired accuracy is obtained.
  
  bound <- NULL
  
  ## create data frame to store results in
  bounds <- data.frame(totInfec=1:N, vaccInfecBound=NA, alphaLevelBound=NA,
                       nextHigherAlphaLevel=NA, cutoff=NA )
  
  for (j in 1:nrow(bounds))
  {
    totInfec <- bounds$totInfec[j]
    
    alphaVal <- semiConstSpending( totInfec, alphaVals=c(0, per.test),
                                   startHarmMonitor = harmMonitorRange)
    
    ## we don't need to do the next few steps unless alphaVal is > 0
    if (alphaVal <= 0) next
    
    ## choose the lowerBound for searching for the next cutoff value.
    ## Under our framework it will always be the same as, or higher than
    ## the cutoff from the previous (smaller) value of 'totInfec'
    if ( is.null(bound) ) {
      lowerBnd <- ceiling( null.p * totInfec )
    } else lowerBnd <- bound
    
    # startVal <- min(totInfec, bound)
    # lowerBnd <- startVal - 2
    valSeq <- totInfec:lowerBnd
    
    upperTailProbs <- cumsum( dbinom(valSeq, totInfec, null.p) )
    signif <- ( upperTailProbs <= alphaVal )
    
    ## if we have at least one significant value then do...
    if ( isTRUE(signif[1]) )
    {
      ## if we have all signif. values, then that's not good...
      #if ( isTRUE( rev(signif)[1] ) )
      #{
      #    print( totInfec )
      #    print( valSeq )
      #    print( signif )
      #    stop("Need to include more values in valSeq")
      #}
      
      ## get "largest" (last) index for which signif == TRUE
      largest.index <- max( which( signif ) )
      
      ## define 'bound' to be the infection count corresponding to
      ## the 'largest.index' (i.e. the smallest infection count for
      ## which we have significance at per-test-level 'alphaVal'
      bound <- valSeq[ largest.index ]
      
      bounds$vaccInfecBound[ j ] <- bound
      bounds$alphaLevelBound[ j ] <- upperTailProbs[ largest.index ]
      bounds$cutoff[ j ] <- alphaVal
    }
    
  }
  
  out <- pNS(Bound=bounds$vaccInfecBound, p=null.p, N=N)
  
  names(out$Bounds)[ names(out$Bounds)=="StoppingBound" ] <- "Nvacc"
  boundOut <- transform(out$Bounds, Nplac= n-Nvacc, RR=round(Nvacc/(n-Nvacc),digits=2))
  boundOut <- cbind( boundOut, stopProb=round(out$Stop,4),
                     cumStopProb=round(cumsum(out$Stop),4),
                     alphaVal = bounds$cutoff )
  
  
  overall.alpha <- out$totalStopProb
  
  summ$totAlpha <- overall.alpha
  
  out <- pNS(Bound=bounds$vaccInfecBound, p=null.p, N=N)
  
  
  ## Add info on stopping probabilities to 'bounds' object
  bounds[, "stoppingProb"] <- out$Stop
  bounds[, "cumStoppingProb"] <- cumsum( bounds[, "stoppingProb"] )
  
  harmBounds =  boundOut
  names(harmBounds)[1:3]=c("N", "V", "P") 
  if (!is.null(dataDir)){
    fileName <- fileNameFunc(round(per.test, 4), N, round(null.p, 2))
    write.csv(harmBounds, file.path(dataDir, fileName), row.names=FALSE)
    if (verbose){ cat("Potential-harm stopping boundaries saved in:\n", file.path(dataDir, fileName), "\n\n") }
  }    
  return(harmBounds)  
}


simTrial <- function(N,
                    aveVE,
                    VEmodel=c("half", "constant"),
                    vePeriods,
                    enrollPeriod,
                    enrollPartial,
                    enrollPartialRelRate,
                    dropoutRate,
                    infecRate,
                    fuTime,
                    visitSchedule,
                    missVaccProb = NULL,
                    VEcutoffWeek,
                    nTrials,
                    blockSize = NULL,
                    stage1,
                    saveDir = NULL,
                    verbose = TRUE,
                    randomSeed = NULL){
VEmodel <- match.arg(VEmodel)
  
## verify whether length of 'N' = length of 'aveVE'
if ( length(aveVE) != length(N) ) 
   stop( "Length of 'aveVE' does not match the number of treatment arms given in 'N'.\n" )

## total number of trial participants
Nppt = sum(N)

## VE for placebo arm
nullVE = aveVE[1]

## VE for vaccine arms
aveVE = aveVE[-1] 

## number of vaccine arms
nVaccArms = length(aveVE)

## total number of arms (assumes a single control arm)
nArms <- nVaccArms + 1

## create vector of treatment group names: 
trtArms <- c("C1", paste("T", 1:nVaccArms, sep=""))

## 'enrollRate' is the required weekly enrollment rate needed to enroll the expected 'Nppt' subjects 
## in 'enrollPeriod' accounting for partial enrollment rate during the initial 'enrollPartial' weeks
enrollRate <- Nppt/(enrollPartialRelRate * enrollPartial + enrollPeriod - enrollPartial)

## 'trtAssgnProbs' contains treatment assignment probabilities
trtAssgnProbs <- structure(N/Nppt, names=trtArms)

## specify block size if not provided by user
if ( is.null(blockSize) ){
  ## gets the minimum valid blockSize within the interval given by 'range' [10, Infinty)
  blockSize <- getBlockSize( N, range=c(10,Inf) )
}

## The rates in 'parSet' are "base rates" that will be 
## modified with "relative rates" for specific groups and/or time periods.

## For example, the infection rate here should be assumed infection rate
## in the population being recruited from.  Participants assigned to the
## control group and/or treatments with no effect on acquisition will
## have this base infection rate. Participants assigned to treatments
## that lower acquistion will have the base infection rate for some
## initial period of time (before the vaccination series is complete)
## and a lower rate later.

## 'parSet' contains weekly rates
parSet <- list(enrollment=enrollRate, dropout=dropoutRate/52, infection=infecRate/52)

## 'enrollSchedule' contains information on enrollment periods and corresponding enrollment rates  
## 'enrollSchedule' is passed as an argument to the 'simFullEDIdata' function
## partial enrollment within 'enrollPartial' weeks; full enrollment thereafter until the end of week 'enrollPeriod'
enrollSchedule <- data.frame(start = c(1, enrollPartial+1), end = c(enrollPartial,NA), relativeRate=c(enrollPartialRelRate, 1))

if(VEmodel!="constant"){
   ## 'vaccEff' is a vector of true *full* VEs for each treatment (defined as a function of 'aveVE')
   vaccEff <- aveVE * (vePeriods[3]-1)/((vePeriods[3]-1)-(vePeriods[2]-1)/2)

   ## If any value in 'vaccEff' is larger than 1, compute the largest possible average VE,
   ## return it as part of the error message, and stop.
   if ( any(vaccEff >= 1) ){

       ## which VE values are too large?
       w.too.large <-which( vaccEff >= 1 )

       ## compute maximum possible average VE such that 'vaccEff' is equal to 1
       ## we need to use a value that is less than 'maxEq'
       maxEq <- ((vePeriods[3]-1)-(vePeriods[2]-1)/2)/(vePeriods[3]-1)

       useAsMax <- floor(1000 * maxEq)/1000
       if (!(useAsMax < (maxEq - sqrt(.Machine$double.eps)))) {
           useAsMax <- useAsMax - 0.001
       }
       stop("The following average vaccine efficacy(ies) is/are not attainable when",
            " using the \n", "halved VE model: ", aveVE[w.too.large], "\n", 
            "The maximum attainable VE is (approximately) ", useAsMax, "\n")
   }
} else {
   vaccEff <- aveVE
}

## 'VEs' are true vaccine efficacies used for data generation
## 'VEs' is a list with one component per treatment
## each component is a vector of vaccine efficacies applied in various time periods of the trial (e.g.,
## partial-VE, full-VE, and waning-VE period)
VEs <- vector("list", nVaccArms)

for (ii in 1:nVaccArms) {
  if (VEmodel=="half"){
    VEs[[ii]] = c(vaccEff[ii]/2, vaccEff[ii], aveVE[ii])    
  }
  if (VEmodel=="constant"){
    VEs[[ii]] = aveVE[ii]
  }
}
names(VEs) <- paste("T", 1:length(vaccEff), sep="")

if ( length(VEs) != nVaccArms ) 
    stop( "VEs specified don't match the number of vaccine arms given in nVaccArms.\n" )

## 'infecRateTbl' contains information on relative infection rates (hazard ratios) for each treatment.
## Please use "Inf" rather than NA to represent intervals that continue indefinitely.

## Each treatment must have a record starting at time 1 and must not have any time gaps in it.
## It does not need to extend to time "Inf" but typically should.
infecRateList <- vector("list", nArms)
names(infecRateList) <- c("C1", names(VEs))

infecRateList[[1]] <- data.frame( trt = "C1", start = 1, end = Inf, relRate = 1)
for (ii in 2:nArms) {
  trtName <- names(VEs)[ii-1]
  if (VEmodel=="half"){
    infecRateList[[ii]] <- data.frame( trt = trtName, start = vePeriods,
                                       end = c(vePeriods[-1]-1, Inf), relRate = 1 - VEs[[trtName]] )
  }
  if (VEmodel=="constant"){
    infecRateList[[ii]] <- data.frame( trt = trtName, start = vePeriods,
                                       end = Inf, relRate = 1 - VEs[[trtName]] )
  }
}
infecRateTbl <- do.call(rbind, infecRateList)

## specify prior distributions for data simulation
simParList <- parSet
simFuncList <- list(enrollment=Constant, dropout=Constant, infection=Constant) # 'Constant' is a point-mass distribution
simPrior <- list( Function=simFuncList, params=simParList )

## set seed of random number generator
if(!is.null(randomSeed)) {
    set.seed( randomSeed )
}

## get rates for use in generating data
## the current implementation considers constant rates, thus no need for inclusion inside the below 'for' loop
## for 'Constant' prior distributions, 'sampleRates' returns a list identical to 'parSet'
rates <- sampleRates(n=1, from=simPrior)

## create lists for storage of trial data
trialList <- vector("list", nTrials)
infecList = vector("list", nTrials)
infecList2 = vector("list", nTrials)
infecListAll = vector("list", nTrials)
trialResult = vector("list", nTrials)

## 1. Generate enrollment times
## the number of enrolled subjects during a specific time interval is Poisson distributed with rate = 'rate' * (end - start + 1), i.e.,
## N <- rpois(1, lambda = rate * (end - (start-1)))
## enrollment times are uniformly distributed in the (start, end) interval, i.e.,
## runif(N, min=start-1, max=end)

## 2. Generate dropout times
## rexp(N, rate= -log( 1 - rate))

## 3. Generate infection times
## changing the 'rate' parameter to a form that gives the prob(event in interval (0,1)) = rate
## rexp(N, rate= -log( 1 - rate))

for ( i in 1:nTrials )
{
    ## generate data
    EDI.i <- simFullEDIdata(
                  rateParams = rates,
                  trtAssignProb = trtAssgnProbs,
                  blockSize = blockSize,
                  infecRates = infecRateTbl,
                  protocolVisits = visitSchedule,
                  enrollPeriod = enrollSchedule,
                  nEnroll= Nppt,
                  maxEnrollment = Nppt,
                  nWeeksFU = fuTime
              )

    ## code treatments as integers: C1=0, T1=1, T2=2...
    trtCode <- match(EDI.i$trt, trtArms) - 1


    entry <- EDI.i$enrollTime
    exit  <- entry + EDI.i$futime


    ## create flag for HIV infection
    infected <- !is.na(EDI.i$infecDxTime)

    ## make EDI.i data into a simpler form
    out <- data.frame( trt   = as.integer(trtCode),
                       entry = entry,
                       exit  = exit,
                       event = as.integer(infected)
                       )
    
    if (!is.null(missVaccProb)){
      # create a set of indicators of belonging to a per-protocol cohort
      u <- runif(NROW(out))
      ppnames <- paste0("pp", 1:length(missVaccProb))
      for (ppIdx in 1:length(missVaccProb)){
        out[[ppnames[ppIdx]]] <- ifelse(out$exit - out$entry > VEcutoffWeek & u >= missVaccProb[ppIdx], 1, 0)
      }
    }    
                         
    ## count number of infections by arm
    ## stg1 infections: infections occurring in first 'stage1" weeks of trial
    stg1_Infec <- infected & EDI.i$futime <= stage1

    ## add infection counts for 24-36 stage 2
    stg2_Infec <- infected & (EDI.i$futime > stage1 & EDI.i$futime <=fuTime)

    cntVec <- as.vector( table( as.factor(trtCode)[ stg1_Infec ] ) )
    cntVec2 <- as.vector( table( as.factor(trtCode)[ stg2_Infec ]))
    cntVecAll <- as.vector( table( as.factor(trtCode)[ infected ]) )
 
    infecList[[ i ]] <- cntVec
    infecList2[[ i ]] <- cntVec2
    infecListAll[[ i ]] <- cntVecAll

    ## store summary data into trialList
    trialList[[ i ]] <- out
}

## summary number of infections
infectCnts <- vector("list", length=nVaccArms)
infectCnts2 <- vector("list", length=nVaccArms)
infectCntsAll <- vector("list", length=nVaccArms)
    
## Put everything into a "trial Object"
trialObj <- list( trialData = trialList,
                  nTrials = nTrials,
                  N = Nppt,
                  nArms = nArms,
                  trtAssgnProbs = trtAssgnProbs,
                  blockSize = blockSize,
                  fuTime = fuTime,
                  rates = parSet,
                  enrollSchedule = enrollSchedule,
                  VEs = VEs,
                  infecRates= infecRateTbl,
                  randomSeed = randomSeed,
                  NinfStage1 = infecList,
                  NinfStage2 = infecList2,
                  NinfAll = infecListAll                  
                 )

  # save trial output and information on used rates
  if (!is.null(saveDir)){
    saveFile <- paste("simTrial_nPlac=", N[1], "_nVacc=", paste(N[-1], collapse="_"), "_aveVE=", paste(round(aveVE,2), collapse="_"), "_infRate=", infecRate,".RData", sep="")
    save(trialObj, file=file.path(saveDir, saveFile))
    if (verbose){ cat("Output saved in:\n", file.path(saveDir, saveFile), "\n\n") }
  } else {
    return( trialObj )
  }
}
########################### End of simTrial function ####################################################

## Function to determine the number (i.e. count) of the first event that
## meets all the following criteria:
##   (a) the cumulate percentage of events after "week1" having some property meets or
##       exceeds the value given in 'minPercent',  and
##   (b) the number/count of the event is at least as large as 'minCount'
##   (c) at least 2 infections after month 12 visit
##
## The arguments to the function are: 
## ----------------------------------
##   'x': an indicator vector (0/1 or FALSE/TRUE) reporting whether each
##        event has some property.  The vector should only contain info
##        on *EVENTS*. 
##
##   'minPercent': the threshold that must be met for the cumulative 
##                 percentage of "positive" entries (1s or TRUEs) in 'x'
##
##   'minCount':   the minimum event total that must be attained along
##                 with the minimum percentage threshold
##
##    week1 = 26:  first time cutoff for 'minPercent'          
##    nInfecAfterwk = 2:   number of infections required after week2
##    week2 = 52:          second time cutoff for 'nInfecAfterwk'
## 
## Example:  
##   We have defined the infection total at which we wish to perform the first
##   non-efficacy analysis for a vaccine as the smallest infection count for
##   which the percentage of infections occurring at least 6-months
##   post-enrollment reaches/exceeds 25% of infections, with at least 2 infections
##   after week 52 and the additional
##   criteria that the total must also be at least 50 (to be of sufficient 
##   size that the "large sample" normal-approximations used in the monitoring
##   criteria are reasonable).
##

getInfecCntFirstNonEff <- 
  function(  x,                 ## 'stage1' placebo and j-vaccine infections ordered by 'exit' time
             minPercent = 0,    ## min. fraction of infections that have occurred after 'week1'
             minCount = 0,      ## min. number of total infection required
             week1 = 26,        ## 
             nInfecAfterwk = 2, ## number of infections required after week2
             week2 = 52)        ## the first interim analysis can only occur after 'week2' 
  {
    ## number of infections at which first non-efficacy monitoring starts
    N1 <- NA
    
    ## time that reach 'minCnt' infections 
    nInfec <- nrow(x)
    
    ## number of infections after 'week2'
    nPostWk2 <- sum( x$futime > week2 )
    
    if ( (nInfec >= minCount) && (nPostWk2 >= nInfecAfterwk) ) {
      
      ## create column containing fraction of infections post-week1
      postWk1_fract <- cumsum(x$futime>week1)/(1:nInfec) 
      
      ## the value of N1 is the smallest 'N' that satisfies all three conditions
      ## SIMULTANEOUSLY
      N1 <- which( (x$nInf >= minCount) & (postWk1_fract >= minPercent) &
                     ( cumsum(x$futime > week2) >= 2 ) )[ 1 ]
      
      ## if nothing meets all three criteria then 'N1' will have length 0
      if ( length(N1) == 0 )
        N1 <- NA
    } 
    return( N1 ) 
  }


monitorTrial = function(dataFile,
                        stage1,
                        stage2,
                        harmMonitorRange,
                        alphaPerTest=NULL,                        
                        minCnt,
                        minPct,
                        week1,
                        minCnt2,
                        week2,
                        nonEffInterval,
                        lowerVEnoneff,    # lower confidence bound to be below 'lowerVEnoneff' to meet criterion 1 for declaring non-efficacy
                        upperVEnoneff,    # upper confidence bound to be below 'upperVEnoneff' to meet criterion 2 for declaring non-efficacy
                        highVE,           # lower confidence bound to be above 'highVE' for declaring high efficacy
                        stage1VE,         # lower confidence bound to be above 'stage1VE' for advancing into Stage 2
                        lowerVEuncPower=NULL,  # lower confidence bound to be above 'lowerVEuncPower' for rejecting H0: VE<=lowerVEuncPower*100%                        
                        alphaNoneff,      # One minus confidence level of 2-sided CI for non-efficacy monitoring
                        alphaHigh,        # One minus confidence level of 2-sided CI for high efficacy monitoring
                        alphaStage1,      # One minus confidence level of 2-sided CI for testing whether an arm advances into Stage 2
                        alphaUncPower=NULL, # One minus confidence level of 2-sided CI for unconditional power to reject H0: VE<=lowerVEuncPower*100%                        
                        estimand=c("combined", "cox", "cuminc"),
                        post6moMonitor = FALSE,
                        VEcutoffWeek,
                        saveDir = NULL,
                        verbose = TRUE){
  estimand <- match.arg(estimand)
  
  upperHRnoneff <- 1 - lowerVEnoneff
  lowerHRnoneff <- 1 - upperVEnoneff
  stage1HR <- 1 - stage1VE
  upperHRuncPower <- NULL
  if (!is.null(lowerVEuncPower)){ upperHRuncPower <- 1 - lowerVEuncPower }
  highHR <- 1 - highVE
  
  if (!is.null(saveDir)){
    ## load in RData object (a list named 'trialObj' )
    load(file.path(saveDir, dataFile))
  } else {
    trialObj <- dataFile
    rm(dataFile)
  }  
  
  ## check contents of 'trialData'
  d= trialObj[["trialData"]][[1]]
  if ( !all( c("entry","exit","event","trt") %in% names(d) ) )
    stop("DataFrame 'd' must contain columns: entry, exit, event, trt\n")
  
  nTrtArms <- as.integer( trialObj$nArms - 1 )
  nTrials <- length(trialObj[["trialData"]])
  
  # calculate 'null.p' the probability of being assigned vaccine considering only 1 vaccine arm
  nPpts <- trialObj$N
  pVaxPla <- trialObj$trtAssgnProbs[1:2]
  nVaxPla <- pVaxPla*nPpts
  null.p <- nVaxPla[2]/sum(nVaxPla)
  names(null.p) <- NULL
  
  # calculate stopping boundaries for harm
  if (is.null(alphaPerTest)){ alphaPerTest <- getAlphaPerTest(harmMonitorRange, null.p) }
  harmBounds <- getHarmBound(N=harmMonitorRange[2], per.test=alphaPerTest, harmMonitorRange=harmMonitorRange,
                             null.p=null.p, dataDir=saveDir, verbose=verbose)
  
  ## creates a list of length 'nTrials' each element of which is a list of 
  ## length 'nTrtArms'
  out <- rep( list(vector("list",nTrtArms)), nTrials )
  
  for (i in 1:nTrials ) {
    
    ## extract data for the i-th trial
    datI <- trialObj[["trialData"]][[ i ]]
    minEnrollTime <- min(datI$entry, na.rm=TRUE)
    datI$futime <- datI$exit - datI$entry
    
    ## restrict data to first stage1 *within each ppt*
    ## first store original (uncensored) exit time into a diff. variable
    datI$exitUncens <- datI$exit
    datI$eventUncens <- datI$event
    
    ## variables for high eff. monitoring,
    ## highEffTime defined earlier in the input parameter
    datI$event.H <- datI$event == 1 & ( datI$futime <= stage2 )
    datI$exit.H  <- pmin( datI$exit, datI$entry + stage2 )
    
    ## variables for all other monitoring, stage 1
    ## stage1 defined earlier in the input parameter
    datI$event <- datI$event == 1 & ( datI$futime <= stage1 )
    datI$exit  <- pmin( datI$exit, datI$entry + stage1 )
    datI$futime <- datI$exit - datI$entry
    
    ## variables for post 6 month VE, infections occur before 6 month are censored at 6 month
    datI$eventPost6m = ifelse (datI$futime<=VEcutoffWeek & datI$event==1, 0, datI$event)
    
    ## create an indicator of post-6 month events - needed for determining
    ## when a certain pct of the infections occurred 6-months-post-trial entry
    datI$post6mo <- (datI$futime > VEcutoffWeek)
    
    ## create separate set of only infection ('events')
    eventDF <- subset(datI, event == 1)
    
    ## order the events by trial time at which they are observed
    eventDF <- eventDF[order(eventDF$exit), ]
    
    ## Now we move to comparing each active trt arm with the placebo arm
    for (j in 1:nTrtArms) {
      
      ## subset *all data* down to the two arms being compared
      datI.j <- subset(datI, trt %in% c(0,j) )
      
      ## convert 'trt' to indicator variable before passing it to
      ## 'applyStopRules' (i.e. convert the non-zero values to 1)
      datI.j$trt <- as.integer(datI.j$trt > 0 )
      
      ## subset events in relevant arms: j-th active trt and placebo(trt=0)
      E.j <- subset(eventDF, trt %in% c(0,j) )
      nInfec <- nrow(E.j) # counts infections through 'stage1'
      E.j$nInf= 1:nInfec
      
      ## data for high efficacy monitoring
      datIH.j = subset(datI.j, select=c("entry", "exit.H", "event.H", "trt"))
      datIH.j$event = datIH.j$event.H
      datIH.j$exit = datIH.j$exit.H
      
      ## Determine the infection total (N1) that we will start monitoring
      ## for futility at:
      if ( nInfec >= minCnt ) {
        
        ## Run function that determines "N1" - the infection total (summing
        ## over *both* trt groups) at which the first futility analysis 
        ## will occur.
        N1 <- getInfecCntFirstNonEff( E.j, minPercent = minPct, 
                                      minCount = minCnt, week1 = week1, nInfecAfterwk = minCnt2, week2 = week2)
        
        if ( is.na(N1) ) {   ## not enough infection to trigger futility monitoring
          N1 <- nInfec
          doFut <- FALSE
        } else {
          doFut <- TRUE
        }
      } else {
        ## we have less than minCnt infections in total-- we have to do testing
        ## *something* though.  Do harm and high-eff testing since don't have
        ## bounds for futility
        N1 <- nInfec
        doFut <- FALSE
      }
      
      # 'highN' are infection counts that trigger high-efficacy monitoring
      highN <- N1 + 4*nonEffInterval
      highN.2 <- (highN + sum(datIH.j$event))/2
      if (highN.2 > highN){ highN <- c(highN, highN.2) }
      
      ## Based on value of 'N1', extract bounds for efficacy, futility and harm
      ## then apply them
      harmBounds.j <- subset(harmBounds, N <= N1, select=c(N,V,P))
      
      if ( nrow(harmBounds.j) == 0 ) 
        stop("Harm bounds have length zero - check your code.\n")
      
      ## Evaluate harm monitoring bounds (just based on infected ppts)
      ## so only uses 'E.j'
      harmRes  <- do_harm_monitoring(E.j, bounds = harmBounds.j ) 
      
      ## if there is indication of "harm" then need to create store output
      ## object then move to next trial 
      if ( harmRes$isHarm ) {
        isHarm <- TRUE
        harmInfecTot <- harmRes$stopInfectCnt
        
        harmList <-
          list( finished = TRUE,
                boundHit = "Harm",
                stopTime = harmRes$stopTime - minEnrollTime,
                stopInfectCnt = harmRes$stopInfectCnt,
                stopInfecSplit = harmRes$stopInfecSplit,
                totInfecCnt = nInfec )
      } else { 
        isHarm <- FALSE
      }
      
      ## do high-efficacy monitoring next, even if we hit a 'harm' boundary.
      ## Just in case we can hit a high eff one first.  Yes it seems crazy,
      ## but I'm not going to rule anything out...  If we hit both then I'll
      ## issue a warning to check things...
      
      ## However, if there was harm, then we only need to check for high eff
      ## at timepoints *before* the harm bound was hit - if there are any
      
      if ( isHarm ) {
        ## if there are high efficacy test times before harm time, 
        ## do high eff monitoring 
        # if ( nrow(highEffBounds)!= 0 ) {  
        if ( any(highN < harmInfecTot) ) {
          highEffRes <- applyStopRules(datIH.j,
                                       infectionTotals = highN[highN < harmInfecTot],
                                       boundLabel = "HighEff",
                                       stage1HR = stage1HR,
                                       upperHRuncPower = upperHRuncPower,
                                       highHR = highHR,
                                       alphaStage1 = alphaStage1,
                                       alphaUncPower = alphaUncPower,
                                       alphaHigh = alphaHigh, 
                                       post6moCut=VEcutoffWeek, 
                                       post6moMonitor = post6moMonitor,
                                       estimand=estimand) 
          
          if ( isTRUE(highEffRes$boundHit == "HighEff") ) {
            stop("We Hit both High Eff and harm bounds! ",
                 "Please check code and bounds.\n",
                 "Trial=",i,", Trt Arm=",j, ", ",
                 "EffStopCnt=", highEffRes$stopInfectCnt, ", ",
                 "HarmStopCnt=", harmInfecTot, "\n",
                 "Harm Infec Split:", harmList$stopInfecSplit, "\n")
          }
        }
        # }
        out[[i]][[j]] <- harmList
        next
      } 
      
      ## We only reach here if no harm bounds were hit
      
      # if ( nrow(highEffBounds)!= 0 ) { 
      highEffRes <- applyStopRules(datIH.j,
                                   infectionTotals = highN,
                                   boundLabel = "HighEff",
                                   stage1HR = stage1HR,
                                   upperHRuncPower = upperHRuncPower,
                                   highHR = highHR, 
                                   alphaStage1 = alphaStage1,
                                   alphaUncPower = alphaUncPower,
                                   alphaHigh=alphaHigh, 
                                   post6moCut=VEcutoffWeek,
                                   post6moMonitor = post6moMonitor,
                                   estimand=estimand)
      isHighEff <- highEffRes$boundHit == "HighEff"          
      if (is.na(isHighEff)){ isHighEff <- FALSE }      
      # } 
      
      if ( isHighEff ) {
        heCnt  <- highEffRes$stopInfectCnt
        
        ## if we declare high efficacy at an infection count that is before
        ## the first futility monitoring total, then we don't need to do
        ## futility and declare high eff 
        if ( heCnt < N1 || !doFut) {
          out[[i]][[j]] <- highEffRes
          next
        } else {
          ## only need to monitor for fut at "times" up until heCnt
          seq(N1, heCnt, by =nonEffInterval)
          futRes <- applyStopRules(datI.j,
                                   infectionTotals = seq(N1, heCnt, by =nonEffInterval),             
                                   boundLabel = "NonEff",
                                   lowerHRnoneff = lowerHRnoneff, 
                                   upperHRnoneff = upperHRnoneff, 
                                   stage1HR = stage1HR,
                                   upperHRuncPower = upperHRuncPower,
                                   alphaNoneff = alphaNoneff, 
                                   alphaStage1 = alphaStage1,
                                   alphaUncPower = alphaUncPower,
                                   post6moCut=VEcutoffWeek, 
                                   post6moMonitor = post6moMonitor,
                                   estimand=estimand)           
          if (futRes$finished) {
            out[[i]][[j]] <- futRes
            next
          } else {
            out[[i]][[j]] <- highEffRes
            next
          }
        }
      }
      
      ## only reach here if 'doFut' is TRUE
      ## get a vector of infections for nonefficacy monitoring
      nonEffInfec = seq(N1, nInfec, by =nonEffInterval)
      if (max(nonEffInfec) < nInfec)
        nonEffInfec = c(nonEffInfec, nInfec)
      
      futRes <- applyStopRules(datI.j,
                               infectionTotals = nonEffInfec,             
                               boundLabel = "NonEff",
                               lowerHRnoneff = lowerHRnoneff, 
                               upperHRnoneff = upperHRnoneff,
                               stage1HR = stage1HR,
                               upperHRuncPower = upperHRuncPower,
                               alphaNoneff = alphaNoneff,
                               alphaStage1 = alphaStage1,
                               alphaUncPower = alphaUncPower,
                               post6moCut=VEcutoffWeek,
                               post6moMonitor = post6moMonitor,
                               estimand=estimand) 
      out[[i]][[j]] <- futRes
      
      out[[i]][[j]]$firstNonEffCnt <- N1
      
      out[[i]][[j]]$highEffCnt <- highN
      
      ## add 'lastExitTime' to the object, in case we have efficacy and
      ## the trial continues until the last person in the arm exits the study
      out[[i]][[j]]$lastExitTime <- max( datI.j$exitUncens, na.rm=TRUE ) - minEnrollTime
      
      out[[i]][[j]]$altTest <- !is.null(lowerVEuncPower)
    }
  }
  
  if (verbose){
    for (i in 1:nTrtArms) {
      cat("Probabilities of reaching each possible conclusion:\n")
      print(round(table(unlist( lapply( out, function(x) x[[i]]$boundHit ) ), useNA="ifany")/nTrials, 2))
      
      if (!is.null(lowerVEuncPower)){
        designPower <- colSums(do.call("rbind",lapply(out, function(oneTrial){ oneTrial[[i]]$altDetected })), na.rm=TRUE)/nTrials
        cat("\nUnconditional power to reject the specified null hypotheses =",round(designPower, 2),"\n\n", sep=" ")
      }    
    }
  }  
  
  ## save monitoring output
  if (!is.null(saveDir)){
    saveFile <- paste("monitorTrial", substr(dataFile, 9, nchar(dataFile)-6), "_", estimand, ".RData", sep="")
    save(out, file = file.path(saveDir, saveFile) )
    if (verbose){ cat("Output saved in:\n", file.path(saveDir, saveFile), "\n\n") }
  } else {
    return(out)
  }  
}

######################################## end of 'monitorTrial'#########################################

##############
## This function does the log rank test and estimate VE

finalLogRankTest = function(datI, stage =78, NullHR, alpha=0.025) {
   ## variables for log rank test, censor at 'stage' 
   ## which can be stage 1 or 2 
   datI$futime <- datI$exit - datI$entry
   datI$event <- datI$event == 1 & ( datI$futime <= stage)
   datI$exit  <- pmin( datI$exit, datI$entry + stage )
   datI$futime <- datI$exit - datI$entry

   ## convert 'trt' to indicator variable before calculating VE 
   ## (i.e. convert the non-zero values to 1)
   datI$trt <- as.integer(datI$trt > 0 )

   ## create the model ('coxFormula') to be fit by coxph()
   coxFormula <- Surv(futime, event) ~ trt
   
   ## run cox model
   coxPH.i <- coxph( coxFormula, data=datI, init=log(NullHR)) 

   ## extract the hazard ratio for 'trt' and store         
   hr.i <- exp( coxPH.i$coef )

   ## 95% CI for HR
   HRci.up = exp(coxPH.i$coef+qnorm(0.975)*sqrt(coxPH.i$var))
   HRci.lw = exp(coxPH.i$coef-qnorm(0.975)*sqrt(coxPH.i$var))

   ## log rank test for stage 1
   ## one-side p value < 0.025, efficacy o.w. nonefficacy
   logRank <- survdiff(coxFormula, data=datI)
   if (logRank$obs[2] <= logRank$exp[2]){ 
     p.logRank <- (1-pchisq(logRank$chisq, 1))/2 
   } else {
     p.logRank <- 1 - (1-pchisq(logRank$chisq, 1))/2 
   }

   if(p.logRank < alpha && hr.i < NullHR){
     bound <- "Eff"
   } else {
     bound <- "nonEff"
   }
   list(bound=bound, VE = 1-hr.i)
}

testVE <- function(datI, lowerVE, stage, alpha){
  upperFR <- 1-lowerVE
  ## variables for cumulative hazard-based Wald test, censor at 'stage' 
  ## which can be stage 1 or 2 
  datI$futime <- datI$exit - datI$entry
  datI$event <- datI$event == 1 & ( datI$futime <= stage)
  datI$exit  <- pmin( datI$exit, datI$entry + stage )
  datI$futime <- datI$exit - datI$entry
  
  ## convert 'trt' to indicator variable before calculating VE 
  ## (i.e. convert the non-zero values to 1)
  datI$trt <- as.integer(datI$trt > 0 )
  
  KM <- survfit(Surv(futime, event) ~ trt, data=datI, error="greenwood")
  KM.sum <- summary(KM)
  # Nelson-Aalen estimates
  na.0 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=0"]/KM.sum$n.risk[KM.sum$strata=="trt=0"])
  varna.0 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=0"]/KM.sum$n.risk[KM.sum$strata=="trt=0"]^2)
  na.0 <- na.0[length(na.0)]
  varna.0 <- varna.0[length(varna.0)]        
  na.1 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=1"]/KM.sum$n.risk[KM.sum$strata=="trt=1"])
  varna.1 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=1"]/KM.sum$n.risk[KM.sum$strata=="trt=1"]^2)
  na.1 <- na.1[length(na.1)]
  varna.1 <- varna.1[length(varna.1)]
  
  # survival estimates
  S.0 <- exp(-na.0)
  varS.0 <- ifelse(is.na(varna.0), NA, exp(-2*na.0) * varna.0)
  S.1 <- exp(-na.1)
  varS.1 <- ifelse(is.na(varna.1), NA, exp(-2*na.1) * varna.1)
  
  # cumulative incidence ratio
  F.0 <- 1 - S.0
  F.1 <- 1 - S.1
  FR <- F.1/F.0
  varlogFR <- ifelse(is.na(varS.0) | is.na(varS.1), NA, varS.1/(F.1^2) + varS.0/(F.0^2))
    
  FRci.up <- ifelse(is.na(varlogFR), NA, exp(log(FR)+qnorm(1-alpha/2)*sqrt(varlogFR)))
  bound <- ifelse(FRci.up < upperFR, "Eff", "nonEff")  
  
  return(list(bound=bound, VE=1-FR))
}

## This function calculates the probability of stopping for each bound for each arm
sumTrialOneTrtArm =  function (nTrtArms = 2,    ## number of vaccine arms
                      percentile.v = c(0, 0.25, 0.50, 0.75, 1),  ## vector of percentiles for infection counts
                      RDataFile,
                      dataDir = NULL  ## dir for simulated trials
                               ) {

## load in RData object (a list named 'out' )
load( file.path(dataDir, RDataFile) )
nTrials = length(out)


## create a matrix to store the prob. for each bound 
## 'power' =  Eff'+'HighEff'
## 'Total' = sum of all bounds to check if it is 1
trialRes = data.frame(matrix(0,nrow=nTrtArms, ncol=7))
names(trialRes) = c("Harm", "nonEffInterim", "nonEffFinal", "Eff", "HighEff", "Power", "Total")

## an object to store the distribution of stop infections, 
stopNRes = data.frame(matrix(NA,nrow=nTrtArms, ncol=length(percentile.v)))

## an object to store the distribution of total infections, 
stage1NRes = data.frame(matrix(NA,nrow=nTrtArms, ncol=length(percentile.v)))

## an object to store the distribution of post 6 month infections
#post6NRes = data.frame(matrix(NA,nrow=nTrtArms, ncol=length(percentile.v)))

# for each vaccine arm, the maximum number of infections across all simulated trials
stopNmax <- numeric(nTrtArms)
  
for (i in 1:nTrtArms) {
  
   ## get the number of infections at where the arm stops
   stopN = unlist( lapply( out, function(x) x[[i]]$stopInfectCnt ))
   stopNmax[i] <- max(stopN, na.rm=TRUE)
   
   ## get the number of total infections at the end of stage1 
   ## (i.e. all subjects complete the stage 1 follow up)
   stage1N = unlist( lapply( out, function(x) x[[i]]$totInfecCnt))
   
   ## get the bounds hit when the arm stops
   bounds <- unlist( lapply( out, function(x) x[[i]]$boundHit) )

   ## get the number of infections from month 6 to end of stage 1
   #post6N = unlist( lapply( out, function(x) x[[i]]$totInfecCntPost6))
   post6N = unlist( lapply( out, function(x) x[[i]]$infectPost6mo))
  
   
   trialRes$Harm[i] = sum(bounds=="Harm")/nTrials
   trialRes$nonEffInterim[i] = sum(bounds=="NonEffInterim")/nTrials
   trialRes$nonEffFinal[i] = sum(bounds=="NonEffFinal")/nTrials
   
   trialRes$Eff[i] = sum(bounds=="Eff")/nTrials
   trialRes$HighEff[i] = sum(bounds=="HighEff")/nTrials

   trialRes$Power[i] = trialRes$Eff[i] + trialRes$HighEff[i]
   trialRes$Total[i] = trialRes$Power[i] + trialRes$Harm[i] + trialRes$nonEffInterim[i] + trialRes$nonEffFinal[i]
   
   ## get the distribution of infections at where the trial stops
   stopNRes [i,] = quantile(stopN, prob=percentile.v)
   
   ## get the distribution of infections at the end of stage 1
   stage1NRes [i,] = quantile(stage1N, prob=percentile.v)
   
   ## get the distribtution of infections from month 6 to the end of stage1
  # post6NRes [i,] = quantile(post6N, prob=percentile.v)
  
   }
   stopNRes$countTypes = "infectionTrialStop"
   stage1NRes$countTypes = "infectionStage1"
   #post6NRes$countTypes = "post6mo"

   #res = rbind(stopNRes, stage1NRes, post6NRes)
   res = rbind(stopNRes, stage1NRes)
   names(res) = c(paste(percentile.v*100, "%", sep=""), "countTypes")

   list(trialRes = trialRes, NRes = res, armMaxNinf = which(stopNmax==max(stopNmax))[1])
}  


censTrial <- function(dataFile,
                      monitorFile,
                      stage1,
                      stage2,
                      saveDir = NULL,
                      verbose = TRUE){
                     
  if (!is.null(saveDir)){
    ## load the trial data in RData object (a list named 'trialObj' )
    load(file.path(saveDir, dataFile))
  } else {
    trialObj <- dataFile
    rm(dataFile)
  }  
  
  nTrials = length(trialObj[["trialData"]])  
  nTrtArms <- as.integer( trialObj$nArm - 1 )   ## one placebo arm
  
  if (!is.null(saveDir)){
    ## load the trial data in RData object (a list named 'trialObj' )
    load(file.path(saveDir, monitorFile))
  } else {
    out <- monitorFile
    rm(monitorFile)
  }    
  
  ## a matrix to store bounds results from monitoring each single arm
  boundsRes = matrix(NA, ncol=nTrtArms, nrow = nTrials)
  
  ## an object to store the arm stop time when the arm stops (since trial starts), 
  stopTime = matrix(NA, ncol=nTrtArms, nrow = nTrials)
  
  ## create a list to store the trial data that all subjects have been correctly
  ## censored based on the monitoring results from all arms
  trialListCensor <- vector("list", nTrials)
  
  for (i in 1:nTrtArms) {  
    ## get the bounds 
    boundsRes [, i]= unlist( lapply( out, function(x) x[[i]]$boundHit) )
    stopTime [, i]= unlist( lapply( out, function(x) x[[i]]$stopTime) )
  }
  
  ## go through each trial
  for (i in 1:nTrials ) {
    
    ## create a list to store censored data for ith trial
    #trialCensorI = vector("list",nTrtArms+1 )
    trialCensorI = NULL
    
    ## extract data for the i-th trial
    datI <- trialObj[["trialData"]][[ i ]]
    datI$futime <- datI$exit - datI$entry
    
    ## maximum possible trial duration  
    trialStop = max(datI$exit)
    
    ## get the placebo arm
    datI.0 = subset(datI, trt == 0)  
    
    ## if none of the arms are efficacy or highefficacy (i.e. all arms are either harm or noneff),
    ## the entire trial stops when the last arm hits noneff/harm
    if (!any(boundsRes[i,]%in%c("Eff", "HighEff"))) {
      
      ## get the time when the trial stops
      trialStop = max(stopTime[i,])  
      
      ## censor the placebo arm        
      datI.0$event = datI.0$event == 1 & (datI.0$exit <= trialStop)
      datI.0$exit = pmin( datI.0$exit, trialStop)
    }
    ## if at least one arm reaches efficacy, placebo arm will continue follow up until stage 2
    ## no need of extra action, store the censored placebo
    trialCensorI = rbind(trialCensorI, datI.0)
    
    ## Now we move to censor the trt arms
    for (j in 1:nTrtArms) {
      
      ## subset jth trt arm
      datI.j <- subset(datI, trt==j )
      
      ## censor the subjects based on bounds results for all arms
      ## first, get the stop time for trial i arm j
      t.j = stopTime[i, j]
      
      ## second, check if hit harm
      if (boundsRes[i, j] %in% "Harm") {  ## if "Harm", stop
        ## censor the trt arm        
        datI.j$event = datI.j$event == 1 & (datI.j$exit <= t.j)
        datI.j$exit = pmin( datI.j$exit, t.j)
        
      } else { 
        
        if (boundsRes[i, j] %in% c("NonEffInterim", "NonEffFinal")) {  ## if hit the non efficacy bound
          ## get the stop time for this arm, 
          ## which is the end of stage 1 or when the last arm hits noneff/harm if no eff. trt arms
          endStage1 = max(datI.j$entry + stage1)
          trialStop = min (c(trialStop, endStage1))
          
          ## censor the trt arm at 'trialStop'      
          datI.j$event = datI.j$event == 1 & (datI.j$exit <= trialStop)
          datI.j$exit = pmin( datI.j$exit, trialStop)                   
          
        } else { ## hit high eff or efficacy, i will remove this later
          ## the arm will follow up to stage 2, no need of action
        }      
      }
      
      ## store the censored trt arm
      #trialCensorI [[j+1]] = datI.j
      trialCensorI = rbind(trialCensorI, datI.j)                
    }      
    trialListCensor [[i]] = trialCensorI      
  }
  
  if (!is.null(saveDir)){
    ## output newly censored data 
    ## filename for censored trial data, it includes:  nTrtArms 1Vac, 2Vac, 3Vac,
    saveFile <- paste("trialDataCens", substr(monitorFile, 13, nchar(monitorFile)), sep="")
    ## save trial output and info on rates used
    save(trialListCensor, file=file.path(saveDir, saveFile) )
    if (verbose){ cat("Trial data with correct censoring saved in:\n", file.path(saveDir, saveFile), "\n\n") }
  } else {
    return(trialListCensor)
  }  
}


rankTrial <- function(censFile,
                      idxHighestVE,
                      headHead=NULL,
                      poolHead=NULL,
                      lowerVE,
                      stage1,
                      stage2,
                      alpha,
                      saveDir = NULL,
                      verbose = TRUE){                 
  
  if (!is.null(saveDir)){
    ## load the trial data in RData object (a list named 'trialObj' )
    load(file.path(saveDir, censFile))
  } else {
    trialListCensor <- censFile
    rm(censFile)
  }  
  
  nTrials = length(trialListCensor)  
  nTrtArms <- length(unique(trialListCensor[[1]]$trt)) - 1    ## one placebo arm
  
  if (!is.null(headHead)){
    if (NCOL(headHead)!=2){ stop("Number of columns in headHead must equal to 2.'\n") }
    
    ## a matrix to store power for head to head comparison of single arm
    ## column 1 = power for VE(0-18), column 2 = power for VE(0-36)
    headHeadPw = matrix(0, nrow=NROW(headHead), ncol=2)    
  }  
  
  if (!is.null(poolHead)){
    if (!(NCOL(poolHead) %in% 3:4)){ stop("Number of columns in poolHead must equal to 3 or 4.'\n") }
    
    ## a matrix to store power for comparison of pooled arms
    ## column 1 = power for VE(0-18), column 2 = power for VE(0-36)
    poolHeadPw = matrix(0, nrow=NROW(poolHead), ncol=2)    
    
    # power to correctly identify the best pooled vaccine arm
    poolRankSelectPw <- numeric(NROW(poolHead))
  }  
  
  # power to correctly identify the best vaccine arm
  rankSelectPw <- 0
  
  ## go through each trial
  for (i in 1:nTrials ) {
    ## extract data for the i-th trial
    datI <- trialListCensor[[i]] 
    
    ## a matrix to store estimated VE for each trt arm 
    ## column 1 = power for VE(0-18), column 2 = power for VE(0-36)
    VE.I = matrix(0, nrow=nTrtArms, ncol=2)
    
    ## cum hazard-based Wald test results for each trt arm
    ## column 1 = power for VE(0-18), column 2 = power for VE(0-36)
    test.I = matrix(0, nrow=nTrtArms, ncol=2)
    
    ## now calculate power for rank and select
    for (j in 1:nTrtArms){
      ## extract data for the j-th trt arm & placebo
      testData = subset(datI, trt %in% c(0, j))
      
      ## power for stage 1
      bnd = testVE(datI=testData, lowerVE=lowerVE, stage =stage1, alpha=alpha)
      VE.I[j, 1] = bnd$VE
      test.I[j, 1] = bnd$bound
      
      ## power for stage 2
      bnd = testVE(datI=testData, lowerVE=lowerVE, stage =stage2, alpha=alpha)
      VE.I[j, 2] = bnd$VE
      test.I[j, 2] = bnd$bound
    }
    
    # identify the vaccine arm with the highest estimated VE(0-36)
    bestInd = which(VE.I[,2] == max(VE.I[,2]))
      
    # to estimate the probability that the vaccine arm with the highest estimated VE(0-36) is the one with the
    # true highest VE(0-36), AND for that regimen the hypothesis H0: VE(0-18)<=0% was rejected    
    if (bestInd==idxHighestVE && test.I[bestInd, 1]=="Eff"){ rankSelectPw = rankSelectPw + 1 }    
    
    if (!is.null(headHead)){  
      # head-to-head comparison of individual vaccine arms at 18 and 36 months
      pw = headTestVE(datI, headHead, stage1, stage2, alpha)
      headHeadPw = headHeadPw + pw  # counts of detected relative efficacy for each comparison            
    }    
    
    if (!is.null(poolHead)) {
      # head-to-head comparison of pooled vaccine arms at 18 and 36 months
      pw = headTestVE(datI, poolHead, stage1, stage2, alpha)
      poolHeadPw = poolHeadPw + pw   
      
      # ranking and selection for pooled vaccine arms
      for (hi in 1:NROW(poolHead)) {
        arm1 = poolHead[hi,1:2]
        arm2 = poolHead[hi,3:NCOL(poolHead)]
        
        ## change trt index to combine arms listed in "arm1" and "arm2"
        ## and store the new data in datI.p
        datI$trt[datI$trt %in% arm1] = arm1[1]
        
        ## if arm2 is placebo, then no changes
        datI$trt[datI$trt %in% arm2] = arm2[1]
        
        ## now update 'bestVE' to arm1[1]
        bestVE = arm1[1]
        
        ## a matrix to store estimated VE for each trt arm 
        ## column 1 = power for VE(0-18), column 2 = power for VE(0-36)
        
        ## make it very small and impossible to be selected
        VE.I = matrix(-10, nrow=nTrtArms, ncol=2)
        test.I = matrix("empty", nrow=nTrtArms, ncol=2)
        
        ## get index for combined arms
        indTrt = sort(unique(datI$trt))
        indTrt = indTrt[-1]    ## remove placebo
        
        for (j in indTrt) {
          ## extract data for the j-th trt arm & placebo
          testData = subset(datI, trt %in% c(0, j)) 
          
          ## power for stage 1
          bnd = testVE(datI=testData, lowerVE=lowerVE, stage=stage1, alpha=alpha)
          VE.I[j, 1] = bnd$VE
          test.I[j, 1] = bnd$bound
          
          ## power for stage 2
          bnd = testVE(datI=testData, lowerVE=lowerVE, stage=stage2, alpha=alpha)
          VE.I[j, 2] = bnd$VE
          test.I[j, 2] = bnd$bound
        }
        
        # identify the pooled vaccine arm with the highest estimated VE(0-36)
        bestInd = which(VE.I[,2]==max(VE.I[,2]))
        
        # to estimate the probability that the pooled vaccine arm with the highest estimated VE(0-36) is the one 
        # with the true highest VE(0-36), AND for that regimen the hypothesis H0: VE(0-18)<=0% was rejected    
        if (bestInd==bestVE && test.I[bestInd, 1]=="Eff"){ poolRankSelectPw[hi] = poolRankSelectPw[hi] + 1 }
      }
    }
  }
  
  out <- list(rankSelectPw = rankSelectPw/nTrials)
  if (!is.null(headHead)){ out$headHeadPw <- headHeadPw/nTrials } 
  if (!is.null(poolHead)){ 
    out$poolRankSelectPw <- poolRankSelectPw/nTrials
    out$poolHeadPw <- poolHeadPw/nTrials    
  }
  
  if (!is.null(saveDir)){
    saveFile <- paste0("rankSelectPwr",substr(censFile, 14, nchar(censFile)))
    save(out, file = file.path(saveDir, saveFile) )
    if (verbose){ cat("Output saved in:\n", file.path(saveDir, saveFile), "\n\n") }
  } else {
    return(out)
  }  
}

  
###################################################
## This function calculates the power for head to head comparison
headHeadTest = function(datI, headHeadInd, stage1, stage2, nullHR, alpha=0.025) {
  ## a matrix to store power for head to head comparison 
  ## column 1 = power for VE(0-18), column 2 = power for VE(0-36)
  headHeadPw = matrix(0, nrow=nrow(headHeadInd), ncol=2)
  
  for (h in 1:nrow(headHeadInd)) {
    if(ncol(headHeadInd)==2) {  ## 2 columns
      arm1 = headHeadInd[h,1]
      arm2 = headHeadInd[h,2]
    } else {                    ## 3 or 4 column
      arm1 = headHeadInd[h,1:2]
      arm2 = headHeadInd[h,3:ncol(headHeadInd)]
    }
       
    testData = subset(datI, trt %in% c(arm1, arm2))
    
    ## need to change the index in arm2 to 0
    testData$trt [testData$trt %in% arm2] = 0 
    
    ## power for stage 1
    bnd = finalLogRankTest(datI=testData, stage =stage1, NullHR=nullHR, alpha=alpha)
    if (bnd$bound =="Eff") headHeadPw[h, 1] = 1           
    
    ## power for stage 2
    bnd = finalLogRankTest(datI=testData, stage =stage2, NullHR=nullHR, alpha=alpha)
    if (bnd$bound =="Eff") headHeadPw[h, 2] = 1 
  }
  headHeadPw
}

headTestVE <- function(datI, headHeadInd, stage1, stage2, alpha){
  ## a matrix to store power for head to head comparison 
  ## column 1 = power for VE(0-18), column 2 = power for VE(0-36)
  headHeadPw = matrix(0, nrow=NROW(headHeadInd), ncol=2)
  
  for (h in 1:NROW(headHeadInd)) {
    if(NCOL(headHeadInd)==2) {  ## 2 columns
      arm1 = headHeadInd[h,1]
      arm2 = headHeadInd[h,2]
    } else {                    ## 3 or 4 column
      arm1 = headHeadInd[h,1:2]
      arm2 = headHeadInd[h,3:NCOL(headHeadInd)]
    }
    
    testDataRel <- subset(datI, trt %in% c(arm1, arm2))
    ## need to change the index in arm2 to 0
    testDataRel$trt [testDataRel$trt %in% arm2] = 0    
    
    testData <- subset(datI, trt %in% c(0, arm1))
        
    # head-to-head power = the 1-sided head-to-head test rejects; AND the superior arm passes 
    # Stage 1 with VE(0-stage1) > 0%
    ## relative VE(0-stage1)    
    bndRel <- testVE(datI=testDataRel, lowerVE=0, stage=stage1, alpha=alpha)
    ## VE(0-stage1) of the superior treatment
    bnd <- testVE(datI=testData, lowerVE=0, stage=stage1, alpha=alpha)    
    if (bndRel$bound =="Eff" & bnd$bound=="Eff"){ headHeadPw[h, 1] = 1 }
    
    ## power for stage 2
    bndRel <- testVE(datI=testDataRel, lowerVE=0, stage=stage2, alpha=alpha)
    bnd <- testVE(datI=testData, lowerVE=0, stage=stage2, alpha=alpha)    
    if (bndRel$bound =="Eff" & bnd$bound=="Eff"){ headHeadPw[h, 2] = 1 }
  }
  return(headHeadPw)
}


#########
buildBounds = function(nInfec, highEffBounds) {
   
   ## high efficacy : approxiamtion 
   idx = which(highEffBounds$N==nInfec)
   
   ## get number of row
   nBounds = nrow (highEffBounds)
   if (length(idx)==0) {  ## no bound for nInfec, add it
      highEffBounds  = rbind(highEffBounds, c(nInfec, highEffBounds[nBounds,2]))   ## approximation  
      
      ## order the bound by infections
      highEffBounds = highEffBounds[order(highEffBounds$N), ]
      
      idx = which(highEffBounds$N==nInfec)
      if(idx>1 && (idx+1)<=nrow (highEffBounds))       ## approximate the bound by averaging over 2 adjacent twos
         highEffBounds[idx,2] = mean(c(highEffBounds[idx-1,2], highEffBounds[idx+1,2]))
    }
    highEffBounds = subset(highEffBounds, N<=nInfec)
 
   list(highEffBounds=highEffBounds)            
}

VEpowerPP <- function(dataList, lowerVEuncPower, alphaUncPower, VEcutoffWeek, stage1, outName=NULL, saveDir=NULL, verbose=TRUE){
  upperFRuncPower <- 1-lowerVEuncPower
  # output list (for each censTrial output object) of lists with components 'VE' and 'VEpwPP'
  pwList <- as.list(NULL)
  for (k in 1:length(dataList)){
    if (!is.null(saveDir)){
      # assumes 'dataList[[k]]' is a character string
      load(file.path(saveDir, dataList[[k]]))
      dataList[[k]] <- trialListCensor
      rm(trialListCensor)
    }
    
    nTrials <- length(dataList[[k]])
    nPPcohorts <- length(grep("pp", colnames(dataList[[k]][[1]])))
    if (nPPcohorts==0){ stop("Missing per-protocol cohort indicator in the data-set.") }
    ppnames <- paste0("pp", 1:nPPcohorts)
    
    VE <- matrix(NA, nrow=nTrials, ncol=nPPcohorts)
    VEpwPP <- matrix(NA, nrow=nTrials, ncol=nPPcohorts)
    for (i in 1:nTrials){
      dataI <- dataList[[k]][[i]]
      dataI$futime <- dataI$exit - dataI$entry
      # count only infections between VEcutoffWeek and the end of Stage 1
      dataI$event <- dataI$event==1 & dataI$futime>VEcutoffWeek & dataI$futime<=stage1
      # censor everyone at the end of Stage 1
      dataI$futime <- pmin(dataI$futime, stage1)

      # AN ALTERNATIVE WAY TO COMPUTE PER-PROTOCOL VE      
#       # restrict to subjects with follow-up time exceeding 'VEcutoffWeek' weeks (per-protocol criterion 1)
#       dataI <- subset(dataI, futime > VEcutoffWeek)
#       # censor all subjects at the Week 'stage1' visit
#       dataI$event <- dataI$event == 1 & (dataI$futime <= stage1)
#       dataI$futime <- pmin(dataI$futime, stage1)
#       # shift the time origin to the Week 'VEcutoffWeek' visit
#       dataI$futime <- dataI$futime - VEcutoffWeek      
      
      for (j in 1:nPPcohorts){
        # restrict to subjects with non-missing vaccinations (per-protocol criterion 2)
        dataIj <- dataI[dataI[[ppnames[j]]]==1,]
        # now we have the final per-protocol data-set at the end of 'stage 1'
        # get the group sizes in the PP cohort
        nAll <- table(dataIj$trt)
        if (length(nAll)<2){
          next
        } else {
          # get the infection counts in the PP cohort
          nInf <- table(dataIj$trt[dataIj$event==1])
          if (length(nInf)<2){
            next
          } else {
            KM <- survfit(Surv(futime, event) ~ trt, data=dataIj, error="greenwood")
            KM.sum <- summary(KM)
            
            # Nelson-Aalen estimates
            na.0 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=0"]/KM.sum$n.risk[KM.sum$strata=="trt=0"])
            varna.0 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=0"]/KM.sum$n.risk[KM.sum$strata=="trt=0"]^2)
            na.0 <- na.0[length(na.0)]
            varna.0 <- varna.0[length(varna.0)]
            
            na.1 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=1"]/KM.sum$n.risk[KM.sum$strata=="trt=1"])
            varna.1 <- cumsum(KM.sum$n.event[KM.sum$strata=="trt=1"]/KM.sum$n.risk[KM.sum$strata=="trt=1"]^2)
            na.1 <- na.1[length(na.1)]
            varna.1 <- varna.1[length(varna.1)]
            
            # survival estimates
            S.0 <- exp(-na.0)
            varS.0 <- ifelse(is.na(varna.0), NA, exp(-2*na.0) * varna.0)
            S.1 <- exp(-na.1)
            varS.1 <- ifelse(is.na(varna.1), NA, exp(-2*na.1) * varna.1)
            
            # cumulative incidence ratio
            F.0 <- 1 - S.0
            F.1 <- 1 - S.1
            FR.i <- F.1/F.0
            varlogFR <- ifelse(is.na(varS.0) | is.na(varS.1), NA, varS.1/(F.1^2) + varS.0/(F.0^2))
            
            # VE(6.5-stage1) estimate
            VE[i,j] <- 1 - FR.i
            
            # 1-sided test of null hypothesis VE(6.5-stage1)<=lowerVEuncPower
            FRci.up <- ifelse(is.na(varlogFR), NA, exp(log(FR.i)+qnorm(1-alphaUncPower/2)*sqrt(varlogFR)))
            # if the upper bound of the CI for the cumulative incidence ratio lies below 'upperFRuncPower', reject H0
            VEpwPP[i,j] <- ifelse(FRci.up < upperFRuncPower, 1, 0)            
          }          
        }        
      }
    }
    VE <- apply(VE, 2, mean, na.rm=TRUE)
    VEpwPP <- apply(VEpwPP, 2, mean, na.rm=TRUE)
    pwList[[k]] <- list(VE=VE, VEpwPP=VEpwPP)    
  }
  if (!is.null(saveDir)){
    saveFile <- ifelse(is.null(outName), "VEpwPP.RData", outName)
    save(pwList, file=file.path(saveDir, saveFile))
    if (verbose){ cat("Output saved in:\n", file.path(saveDir, saveFile), "\n\n") }
  } else {
    return(pwList)
  }  
}
