################################################################################
### Initialization and other basic methods for the S4 class "sts"
###
### Copyright (C) 2007-2014 Michael Hoehle, 2012-2016 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################


######################################################################
# initialize-method -- see ../man/sts-class.Rd for class information
######################################################################

#Ensure that all matrix slots have the same dimnames, which are
#always taken from the observed matrix
fix.dimnames <- function(x) {
  dn <- dimnames(x@observed)
  #Make sure all arrays have the same dimnames
  dimnames(x@alarm) <- dimnames(x@state) <- dimnames(x@upperbound) <-
    dimnames(x@populationFrac) <- dn
  #Special for neighbourhood
  dimnames(x@neighbourhood) <- dn[c(2L,2L)]
  return(x)
}

## a user-level constructor function,
## which calls the standard generator function .sts(),
## which calls initialize() on the "sts" prototype - see init.sts() below
## NOTE: using sts() is the preferred approach since surveillance 1.10-0
## NOTE: NULL arguments are ignored => default slot values
sts <- function (observed,
                 start = c(2000, 1), frequency = 52, # prototype values
                 population = NULL, # an alias for "populationFrac"
                 ...) # further named arguments representing "sts" slots
{
    slots <- list(observed = observed, start = start, freq = frequency, ...)
    if (!is.null(population)) {
        if ("populationFrac" %in% names(slots))
            warning("'population' takes precedence over 'populationFrac'")
        slots$populationFrac <- population
    } # else "populationFrac" is a possible element of ...
    
    ## call the standard generator function with explicitly set slots
    isNULL <- vapply(X = slots, FUN = is.null,
                     FUN.VALUE = FALSE, USE.NAMES = FALSE)
    do.call(.sts, slots[!isNULL])
}

## initialize-method called by new("sts", ...),
## the long-standing default way of creating "sts" objects.
## For backward-compatibility, we keep this customized initialize-method,
## although it would be cleaner to put things into the generator function
## and use the default initialize-method.
init.sts <- function(.Object, ..., # also for slots of classes extending "sts"
                     observed, # use copy constructor if missing(observed)
                     ## the following default arguments depend on dim(observed)
                     epoch = seq_len(nTime),
                     state = matrix(FALSE, nTime, nUnit),
                     alarm = matrix(NA, nTime, nUnit),
                     upperbound = matrix(NA_real_, nTime, nUnit),
                     neighbourhood = matrix(NA, nUnit, nUnit),
                     populationFrac = matrix(1/nUnit, nTime, nUnit),
                     ## FIXME: change default to a matrix of NA_real_ ?
                     ## the map slot needs special treatment (see below)
                     map = .Object@map # old/prototype value
                     ## the remaining slots have useful prototype values
                     ## and are handled as part of ...
                     ##start = c(2000, 1), freq = 52,
                     ##epochAsDate = FALSE, multinomialTS = FALSE,
                     ##control = .Object@control
                     )
{
    if (nargs() < 2) # nothing to do
        return(.Object)
    
    if (missing(observed)) { # use default initialize-method
        ## such that, e.g., initialize(stsObj, map=newMap) will set a new map
        ## and copy other slots from stsObj instead of (re-)setting to defaults,
        ## as well as to support new("stsBP", stsObj, ci=ci, lambda=lambda).
        ## CAVE: automatic dimension correction of matrix slots is not done.
        .Object <- callNextMethod()
        ## Drawback: .Object@map has been forced to "SpatialPolygons"
        if (!missing(map)) # restore the supplied map
            .Object@map <- map
        ## If missing(map), .Object@map = as(stsObj@map, "SpatialPolygons"),
        ## i.e., data will be lost => map=stsObj@map must be passed explicitly
        .Object <- fix.dimnames(.Object)
        return(.Object)
    }
    
    ## Ensure matrix form (auto-conversion is useful for single time series)
    observed <- as.matrix(observed)
    nUnit <- ncol(observed)
    nTime <- nrow(observed)
    state <- as.matrix(state)
    alarm <- as.matrix(alarm)
    upperbound <- as.matrix(upperbound)
    
    ## clear rownames and set colnames for the matrix of observed counts
    if (is.null(namesObs <- colnames(observed))){
        namesObs <- paste0("observed", seq_len(nUnit))
    }
    dimnames(observed) <- list(NULL, namesObs)
    
    ## if there is only one state-vector for more than one area, repeat it
    if (nUnit > 1 && ncol(state) == 1 && length(state) == nTime) {
        state <- rep.int(state, nUnit)
        dim(state) <- c(nTime, nUnit)
    }

    ## time-constant population fractions can be provided as a single vector
    if (is.vector(populationFrac, mode="numeric") &&
        length(populationFrac) == nUnit) {
        populationFrac <- matrix(populationFrac, nTime, nUnit, byrow=TRUE)
    }

    ## we need to set the map manually since the initialize,ANY-method called
    ## next would coerce a "SpatialPolygonsDataFrame" to "SpatialPolygons"
    if (!missing(map))
        .Object@map <- map
    
    ## set all other slots (including for classes extending this class)
    ## using the default initialize-method
    .Object <- callNextMethod(.Object, ...,
                              observed=observed, epoch=epoch,
                              state=state, alarm=alarm, upperbound=upperbound,
                              neighbourhood=neighbourhood,
                              populationFrac=populationFrac)
    ## this also checks validObject(.Object)
    
    ## make sure all arrays have the same dimnames
    .Object <- fix.dimnames(.Object)
    
    return(.Object)
}

setMethod("initialize", "sts", init.sts)


###########################################################################
# Conversion between old "disProg" and new "sts" classes
###########################################################################

## transform a "disProg" object to the new "sts" class
disProg2sts <- function(disProgObj, map=NULL) {
  disProgObj$map <- map
  ## NOTE: we cannot trust disProgObj$week to be a valid "epoch" specification,
  ## e.g., the week in data("ha") refers to the week number _within_ a year.
  ## CAVE: in "disProg" objects, several elements may be missing or NULL,
  ## and there could be further elements not matching any "sts" slot,
  ## e.g., in "disProg" objects generated by sim.pointSource()
  validElements <- names(disProgObj) %in% slotNames("sts") &
      !vapply(X=disProgObj, FUN=is.null, FUN.VALUE=FALSE, USE.NAMES=FALSE)
  ## initialize an "sts" object using the valid "disProg" elements
  stsObj <- do.call(.sts, disProgObj[validElements])
  return(stsObj)
}

## The reverse action
sts2disProg <- function(sts) {
  disProgObj <- create.disProg(week=sts@epoch, start=sts@start, freq=sts@freq,
                               observed=sts@observed, state=sts@state, neighbourhood=sts@neighbourhood,
                               populationFrac=sts@populationFrac, epochAsDate=sts@epochAsDate)
  #For survRes: alarm=sts@alarm, upperbound=sts@upperbound)
  return(disProgObj)
}


###########################################################################
#Method to aggregate over all units, either the time series is aggregated
#so a new sampling frequency of nfreq units per time slot is obtained.
#The other alternative is to aggregate all units.
#
# Note: The function is not 100% consistent with what the generic
#       aggregate does. 
#
# Warning: In case the aggregation is by unit the upperbound slot is set
#          to NA. Furthermore the MAP object is left as.is, but
#          the object cannot be plotted anymore.
#
# Params:
#   by - a string being either "time" or "unit"
#   nfreq - new sampling frequency if by=="time". If "all" then all
#           time instances are summed.
###########################################################################

setMethod("aggregate", signature(x="sts"), function(x,by="time",nfreq="all",...) {
  
 #Action of aggregation for populationFrac depends on the type 
 binaryTS <- sum( x@populationFrac > 1 ) > 1  # FIXME: x@multinomialTS?

  #Aggregate time
  if (by == "time") {
    if (nfreq == "all") {
      howmany <- dim(x@observed)[1]
    } else if (nfreq == x@freq) { # nothing to do
      return(x)
    } else { # nfreq != x@freq
      howmany <- x@freq / nfreq
      if (howmany - ceiling(howmany) != 0)
          stop("nfreq has to be a multiple of x@freq.")
    }
    
    n <- dim(x@observed)[1]
    m <- ceiling(n/howmany)
    new <- rep(1:m,each=howmany)[1:n]
    x@freq <- ifelse(nfreq == "all", howmany, nfreq)
    x@epoch <- 1:m
    
    x@observed <- as.matrix(aggregate(x@observed,by=list(new),sum)[,-1])
    x@state <- as.matrix(aggregate(x@state,by=list(new),sum)[,-1])>0
    x@alarm <- as.matrix(aggregate(x@alarm,by=list(new),sum)[,-1]) # number of alarms
    x@upperbound <- as.matrix(aggregate(x@upperbound,by=list(new),sum)[,-1])
    x@populationFrac <- as.matrix(aggregate(x@populationFrac,by=list(new),sum)[,-1])
    ## FIXME: should make clear (warn?) that population is summed over time

    #the population fractions need to be recomputed if not a binary ts
    if (!binaryTS) {
      sums <- matrix(rep(apply(x@populationFrac,1,sum),times=ncol(x)),ncol=ncol(x))
      x@populationFrac <-x@populationFrac/sums
    }
  }
  if (by == "unit") {
    #Aggregate units
    x@observed <- as.matrix(apply(x@observed, MARGIN=1, sum))
    x@state <- as.matrix(apply(x@state, MARGIN=1, sum))>0
    x@alarm <- as.matrix(apply(x@alarm, MARGIN=1, sum))>0 # contrary to counting for by="time"!
    #There is no clever way to aggregate the upperbounds
    x@upperbound <- matrix(NA_real_,ncol=ncol(x@alarm),nrow=nrow(x@alarm))
    x@populationFrac <- as.matrix(apply(x@populationFrac, MARGIN=1, sum))#>0
    x@neighbourhood <- matrix(NA, 1, 1) # consistent with default for new("sts")
    ## FIXME: x@map will be invalid, remove or unionSpatialPolygons()?
  }

  #validObject(x) #just a check

  return(x)
})
  

#####################################################################
# Miscellaneous access methods
####################################################################

setMethod("dim", "sts", function (x) dim(x@observed))
setMethod("dimnames", "sts", function (x) dimnames(x@observed))

#Extract which observation within year we have
setMethod("epochInYear", "sts", function(x,...) {
  #Strptime format strings available as:
  #http://www.opengroup.org/onlinepubs/009695399/functions/strptime.html
  if (x@epochAsDate) {
    epochStr <- switch( as.character(x@freq), "12" = "%m","52" =  "%V","365" = "%j")
    return(as.numeric(formatDate(epoch(x),epochStr)))
  } else {
    return( (x@epoch-1 + x@start[2]-1) %% x@freq + 1)
  }
})

#Extract the corresponding year for each observation using
setMethod("year", "sts", function(x,...) {
  if (x@epochAsDate) {
    return(as.numeric(formatDate(epoch(x),"%G")))
  } else {
    ((x@epoch-1 + x@start[2]-1) + (x@freq*x@start[1])) %/% x@freq 
  }
})


#####################################################################
#[-method for accessing the observed, alarm, etc. objects
#####################################################################

setMethod("[", "sts", function(x, i, j, ..., drop) {
  #default value for i and j
  if(missing(i)) {i <- min(1,nrow(x@observed)):nrow(x@observed)}
  if(missing(j)) {j <- min(1,ncol(x@observed)):ncol(x@observed)}

  x@epoch <- x@epoch[i]
  x@observed <- x@observed[i,j,drop=FALSE]
  x@state <- x@state[i,j,drop=FALSE]
  x@alarm <- x@alarm[i,j,drop=FALSE]

  x@populationFrac <- x@populationFrac[i,j,drop=FALSE]
  #If not binary TS the populationFrac is normed
  binaryTS <- sum( x@populationFrac > 1 ) > 1 # FIXME @ Michael: x@multinomialTS
  if (!binaryTS) {
    x@populationFrac <- x@populationFrac / apply(x@populationFrac,MARGIN=1,sum)
   }  
  x@upperbound <- x@upperbound[i,j,drop=FALSE]

  #Neighbourhood matrix
  if (ncol(x@observed) != ncol(x@neighbourhood) &&  # selected units
      !all(x@neighbourhood %in% c(NA,0,1))) { # no adjacency matrix
      message("Note: selection of units could invalidate the 'neighbourhood'")
      ## e.g., if 'neighbourhood' specifies neighbourhood orders
  }
  x@neighbourhood <- x@neighbourhood[j,j,drop=FALSE]

  #Fix the corresponding start entry. it can either be a vector of
  #logicals or a specific index. Needs to work in both cases.
  #Note: This code does not work if we have week 53s!
  if (is.logical(i)) {
    i.min <- which.max(i) #first TRUE entry
  } else {
    i.min <- min(i)
  }
  start <- x@start
  new.sampleNo <- start[2] + i.min - 1
  start.year <- start[1] + (new.sampleNo - 1) %/% x@freq 
  start.sampleNo <- (new.sampleNo - 1) %% x@freq + 1
  x@start <- c(start.year,start.sampleNo)
  ## If !epochAsDate, we also have to update epoch since it is relative to start
  if (!x@epochAsDate) x@epoch <- x@epoch - i.min + 1
                                                     
  ## Note: We do not automatically subset the map according to j, since
  ##       identical(row.names(map), colnames(observed))
  ##       is not a property of the sts-class; Unmonitored regions are allowed.
  
  #Done
  return(x)
})


#########################################################################
## Plot method ... the type argument specifies what type of plot to make
##
## plot as multivariate time series:  type = observed ~ time | unit 
## plot as map object aggregated over time: type = observed ~ 1 | unit
## new map implementation via: type = observed ~ unit
## the specific plot functions are in separate files (stsplot_*.R)
########################################################################

setMethod("plot", signature(x="sts", y="missing"),
          function (x, type = observed ~ time | unit, ...) {
  
  # catch new implementation of time-aggregate map plot
  if (isTRUE(all.equal(observed ~ unit, type)))
      return(stsplot_space(x, ...))

  #Valid formula?
  valid <- lapply(as.list(type[[3]]), function(i)
                  is.na(pmatch(i,c("1","unit","|","time","*","+"))))
  valid <- all(!unlist(valid))
  obsOk <- (type[[2]] == "observed")
  alarmOk <- (type[[2]] == "alarm")
  if (!valid || !(obsOk | alarmOk))
      stop("Not a valid plot type")

  #Parse the formula, i.e. extract components
  map   <- (length(type[[3]])==3) && (type[[3]][[1]] == "|") && (type[[3]][[2]] == "1")
  time  <- pmatch("time",type[[3]]) > 0
  #All-in-one if type=time+unit -> no, use argument "as.one" for stsplot_time
  #as.one <- all(!is.na(pmatch(c("time","unit"),type[[3]] ))) && is.na(pmatch("|",type[[3]]))

  #No unit dimension?
  justTime <- type[[3]] == "time"

  #space-time plots
  if (map) {
    stsplot_spacetime(x, type, ...)
    return(invisible())
  }
  #time plots
  if (time) {
    if (obsOk) {
      #In case observed ~ time, the units are aggregated
      stsplot_time(if(justTime) aggregate(x,by="unit") else x, ...)
      return(invisible())
    }
    if (alarmOk) {
      stsplot_alarm(x, ...)
      return(invisible())
    }
  }
})


## define how "sts" objects get printed

setMethod( "show", "sts", function( object ){
  cat( "-- An object of class ", class(object), " -- \n", sep = "" )
  if (!object@epochAsDate) {
    cat( "freq:\t\t", object@freq,"\n" )
  } else {
    epochStr <- switch( as.character(object@freq), "12" = "%m","52" =  "%V","365" = "%j")
    cat( "freq:\t\t", paste(object@freq," with strptime format string ",epochStr,"\n",sep=""))
  }
  if (!object@epochAsDate) {
    cat( "start:\t\t",object@start,"\n" )
  } else {
    cat( "start:\t\t",paste(epoch(object)[1]),"\n" )
  }
  cat( "dim(observed):\t", dim(object@observed), "\n\n")

  n <- 1
  cat("Head of observed:\n")
  print(head(object@observed,n))

  if (npoly <- length(object@map)) {
      cat("\nmap:\n")
      print(modifyList(summary(object@map), list(data=NULL))) # no data summary
      cat("Features    :", npoly, "\n")
      if (inherits(object@map, "SpatialPolygonsDataFrame"))
          cat("Data slot   :", ncol(object@map), "variables\n")
  }

  if (ncol(object@observed) > 1) {
      cat("\nhead of neighbourhood:\n")
      print( head(object@neighbourhood,n))
  }
} )
