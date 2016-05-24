
# Enhancing the Documentation Slot

# Series:
#    @.Data
#    @ positions
#    @ format
#    @ FinCenter
#    @ units  
#    @ recordIDs
#    @ title
#    @ documentation
#       attributes(@documentation, "Attributes")


# inspect the Information use
#   slotNames(object)
#   slot(object, name)


# Load Library:
require(timeSeries)


###############################################################################
# Data

obj1 <- timeSeries(rnorm(12), timeCalendar())
getAttributes(obj1)
setAttributes(obj1) <- list(series=series(obj1)[1:6, , drop=FALSE])
getAttributes(obj1)


obj2 <- timeSeries(rnorm(12), timeCalendar())
getAttributes(obj2)
setAttributes(obj2) <- list(series=as.matrix(obj2)[7:12, , drop=FALSE])
getAttributes(obj2)


###############################################################################
# Base Functions:


# base-apply.R
getAttributes( apply(obj1, 1, mean) )                                # ok


# base-applySeries.R
# should be deprecated, use generic apply() and aggregate() functions


# base-cbind.R
cbind(obj1, obj2)
getAttributes( cbind(obj1, obj2) )                                   # ok 
getAttributes( rbind(obj1, obj2) )                                   # ok 


# ... more base-cbind.R
# cbind ...
getAttributes(cbind(obj1, obj2))                                     # ok
getAttributes(cbind(obj1, as.matrix(obj2)))                          # ok
getAttributes(cbind(as.matrix(obj1), obj2))                          # ok
getAttributes(cbind(obj1))                                           # ok
# rbind ...
getAttributes(rbind(obj1, obj2))                                     # ok
getAttributes(rbind(obj1, as.matrix(obj2)))                          # ok
getAttributes(rbind(as.matrix(obj1), obj2))                          # ok
getAttributes(rbind(obj1))                                           # ok


# base-diff.R
getAttributes( diff(obj1) )                                          # ok


# base-merge.R
getAttributes(  merge(obj1, obj2) )                                  # ok


# base-rank.R
getAttributes( rank(obj1) )                                          # ok


# base-rev.R
getAttributes( rev(obj1) )                                           # ok


# base-sample.R
getAttributes( sample(obj1) )                                        # ok


# base-scale.R
getAttributes( scale(obj1) )                                         # ok


# base-sort.R
getAttributes( sort(obj1) )                                          # ok


################################################################################
# Subsetting:


# base-subsetting.R
#  .subset_timeSeries
#  .findIndex
#  $,timeSeries              Subsets a time series by column names
#  $<-,timeSeries            Replaces subset by column names
#  [,timeSeries              Subsets a time series object
#  [<-,timeSeries            Assigns value to subsets of a time series


# Should work by dafault ...
getAttributes( obj1[3:4, 1] )                                        # ok
getAttributes( head(obj1) )                                          # ok 
getAttributes( tail(obj1) )                                          # ok
                                                 

################################################################################
# Methods:


# methods-mathOps.R
# here the multiplications "*", works also with "+", "=", "/". ...
getAttributes( obj1 * 2)                                             # ok
getAttributes( obj1 * (1:12) )                                       # ok
getAttributes( obj1 * matrix(1:12, ncol=1) )                         # ok
getAttributes( obj1 * as.ts(1:12) )                                  # ok
getAttributes( obj1 * obj2 )                                         # ok ???
getAttributes( 2 * obj2 )                                            # ok
getAttributes( (1:12) * obj2 )                                       # ok
getAttributes( matrix(1:12, ncol=1) * obj2)                          # ok
getAttributes( as.ts(1:12) * obj2)                                   # ok
getAttributes( obj2 * obj1)                                          # ok ???


# More Math Functions
getAttributes( abs(obj1) )                                           # ok  
getAttributes( exp(obj1) )                                           # ok
getAttributes( obj1^2 )                                              # ok
# ...


# Round and Truncate:
getAttributes( round(obj1, digits=2) )                               # ok
getAttributes( trunc(obj1, digits=2) )                               # ok
getAttributes( signif(obj1, digits=3) )                              # ok
getAttributes( ceiling(100*obj1) )                                   # ok
getAttributes( floor(100*obj1) )                                     # ok


################################################################################
# Financial 'timeSeries' Functions


# fin-align.R
getAttributes( align(obj1) )                                         # ok


# fin-cumulated.R
getAttributes( cumulated(obj1) )                                     # ok


# fin-daily.R
# align(obj1) and alignDailySeries(obj1) are the same
# deprecate align Daily Series
getAttributes( alignDailySeries(obj1) )                              # ok
# Can we use the generic function aggregate ?
getAttributes( rollDailySeries(obj1, FUN=mean) )                     # ok


# fin-drawdowns.R
getAttributes( drawdowns(obj1) )                                     # ok


# fin-durations.R
getAttributes( durations(obj1) )                                     # ok


# fin-monthly.R
getAttributes( rollMonthlySeries(obj1, "3m", FUN=mean) )             # ok
getAttributes( countMonthlyRecords(obj1) )                           # ok


# fin-periodical.R                                                   # todo
.endOfPeriodSeries
.endOfPeriodStats
.endOfPeriodBenchmarks


# fin-returns.R
OBJ1 <- cumulated(obj1)
getAttributes( returns(OBJ1) )                                       # ok


# fin-runlengths.R
getAttributes( runlengths(obj1) )                                    # ok


# fin-splits.R
getAttributes( outlier(obj1) )                                       # ok


# fin-spreads.R
SPREADS <- spreads(obj3, which=c(1, 2))
getAttributes( SPREADS)                                              # ok fails
MIDQUOTES <- midquotes(obj3, which=c(1,2))
getAttributes( MIDQUOTES )                                           # ok fails


# fin-turns.R
INDEX <- cumulated(obj1) 
getAttributes( turns(INDEX) )                                        # ok


################################################################################
# Statistics timeSeries  Functions                                 


# statistics-colCumsums.R
getAttributes( colCumsums(obj1) )                                    # ok


# statistics-colSums.R
# returns no timeSeries objects                  


# statistics-orderColnames.R
# returns no timeSeries objects


# statistics-orderStatistics.R
# returns no timeSeries objects


# statistics-rollMean.R
getAttributes( rollStats(obj1, k=1, FUN=mean) )                      # ok
getAttributes( rollMean(obj1, k=1) )                                 # ok
getAttributes( rollMin(obj1, k=1) )                        # FAILS
getAttributes( rollMax(obj1, k=1) )                        # FAILS
getAttributes( rollMedian(obj1, k=1) )                               # ok


# statistics-rowCumsums.R


# statistics-smoothLowess.R
getAttributes( smoothLowess(obj1) )                                  # ok
getAttributes( smoothSupsmu(obj1) )                                  # ok
getAttributes( smoothSpline(obj1) )                                  # ok


################################################################################
# stats


# stats-aggregate.R
by1 <- time(obj1[3*(1:4),])
getAttributes( aggregate(obj1, by=by1, FUN=mean) )                   # ok


# stats-filter.R
getAttributes( filter(obj1, filter=c(1,1)) )                         # ok


# stats-lag.R
getAttributes( lag(obj1) )                                           # ok


# stats-na.contiguous.R
# returns no timeSeries objects


# stats-na.omit.R
obj3 <- obj1; obj3[4, 1] <- NA; obj3
getAttributes( na.omit(obj3) )                                       # ok
# What about? - They should be deprecated.
#   removeNA
#   substituteNA
#   interp NA 


# stats-window.R
Time <- time(obj1)
getAttributes( window(obj1, Time[3], Time[6]) )                      # ok


################################################################################
# Attributes Functions


getAttributes <- 
  function (obj) 
{
  # FUNCTION:
    
  # Check Argument:
  stopifnot(class(obj) == "timeSeries")
  
  # Extract Attributes:
  ans <- attr(obj@documentation, "Attributes")
    
  # return Value:
  ans
}


# -----------------------------------------------------------------------------


`setAttributes<-` <- 
  function(obj, value)
{
  # Example:
  #   obj <- dummySeries(); getAttributes(obj)
  #   setAttributes(obj) <- list(mat=matrix(1:4, ncol=2)); getAttributes(obj)
  #   getAttributes(obj)$mat[[1]]
    
  # FUNCTION:
    
  # Check Arguments:
  stopifnot(class(obj) == "timeSeries")
  stopifnot(is.list(value))
  stopifnot(length(value) == 1)
  stopifnot(!is.null(value))
  
  # Compose New Attribute:
  name <- names(value)
  names(value) <- NULL
  A <- list(value)
  names(A) <- name
  # print(A)

  # Get Already Existing Attribute
  B <- getAttributes(obj)
  if(is.null(B)) B <- list()
  # print(B)

  # Join Attributes:
  JOINED <- sapply(unique(c(names(A), names(B))), 
    function(x) list(c(A[[x]], B[[x]])))
  # print(JOINED)

  # Assign Attribute:
  attr(obj@documentation, "Attributes") <- JOINED
  
  # Return Value:
  obj 
}


###############################################################################

