### ============================================================================
### Check events data set 
### Changes version 1.11: event can be an R-function, even if DLL model
###                       continueeroot: to continue even if a root is found
### ============================================================================

checkevents <- function (events, times, vars, dllname, root = FALSE) {

  if (is.null(events)) return(list())
  if (is.null(events$data) && is.null(events$func) && 
      is.null(events$terminalroot)) return(list())

  funevent <- events$func

  if (root) {  # check if root should trigger an event...
    Root <- events$root
    if (is.null(Root)) Root <- 0
    Root <- as.integer(Root)
  } else Root <- 0L

  maxroot <- events$maxroot
  if (is.null(maxroot)) maxroot <- 100  # number of roots to save.
  if (maxroot < 0)
    stop("events$maxroot should be > 0 in events")
  Terminalroot <- events$terminalroot

  if (! is.null(Terminalroot) && is.null(funevent))
    funevent <- function(t,y,p) return(y)  # dummy event function 

  if (is.null(Terminalroot)) 
    Terminalroot <- 0  # at which roots simulation should continue


## ----------------------
## event in a function
## ----------------------
  if (!is.null(funevent)) {
    if (class (funevent) == "CFunc") {
      funevent <- body(funevent)[[2]]
      Type <- 3
    } else if (is.character(funevent)){ 
      if (is.null(dllname))
        stop("'dllname' should be given if 'events$func' is a string")
      if (is.loaded(funevent, PACKAGE = dllname, type = "") ||
      is.loaded(funevent, PACKAGE = dllname, type = "Fortran")) {
        funevent <- getNativeSymbolInfo(funevent, PACKAGE = dllname)$address
      } else
        stop(paste("'events$func' should be loaded ",funevent))
      Type <- 3  
    } else {
      Type <- 2  # SHOULD ALSO CHECK THE FUNCTION if R-function....
#      if (!is.null(dllname))      KARLINE: removed that 02/07/2011
#       stop("'events$func' should be a string, events specified in compiled code if 'dllname' is not NULL")
    }
    if (Root == 0) {
      if (is.null(events$time)) 
        stop("either 'events$time' should be given and contain the times of the events, if 'events$func' is specified and no root function or your solver does not support root functions")
      eventtime <- sort(as.double(events$time)) # Karline: sorted that 4-01-2016

      if (any(!(eventtime %in% times))) {
        warning("Not all event times 'events$time' are in output 'times' so they are automatically included.")
        uniqueTimes <- cleanEventTimes(times, eventtime)
        if (length(uniqueTimes) < length(times))
          warning("Some time steps were very close to events - only the event times are used in these cases.")
        times <- sort(c(uniqueTimes, eventtime))
      }
    } else eventtime <- min(times) - 1  # never reached....
      return (list (Time = eventtime, SVar = NULL, Value = NULL,
        Method = NULL, Type = as.integer(Type), func = funevent,
        Rootsave = as.integer(maxroot), Root = Root,
        Terminalroot = as.integer(Terminalroot), newTimes = times))    # added newTimes - Karline 4-01-2016

  }
## ----------------------
## event as a data series
## ----------------------
  eventdata <- events$data
  if (is.matrix(eventdata)) eventdata <- as.data.frame(eventdata)

  if (ncol(eventdata) < 3)
    stop("'event' should have at least 3 columns: state variable, time, value")

  if (!is.data.frame(eventdata))
    stop("'event' should be a data.frame with 3(4) columns: state variable, time, value, (method)")
    
  ## this should make check < 3 columns obsolete
  evtcols <-  c("var", "time", "value", "method")
  if (!all(evtcols %in% names(eventdata)))
    stop("structure of events does not match specification, see help('events')")
  
  ## make sure that event data frame has correct order
  eventdata <- eventdata[evtcols]

## variables, 1st column should be present
  if (is.factor(eventdata[,1]))
    eventdata[,1] <- as.character(eventdata[,1])

  if (is.character(eventdata[,1]))  {
    vv <- match(eventdata[,1], vars)
  if (is.character(eventdata[,1]))  {
    vv <- match(eventdata[,1],vars)
    if (any(is.na(vv)))
      stop("unknown state variable in 'event': ", paste(eventdata[,1][which(is.na(vv))], ","))
    eventdata[,1] <- vv
  } else if (max(eventdata[,1]) > length(vars))
      stop("unknown state variable in 'event': ", paste(eventdata[,1][which(is.na(vv))],","))
    eventdata[,1] <- vv
  } else if (max(eventdata[,1])>length(vars))
      stop("too many state variables in 'event'; should be < ", paste(length(vars)))

## 2nd and 3rd columns should be numeric
  if (!is.numeric(eventdata[,2]))
      stop("times in 'event', 2nd column should be numeric")

  if (!is.numeric(eventdata[,3]))
      stop("values in 'event', 3rd column should be numeric")

## Times in 'event' should be embraced by 'times'
  rt <- range(times)
  ii <- c(which(eventdata[,2] < rt[1]), which(eventdata[,2] > rt[2]))
  if (length(ii) > 0) 
    eventdata <- eventdata [-ii,]
  if (any(!(eventdata[,2] %in% times))) {
        warning("Not all event times 'events$times' were in output 'times' so they are automatically included.")
        uniqueTimes <- cleanEventTimes(times, eventdata[,2])
        if (length(uniqueTimes) < length(times))
          warning("Some time steps were very close to events - only the event times are used in these cases.")
        times <- sort(c(uniqueTimes, eventdata[,2]))
      }  

  if (any(!(eventdata[,2] %in% times))) {
    warning("Not all event times 'events$times' where in output 'times' so they are automatically included.")
    uniqueTimes <- cleanEventTimes(times, eventdata[,2])
    if (length(uniqueTimes) < length(times))
      warning("Some time steps were very close to events - only the event times are used in these cases.")
    times <- sort(c(uniqueTimes, eventdata[,2]))
  }  


## 4th column: method; if not available: "replace" = method 1 - to date: 3 methods
  if (ncol(eventdata) ==3)
    eventdata$method <- rep(1,nrow(eventdata))
  else if (is.numeric(eventdata[,4])) {
    if (max(eventdata[,4]) > 3 | min(eventdata[,4]) < 1)
      stop("unknown method in 'event': should be >0 and < 4") 
  } else {
    vv <- charmatch(eventdata[,4],c("replace","add","multiply"))
    if (any(is.na(vv)))
      stop("unknown method in 'event': ", paste(eventdata[,3][which(is.na(vv))],","),
        " should be one of 'replace', 'add', 'multiply'")
    eventdata$method <- vv
  }

## Check the other events elements (see optim code)
  con <- list(ties = "notordered", time = NULL, data = NULL, func = NULL, root = NULL)
  nmsC <- names(con)
  con[(namc <- names(events))] <- events
  if (length(noNms <- namc[!namc %in% nmsC]) > 0)
     warning("unknown names in events: ", paste(noNms, collapse = ", "))

## Check what needs to be done in case the time series is not "ordered"

  if (!identical(con$ties, "ordered")) { # see approx code

## first order with respect to time (2nd col), then to variable (1st col)
    if(length(x <- unique(eventdata[,1:2])) < nrow(eventdata)){
      ties <- mean
      if (missing(ties))
        warning("collapsing to unique 'x' values")
      eventdata <- aggregate(eventdata[,c(3, 4)], eventdata[,c(1, 2)], ties)
         ties <- mean
         if (missing(ties))
           warning("collapsing to unique 'x' values")
          eventdata <- aggregate(eventdata[,c(3,4)], eventdata[,c(1,2)], ties)
    }
  }

  return (list (Time = as.double(eventdata[,2]), SVar = as.integer(eventdata[,1]),
    Value = as.double(eventdata[,3]), Method = as.integer(eventdata[,4]),
    Rootsave = as.integer(maxroot),
    Type = 1L, Root = Root,
    Terminalroot = as.integer(Terminalroot),
    newTimes = times))
}


