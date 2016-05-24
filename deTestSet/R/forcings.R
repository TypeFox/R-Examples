### ============================================================================
### Check forcing function data set, event inputs and time-lag input
### ============================================================================


checkforcings <- function (forcings, times, dllname,
                           initforc, verbose, fcontrol = list()) {


## Check the names of the initialiser function

 if (is.null(initforc))
   stop(paste("initforc should be loaded if there are forcing functions ",initforc))
 if (is.loaded(initforc, PACKAGE = dllname, type = "") ||
     is.loaded(initforc, PACKAGE = dllname, type = "Fortran")) {
 ModelForc <- getNativeSymbolInfo(initforc, PACKAGE = dllname)$address
 } else
   stop(paste("initforc should be loaded if there are forcing functions ",initforc))

## Check the type of the forcing function data series

  if (is.data.frame(forcings)) forcings <- list(a=forcings)
  if (! is.list(forcings)) forcings <- list(a=forcings)
  nf <- length(forcings)
  #1 check if each forcing function consists of a 2-columned matrix
  for (i in 1:nf) {
    if (ncol(forcings[[i]]) != 2)
      stop("forcing function data sets should consist of two-colum matrix")
  }

## Check the control elements (see optim code)

  con <- list(method="linear", rule = 2, f = 0, ties = "ordered")
  nmsC <- names(con)
  con[(namc <- names(fcontrol))] <- fcontrol
  if (length(noNms <- namc[!namc %in% nmsC]) > 0)
     warning("unknown names in fcontrol: ", paste(noNms, collapse = ", "))

  method <- pmatch(con$method, c("linear", "constant"))
    if (is.na(method))
        stop("invalid interpolation method for forcing functions")
  # 1 if linear, 2 if constant...

## Check the timespan of the forcing function data series

  # time span of forcing function data sets should embrace simulation time...
  # although extrapolation is allowed if con$rule = 2 (the default)
  r_t <- range(times)

  for (i in 1:nf) {
    r_f <- range(forcings[[i]][,1])   # time range of this forcing function

    if (r_f[1] > r_t[1]) {
      if (con$rule == 2) {
        mint <- c(r_t[1],forcings[[i]][1,2] )
        forcings[[i]] <- rbind(mint,forcings[[i]])
        if(verbose)
          warning(paste("extrapolating forcing function data sets to first timepoint",i))
       } else
         stop(paste("extrapolating forcing function data sets to first timepoint",i))
    }

    nr   <- nrow(forcings[[i]])
    if (r_f[2] < r_t[2]) {
      if (con$rule == 2) {
        maxt <- c(r_t[2],forcings[[i]][nr,2] )
        forcings[[i]] <- rbind(forcings[[i]],maxt)
        if(verbose)
          warning(paste("extrapolating forcing function data sets to last timepoint",i))
      } else
        stop(paste("extrapolating forcing function data sets to last timepoint",i))
    }

  }

## Check what needs to be done in case the time series is not "ordered"

   if (!identical(con$ties, "ordered")) { # see approx code

     for (i in 1:nf) {

       x <- forcings[[i]][,1]
       nx <- length(x)
       if (length(ux <- unique(x)) < nx) {  # there are non-unique values
         y <- forcings[[i]][,2]
         ties <- con$tiesn
         if (missing(ties))
           warning("collapsing to unique 'x' values")
          y <- as.vector(tapply(y, x, ties))
          x <- sort(ux)
          forcings[[i]] <- cbind(x, y)

       } else {                             # values are unique, but need sorting
          y <- forcings[[i]][,2]
          o <- order(x)
          x <- x[o]
          y <- y[o]
          forcings[[i]] <-  cbind(x,y)
       }
    } # i
  }

## In case the interpolation is of type "constant" and f not equal to 0
## convert y-series, so that always the left value is taken
  if (method == 2 & con$f != 0) {
     for (i in 1:nf) {
       y <- forcings[[i]][,2]
       YY <- c(y,y[length(y)])[-1]
       forcings[[i]][,2] <- (1-con$f)*y + con$f*YY
     }
  }
## all forcings in one vector; adding index to start/end

  fmat <- tmat <- NULL
  imat <- rep(1,nf+1)

  for (i in 1:nf) {
    tmat <- c(tmat, forcings[[i]][,1])
    fmat <- c(fmat, forcings[[i]][,2])
    imat[i+1]<-imat[i]+nrow(forcings[[i]])
  }

  storage.mode(tmat) <- storage.mode(fmat) <- "double"
  storage.mode(imat) <- "integer"

  # DIRTY trick not to inflate the number of arguments:
  # add method (linear/constant) to imat
  return(list(tmat=tmat,fmat=fmat,imat=c(imat,method),ModelForc=ModelForc))
}

### ============================================================================
### Check events data set
### ============================================================================

checkevents <- function (events, times, vars, dllname, root = FALSE) {

  if (is.null(events)) return(list())
  if (is.null(events$data) && is.null(events$func)) return(list())
  # only effective if lsodar, lsode,... "root" triggers an event, does not stop
  if (root) {  # check if root should trigger an event...
    Root <- events$root
    if (is.null(Root)) Root <- 0
    Root <- as.integer(Root)
  } else Root <- as.integer(0)

  Rootsave <- events$maxroot
  if (is.null(Rootsave)) Rootsave <- 100  # number of roots to save.
  if (Rootsave < 0)
    stop("events$Rootsave should be > 0 in events")

  funevent <- events$func
  if (!is.null(funevent)) {
    if (is.character(funevent)){
     if (is.null(dllname))
       stop("'dllname' should be given if 'events$func' is a string")
     if (is.loaded(funevent, PACKAGE = dllname, type = "") ||
     is.loaded(funevent, PACKAGE = dllname, type = "Fortran")) {
       funevent <- getNativeSymbolInfo(funevent, PACKAGE = dllname)$address
     } else
       stop(paste("'events$func' should be loaded ",funevent))
     Type = 3
    } else {
      Type = 2  # SHOULD ALSO CHECK THE FUNCTION if R-function....
      if (!is.null(dllname))
       stop("'events$func' should be a string, events specified in compiled code if 'dllname' is not NULL")
    }
    if (Root == 0) {
      if (is.null(events$time))
        stop("'events$time' should be given and contain the times of the events, if 'events$func' is specified and no root function")
      Time <- as.double(events$time)
      # Karline: added this extra check ....
      if (prod(Time %in% times) != 1)
        stop ("Not all event times 'events$times' are in output 'times'; include event$times in 'times'")
    } else Time <- min(times) - 1  # never reached....
      return (list (Time = Time, SVar = NULL, Value = NULL,
        Method = NULL, Type = as.integer(Type), func = funevent,
        Rootsave = as.integer(Rootsave), Root = Root))

  }  ## Check the event data series
  event <- events$data
  if (is.matrix(event)) event <- as.data.frame(event)

  if (ncol(event) < 3)
    stop("'event' should have at least 3 columns: state variable, time, value")

  if (!is.data.frame(event))
    stop("'event' should be a data.frame with 3(4) columns: state variable, time, value, (method)")

## variables, 1st column should be present
  if (is.factor(event[,1]))
    event[,1] <- as.character(event[,1])

  if (is.character(event[,1]))  {
    vv <- match(event[,1],vars)
    if (any(is.na(vv)))
      stop("unknown state variable in 'event': ", paste(event[,1][which(is.na(vv))],","))
    event[,1] <- vv
  } else if (max(event[,1])>length(vars))
      stop("too many state variables in 'event'; should be < ", paste(length(vars)))

## 2nd and 3rd columns should be numeric
  if (!is.numeric(event[,2]))
      stop("times in 'event', 2nd column should be numeric")

  if (!is.numeric(event[,3]))
      stop("values in 'event', 3rd column should be numeric")

## Times in 'event' should be embraced by 'times'
  rt <- range(times)
  ii <- c(which(event[,2]<rt[1]),which(event[,2]>rt[2]))
  if (length(ii) > 0)
    event <- event [-ii,]
  if (prod(event[,2] %in% times) != 1)
    stop ("Not all event times 'events$times' are in output 'times'; include event$times in 'times'")


## 4th column: method; if not available: "replace" = method 1 - to date: 3 methods
  if (ncol(event) ==3)
    event$method <- rep(1,nrow(event))
  else if (is.numeric(event[,4])) {
      if (max(event[,4]) > 3 | min(event[,4]) < 1)
        stop("unknown method in 'event': should be >0 and < 4")
  } else {
    vv <- charmatch(event[,4],c("replace","add","multiply"))
    if (any(is.na(vv)))
      stop("unknown method in 'event': ", paste(event[,3][which(is.na(vv))],","),
        " should be one of 'replace', 'add', 'multiply'")
    event$method <- vv
  }

## Check the other events elements (see optim code)
  con <- list(ties = "notordered", time = NULL, data=NULL, func = NULL, root = NULL)
  nmsC <- names(con)
  con[(namc <- names(events))] <- events
  if (length(noNms <- namc[!namc %in% nmsC]) > 0)
     warning("unknown names in events: ", paste(noNms, collapse = ", "))

## Check what needs to be done in case the time series is not "ordered"

  if (!identical(con$ties, "ordered")) { # see approx code

## first order with respect to time (2nd col), then to variable (1st col)
    if(length(x<-unique(event[,1:2])) < nrow(event)){
         ties <- mean
         if (missing(ties))
           warning("collapsing to unique 'x' values")
          event <- aggregate(event[,c(3,4)], event[,c(1,2)], ties)
    }
  }

  return (list (Time = as.double(event[,2]), SVar = as.integer(event[,1]),
    Value = as.double(event[,3]), Method = as.integer(event[,4]),
    Rootsave = as.integer(Rootsave),
    Type = as.integer(1), Root = Root))
}

### ============================================================================
### Check timelags data set - also passes "dllname" now  (not yet used)
### ============================================================================

checklags <- function (lags, dllname) {
  if (!is.null(lags)) {
    lags$islag = as.integer(1)
    if (is.null(lags$mxhist))
       lags$mxhist <- 1e4
    if (lags$mxhist <1)
      lags$mxhist <- 1e4
    lags$mxhist<-as.integer(lags$mxhist)
    if (is.null(lags$interpol))   # 1= hermitian, 2 = higher order interpolation
       lags$interpol <- 1
    lags$interpol<-as.integer(lags$interpol)
    lags$isfun <- as.integer(0)
  } else
    lags$islag=as.integer(0)
  return(lags)
}




lagvalue <- function (t, nr = NULL)
{
    if (is.null(nr))
        nr <- 0
    out <- .Call("getLagValue", t = t, PACKAGE = "deTestSet", as.integer(nr))
    return(out)
}
lagderiv <- function (t, nr = NULL)
{
    if (is.null(nr))
        nr <- 0
    out <- .Call("getLagDeriv", t = t, PACKAGE = "deTestSet", as.integer(nr))
    return(out)
}
