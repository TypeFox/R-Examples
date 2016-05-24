### ============================================================================
### Chekc forcing function data set
### ============================================================================


checkforcings <- function (forcings,times,dllname,
                           initforc,verbose,fcontrol=list()) {


## Check the names of the initialiser function

 if (is.null(initforc))
   stop(paste("initforc should be loaded if there are forcing functions ",initforc))

 if (is.compiled(initforc)) { 
      if (class(initforc) == "CFunc")
        ModelForc <- body(initforc)[[2]]
 else if (is.loaded(initforc, PACKAGE = dllname, type = "") ||
     is.loaded(initforc, PACKAGE = dllname, type = "Fortran")) 
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

  # add method (linear/constant) to imat
  return(list(tmat=tmat,fmat=fmat,imat=c(imat,method),ModelForc=ModelForc))
}

