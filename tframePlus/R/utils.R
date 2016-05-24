tfI <- function (x, I=rep(TRUE, NCOL(x)), t0=rep(0, NCOL(x))) { 
   nx <- as.matrix(x)
   for (i in seq(NCOL(nx))) if(I[i]) nx[,i] <- cumsum(c(t0[i], nx[,i]))[-1] 
   tframed(nx,  tf=tframe(x))
   }

tfpersp <- function (x, tf=tfspan(x), start=tfstart(tf), end=tfend(tf),
       theta = -30, phi = 15, scale = FALSE, 
       xlab = "Time", ylab = "", zlab = "", 
       aspect= c(0.5, 0.5), #y/time, z/time,
       ticktype="detailed",ltheta = -120, lphi = 15,
       ...) {
    # 
    if (!is.null(start)) x <- tframe::tfwindow(x, start = start, warn = FALSE)
    if (!is.null(end))   x <- tframe::tfwindow(x, end = end, warn = FALSE)    
    tline <- time(x)
    rx <- max(tline, na.rm=TRUE ) - min(tline, na.rm=TRUE)
    rz <- max(x, na.rm=TRUE) - min(x, na.rm=TRUE)
    ry <- NCOL(x) 
    persp(x=tline, y=(1:ncol(x))* aspect[1]*rx/ry, 
       z=x, expand=aspect[2]*rx/rz,
       theta = theta, phi = phi, scale = scale, 
       xlab=xlab, ylab=ylab,
       ticktype=ticktype,
       ltheta=ltheta, lphi=lphi) #, col = fcol, shade = 0.4,...)
    }

changeTSrepresentation <- function(x, newRepresentation){
   if (is.null(newRepresentation)) newRepresentation <- "default"
   if (is.function(newRepresentation))
       return(newRepresentation(x))
   else if (is.character(newRepresentation))  
       if (newRepresentation  == "default"){
          return(x)
	  # above assumes that x is already in default! 
	  # and also that the user or calling routine required zoo if needed. 
	  # Next fails by converting daily from zoo to ts.
          #if(tffrequency(x) %in% c( 1, 4, 12, 2))  return(as.ts(x))
	  #else {
          #   requireNamespace("zoo")
	  #   return(zoo::as.zoo(x))
	  #   }
	  }
       else if (newRepresentation  == "ts") 
          return(as.ts(x))
       else if (newRepresentation  == "zoo"){
          requireNamespace("zoo")
          #require("zoo") # because of bug in zoo <= 1.7-11
	  return(zoo::as.zoo(x))
	  }
       else if (newRepresentation == "timeSeries") {
          if (! "timeSeries" %in% loadedNamespaces())
	      stop("timeSeries package must be attached.")
	  # zoo and base have as.Date(). 
	  # zoo has index() while timeSeries and stats have time().
	  # note that timeSeries does not seem to support start()
	  return(timeSeries::as.timeSeries(x))
          }
       else if (newRepresentation == "tis") {
          if (! "tis" %in% loadedNamespaces())
	      stop("tis package must be attached.")
          return(tis::as.tis(x))
          }
       else {
  	  return(get(newRepresentation)(x))
          }
    else  stop("mode of time series representation ", 
                mode(newRepresentation ), "not recognized.")
   "should not get here"
   }

TSwriteXLS <- function(x, ..., FileName="R.xls", SheetNames=NULL,
               dateHeader="date", verbose = FALSE){
  # consider tempfile() in overwrite case 
  xx <- list(x, ...)
  env <- sys.frame(sys.nframe())
   genSheetNames <- if (is.null(SheetNames)) TRUE else FALSE
  frNames <- paste("seriesData", seq(length(xx)), sep="")
  if (1 == length(dateHeader)) dateHeader <- rep(dateHeader, length(xx))
  for (i in seq(length(xx))) {
    x <- xx[[i]]
    if (genSheetNames) SheetNames <- c(SheetNames, seriesNames(x)[1])
    tm <- time(x)
    y <- floor(time(tm))
    if (frequency(x) == 4) {
       p <- cycle(tm)
       pp <- c("Q1", "Q2", "Q3", "Q4")[p]
       dt <- paste(pp, y)
       seriesData <- data.frame(date=dt, year=y, period=p, quarter=pp, x)
       names(seriesData) <- c(dateHeader[i], "year", "period", "quarter", seriesNames(x))
       }
    else if (frequency(x) == 12) {
       p <- cycle(tm)
       pp <- c("Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")[p]
       dt <- paste(pp, y)
       seriesData <- data.frame(date=dt, year=y, period=p, month=pp, x)
       names(seriesData) <- c(dateHeader[i], "year", "period", "month", seriesNames(x))
       }
    else  { # annual, daily and weekly are freq 1
       seriesData <- data.frame(date=tm, x ) 
       names(seriesData) <- c(dateHeader[i], seriesNames(x))
       }
    assign(frNames[i], seriesData, envir=env)
    }
  if (requireNamespace("WriteXLS", quietly = TRUE) &&
           WriteXLS::testPerl(verbose=FALSE)) {
     rr <- WriteXLS::WriteXLS(frNames, ExcelFileName=FileName,
               SheetNames=SheetNames, verbose=verbose, envir=env) 
     }
  else { # work around with save and transfer ...
     warning("WriteXLS not usable. Writing txt file.")
     r <- try(save(seriesData, file = paste(FileName, ".txt", sep="")),
               silent = TRUE)
     rr <- ! inherits(r, "try-error")
     }
  rr
  }

TSwriteCSV <- function(x, FileName="R.csv",  dateFormat=1, dateHeader="date"){
  #dateFormat  0=none 1="Jan 1969"   2=1969,1   3=1969,"Jan"
  tm <- time(x)
  y <- floor(time(tm))
  if (frequency(x) == 4) {
     p <- cycle(tm)
     pp <- c("Q1", "Q2", "Q3", "Q4")[p]
     dt <- paste(pp, y)
     if ( dateFormat==0){
       seriesData <- data.frame(x)
       names(seriesData) <- c(seriesNames(x))
       }
     else if ( dateFormat==1){
       seriesData <- data.frame(date=dt, x)
       names(seriesData) <- c(dateHeader, seriesNames(x))
       }
     else if ( dateFormat==2){
       seriesData <- data.frame(year=y, period=p, x)
       names(seriesData) <- c("year", "period", seriesNames(x))
       }
     else if ( dateFormat==3){
       seriesData <- data.frame(year=y, quarter=pp, x)
       names(seriesData) <- c("year", "quarter", seriesNames(x))
       }
     else stop("dateFormat not supported.")
     }
  else if (frequency(x) == 12) {
     p <- cycle(tm)
     pp <- c("Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")[p]
     dt <- paste(pp, y)
     if ( dateFormat==0){
       seriesData <- data.frame(x)
       names(seriesData) <- c(seriesNames(x))
       }
     else if ( dateFormat==1){
       seriesData <- data.frame(date=dt, x)
       names(seriesData) <- c(dateHeader, seriesNames(x))
       }
     else if ( dateFormat==2){
       seriesData <- data.frame(year=y, period=p, x)
       names(seriesData) <- c("year",  "period",  seriesNames(x))
       }
     else if ( dateFormat==3){
       seriesData <- data.frame(year=y, month=pp, x)
       names(seriesData) <- c("year",  "month",   seriesNames(x))
       }
     else stop("dateFormat not supported.")
     }
  else  { # annual, daily and weekly are freq 1
     if ( dateFormat==0){
       seriesData <- data.frame(x ) 
       names(seriesData) <- c(seriesNames(x))
       }
     else if ( dateFormat==1){
       seriesData <- data.frame(date=tm, x ) 
       names(seriesData) <- c(dateHeader, seriesNames(x))
       }
     else stop("dateFormat not supported for this data's dates.")
     }
  write.csv(seriesData, file=FileName, row.names = FALSE) 
  }

as.weekly <- function(x, FUN=sum, na.rm=FALSE, foldFrom=end(x), periodicity = 7){
   # To weekly by periodicity groupings backward from foldFrom
   # and  drop any partial week from the beginning
   # 7 for 7 day weeks.  tested only with daily to weekly, 
   #  periodicity <- 1/frequency(x) unfortunately does not work for daily
   #drop <- length(x) %% periodicity
   addst <- periodicity - (Tobs(tframe::tfwindow(x, end= foldFrom)) %% periodicity) 
   adden <- (Tobs(tframe::tfwindow(x, start= foldFrom))-1) %% periodicity 
   x <- tfExpand.zoo(x, add.start=addst, add.end  =adden)
   r <- as.matrix(x)
   #if (drop > 0) r <- r[ -(1:drop),, drop=FALSE]
   C <- NCOL(r)
   R <- NROW(r)/periodicity
   rr <- matrix(NA, R, C)
   for (i in 1:C) rr[,i ] <- apply(matrix(r[,i],periodicity, R),2, FUN=FUN)
   r <- zoo::zoo(rr, foldFrom - (NROW(rr)-1):0 * periodicity)
   if(na.rm) trimNA(r) else r
   }
   
as.quarterly <- function (x, FUN=sum, na.rm=FALSE, ...){
    # convert to quarterly (from monthly only, so far
    #  use aggregate, but shift to match quarters
    if (4 == frequency(x)) return(if(na.rm) trimNA(x) else x) 
    if (12 != frequency(x)) stop("only monthly conversion supported for now.")
    tf <- tframe(x)
    nm <- seriesNames(x)
    x <- tfExpand(x, add.start=(tfstart(tf)[2] %% 3)-1,
                     add.end  =(3 - tfend(tf)[2]) %% 3)
    r <- aggregate(x, nfrequency=4, FUN=FUN, 
        ndeltat=1, ts.eps=getOption("ts.eps"), ...) 
    if(na.rm) trimNA(r) else r
    }

as.annually <- function (x, FUN=sum, na.rm=FALSE, ...){
    # convert to annual (from quarterly or monthly only, so far)
    #  use aggregate, but shift to match years
    # Monthly to annual gives the aggregate through quarterly which
    #  is not exactly correct.
    if (1 == frequency(x)) return(if(na.rm) trimNA(x) else x) 
    if (12 == frequency(x)) x <- as.quarterly(x, FUN=FUN, na.rm=na.rm, ...)
    if (4 != frequency(x))
       stop("currently only quarterly and monthly series are supported." )
    tf <- tframe(x)
    nm <- seriesNames(x)
    x <- tfExpand(x, add.start=(tfstart(tf)[2] %% 4)-1,
                     add.end  =(4 - tfend(tf)[2]) %% 4)
    r <- aggregate(x, nfrequency=1, FUN=FUN, 
        ndeltat=1, ts.eps=getOption("ts.eps"), ...) 
    if(na.rm) trimNA(r) else r
    }


rollAggregate <- function(x, FUN=sum, na.rm=FALSE, aggPeriods=4, ...){
   r <- as.matrix(x)
   N <- Tobs(x)
   rr <- array(NA, c(N, NCOL(r), aggPeriods))
   rr[,,1] <- r
   for (i in 1:(aggPeriods-1)) rr[(1+i):N,,(1+i)] <- r[1:(N-i),]
   rr <- apply(rr,1:2, FUN=FUN, ...)
   rr <- tframed(rr, tf=tframe(x))
   if(na.rm) trimNA(rr) else rr
   }
