


.RandomFields.env <- new.env()

sleep.milli <- function(milli) {
  .C("sleepMilli", as.integer(milli))
  invisible(NULL)
}
	
hostname<-function(){.C("hostname", h=paste(seq(0,0,l=100), collapse=""),
                        as.integer(100), PACKAGE="RandomFields")$h}

pid <- function() {.C("pid", i=integer(1), PACKAGE="RandomFields")$i}


FileExists <- function(file, printlevel=RFoptions()$general$printlevel) {
    ## for parallel simulation studies: the same data output file should not
  ## be created twice. So:
  ## 1. if file exists then assume another process has done the work already
  ## 2. if file.lock existss then assume another process is doing the work
  ## 3.a. otherwise create file.lock to show other processes that the process
  ##      will do the work
  ## 3.b. check if another process has started with the same work at the same
  ##      time it may happen that in case of simulatenous creation of file.lock
  ##      no process will do the work...(then the lock file will rest.)
  lock.ext <- ".lock";
  if (file.exists(file)) { #1.
    if (printlevel>=PL_ERRORS ) cat("'", file, "' already exists.\n");
    return(1)
  } else { 
    LockFile <- paste(file, lock.ext, sep="")
    if (file.exists(LockFile)) { #2.
      if (printlevel>=PL_ERRORS ) cat("'",file,"' is locked.\n");
      return(2);
    }
    PID <- pid();
    write(file=LockFile,c(PID,hostname()),ncolumns=2,append=TRUE); #3.a.
    Pid <- matrix(scan(LockFile,what=character(0), quiet=TRUE),nrow=2)
    if ((sum(Pid[1,]==PID)!=1) || (sum(Pid[1,]>PID)>0)){ #3.b.
      if (printlevel>PL_ERRORS )
        cat("Lock file of '", file, "' is knocked out.\n");
      return(3);
    }
  }
  return(0);
}

LockRemove <- function(file) {
  ## removes auxiliary files created by FileExists
  lock.ext <- ".lock";
  file.remove(paste(file, lock.ext, sep=""))
}




vectordist <- function(x, diag=FALSE) {
  storage.mode(x) <- "double"
  res <- .Call("vectordist", t(x), diag)
  dimnames(res) <- list(dimnames(x)[[2]], NULL)
  return(t(res));
}



my.legend <- function(lu.x, lu.y, zlim, col, cex=1, ...) {
  ## uses already the legend code of R-1.3.0
  cn <- length(col)
  if (cn < 43) {
    col <- rep(col, each=ceiling(43 / cn))
    cn <- length(col)
  }
  filler <- vector("character", length=(cn-3)/2)
  legend(lu.x, lu.y, y.intersp=0.03, x.intersp=0.1, 
         legend=c(format(zlim[2], dig=2), filler,
             format(mean(zlim), dig=2), filler,
             format(zlim[1], dig=2)),
         lty=1, col=rev(col),cex=cex, ...)
}
