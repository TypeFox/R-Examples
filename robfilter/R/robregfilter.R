# robreg.filter - robust time series filters based on robust regression estimators                                         #
#                                                                                                                          #
# Authors:                                                                                                                 #
#   Prof. Dr. Roland Fried         <fried@statistik.uni-dortmund.de>                                                       #
#   Dipl.-Stat. Karen Schettlinger <schettlinger@statistik.uni-dortmund.de>                                                #

# FIXME: h= also subset size for the LQD???
# Bug: If one window contains only missing values the programme crashes!

robreg.filter <- function(y,                    # input time series data (numeric vector or ts object)                                
                          width,                # window width                                                                        
                          method="all",          # filtering (regression) methods                                                      
                          h = floor(width/2)+1, # subset size / trimming fraction for LTS                                             
                          minNonNAs=5,          # required minimal number of non-missing values in each window                        
                          online=FALSE,         # indicator for estimation at the rightmost point or in the centre of each time window
                          extrapolate=TRUE      # indicator for extrapolation of the estimations to the edges                         
                          ) {

  ## save the name of the input time series
  ts.name <- deparse(substitute(y))

  ## Stopping rules and error messages

  # Errors concerning the input time series
  if(is.ts(y)){y <- as.numeric(y)} # coerce a time series object into a numeric vector
  if (!is.null(dim(y))){
    stop("\n        argument 'y' must be a one dimensional numeric vector or time series object")
  }
  if ( any(is.infinite(y), na.rm=TRUE) ){
    stop("argument 'y' contains Inf or -Inf values")
  }
  # length of the time series
  N <- length(y)

  # Errors concerning the filter methods
  all.methods <- c("LQD", "RM", "LMS", "LTS", "DR", "MED") # character string of all possible methods  
  # character string of chosen methods (method.names):
  if( any(method == "all") ){
    method.names <- all.methods
  } else {
    method.names <- method
  }
  n.met  <- length(method.names)
  if (!all(method.names %in% all.methods)){
    stop("invalid specification of 'method': possible values are ", paste(all.methods, collapse=", "))
  }

  # Errors concerning the window widths and minimum number of non-missing values
  if (missing(width)){
    stop("argument 'width' is missing with no default")
  }
  #   [ (          no integer                  ) | (width < 5) | (      even number                    ) ] & [!online] | (identical(all.equal(width%%2,0),TRUE)) ) 
#  if( ((!identical(all.equal(width%%1,0),TRUE)) | (width < 5)) | (identical(all.equal(width%%2,0),TRUE)) & (!online) ){ 
#    stop("argument 'width' must be an odd positive integer >= 5 for offline estimation ('online=FALSE')")
#  }

#  #   [ (          no integer                  ) | (width < 5) ] & [online]
#  if( ( (!identical(all.equal(width%%1,0),TRUE)) | (width < 5) ) & (online) ){
#      stop("argument 'width' must be a positive integer >= 5")
#  }
#  #   [ (          no integer                  ) | (width < 5) ] & [online]
  if( ( (!identical(all.equal(width%%1,0),TRUE)) | (width < 5) ) ){
      stop("argument 'width' must be a positive integer >= 5")
  }
							  
							  
### FIXME: minNonNAS >= 0 !!!  (for all methods!)
  if ( (!identical(all.equal(minNonNAs%%1,0),TRUE)) | (minNonNAs < 5) ){
    stop("argument 'minNonNAs' must be an integer >= 5")
  }
  if (width < minNonNAs)  { 
    stop(paste("argument 'width' cannot be smaller than 'minNonNAs' (=",minNonNAs,")",sep=""))
  }
  if (N < width){
    stop("argument 'y' must contain at least 'width' elements")
  }


#### FIXME (in C++ code): Too many subsequent missing values in the time series cause a crash
#  wNAs <- vector(N-width+1)
#  for( t in 1:(N-width+1) ){
#    wNAs[t] <- sum(is.na(y[t:(t+width-1)])) 
#  }
#  if( any(wNAs > width-minNonNAs) ) {
#    num.toomanyNAs <- sum(wNAs > width-minNonNAs)
#    if(num.toomanyNAs == 1){
#      stop(paste("1 window of width ",width," contains less than minNonNAs=",minNonNAs," elements\n",sep=""))
#    } else {
#      stop(paste(num.toomanyNAs," windows of width ",width," contain less than minNonNAs=",minNonNAs," elements\n",sep=""))
#    }
#  }
#### 

#### Andersch:
#  longestNAspan <- 0
#  for ( t in 1:N ){
#    s <- t
#    while ( is.na(y[s]) ){
#      s <- s+1          
#    }
#    if (s - t > longestNAspan) {
#      longestNAspan <- s - t
#    }
#  }
#  if (longestNAspan >= width)
#    stop(paste("There is at least one passage of ", longestNAspan, " subsequent NAs, which is longer than width=", width, "\n",sep=""))
####

  
  ## Information message + error message in case of missing values for methods which cannot handle them properly
  if( any(is.na(y)) ){ 
   cat(paste( sum(is.na(y)), "out of", N, "time series values in", ts.name, "are missing. \n"))
### FIXME: check C++ code of other methods than RM and MED for behaviour in the presence of missing values
#  if( any(method.names %in% c("DR")) ) {
#     stop(paste("method DR (dr.filter) can not handle missing values yet"))
#   }
###
  }

  # function call
  as.i <- as.integer
  rr <- .Call("robustRegression", y, as.i(width), method.names, !online, as.i(h), as.logical(extrapolate), as.i(minNonNAs), NAOK=TRUE)
  if (rr[[1]] == -1){
  	print("An error occured during calculation. The resulting data could be inconsistent");
  }
  r <- rr[[2]];
  # output results
  level <- matrix(ncol=n.met, nrow=N, dimnames=list(1:N, all.methods[which(all.methods %in% method.names)]))
  slope <- matrix(ncol=n.met, nrow=N, dimnames=list(1:N, all.methods[which(all.methods %in% method.names)]))
  for (i in 1:n.met) {
    level[,i] <- r[[i]]$level
    slope[,i] <- r[[i]]$slope
  }

   if ((!all(is.na(level)))&(!all(is.na(slope)))){# die NA's waren in der Eingabe so unguenstig verteilt, dass nur NA geschätzt wurden
        return( structure( list( level=as.data.frame(level), slope=as.data.frame(slope), 
                           y=y, width=width, method=method.names, 
                           minNonNAs=minNonNAs, online=online, extrapolate=extrapolate, ts.name=ts.name),
          class="robreg.filter")
        )
   }
   return (NA)
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default output
print.robreg.filter <- function(x, ...) {

    # Number of filter methods
    n.met   <- dim(x$level)[2]

    S <- summary(x)
    N <- length(x$y)
    start <- ifelse(x$online, x$width, (x$width+1)/2)
    if(all(is.na(x$level[start+c(0:2),]))){ 
        start <- min(which(apply(x$level, 1, function(L){any(!is.na(L))}) )) 
    }

    # Print level and slope output

    if(n.met ==1) {
      level.output <- as.data.frame(c("...",round(x$level[start+c(0:2),],6),"..."))
      slope.output <- as.data.frame(c("...",round(x$slope[start+c(0:2),],6),"..."))
      dimnames(level.output)[[2]] <- names(x$level)
      dimnames(level.output)[[1]] <- c(" ",start+c(0:2),"  ")
      dimnames(slope.output)[[2]] <- names(x$slope)
      dimnames(slope.output)[[1]] <- c(" ",start+c(0:2),"  ")
      vars <- "variable"
    } else {
      level.output <- rbind(rep("...",n.met),round(x$level[start+c(0:2),],6),rep("...",n.met))
      rownames(level.output)[c(1,5)] <- c(" ","  ")
      slope.output <- rbind(rep("...",n.met),round(x$slope[start+c(0:2),],6),rep("...",n.met))
      rownames(slope.output)[c(1,5)] <- c(" ","  ")
      vars <- "variables"
    }
    cat("$level \n") # print level
    print(level.output)
    cat("('",S["level","Class"],"' with ",N," obs. of", S["level","Length"], " ", vars,")\n \n",sep="")
    cat("$slope \n")
    print(slope.output)
    cat("('",S["slope","Class"],"' with ",N," obs. of ", S["slope","Length"], " ", vars,")\n \n",sep="")

}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default plot
plot.robreg.filter <- function(x, ...) {
    # Length of the time series
    N <- length(x$y)
    
    # Names & number of filter methods
    m.names <- names(x$level)
    n.met   <- length(m.names)
    
    # Setting the y-limits
    ylims <- c(min(x$y,min(x$level,na.rm=TRUE),na.rm=TRUE),max(x$y,max(x$level,na.rm=TRUE),na.rm=TRUE))
    xlims <- c(1,N)
    
    # Defining the title
    if(n.met > 1) { 
        t1 <- "Regression Filters"
    } else {
        t1 <- "Regression Filter"
    }
    titel <- ifelse(x$online,paste("Robust Online ",t1,sep=""),paste("Robust ",t1,sep=""))

    # Possible colors
    mcols <- c("red","green3","blue","skyblue2","orange","yellow")

    # Plot
    par(mar=c(4,4,4,7),oma=rep(0,4),mgp=c(2.5,1,0))
    plot(x$y,xlim=xlims,ylim=ylims,type="l",xlab="Time", ylab=x$ts.name, main=titel)
    
    for(i in seq(along=m.names)){
        lines(x$level[[i]],col=mcols[i],lwd=2)
    }

    # Legend
    par(xpd=TRUE,cex=0.8)
#    legend(par("usr")[2],mean(c(par("usr")[3],par("usr")[4])),c(x$ts.name,m.names),xjust=0,yjust=0.5,lty=rep(1,n.met+1),lwd=c(1,rep(2,n.met)),col=c("black",mcols[1:n.met]),bty="n")
    legend(par("usr")[2],mean(c(par("usr")[3],par("usr")[4])),c("Time Series",m.names),xjust=0,yjust=0.5,lty=rep(1,n.met+1),lwd=c(1,rep(2,n.met)),col=c("black",mcols[1:n.met]),bty="n")
    par(xpd=FALSE,cex=1)
}
